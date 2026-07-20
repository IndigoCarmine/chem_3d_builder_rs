//! Holds the molecule currently shown in the 3D panel and drives minimization.
//!
//! OpenBabel stops a minimization for two very different reasons — the geometry
//! converged, or the step budget ran out — and its own API reports both the same
//! way. `openbabel_rs` now separates them (`StopReason`), which is what lets this
//! module keep going until the geometry is *actually* minimized instead of
//! stopping at an arbitrary step count and presenting a half-relaxed structure as
//! finished.
//!
//! A run is therefore a *sequence* of bounded segments rather than one call: the
//! worker minimizes `SEGMENT_STEPS` at a time and repeats while the answer keeps
//! coming back `MaxSteps`. Each segment's frames are streamed to the UI as they
//! are produced, so the animation plays during the run rather than after it, and
//! the segment boundary doubles as the cancellation point — no cancel hook is
//! needed inside OpenBabel itself.
//!
//! While a run is in flight the worker holds OpenBabel's process-global lock for
//! its whole duration, so the UI thread must not touch any OpenBabel API until
//! the molecule comes home. `Running` is the state that encodes that: it owns no
//! `Molecule`, and callers gate their buttons on it. Playback is safe because a
//! frame is plain coordinate data.

use std::collections::VecDeque;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError};

use openbabel::{Algorithm, Minimizer, Molecule, OptStep, StopReason};

use crate::forcefield_kind::FfKind;

/// Steps per segment. Also the cancellation granularity: the worker only notices
/// a stop request between segments, so this trades responsiveness against the
/// cost of restarting conjugate gradients each time.
const SEGMENT_STEPS: u32 = 250;

/// Absolute ceiling across all segments of one run, so a molecule that will not
/// settle cannot spin forever.
const MAX_TOTAL_STEPS: u32 = 20_000;

/// Energy-convergence threshold handed to OpenBabel, in the force field's unit.
const ENERGY_CONVERGENCE: f64 = 1e-6;

/// How a run ended, for the read-out badge.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Outcome {
    /// Still descending — the state every frame mid-run carries.
    Running,
    Converged,
    /// Hit [`MAX_TOTAL_STEPS`] with the geometry still moving.
    StepLimit,
    /// Stopped by the user.
    Cancelled,
}

/// The read-out badge: energy, its unit, cumulative steps, and how the run ended.
///
/// The unit travels with the number because it varies by force field — UFF/GAFF/
/// Ghemical report kJ/mol, MMFF94/MMFF94s kcal/mol — so a bare energy is not
/// interpretable. The outcome travels with it for the same reason: an energy
/// alone cannot say whether the geometry is actually minimized.
#[derive(Clone, Copy)]
pub struct Badge {
    pub energy: f64,
    pub unit: &'static str,
    pub steps: usize,
    pub outcome: Outcome,
}

/// What the worker sends back as it goes.
///
/// `pub(crate)` only because it appears in `MmState::Running`'s field types;
/// nothing outside this module constructs or matches on it.
pub(crate) enum WorkerMsg {
    /// One segment's frames, ready to animate.
    Chunk(Vec<OptStep>),
    /// The run is over; the molecule comes home whether or not it produced anything.
    Done {
        mol: Molecule,
        outcome: Outcome,
        /// True if no segment ever produced a frame — the only shape a
        /// force-field setup failure takes.
        failed: bool,
    },
}

#[derive(Default)]
pub enum MmState {
    #[default]
    Empty,
    Idle {
        mol: Molecule,
        badge: Option<Badge>,
    },
    /// A worker owns the molecule *and* OpenBabel's global lock. Frames it has
    /// sent so far wait in `queue` and play one per repaint.
    Running {
        rx: Receiver<WorkerMsg>,
        cancel: Arc<AtomicBool>,
        queue: VecDeque<OptStep>,
        /// Set once `Done` arrives; the queue still drains before we go idle.
        finished: Option<(Molecule, Outcome, bool)>,
        ff: FfKind,
        /// Steps accumulated by earlier runs — OpenBabel restarts its step count
        /// from zero each run, but the badge counts across runs.
        base: usize,
        /// Steps counted within this run, across segments.
        run_steps: usize,
        /// Last frame step seen; resets each segment, so deltas drive run_steps.
        last_step: u32,
        /// Set by [`MmState::cancel`]. Playback then keeps only the newest frame
        /// of whatever is queued, so stopping is immediate even when the worker
        /// has already run far ahead of the animation.
        cancelled: bool,
        badge: Option<Badge>,
    },
}

/// What one `poll` produced, for the caller to act on.
pub enum PollOutcome {
    /// Nothing in flight.
    Idle,
    /// A run is going; keep repainting so we come back to `try_recv`.
    Waiting,
    /// New geometry to push to the viewport (Ångström).
    Frame(Vec<[f32; 3]>),
    Failed(String),
}

/// Minimize `mol` until it converges, handing each segment's frames to `on_chunk`.
///
/// Runs on the worker thread; tests call it directly. Returns the molecule at its
/// final geometry, how the run ended, and whether it ever produced a frame.
fn run_minimize(
    mut mol: Molecule,
    ff: FfKind,
    steps_per_frame: u32,
    cancel: &AtomicBool,
    mut on_chunk: impl FnMut(Vec<OptStep>),
) -> (Molecule, Outcome, bool) {
    let mut total = 0u32;
    let mut produced_any = false;

    loop {
        if cancel.load(Ordering::Relaxed) {
            return (mol, Outcome::Cancelled, produced_any);
        }

        // A `Minimizer` owns a `Constraints`, which wraps a cxx opaque type and
        // is therefore `!Send` — it cannot be built by the caller and moved here.
        let mut cfg = Minimizer::new(ff.ob_id());
        cfg.algorithm(Algorithm::ConjugateGradients)
            // Never offer L-BFGS: paired with UFF it corrupts the heap in
            // OpenBabel 3.2.1.
            .max_steps(SEGMENT_STEPS.min(MAX_TOTAL_STEPS - total))
            .energy_convergence(ENERGY_CONVERGENCE)
            .steps_per_frame(steps_per_frame);

        let run = mol.minimize(&cfg);
        let stop = run.stop_reason();
        let taken = run.steps_taken();
        let frames: Vec<OptStep> = run.collect();

        produced_any |= !frames.is_empty();
        if !frames.is_empty() {
            on_chunk(frames);
        }
        total += taken;

        match stop {
            StopReason::Converged => return (mol, Outcome::Converged, produced_any),
            // No trajectory at all: unknown or un-setupable force field. The
            // caller turns `produced_any == false` into the user-facing error.
            StopReason::Failed => return (mol, Outcome::StepLimit, produced_any),
            StopReason::MaxSteps if total >= MAX_TOTAL_STEPS => {
                return (mol, Outcome::StepLimit, produced_any);
            }
            // Still descending. The next segment restarts conjugate gradients
            // from the current coordinates — the molecule carries them forward,
            // since the shim writes each frame's geometry back into it — which is
            // a plain CG restart, costing only the lost direction history.
            StopReason::MaxSteps => {}
        }
    }
}

/// Drop every queued frame but the newest, so playback jumps to the latest
/// geometry instead of replaying the backlog.
fn keep_last(queue: &mut VecDeque<OptStep>) {
    if let Some(last) = queue.pop_back() {
        queue.clear();
        queue.push_back(last);
    }
}

impl MmState {
    /// A freshly generated 3D structure, with its starting energy read once.
    pub fn ready(mol: Molecule, ff: FfKind) -> Self {
        let badge = mol.energy(ff.ob_id()).map(|energy| Badge {
            energy,
            unit: ff.energy_unit(),
            steps: 0,
            outcome: Outcome::Running,
        });
        MmState::Idle { mol, badge }
    }

    /// The molecule, when the UI thread is allowed to touch it.
    pub fn mol(&self) -> Option<&Molecule> {
        match self {
            MmState::Idle { mol, .. } => Some(mol),
            MmState::Empty | MmState::Running { .. } => None,
        }
    }

    pub fn badge(&self) -> Option<Badge> {
        match self {
            MmState::Empty => None,
            MmState::Idle { badge, .. } | MmState::Running { badge, .. } => *badge,
        }
    }

    /// True while a worker holds the molecule and the OpenBabel lock.
    pub fn is_running(&self) -> bool {
        matches!(self, MmState::Running { .. })
    }

    /// Hand the molecule to a worker and start minimizing. A no-op unless idle.
    pub fn start(&mut self, ff: FfKind, steps_per_frame: usize, ctx: &egui::Context) {
        let (mol, badge) = match std::mem::replace(self, MmState::Empty) {
            MmState::Idle { mol, badge } => (mol, badge),
            other => {
                *self = other;
                return;
            }
        };

        let (tx, rx) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let spf = steps_per_frame.max(1) as u32;
        let base = badge.map(|b| b.steps).unwrap_or(0);
        let worker_cancel = Arc::clone(&cancel);
        let worker_ctx = ctx.clone();
        std::thread::spawn(move || {
            let chunk_tx = tx.clone();
            let chunk_ctx = worker_ctx.clone();
            let (mol, outcome, produced_any) =
                run_minimize(mol, ff, spf, &worker_cancel, move |frames| {
                    // Err only if the UI dropped the receiver, which is harmless.
                    let _ = chunk_tx.send(WorkerMsg::Chunk(frames));
                    // Wake the UI so it animates during the run, not after it.
                    chunk_ctx.request_repaint();
                });
            let _ = tx.send(WorkerMsg::Done {
                mol,
                outcome,
                failed: !produced_any,
            });
            worker_ctx.request_repaint();
        });

        *self = MmState::Running {
            rx,
            cancel,
            queue: VecDeque::new(),
            finished: None,
            ff,
            base,
            run_steps: 0,
            last_step: 0,
            cancelled: false,
            badge,
        };
    }

    /// Stop a run: signal the worker and cut the animation short. A no-op unless
    /// running.
    ///
    /// Signalling the worker alone is not enough to make the button feel like it
    /// did anything. The worker runs far ahead of playback — a 20 000-step run
    /// finishes computing long before its ~2 500 frames have been drawn one per
    /// repaint — so a stop that only sets the flag leaves the user watching tens
    /// of seconds of animation they just asked to end. Collapsing the queue to
    /// its newest frame ends it now, and keeps the displayed geometry equal to
    /// the molecule's, which is what makes an export match what is on screen.
    pub fn cancel(&mut self) {
        if let MmState::Running {
            cancel,
            queue,
            cancelled,
            ..
        } = self
        {
            cancel.store(true, Ordering::Relaxed);
            *cancelled = true;
            keep_last(queue);
        }
    }

    /// Advance one frame: drain whatever the worker has sent, then emit the next
    /// queued frame.
    pub fn poll(&mut self) -> PollOutcome {
        // Check before taking ownership: a `replace` whose value then fails to
        // match would *drop* the state it pulled out — and this runs on every
        // repaint, so an `Idle` molecule would vanish the moment it went idle.
        if !self.is_running() {
            return PollOutcome::Idle;
        }
        let MmState::Running {
            rx,
            cancel,
            mut queue,
            mut finished,
            ff,
            base,
            run_steps,
            last_step,
            cancelled,
            badge,
        } = std::mem::replace(self, MmState::Empty)
        else {
            unreachable!("just checked is_running");
        };

        // Take everything that has arrived since the last repaint.
        loop {
            match rx.try_recv() {
                Ok(WorkerMsg::Chunk(frames)) => queue.extend(frames),
                Ok(WorkerMsg::Done {
                    mol,
                    outcome,
                    failed,
                }) => finished = Some((mol, outcome, failed)),
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    if finished.is_none() {
                        return PollOutcome::Failed(
                            "最小化スレッドが異常終了しました。".to_string(),
                        );
                    }
                    break;
                }
            }
        }

        // A cancelled run keeps only the newest frame of whatever the worker
        // sent while we were away, so late chunks cannot restart the backlog.
        if cancelled {
            keep_last(&mut queue);
        }

        // Frames first: play out everything received before settling the run, so
        // the animation never skips its tail.
        if let Some(f) = queue.pop_front() {
            let coords = f
                .coordinates
                .iter()
                .map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
                .collect();
            // `f.step` counts within its own segment and restarts at the next
            // one, so accumulate deltas: a step that failed to advance means a
            // fresh segment, whose delta is the whole of `f.step`.
            let delta = if f.step > last_step {
                f.step - last_step
            } else {
                f.step
            };
            let run_steps = run_steps + delta as usize;
            let last_step = f.step;
            let badge = Some(Badge {
                energy: f.energy,
                unit: ff.energy_unit(),
                steps: base + run_steps,
                outcome: Outcome::Running,
            });
            *self = MmState::Running {
                rx,
                cancel,
                queue,
                finished,
                ff,
                base,
                run_steps,
                last_step,
                cancelled,
                badge,
            };
            return PollOutcome::Frame(coords);
        }

        match finished {
            // Worker done and the queue is empty: settle.
            Some((mol, outcome, failed)) => {
                if failed {
                    *self = MmState::Idle { mol, badge };
                    return PollOutcome::Failed(format!(
                        "{} 力場での最小化に失敗しました。",
                        ff.label()
                    ));
                }
                // Stamp how the run ended onto the last frame's readings.
                let badge = badge.map(|b| Badge { outcome, ..b });
                *self = MmState::Idle { mol, badge };
                PollOutcome::Idle
            }
            // Still working, nothing queued yet.
            None => {
                *self = MmState::Running {
                    rx,
                    cancel,
                    queue,
                    finished,
                    ff,
                    base,
                    run_steps,
                    last_step,
                    cancelled,
                    badge,
                };
                PollOutcome::Waiting
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ethanol() -> Molecule {
        let mut mol = Molecule::parse("CCO", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        mol
    }

    /// Collect a whole run synchronously, the way the worker thread would.
    fn drive(mol: Molecule, ff: FfKind, spf: u32) -> (Vec<OptStep>, Outcome, bool) {
        let cancel = AtomicBool::new(false);
        let mut frames = Vec::new();
        let (_mol, outcome, produced) =
            run_minimize(mol, ff, spf, &cancel, |chunk| frames.extend(chunk));
        (frames, outcome, produced)
    }

    #[test]
    fn uff_trajectory_does_not_raise_energy() {
        let mol = ethanol();
        let e0 = mol.energy(FfKind::Uff.ob_id()).expect("initial energy");
        let (frames, _outcome, produced) = drive(mol, FfKind::Uff, 8);
        assert!(produced, "expected some minimization frames");
        let last = frames.last().unwrap();
        assert!(
            last.energy <= e0 + 1e-6,
            "energy should not increase: {e0} -> {}",
            last.energy
        );
        assert!(last.step > 0, "expected some minimization steps");
    }

    /// The regression this whole change exists for.
    ///
    /// This molecule needs several thousand steps to settle — far past the old
    /// hard 500-step budget, which would have stopped here around a tenth of the
    /// way down and reported the result as if it were minimized. The run must
    /// now cross segment boundaries and keep going until it genuinely converges.
    #[test]
    fn a_long_run_continues_past_the_first_segment() {
        let mut mol = Molecule::parse("c1ccccc1C(=O)NC2CCCCC2", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");

        // Each segment's frame steps restart at zero, so a step that fails to
        // advance marks a boundary.
        let cancel = AtomicBool::new(false);
        let mut segments = 1usize;
        let mut last = 0u32;
        let (_mol, outcome, _) = run_minimize(mol, FfKind::Uff, 8, &cancel, |chunk| {
            for f in chunk {
                if f.step <= last {
                    segments += 1;
                }
                last = f.step;
            }
        });

        assert_eq!(outcome, Outcome::Converged, "run ended without converging");
        assert!(
            segments > 1,
            "converged inside one segment — this molecule no longer exercises \
             multi-segment continuation, so the regression is untested"
        );
    }

    /// Cancellation is checked at the segment boundary, so a flag set before the
    /// first segment stops the run before any work happens.
    #[test]
    fn cancel_before_the_first_segment_stops_immediately() {
        let cancel = AtomicBool::new(true);
        let mut frames = Vec::new();
        let (_mol, outcome, produced) =
            run_minimize(ethanol(), FfKind::Uff, 8, &cancel, |c| frames.extend(c));
        assert_eq!(outcome, Outcome::Cancelled);
        assert!(!produced && frames.is_empty(), "cancelled run did work");
    }

    /// Grounds the FfKind line-up: every id we offer must resolve in OpenBabel.
    #[test]
    fn every_forcefield_is_available() {
        let mol = ethanol();
        for k in FfKind::ALL {
            assert!(mol.energy(k.ob_id()).is_some(), "{} unavailable", k.label());
        }
    }
}
