//! Holds the molecule currently shown in the 3D panel and drives minimization.
//!
//! OpenBabel's minimizer runs to completion in one call and hands back the whole
//! trajectory, so the run happens on a worker thread (`Molecule` is `Send`) and
//! the UI plays the recorded frames back one per repaint — which is what keeps
//! the animation this app has always had.
//!
//! While a run is in flight the worker holds OpenBabel's process-global lock for
//! its whole duration, so the UI thread must not touch any OpenBabel API until
//! the molecule comes home. `Running` is the state that encodes that: it owns no
//! `Molecule`, and callers gate their buttons on it. Playback is safe because a
//! frame is plain coordinate data.

use std::sync::mpsc::{self, Receiver, TryRecvError};

use openbabel::{Algorithm, Minimizer, Molecule, OptStep};

use crate::forcefield_kind::FfKind;

/// Step budget for one press of the minimize button.
const MAX_STEPS: u32 = 500;

/// `(energy, unit, cumulative steps)` for the read-out badge. The unit travels
/// with the number because it varies by force field — UFF/GAFF/Ghemical report
/// kJ/mol, MMFF94/MMFF94s kcal/mol — so a bare energy is not interpretable.
type Badge = Option<(f64, &'static str, usize)>;
/// The molecule always comes back, whether or not the run produced anything.
type WorkerMsg = (Molecule, Vec<OptStep>);

#[derive(Default)]
pub enum MmState {
    #[default]
    Empty,
    Idle {
        mol: Molecule,
        badge: Badge,
    },
    /// A worker owns the molecule *and* OpenBabel's global lock.
    Running {
        rx: Receiver<WorkerMsg>,
        ff: FfKind,
        badge: Badge,
    },
    Playing {
        mol: Molecule,
        traj: Vec<OptStep>,
        frame: usize,
        /// Steps accumulated by earlier runs — OpenBabel restarts `OptStep.step`
        /// from zero each time, but the badge counts across runs.
        base: usize,
        /// Unit of the force field this trajectory was produced with.
        unit: &'static str,
        badge: Badge,
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

/// Minimize `mol` and return it with its trajectory. Runs on the worker thread;
/// tests call it directly.
fn run_minimize(mut mol: Molecule, ff: FfKind, steps_per_frame: u32) -> WorkerMsg {
    // A `Minimizer` owns a `Constraints`, which wraps a cxx opaque type and is
    // therefore `!Send` — it cannot be built by the caller and moved here.
    let mut cfg = Minimizer::new(ff.ob_id());
    cfg.algorithm(Algorithm::ConjugateGradients)
        .max_steps(MAX_STEPS)
        .steps_per_frame(steps_per_frame);
    // Never offer L-BFGS: paired with UFF it corrupts the heap in OpenBabel 3.2.1.
    let traj = mol.minimize(&cfg).collect();
    (mol, traj)
}

impl MmState {
    /// A freshly generated 3D structure, with its starting energy read once.
    pub fn ready(mol: Molecule, ff: FfKind) -> Self {
        let badge = mol.energy(ff.ob_id()).map(|e| (e, ff.energy_unit(), 0));
        MmState::Idle { mol, badge }
    }

    /// The molecule, when the UI thread is allowed to touch it.
    pub fn mol(&self) -> Option<&Molecule> {
        match self {
            MmState::Idle { mol, .. } | MmState::Playing { mol, .. } => Some(mol),
            MmState::Empty | MmState::Running { .. } => None,
        }
    }

    pub fn badge(&self) -> Badge {
        match self {
            MmState::Empty => None,
            MmState::Idle { badge, .. }
            | MmState::Running { badge, .. }
            | MmState::Playing { badge, .. } => *badge,
        }
    }

    /// True while a worker holds the molecule and the OpenBabel lock.
    pub fn is_running(&self) -> bool {
        matches!(self, MmState::Running { .. })
    }

    pub fn is_playing(&self) -> bool {
        matches!(self, MmState::Playing { .. })
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
        let spf = steps_per_frame.max(1) as u32;
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let msg = run_minimize(mol, ff, spf);
            // Err only if the UI dropped the receiver, which is harmless.
            let _ = tx.send(msg);
            // Wake the UI even if it went to sleep waiting.
            ctx.request_repaint();
        });

        *self = MmState::Running { rx, ff, badge };
    }

    /// Stop playback where it stands. A no-op unless playing.
    pub fn stop(&mut self) {
        *self = match std::mem::replace(self, MmState::Empty) {
            MmState::Playing {
                mut mol,
                traj,
                frame,
                badge,
                ..
            } => {
                // `minimize` left the molecule at the *final* geometry, but the
                // screen shows `frame`. Rewind so what you export is what you see.
                if let Some(f) = traj.get(frame.saturating_sub(1)) {
                    let flat: Vec<f64> = f.coordinates.iter().flatten().copied().collect();
                    mol.set_coordinates(&flat);
                }
                MmState::Idle { mol, badge }
            }
            other => other,
        };
    }

    /// Advance one frame: collect a finished run, or emit the next recorded frame.
    pub fn poll(&mut self) -> PollOutcome {
        match std::mem::replace(self, MmState::Empty) {
            MmState::Empty => PollOutcome::Idle,
            idle @ MmState::Idle { .. } => {
                *self = idle;
                PollOutcome::Idle
            }

            MmState::Running { rx, ff, badge } => match rx.try_recv() {
                Err(TryRecvError::Empty) => {
                    *self = MmState::Running { rx, ff, badge };
                    PollOutcome::Waiting
                }
                Err(TryRecvError::Disconnected) => {
                    PollOutcome::Failed("最小化スレッドが異常終了しました。".to_string())
                }
                // An empty trajectory is the only failure signal OpenBabel gives
                // us — `minimize` reports an unknown or unsetupable force field
                // by producing no frames.
                Ok((mol, traj)) if traj.is_empty() => {
                    *self = MmState::Idle { mol, badge };
                    PollOutcome::Failed(format!("{} 力場での最小化に失敗しました。", ff.label()))
                }
                Ok((mol, traj)) => {
                    let base = badge.map(|(_, _, steps)| steps).unwrap_or(0);
                    *self = MmState::Playing {
                        mol,
                        traj,
                        frame: 0,
                        base,
                        unit: ff.energy_unit(),
                        badge,
                    };
                    self.poll() // emit the first frame now rather than next repaint
                }
            },

            MmState::Playing {
                mol,
                traj,
                frame,
                base,
                unit,
                badge,
            } => match traj.get(frame) {
                Some(f) => {
                    let coords = f
                        .coordinates
                        .iter()
                        .map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
                        .collect();
                    let badge = Some((f.energy, unit, base + f.step as usize));
                    *self = MmState::Playing {
                        mol,
                        traj,
                        frame: frame + 1,
                        base,
                        unit,
                        badge,
                    };
                    PollOutcome::Frame(coords)
                }
                None => {
                    // Played out. The molecule already holds this final geometry.
                    *self = MmState::Idle { mol, badge };
                    PollOutcome::Idle
                }
            },
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

    #[test]
    fn uff_trajectory_does_not_raise_energy() {
        let mol = ethanol();
        let e0 = mol.energy(FfKind::Uff.ob_id()).expect("initial energy");
        let (_mol, traj) = run_minimize(mol, FfKind::Uff, 8);
        assert!(!traj.is_empty(), "expected some minimization frames");
        let last = traj.last().unwrap();
        assert!(
            last.energy <= e0 + 1e-6,
            "energy should not increase: {e0} -> {}",
            last.energy
        );
        assert!(last.step > 0, "expected some minimization steps");
    }

    /// Grounds the FfKind line-up: every id we offer must resolve in OpenBabel.
    #[test]
    fn every_forcefield_is_available() {
        let mol = ethanol();
        for k in FfKind::ALL {
            assert!(
                mol.energy(k.ob_id()).is_some(),
                "{} unavailable",
                k.label()
            );
        }
    }
}
