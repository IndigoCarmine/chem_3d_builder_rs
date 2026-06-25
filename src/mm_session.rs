//! Holds the live MM state for the molecule currently shown in the 3D panel,
//! and advances minimization a few steps at a time so the UI can animate it.

use crate::bridge::{self, Mol3D};
use crate::forcefield_kind::FfKind;
use zunda_rs::molecule::Vec3;
use zunda_rs::ForceField;

// Minimizer tuning (steepest descent + backtracking line search), mirroring the
// approach in zunda_rs/examples/benzene_uff.rs.
const GRAD_TOL: f64 = 1.0e-4;
const ENERGY_TOL: f64 = 1.0e-7;
const STEP_SIZE: f64 = 0.05;
const TRUST_RADIUS: f64 = 0.3;

pub struct MmSession {
    ff: Box<dyn ForceField>,
    mol: Mol3D,
    positions: Vec<Vec3>,
    /// Whether minimization is actively running (advanced each frame).
    pub minimizing: bool,
    pub total_steps: usize,
    pub energy: f64,
    pub converged: bool,
}

impl MmSession {
    /// Build a session for `mol`: set up the force field and read the initial energy.
    pub fn new(mol: Mol3D, ff_kind: FfKind) -> Result<Self, String> {
        if mol.is_empty() {
            return Err("構造が空です。まず左側で原子を描いてください。".to_string());
        }
        let mm = bridge::build_mm_molecule(&mol, &mol.positions);
        let mut ff = ff_kind.make();
        ff.setup(&mm)
            .map_err(|e| format!("力場のセットアップに失敗しました: {e}"))?;
        let positions = mm.atom_positions();
        let energy = ff
            .energy(&positions)
            .map_err(|e| format!("エネルギー計算に失敗しました: {e}"))?;
        Ok(Self {
            ff,
            mol,
            positions,
            minimizing: false,
            total_steps: 0,
            energy,
            converged: false,
        })
    }

    pub fn atom_count(&self) -> usize {
        self.mol.elements.len()
    }

    /// Run up to `max_steps` steepest-descent iterations (with a backtracking line
    /// search). Calls the force field through `&dyn ForceField`, so it works for any
    /// of the boxed force fields. Stops `minimizing` on convergence or when stuck.
    pub fn step(&mut self, max_steps: usize) {
        let mut prev_e = match self.ff.energy(&self.positions) {
            Ok(e) => e,
            Err(_) => {
                self.minimizing = false;
                return;
            }
        };

        for _ in 0..max_steps {
            let grad = match self.ff.gradient(&self.positions) {
                Ok(g) => g,
                Err(_) => {
                    self.minimizing = false;
                    break;
                }
            };

            let max_g = grad.iter().map(|g| g.norm()).fold(0.0_f64, f64::max);
            if max_g <= GRAD_TOL {
                self.converged = true;
                self.minimizing = false;
                break;
            }

            // Step downhill along the negative gradient with a backtracking line search.
            let dir: Vec<Vec3> = grad.iter().map(|g| -*g).collect();
            let mut ls = STEP_SIZE;
            let mut e_best = prev_e;
            let mut pos_best = self.positions.clone();
            let mut progress = false;

            for _ in 0..10 {
                let mut trial = pos_best.clone();
                for (p, d) in trial.iter_mut().zip(dir.iter()) {
                    let mut delta = *d * ls;
                    delta.x = delta.x.clamp(-TRUST_RADIUS, TRUST_RADIUS);
                    delta.y = delta.y.clamp(-TRUST_RADIUS, TRUST_RADIUS);
                    delta.z = delta.z.clamp(-TRUST_RADIUS, TRUST_RADIUS);
                    *p += delta;
                }
                match self.ff.energy(&trial) {
                    Ok(e_trial) if e_trial < e_best => {
                        e_best = e_trial;
                        pos_best = trial;
                        progress = true;
                        ls = (ls * 2.15).min(1.0);
                    }
                    Ok(_) => {
                        ls *= 0.1;
                        if ls < 1.0e-12 {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }

            self.positions = pos_best;
            self.total_steps += 1;
            self.energy = e_best;

            if !progress {
                // Line search made no progress → at (or near) a minimum.
                self.minimizing = false;
                break;
            }
            if (prev_e - e_best).abs() <= ENERGY_TOL {
                self.converged = true;
                self.minimizing = false;
                break;
            }
            prev_e = e_best;
        }
    }

    /// Snapshot of the current geometry as a viewer molecule.
    pub fn viewer_molecule(&self) -> moleucle_3dview_rs::Molecule {
        bridge::to_viewer_molecule(&self.mol, &self.positions)
    }
}
