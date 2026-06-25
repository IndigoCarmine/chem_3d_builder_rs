//! Bridge layer between the three crates.
//!
//! The 2D editor (`chembuider_rs`) produces a flat structure with **no explicit
//! hydrogens** and z = 0. The MM engine (`zunda_rs`) only *minimizes* existing 3D
//! coordinates — it neither adds hydrogens nor generates 3D from 2D. The 3D viewer
//! (`moleucle_3dview_rs`) wants nm units (Å × 0.1).
//!
//! So this module does the missing 2D → 3D work:
//!   1. add implicit hydrogens from a simple valence model,
//!   2. lift the 2D layout into 3D with a small random z-jitter on every atom
//!      (a flat, symmetric start is a stationary point for the force field — the
//!      jitter breaks planarity so MM can relax into a real 3D geometry),
//! then converts to/from the zunda and viewer molecule types.
//!
//! Limitations (acceptable for a first version, flagged for later work):
//!   - single conformer from a planar layout + MM relaxation; no distance-geometry
//!     embedding (ETKDG-style), so flexible/large molecules may land in a poor local
//!     minimum.
//!   - formal-charge handling in the valence model is simplified.

use std::collections::HashMap;

use chembuider_rs::BondOrder;

/// Ångström. Initial bond length used when placing added hydrogens.
const H_BOND_LEN: f64 = 1.0;
/// Ångström. Magnitude of the z-jitter applied to break planarity.
const Z_JITTER: f64 = 0.3;

/// A flat, crate-agnostic 3D molecule used as the common currency between
/// `chembuider_rs`, `zunda_rs` and `moleucle_3dview_rs`.
#[derive(Debug, Clone, Default)]
pub struct Mol3D {
    pub elements: Vec<String>,
    pub positions: Vec<[f64; 3]>,
    /// `(atom_a, atom_b, order)` with 0-based, contiguous indices.
    pub bonds: Vec<(usize, usize, u8)>,
}

impl Mol3D {
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
}

fn order_num(o: &BondOrder) -> u8 {
    match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
    }
}

/// Default valence per element (number of bonds to fill with hydrogens).
/// Returns 0 for elements we don't auto-hydrogenate.
fn base_valence(element: &str) -> i32 {
    match element {
        "H" => 1,
        "B" => 3,
        "C" => 4,
        "N" => 3,
        "O" => 2,
        "F" => 1,
        "Si" => 4,
        "P" => 3,
        "S" => 2,
        "Cl" => 1,
        "Br" => 1,
        "I" => 1,
        _ => 0,
    }
}

/// Tiny deterministic PRNG (splitmix64-style) → value in `[-1.0, 1.0]`.
/// Deterministic so the generated geometry is reproducible across runs.
fn jitter(seed: u64) -> f64 {
    let mut x = seed
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(0x1234_5678_9ABC_DEF1);
    x ^= x >> 30;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    ((x as f64) / (u64::MAX as f64)) * 2.0 - 1.0
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Convert a 2D structure from the editor into an initial 3D structure with
/// explicit hydrogens added.
pub fn to_mol_3d(mol: &chembuider_rs::Molecule) -> Mol3D {
    let mut elements: Vec<String> = Vec::new();
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut bonds: Vec<(usize, usize, u8)> = Vec::new();
    let mut id_to_idx: HashMap<u32, usize> = HashMap::new();

    // Heavy atoms (everything the user drew), lifted into 3D with a z-jitter.
    for atom in &mol.atoms {
        let idx = elements.len();
        id_to_idx.insert(atom.id, idx);
        elements.push(atom.element.clone());
        let z = Z_JITTER * jitter(atom.id as u64 ^ 0xABCD);
        positions.push([atom.pos[0] as f64, atom.pos[1] as f64, z]);
    }

    let heavy_count = elements.len();
    let mut bond_sum = vec![0i32; heavy_count];
    for b in &mol.bonds {
        if let (Some(&a), Some(&c)) = (id_to_idx.get(&b.begin), id_to_idx.get(&b.end)) {
            let o = order_num(&b.order);
            bonds.push((a, c, o));
            bond_sum[a] += o as i32;
            bond_sum[c] += o as i32;
        }
    }

    // Add implicit hydrogens to each heavy atom.
    for hi in 0..heavy_count {
        let h_needed = (base_valence(&elements[hi]) - bond_sum[hi]).max(0);
        if h_needed == 0 {
            continue;
        }

        let p = positions[hi];
        // Average direction toward existing (heavy) neighbors.
        let mut toward = [0.0f64; 3];
        for &(a, c, _) in &bonds {
            let other = if a == hi {
                Some(c)
            } else if c == hi {
                Some(a)
            } else {
                None
            };
            if let Some(o) = other {
                let q = positions[o];
                let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
                let n = norm3(d).max(1e-6);
                for k in 0..3 {
                    toward[k] += d[k] / n;
                }
            }
        }

        for j in 0..h_needed {
            let seed = ((hi as u64) << 20) ^ (j as u64).wrapping_mul(0x0100_0193);
            // Point away from the existing neighbors, then jitter so multiple
            // hydrogens on the same atom don't overlap. MM fixes the geometry.
            let mut dir = [-toward[0], -toward[1], -toward[2]];
            for k in 0..3 {
                dir[k] += 0.8 * jitter(seed.wrapping_add(k as u64));
            }
            let n = norm3(dir);
            let dir = if n < 1e-6 {
                [jitter(seed), jitter(seed + 7), jitter(seed + 13)]
            } else {
                dir
            };
            let n = norm3(dir).max(1e-6);

            let h_idx = elements.len();
            elements.push("H".to_string());
            positions.push([
                p[0] + dir[0] / n * H_BOND_LEN,
                p[1] + dir[1] / n * H_BOND_LEN,
                p[2] + dir[2] / n * H_BOND_LEN,
            ]);
            bonds.push((hi, h_idx, 1));
        }
    }

    Mol3D {
        elements,
        positions,
        bonds,
    }
}

/// Build a `zunda_rs::MmMolecule` from a `Mol3D` at the given positions.
///
/// `atom_class` is a placeholder: every force field re-assigns atom types from
/// topology inside `setup()`, so the value here is ignored.
pub fn build_mm_molecule(mol: &Mol3D, positions: &[[f64; 3]]) -> zunda_rs::MmMolecule {
    use zunda_rs::molecule::Vec3;
    use zunda_rs::{MmAtom, MmBond, MmMolecule};

    let atoms = mol
        .elements
        .iter()
        .enumerate()
        .map(|(i, el)| MmAtom {
            index: i,
            symbol: el.clone(),
            atom_class: el.clone(),
            atomic_number: zunda_rs::atomic_number(el),
            charge: 0.0,
            partial_charge: 0.0,
            position: Vec3::new(positions[i][0], positions[i][1], positions[i][2]),
        })
        .collect();

    let bonds = mol
        .bonds
        .iter()
        .map(|&(a, b, o)| MmBond::new(a, b, o))
        .collect();

    MmMolecule {
        atoms,
        bonds,
        residue: None,
    }
}

/// Build the viewer molecule from a `Mol3D` and the current MM positions (Å).
/// Positions are scaled Å → nm (× 0.1) to match the viewer's loaders.
pub fn to_viewer_molecule(
    mol: &Mol3D,
    positions: &[zunda_rs::molecule::Vec3],
) -> moleucle_3dview_rs::Molecule {
    use lin_alg::f32::Vec3 as FVec3;
    use moleucle_3dview_rs::molecule::{Atom, Bond};

    const ANGSTROM_TO_NM: f32 = 0.1;

    let atoms = mol
        .elements
        .iter()
        .enumerate()
        .map(|(i, el)| Atom {
            position: FVec3::new(
                positions[i].x as f32 * ANGSTROM_TO_NM,
                positions[i].y as f32 * ANGSTROM_TO_NM,
                positions[i].z as f32 * ANGSTROM_TO_NM,
            ),
            element: el.to_uppercase(),
            id: i,
            name: None,
            res_name: None,
            chain_id: None,
            res_seq: None,
            occupancy: None,
            temp_factor: None,
            charge: None,
        })
        .collect();

    let bonds = mol
        .bonds
        .iter()
        .map(|&(a, b, o)| Bond {
            atom_a: a,
            atom_b: b,
            order: o,
        })
        .collect();

    moleucle_3dview_rs::Molecule { atoms, bonds }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbon(mol: &mut chembuider_rs::Molecule, x: f32, y: f32) -> u32 {
        mol.add_atom("C".to_string(), [x, y], 0)
    }

    fn count_h(m: &Mol3D) -> usize {
        m.elements.iter().filter(|e| e.as_str() == "H").count()
    }

    #[test]
    fn methane_gets_four_hydrogens() {
        let mut mol = chembuider_rs::Molecule::default();
        carbon(&mut mol, 0.0, 0.0);
        let m3 = to_mol_3d(&mol);
        assert_eq!(m3.elements.len(), 5, "C + 4H");
        assert_eq!(count_h(&m3), 4);
        assert_eq!(m3.bonds.len(), 4);
    }

    #[test]
    fn ethene_double_bond_leaves_two_h_per_carbon() {
        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, chembuider_rs::BondOrder::Double);
        let m3 = to_mol_3d(&mol);
        assert_eq!(count_h(&m3), 4, "sp2 C contributes 2 H each");
        assert_eq!(m3.elements.len(), 6);
    }

    #[test]
    fn z_jitter_breaks_planarity() {
        let mut mol = chembuider_rs::Molecule::default();
        carbon(&mut mol, 0.0, 0.0);
        carbon(&mut mol, 1.4, 0.0);
        let m3 = to_mol_3d(&mol);
        assert!(
            m3.positions.iter().any(|p| p[2].abs() > 1e-9),
            "expected non-zero z somewhere"
        );
    }

    #[test]
    fn uff_minimization_reduces_energy() {
        use crate::forcefield_kind::FfKind;
        use crate::mm_session::MmSession;

        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, chembuider_rs::BondOrder::Single);

        let m3 = to_mol_3d(&mol);
        let mut session = MmSession::new(m3, FfKind::Uff).expect("session setup");
        let e0 = session.energy;
        session.minimizing = true;
        session.step(300);
        assert!(
            session.energy <= e0 + 1e-6,
            "energy should not increase: {e0} -> {}",
            session.energy
        );
        assert!(session.total_steps > 0, "expected some minimization steps");
    }
}
