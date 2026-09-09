//! Interactive dihedral rotation: pick an anchor atom, pick a bond, turn the
//! wheel, and the side of the bond the anchor is *not* on rotates.
//!
//! The anchor is what makes the gesture unambiguous. A bond splits the molecule
//! in two and either half could be the one that moves; naming an atom that must
//! stay put picks the half.

use crate::geom3d::{self, Dof};
use openbabel::Molecule;

/// What the user has picked in the 3D view for a rotation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Selection {
    /// The atom that must not move.
    pub anchor: Option<u32>,
    /// The bond to turn about, as a 0-based index into `Molecule::bonds()`.
    pub bond: Option<u32>,
}

impl Selection {
    pub fn is_armed(&self) -> bool {
        self.anchor.is_some() && self.bond.is_some()
    }

    pub fn clear(&mut self) {
        *self = Self::default();
    }
}

/// A bond resolved into "which end holds still, which end swings", plus the
/// atoms that swing with it.
#[derive(Debug)]
pub struct Rotation {
    dof: Dof,
    /// The two atoms flanking the bond that define the reported torsion angle,
    /// on the fixed and moving sides respectively.
    reference: Option<(u32, u32)>,
}

/// Why a bond cannot be rotated. Rendered to the user, so the variants are the
/// distinctions worth telling apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reject {
    NoSuchBond,
    /// Ring, multiple or terminal bond — `Bond::is_rotor()` said no.
    NotRotatable,
    /// The anchor is in a different fragment, so the bond does not separate it
    /// from anything.
    AnchorNotAttached,
}

impl Reject {
    pub fn message(self) -> &'static str {
        match self {
            Reject::NoSuchBond => "その結合は存在しません。",
            Reject::NotRotatable => "この結合は回転できません（環内・多重結合・末端のいずれか）。",
            Reject::AnchorNotAttached => "アンカー原子がこの結合と繋がっていません。",
        }
    }
}

impl Rotation {
    /// Resolve `bond` into a rotation that holds `anchor` still.
    ///
    /// `Molecule::find_children(b, c)` returns the atoms on `c`'s side. If the
    /// anchor is over there, the two ends are swapped so the anchor's side is
    /// the one that stays.
    pub fn resolve(mol: &Molecule, bond: u32, anchor: u32) -> Result<Self, Reject> {
        let Some(b) = mol.bond(bond) else {
            return Err(Reject::NoSuchBond);
        };
        if !b.is_rotor() {
            return Err(Reject::NotRotatable);
        }
        let (mut begin, mut end) = (b.begin_atom_index(), b.end_atom_index());

        // The anchor has to end up on the fixed side. Being the far endpoint
        // itself counts: `end` sits on the axis and would not move, but the
        // fragment hanging off it would, which is not what "anchor" means.
        if anchor == end || mol.find_children(begin, end).contains(&anchor) {
            std::mem::swap(&mut begin, &mut end);
        } else if anchor != begin && !mol.find_children(end, begin).contains(&anchor) {
            return Err(Reject::AnchorNotAttached);
        }

        let dof = Dof::new(mol, begin, end).ok_or(Reject::NotRotatable)?;
        Ok(Self {
            reference: reference_atoms(mol, begin, end),
            dof,
        })
    }

    /// The `a`–`b`–`c`–`d` quadruple whose angle this rotation changes, for
    /// reading the torsion back out of OpenBabel.
    pub fn torsion_atoms(&self) -> Option<[u32; 4]> {
        let (a, d) = self.reference?;
        Some([a, self.dof.b, self.dof.c, d])
    }

    /// The current torsion in degrees, or `None` for a bond with no substituent
    /// on one side to measure against.
    pub fn angle_degrees(&self, mol: &Molecule) -> Option<f64> {
        let [a, b, c, d] = self.torsion_atoms()?;
        Some(mol.torsion(a, b, c, d))
    }

    pub fn moving_atoms(&self) -> &[u32] {
        &self.dof.moving
    }

    /// Turn by `degrees` and hand back the new coordinates. The molecule is not
    /// touched — the caller decides whether to commit.
    pub fn rotated(&self, coords: &[[f64; 3]], degrees: f64) -> Vec<[f64; 3]> {
        let mut out = coords.to_vec();
        self.dof.apply(&mut out, degrees.to_radians());
        out
    }
}

/// Pick a neighbour of each endpoint, away from the bond itself, to define the
/// torsion angle. Heavy atoms win over hydrogens so the reported angle matches
/// the one a chemist would quote.
fn reference_atoms(mol: &Molecule, b: u32, c: u32) -> Option<(u32, u32)> {
    let pick = |centre: u32, exclude: u32| -> Option<u32> {
        let atom = mol.atom(centre)?;
        let mut best: Option<(u32, u32)> = None; // (atomic number, index)
        for n in atom.neighbors() {
            let idx = n.index();
            if idx == exclude {
                continue;
            }
            let z = n.atomic_number();
            if best.is_none_or(|(bz, _)| z > bz) {
                best = Some((z, idx));
            }
        }
        best.map(|(_, idx)| idx)
    };
    Some((pick(b, c)?, pick(c, b)?))
}

/// Commit `coords` to the molecule. Separate from `Rotation` because `ik.rs`
/// writes back the same way.
pub fn commit(mol: &mut Molecule, coords: &[[f64; 3]]) -> Result<(), String> {
    if mol.set_coordinates(&geom3d::flatten(coords)) {
        Ok(())
    } else {
        Err("座標の書き戻しに失敗しました。".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::coordinates;

    fn butane() -> Molecule {
        let mut mol = Molecule::parse("CCCC", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        mol
    }

    /// The central C–C bond of the butane above, found by scanning rather than
    /// hard-coded: SMILES parsing order is not a promise, and `add_hydrogens`
    /// inside `generate_3d` appends atoms whose indices are not either.
    fn central_bond(mol: &Molecule) -> u32 {
        mol.bonds()
            .find(|b| b.is_rotor())
            .map(|b| b.index())
            .expect("butane has a rotatable bond")
    }

    /// The whole point of the anchor: whichever side the user named must come
    /// out of the rotation bit-for-bit unmoved. A rotation applied about the
    /// wrong end still changes the torsion by the right amount, so a test that
    /// only checked the angle would pass on a molecule that had swung the wrong
    /// half around the screen.
    #[test]
    fn anchor_side_does_not_move() {
        let _ob = crate::test_support::ob_guard();
        let mol = butane();
        let bond = central_bond(&mol);
        let anchor = mol.bond(bond).expect("bond").begin_atom_index();

        let rot = Rotation::resolve(&mol, bond, anchor).expect("resolve");
        let before = coordinates(&mol);
        let after = rot.rotated(&before, 60.0);

        let moving = rot.moving_atoms();
        assert!(
            !moving.contains(&anchor),
            "the anchor itself must never be in the moving set"
        );
        for i in 0..before.len() as u32 {
            if moving.contains(&i) {
                continue;
            }
            assert_eq!(
                before[i as usize], after[i as usize],
                "atom {i} is on the anchor's side but moved"
            );
        }
    }

    /// A rotation of n degrees has to land n degrees away, measured by
    /// OpenBabel rather than by our own arithmetic — that is what pins our
    /// rotation to the same sign convention as `Molecule::torsion`.
    #[test]
    fn rotating_changes_the_torsion_by_that_much() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();
        let bond = central_bond(&mol);
        let anchor = mol.bond(bond).expect("bond").begin_atom_index();
        let rot = Rotation::resolve(&mol, bond, anchor).expect("resolve");

        let before = rot.angle_degrees(&mol).expect("torsion");
        let coords = rot.rotated(&coordinates(&mol), 60.0);
        commit(&mut mol, &coords).expect("commit");
        let after = rot.angle_degrees(&mol).expect("torsion");

        // Signed, not just in magnitude: the panel's angle box works out how far
        // to turn as `target - current`, so a rotation that moved the torsion
        // the other way would drive it away from whatever the user typed.
        let delta = (after - before + 540.0).rem_euclid(360.0) - 180.0;
        assert!(
            (delta - 60.0).abs() < 1e-6,
            "expected +60 degrees, got {delta}"
        );
    }

    /// Rotating a ring bond would tear the ring open. `is_rotor()` is the guard,
    /// and this pins that we actually consult it.
    #[test]
    fn ring_bonds_are_rejected() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = Molecule::parse("c1ccccc1", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");

        let ring_bond = mol
            .bonds()
            .find(|b| b.is_in_ring())
            .map(|b| b.index())
            .expect("benzene has ring bonds");
        assert_eq!(
            Rotation::resolve(&mol, ring_bond, 0).unwrap_err(),
            Reject::NotRotatable
        );
    }

    /// Bond lengths belong to the force field, not to this gesture. A rotation
    /// that quietly stretched one would show up as an energy jump the user
    /// could not explain.
    #[test]
    fn rotation_preserves_bond_lengths() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();
        let bond = central_bond(&mol);
        let anchor = mol.bond(bond).expect("bond").begin_atom_index();
        let rot = Rotation::resolve(&mol, bond, anchor).expect("resolve");

        let before: Vec<f64> = mol.bonds().map(|b| b.length()).collect();
        let coords = rot.rotated(&coordinates(&mol), 137.0);
        commit(&mut mol, &coords).expect("commit");
        let after: Vec<f64> = mol.bonds().map(|b| b.length()).collect();

        for (i, (x, y)) in before.iter().zip(&after).enumerate() {
            assert!((x - y).abs() < 1e-9, "bond {i}: {x} -> {y}");
        }
    }
}
