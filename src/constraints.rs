//! Geometry restraints handed to the MM minimizer.
//!
//! OpenBabel's `OBFFConstraints` is a cxx opaque type and therefore `!Send`, so
//! it cannot be built here and moved to the worker thread that runs the
//! minimization (see `mm_session::run_minimize`). What crosses the thread
//! boundary is this module's plain [`Constraint`] — indices and a target value,
//! nothing borrowed and nothing from the FFI — which the worker turns into an
//! `openbabel::Constraints` on its own side, once per segment.

/// One restraint on the geometry.
///
/// Atom indices are 0-based, matching `openbabel::Atom::index` and the indices
/// the 3D viewport reports from a pick — `constraint_indices_are_zero_based`
/// pins that the binding agrees, because an off-by-one here would silently
/// restrain a neighbouring atom instead of failing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    /// Hold the `a`–`b` separation at `length` Å.
    Distance { a: u32, b: u32, length: f64 },
    /// Hold the `a`–`b`–`c` valence angle at `degrees`.
    Angle {
        a: u32,
        b: u32,
        c: u32,
        degrees: f64,
    },
    /// Hold the `a`–`b`–`c`–`d` torsion at `degrees`.
    Torsion {
        a: u32,
        b: u32,
        c: u32,
        d: u32,
        degrees: f64,
    },
}

/// How stiff the restraints are.
///
/// OpenBabel applies one factor to the whole constraint set rather than per
/// restraint, so this has to satisfy the most demanding kind -- the torsion.
/// Distance and angle restraints are harmonic and unbounded, and hold at
/// OpenBabel's own default of 50 000; a torsion restraint is a cosine well only
/// `0.002 * factor` deep, which the force field's own torsional term can simply
/// outweigh. Measured on butane's C-C-C-C dihedral pinned to 90 deg: it does
/// not move at 50 000, lands within 5 deg at 1e6, and on the target at 1e8.
/// Distance and angle restraints still hold at 1e8, so one factor serves all
/// three.
pub const FORCE_FACTOR: f64 = 1.0e8;

impl Constraint {
    /// How many atoms this kind of restraint needs.
    pub fn arity(kind: Kind) -> usize {
        match kind {
            Kind::Distance => 2,
            Kind::Angle => 3,
            Kind::Torsion => 4,
        }
    }

    pub fn kind(&self) -> Kind {
        match self {
            Constraint::Distance { .. } => Kind::Distance,
            Constraint::Angle { .. } => Kind::Angle,
            Constraint::Torsion { .. } => Kind::Torsion,
        }
    }

    /// The atoms this restraint acts on, in order.
    pub fn atoms(&self) -> Vec<u32> {
        match *self {
            Constraint::Distance { a, b, .. } => vec![a, b],
            Constraint::Angle { a, b, c, .. } => vec![a, b, c],
            Constraint::Torsion { a, b, c, d, .. } => vec![a, b, c, d],
        }
    }

    /// The target, in Å for a distance and degrees for the two angles.
    pub fn target(&self) -> f64 {
        match *self {
            Constraint::Distance { length, .. } => length,
            Constraint::Angle { degrees, .. } => degrees,
            Constraint::Torsion { degrees, .. } => degrees,
        }
    }

    pub fn set_target(&mut self, value: f64) {
        match self {
            Constraint::Distance { length, .. } => *length = value,
            Constraint::Angle { degrees, .. } => *degrees = value,
            Constraint::Torsion { degrees, .. } => *degrees = value,
        }
    }

    /// Build a restraint of `kind` over `atoms` aimed at `target`. `None` if the
    /// wrong number of atoms is given.
    pub fn new(kind: Kind, atoms: &[u32], target: f64) -> Option<Self> {
        if atoms.len() != Self::arity(kind) {
            return None;
        }
        Some(match kind {
            Kind::Distance => Constraint::Distance {
                a: atoms[0],
                b: atoms[1],
                length: target,
            },
            Kind::Angle => Constraint::Angle {
                a: atoms[0],
                b: atoms[1],
                c: atoms[2],
                degrees: target,
            },
            Kind::Torsion => Constraint::Torsion {
                a: atoms[0],
                b: atoms[1],
                c: atoms[2],
                d: atoms[3],
                degrees: target,
            },
        })
    }

    /// Build a restraint of `kind` over `atoms`, targeting whatever the
    /// geometry currently reads. `None` if the wrong number of atoms is given.
    pub fn from_current(kind: Kind, atoms: &[u32], mol: &openbabel::Molecule) -> Option<Self> {
        if atoms.len() != Self::arity(kind) {
            return None;
        }
        Some(match kind {
            Kind::Distance => Constraint::Distance {
                a: atoms[0],
                b: atoms[1],
                length: mol.distance(atoms[0], atoms[1]),
            },
            Kind::Angle => Constraint::Angle {
                a: atoms[0],
                b: atoms[1],
                c: atoms[2],
                degrees: mol.angle(atoms[0], atoms[1], atoms[2]),
            },
            Kind::Torsion => Constraint::Torsion {
                a: atoms[0],
                b: atoms[1],
                c: atoms[2],
                d: atoms[3],
                degrees: mol.torsion(atoms[0], atoms[1], atoms[2], atoms[3]),
            },
        })
    }

    /// What the geometry currently reads for this restraint's atoms, for showing
    /// the user how far the structure is from the target.
    pub fn measure(&self, mol: &openbabel::Molecule) -> f64 {
        match *self {
            Constraint::Distance { a, b, .. } => mol.distance(a, b),
            Constraint::Angle { a, b, c, .. } => mol.angle(a, b, c),
            Constraint::Torsion { a, b, c, d, .. } => mol.torsion(a, b, c, d),
        }
    }

    /// Every index this restraint names must exist in the molecule, and no atom
    /// may appear twice — OpenBabel answers a degenerate restraint with a NaN
    /// energy rather than an error.
    pub fn is_valid_for(&self, num_atoms: u32) -> bool {
        let atoms = self.atoms();
        atoms.iter().all(|&i| i < num_atoms)
            && atoms
                .iter()
                .enumerate()
                .all(|(n, a)| !atoms[..n].contains(a))
    }
}

/// Which restraint the panel is currently collecting atoms for.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Kind {
    #[default]
    Distance,
    Angle,
    Torsion,
}

impl Kind {
    pub const ALL: &'static [Kind] = &[Kind::Distance, Kind::Angle, Kind::Torsion];

    pub fn label(self) -> &'static str {
        match self {
            Kind::Distance => "距離 (2原子)",
            Kind::Angle => "角度 (3原子)",
            Kind::Torsion => "二面角 (4原子)",
        }
    }

    pub fn unit(self) -> &'static str {
        match self {
            Kind::Distance => "Å",
            Kind::Angle | Kind::Torsion => "°",
        }
    }

    /// The range the target box clamps to. A distance shorter than any real
    /// contact, or an angle outside 0–180, is not a geometry to restrain towards.
    pub fn range(self) -> std::ops::RangeInclusive<f64> {
        match self {
            Kind::Distance => 0.5..=20.0,
            Kind::Angle => 0.0..=180.0,
            // Torsions are signed: OpenBabel reports them in -180..180.
            Kind::Torsion => -180.0..=180.0,
        }
    }

    /// Chemically meaningful targets, offered beside the free-text box so the
    /// common intentions — "make this a hydrogen bond", "flatten this" — do not
    /// have to be recalled as numbers.
    pub fn presets(self) -> &'static [(&'static str, f64)] {
        match self {
            // Hydrogen bonds are quoted two ways and the difference is a whole
            // Ångström, so each says which atoms it means rather than leaving
            // the user to guess which convention this box wants.
            Kind::Distance => &[
                ("水素結合 H···受容体", 1.9),
                ("水素結合 供与体···受容体", 2.9),
                ("C–C 単結合", 1.54),
                ("C=C 二重結合", 1.34),
            ],
            Kind::Angle => &[
                ("四面体 sp³", 109.5),
                ("三方平面 sp²", 120.0),
                ("直線 sp", 180.0),
            ],
            Kind::Torsion => &[
                ("平面 シス", 0.0),
                ("平面 トランス", 180.0),
                ("ゴーシュ +", 60.0),
                ("ゴーシュ −", -60.0),
                ("直交", 90.0),
            ],
        }
    }
}

/// The combination OpenBabel cannot be trusted with: a dihedral restraint in a
/// set that also holds an angle restraint.
///
/// A dihedral on its own is honoured exactly, and so is a dihedral beside a
/// distance restraint. Add an angle restraint and the result stops being
/// predictable — measured on butane asking for a 90 deg C-C-C-C dihedral at the
/// shipping force factor:
///
/// ```text
/// dihedral only                  89.9   ok
/// dihedral + distance            90.2   ok
/// dihedral + angle 1-2-3         92.9   ok
/// dihedral + angle 0-1-2        118.9   wrong by 29 deg
/// dihedral + distance + angle   180.0   ignored outright
/// ```
///
/// Which angle spoils it depends on the angle, not on whether it shares atoms
/// with the dihedral, and no amount of stiffness, step budget or convergence
/// tightening changes the failures: from 1e3 to 1e10, to 50 000 steps, to a
/// convergence of 1e-12, the bad cases stay bad. There is no rule here worth
/// relying on, so the whole family is refused rather than shipped as a control
/// that works four times in six.
///
/// What matters is that the app never runs such a set and hands back a
/// structure whose dihedral was quietly ignored.
pub fn conflict(constraints: &[Constraint]) -> Option<&'static str> {
    let has = |k: Kind| constraints.iter().any(|c| c.kind() == k);
    (has(Kind::Torsion) && has(Kind::Angle)).then_some(
        "二面角と角度の制約は同時に使えません（OpenBabel の制限で、二面角が無視されます）。         どちらかを削除してください。",
    )
}

/// Turn the plain description into OpenBabel's own constraint set.
///
/// Called on the worker thread — `openbabel::Constraints` is `!Send`, so this
/// cannot be hoisted to the caller.
pub fn to_ob(constraints: &[Constraint], num_atoms: u32) -> openbabel::Constraints {
    let mut ob = openbabel::Constraints::new();
    for c in constraints.iter().filter(|c| c.is_valid_for(num_atoms)) {
        match *c {
            Constraint::Distance { a, b, length } => {
                ob.distance(a, b, length);
            }
            Constraint::Angle { a, b, c, degrees } => {
                ob.angle(a, b, c, degrees);
            }
            Constraint::Torsion {
                a,
                b,
                c,
                d,
                degrees,
            } => {
                ob.torsion(a, b, c, d, degrees);
            }
        }
    }
    ob.force_factor(FORCE_FACTOR);
    ob
}

#[cfg(test)]
mod tests {
    use super::*;
    use openbabel::{Algorithm, Minimizer, Molecule};

    fn butane() -> Molecule {
        let mut mol = Molecule::parse("CCCC", "smi").expect("parse");
        mol.add_hydrogens();
        assert!(mol.generate_3d(), "gen3d");
        mol
    }

    fn minimize_with(mol: &mut Molecule, constraints: &[Constraint]) {
        let n = mol.num_atoms();
        let mut cfg = Minimizer::new("UFF");
        cfg.algorithm(Algorithm::ConjugateGradients)
            .max_steps(2000)
            .constraints(to_ob(constraints, n));
        let run = mol.minimize(&cfg);
        let _: Vec<_> = run.collect();
    }

    /// The test the rest of this module rests on.
    ///
    /// OpenBabel's own `OBFFConstraints` numbers atoms from 1 (`OBAtom::GetIdx`),
    /// while this crate's `Atom::index` and everything the app passes around
    /// numbers them from 0. If the binding did not absorb that difference, every
    /// restraint would quietly land on the neighbouring atom — a wrong answer
    /// with no error anywhere, which is the worst way for this to be wrong.
    ///
    /// So: stretch the 0–1 bond well past any real C–C length and check that
    /// the pair that moved is 0–1 and not 1–2.
    #[test]
    fn constraint_indices_are_zero_based() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        const TARGET: f64 = 2.50;
        let before_01 = mol.distance(0, 1);
        let before_12 = mol.distance(1, 2);
        assert!(
            (before_01 - TARGET).abs() > 0.5 && (before_12 - TARGET).abs() > 0.5,
            "both bonds must start far from the target for this to discriminate"
        );

        minimize_with(
            &mut mol,
            &[Constraint::Distance {
                a: 0,
                b: 1,
                length: TARGET,
            }],
        );

        let after_01 = mol.distance(0, 1);
        let after_12 = mol.distance(1, 2);
        assert!(
            (after_01 - TARGET).abs() < 0.15,
            "the 0-1 bond should have been pulled to {TARGET} Å, but reads {after_01:.3} \
             (1-2 reads {after_12:.3}) — the binding is not 0-based"
        );
        assert!(
            (after_12 - TARGET).abs() > 0.3,
            "the 1-2 bond moved to {after_12:.3} instead: the indices are off by one"
        );
    }

    /// An angle restraint holds, and holds the angle it names.
    #[test]
    fn an_angle_constraint_holds() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        const TARGET: f64 = 95.0;
        assert!(
            (mol.angle(0, 1, 2) - TARGET).abs() > 5.0,
            "the angle must start away from the target"
        );

        minimize_with(
            &mut mol,
            &[Constraint::Angle {
                a: 0,
                b: 1,
                c: 2,
                degrees: TARGET,
            }],
        );

        let got = mol.angle(0, 1, 2);
        assert!(
            (got - TARGET).abs() < 5.0,
            "the 0-1-2 angle should sit near {TARGET}°, reads {got:.1}°"
        );
    }

    /// A torsion restraint holds, and holds the angle it names.
    ///
    /// Butane's C-C-C-C dihedral is the case worth pinning: it relaxes to anti
    /// (180 deg) on its own, so holding it at 90 deg is a real demand on the
    /// restraint rather than something the force field would have done anyway.
    ///
    /// This is also the regression test for the bug that made dihedral
    /// restraints look impossible. OpenBabel copies its force factor into each
    /// restraint as that restraint is added, so a factor set afterwards reached
    /// nothing already there, and every torsion ran at the default 50 000 --
    /// far too soft for a bounded cosine well. It looked exactly like a broken
    /// torsion term, which is what it was mistaken for. `openbabel_rs` now
    /// rebuilds the set so the factor applies whatever the order.
    #[test]
    fn a_torsion_constraint_holds() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        const TARGET: f64 = 90.0;
        let before = mol.torsion(0, 1, 2, 3).abs();
        assert!(
            (before - TARGET).abs() > 30.0,
            "the dihedral must start far from the target: {before:.1} deg"
        );

        minimize_with(
            &mut mol,
            &[Constraint::Torsion {
                a: 0,
                b: 1,
                c: 2,
                d: 3,
                degrees: TARGET,
            }],
        );

        let got = mol.torsion(0, 1, 2, 3).abs();
        assert!(
            (got - TARGET).abs() < 10.0,
            "the 0-1-2-3 torsion should sit near {TARGET} deg, reads {got:.1}              (an unrestrained butane relaxes to ~180)"
        );
    }

    /// A dihedral holds alongside a distance restraint — the combination that
    /// does work — with each checked separately, because a set that silently
    /// dropped one would still pass a test that only looked at the other.
    #[test]
    fn a_torsion_holds_alongside_a_distance() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        let distance = Constraint::Distance {
            a: 0,
            b: 1,
            length: 1.7,
        };
        let torsion = Constraint::Torsion {
            a: 0,
            b: 1,
            c: 2,
            d: 3,
            degrees: 90.0,
        };
        minimize_with(&mut mol, &[distance, torsion]);

        let (d, t) = (distance.measure(&mol), torsion.measure(&mol).abs());
        assert!(
            (d - 1.7).abs() < 0.2,
            "distance restraint dropped: {d:.2} A"
        );
        assert!(
            (t - 90.0).abs() < 10.0,
            "torsion restraint dropped: {t:.1} deg"
        );
    }

    /// Distance and angle together, the other pair that works.
    #[test]
    fn a_distance_and_an_angle_hold_together() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        let distance = Constraint::Distance {
            a: 0,
            b: 3,
            length: 3.6,
        };
        let angle = Constraint::Angle {
            a: 0,
            b: 1,
            c: 2,
            degrees: 100.0,
        };
        minimize_with(&mut mol, &[distance, angle]);

        let (d, a) = (distance.measure(&mol), angle.measure(&mol));
        assert!(
            (d - 3.6).abs() < 0.3,
            "distance restraint dropped: {d:.2} A"
        );
        assert!(
            (a - 100.0).abs() < 8.0,
            "angle restraint dropped: {a:.1} deg"
        );
    }

    /// The measurement behind [`conflict`]: with an angle restraint in the set
    /// the dihedral is ignored outright. Pinned so that the day OpenBabel
    /// honours both, this test fails and the refusal can go.
    #[test]
    fn an_angle_restraint_defeats_a_dihedral() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = butane();

        let torsion = Constraint::Torsion {
            a: 0,
            b: 1,
            c: 2,
            d: 3,
            degrees: 90.0,
        };
        let distance = Constraint::Distance {
            a: 0,
            b: 1,
            length: 1.7,
        };
        let angle = Constraint::Angle {
            a: 0,
            b: 1,
            c: 2,
            degrees: 100.0,
        };
        minimize_with(&mut mol, &[torsion, distance, angle]);

        let t = torsion.measure(&mol).abs();
        assert!(
            (t - 90.0).abs() > 45.0,
            "OpenBabel honoured the dihedral ({t:.1} deg) alongside an angle restraint —              `conflict()` and the panel's refusal are no longer needed"
        );
    }

    /// The refusal fires on exactly that pair and nothing else.
    #[test]
    fn conflict_names_only_the_dihedral_and_angle_pair() {
        let d = Constraint::Distance {
            a: 0,
            b: 1,
            length: 1.5,
        };
        let a = Constraint::Angle {
            a: 0,
            b: 1,
            c: 2,
            degrees: 109.5,
        };
        let t = Constraint::Torsion {
            a: 0,
            b: 1,
            c: 2,
            d: 3,
            degrees: 90.0,
        };
        assert!(conflict(&[]).is_none());
        assert!(conflict(&[d, a]).is_none());
        assert!(conflict(&[d, t]).is_none());
        assert!(conflict(&[t]).is_none());
        assert!(conflict(&[a, t]).is_some());
        assert!(conflict(&[d, a, t]).is_some());
    }
}
