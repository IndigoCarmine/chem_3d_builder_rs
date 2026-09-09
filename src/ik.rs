//! Inverse kinematics for a pair of atoms: name two atoms and a target
//! relationship, and the dihedrals of the rotatable bonds between them are
//! solved so the pair ends up that way.
//!
//! This is the inverse of what `dihedral.rs` does. There the user turns a bond
//! and watches where the atoms land; here the user states where the atoms
//! should land and the bond angles are worked out.
//!
//! ## Method
//!
//! Cyclic coordinate descent, with each degree of freedom solved in closed
//! form. Rotating a bond moves the target atom `b` on a circle about that
//! bond's axis, so with the axis unit vector `û`, a point `p` on the axis, and
//! the fixed atom `a`:
//!
//! ```text
//! v  = b − p            v∥ = (v·û)û         v⊥ = v − v∥
//! w  = a − p − v∥
//! d²(φ) = C − 2R·cos(φ − α)
//!   C = |w|² + |v⊥|²    P = w·v⊥            Q = w·(û × v⊥)
//!   R = √(P² + Q²)      α = atan2(Q, P)
//! ```
//!
//! So one rotation sweeps the distance over exactly `[√(C−2R), √(C+2R)]`, the
//! extremes sit at `φ = α` and `φ = α + π`, and a requested distance inside the
//! range is one `acos` away. No line search, no numerical gradient.

use crate::geom3d::{self, Dof};
use openbabel::Molecule;

/// Stop once the distance is this close to the target, in Å. Well below what
/// any of this is accurate to chemically, and reached in a handful of sweeps.
const TOLERANCE: f64 = 1e-4;
/// Give up on a sweep that improves the objective by less than this.
const MIN_IMPROVEMENT: f64 = 1e-9;
const MAX_SWEEPS: usize = 200;

/// What the two atoms should end up as.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Goal {
    /// A specific separation in Å.
    Distance(f64),
    AsFarAsPossible,
    AsCloseAsPossible,
}

/// Why a request could not be set up at all. Distinct from "solved, but could
/// not reach the number you asked for", which is a successful solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    SameAtom,
    /// The two atoms are in different fragments.
    NotConnected,
    /// Every bond between them is a ring bond, a multiple bond or terminal.
    NoRotatableBond,
}

impl Error {
    pub fn message(self) -> &'static str {
        match self {
            Error::SameAtom => "同じ原子が2回選ばれています。",
            Error::NotConnected => "2つの原子は繋がっていません。",
            Error::NoRotatableBond => {
                "この2原子の間に回転できる結合がありません（環内・多重結合のみ）。"
            }
        }
    }
}

/// The outcome of a solve.
#[derive(Debug)]
pub struct Solution {
    /// The coordinates that achieve `achieved`. Not yet written to the molecule.
    pub coords: Vec<[f64; 3]>,
    pub start: f64,
    pub achieved: f64,
    /// Set for `Goal::Distance` when the target lies outside what the rotatable
    /// bonds can reach. `achieved` is then the closest the solver could get.
    pub unreachable: bool,
    /// How far apart the pair can be driven with these degrees of freedom, from
    /// the geometry the solve ended on. Reported alongside `unreachable` so the
    /// user is told what to ask for instead.
    pub reach: (f64, f64),
    pub sweeps: usize,
}

/// The rotatable bonds on the shortest path between `a` and `b`.
///
/// Restricting the degrees of freedom to the path is what keeps the result
/// legible: every bond in this set demonstrably lies between the two atoms, so
/// nothing off to the side of the molecule swings around for reasons the user
/// cannot see. Ring and multiple bonds drop out on their own — `is_rotor()` is
/// false for both.
fn degrees_of_freedom(mol: &Molecule, a: u32, b: u32) -> Result<Vec<Dof>, Error> {
    if a == b {
        return Err(Error::SameAtom);
    }
    let adj = geom3d::adjacency(mol);
    let path = geom3d::shortest_path(&adj, a, b).ok_or(Error::NotConnected)?;

    let mut dofs = Vec::new();
    for pair in path.windows(2) {
        let (u, v) = (pair[0], pair[1]);
        let Some(bond) = mol.bond_between(u, v) else {
            continue;
        };
        if !bond.is_rotor() {
            continue;
        }
        // Oriented so the moving side is the one holding `b`; the path runs
        // from `a` to `b`, so that is always the `v` side.
        if let Some(dof) = Dof::new(mol, u, v)
            && dof.moving.contains(&b)
        {
            dofs.push(dof);
        }
    }
    if dofs.is_empty() {
        return Err(Error::NoRotatableBond);
    }
    Ok(dofs)
}

/// How rotating `dof` changes the `a`–`b` distance, as the `C`, `R`, `α` of the
/// module docs. `None` when `b` sits on the axis, where the rotation cannot
/// change the distance at all.
fn response(coords: &[[f64; 3]], dof: &Dof, a: u32, b: u32) -> Option<(f64, f64, f64)> {
    let p = coords[dof.b as usize];
    let axis = geom3d::unit(geom3d::sub(coords[dof.c as usize], p))?;

    let v = geom3d::sub(coords[b as usize], p);
    let par = geom3d::dot(v, axis);
    let perp = [
        v[0] - par * axis[0],
        v[1] - par * axis[1],
        v[2] - par * axis[2],
    ];
    let w = {
        let t = geom3d::sub(coords[a as usize], p);
        [
            t[0] - par * axis[0],
            t[1] - par * axis[1],
            t[2] - par * axis[2],
        ]
    };

    let c = geom3d::dot(w, w) + geom3d::dot(perp, perp);
    let p_term = geom3d::dot(w, perp);
    let q_term = geom3d::dot(w, geom3d::cross(axis, perp));
    let r = (p_term * p_term + q_term * q_term).sqrt();
    if r < 1e-12 {
        return None;
    }
    Some((c, r, q_term.atan2(p_term)))
}

/// Fold an angle into `(-π, π]`, so "rotate to α" and "rotate to α − 2π" pick
/// the smaller of the two turns.
fn wrap(radians: f64) -> f64 {
    let x = (radians + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU);
    x - std::f64::consts::PI
}

/// The rotation this degree of freedom should take.
///
/// When the target is outside what this one bond can reach, it goes to whichever
/// extreme is on the right side and the next degree of freedom takes up the
/// rest — that is the point of chaining them, so this is not a failure and does
/// not get reported as one.
fn best_angle(goal: Goal, c: f64, r: f64, alpha: f64) -> f64 {
    match goal {
        Goal::AsCloseAsPossible => wrap(alpha),
        Goal::AsFarAsPossible => wrap(alpha + std::f64::consts::PI),
        Goal::Distance(target) => {
            let cos_theta = (c - target * target) / (2.0 * r);
            if cos_theta.abs() <= 1.0 {
                let theta = cos_theta.acos();
                // Two solutions, α ± θ. Take the smaller turn: the other one
                // reaches the same distance by swinging the fragment much
                // further round, which needlessly rearranges everything else.
                let (x, y) = (wrap(alpha + theta), wrap(alpha - theta));
                if x.abs() <= y.abs() { x } else { y }
            } else if cos_theta > 1.0 {
                // Target shorter than this bond alone can make it.
                wrap(alpha)
            } else {
                wrap(alpha + std::f64::consts::PI)
            }
        }
    }
}

fn objective(goal: Goal, d: f64) -> f64 {
    match goal {
        Goal::Distance(target) => (d - target).abs(),
        Goal::AsFarAsPossible => -d,
        Goal::AsCloseAsPossible => d,
    }
}

/// Solve for `goal` between atoms `a` and `b`, starting from `coords`.
///
/// The molecule is read for connectivity only; the coordinates handed back are
/// the caller's to commit or discard.
pub fn solve(
    mol: &Molecule,
    coords: &[[f64; 3]],
    a: u32,
    b: u32,
    goal: Goal,
) -> Result<Solution, Error> {
    let dofs = degrees_of_freedom(mol, a, b)?;
    let mut coords = coords.to_vec();
    let start = geom3d::distance(coords[a as usize], coords[b as usize]);

    let mut sweeps = 0;
    let mut best = objective(goal, start);

    while sweeps < MAX_SWEEPS {
        sweeps += 1;

        for dof in &dofs {
            let Some((c, r, alpha)) = response(&coords, dof, a, b) else {
                continue;
            };
            let angle = best_angle(goal, c, r, alpha);
            if angle.abs() > 1e-12 {
                dof.apply(&mut coords, angle);
            }
        }

        let d = geom3d::distance(coords[a as usize], coords[b as usize]);
        if let Goal::Distance(target) = goal
            && (d - target).abs() < TOLERANCE
        {
            break;
        }
        // A sweep that no longer moves the objective has converged; without
        // this the loop would spend all 200 sweeps re-solving the same angles.
        let score = objective(goal, d);
        let improved = best - score;
        best = score;
        if improved < MIN_IMPROVEMENT {
            break;
        }
    }

    // Reaching the number that was asked for is the only thing that decides
    // this. The two open-ended goals always succeed by construction.
    let achieved = geom3d::distance(coords[a as usize], coords[b as usize]);
    let unreachable = matches!(goal, Goal::Distance(t) if (achieved - t).abs() >= TOLERANCE);

    Ok(Solution {
        reach: reach(&coords, &dofs, a, b),
        coords,
        start,
        achieved,
        unreachable,
        sweeps,
    })
}

/// The separations the current degrees of freedom can still produce, as a
/// single-sweep estimate from the geometry the solve ended on.
///
/// It is an estimate, not a bound: each bond's own `[√(C−2R), √(C+2R)]` is
/// exact, but the bonds interact, so the combined reach is wider than any one
/// of them. Chaining the per-bond extremes in sequence gets close enough to
/// tell the user roughly what they could have asked for.
fn reach(coords: &[[f64; 3]], dofs: &[Dof], a: u32, b: u32) -> (f64, f64) {
    let sweep = |goal: Goal| {
        let mut work = coords.to_vec();
        for dof in dofs {
            if let Some((c, r, alpha)) = response(&work, dof, a, b) {
                dof.apply(&mut work, best_angle(goal, c, r, alpha));
            }
        }
        geom3d::distance(work[a as usize], work[b as usize])
    };
    let lo = sweep(Goal::AsCloseAsPossible);
    let hi = sweep(Goal::AsFarAsPossible);
    (lo.min(hi), lo.max(hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::coordinates;
    use crate::dihedral::commit;

    /// n-pentane: four heavy atoms of chain between the two ends, so there are
    /// real degrees of freedom and the extremes are far apart.
    fn pentane() -> Molecule {
        let mut mol = Molecule::parse("CCCCC", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        mol
    }

    /// The two terminal carbons of a straight chain, as 0-based indices.
    fn chain_ends(mol: &Molecule) -> (u32, u32) {
        let carbons: Vec<u32> = mol
            .atoms()
            .filter(|a| a.atomic_number() == 6)
            .map(|a| a.index())
            .collect();
        (carbons[0], *carbons.last().expect("carbons"))
    }

    /// The two directions have to actually differ, and in the right order.
    /// A solver that quietly did nothing would still report a "solution" whose
    /// distance matched the starting geometry, so both ends are checked against
    /// each other rather than against a hard-coded number.
    #[test]
    fn far_is_farther_than_close() {
        let _ob = crate::test_support::ob_guard();
        let mol = pentane();
        let (a, b) = chain_ends(&mol);
        let coords = coordinates(&mol);

        let far = solve(&mol, &coords, a, b, Goal::AsFarAsPossible).expect("far");
        let close = solve(&mol, &coords, a, b, Goal::AsCloseAsPossible).expect("close");

        assert!(
            far.achieved > close.achieved + 1.0,
            "far {} vs close {}",
            far.achieved,
            close.achieved
        );
        assert!(far.achieved >= far.start, "far should not shorten the pair");
        assert!(
            close.achieved <= close.start,
            "close should not lengthen the pair"
        );
    }

    /// The point of the closed form: a distance inside the reachable range is
    /// hit exactly, not approached. If this ever degrades into "close enough",
    /// the per-bond `acos` branch has stopped firing.
    #[test]
    fn a_reachable_distance_is_hit_exactly() {
        let _ob = crate::test_support::ob_guard();
        let mol = pentane();
        let (a, b) = chain_ends(&mol);
        let coords = coordinates(&mol);

        let lo = solve(&mol, &coords, a, b, Goal::AsCloseAsPossible)
            .expect("close")
            .achieved;
        let hi = solve(&mol, &coords, a, b, Goal::AsFarAsPossible)
            .expect("far")
            .achieved;
        let target = 0.5 * (lo + hi);

        let got = solve(&mol, &coords, a, b, Goal::Distance(target)).expect("target");
        assert!(
            !got.unreachable,
            "midpoint of the range should be reachable"
        );
        assert!(
            (got.achieved - target).abs() < 1e-3,
            "asked {target}, got {}",
            got.achieved
        );
    }

    /// An impossible request must say so rather than silently returning the
    /// nearest geometry as if it had succeeded — that is the difference between
    /// "your molecule is now 8 Å" and "8 Å is what you asked for".
    #[test]
    fn an_impossible_distance_is_reported() {
        let _ob = crate::test_support::ob_guard();
        let mol = pentane();
        let (a, b) = chain_ends(&mol);
        let coords = coordinates(&mol);

        let got = solve(&mol, &coords, a, b, Goal::Distance(50.0)).expect("solve");
        assert!(got.unreachable, "50 Å is out of reach for pentane");
        assert!(got.achieved < 10.0, "achieved {}", got.achieved);
        assert!(got.reach.0 < got.reach.1);
    }

    /// The solver is only allowed to turn bonds. If it ever moved an atom in a
    /// way that changed a bond length, it would be inventing chemistry the
    /// force field never agreed to.
    #[test]
    fn solving_only_rotates() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = pentane();
        let (a, b) = chain_ends(&mol);
        let before: Vec<f64> = mol.bonds().map(|x| x.length()).collect();

        let solution =
            solve(&mol, &coordinates(&mol), a, b, Goal::AsCloseAsPossible).expect("solve");
        commit(&mut mol, &solution.coords).expect("commit");

        let after: Vec<f64> = mol.bonds().map(|x| x.length()).collect();
        for (i, (x, y)) in before.iter().zip(&after).enumerate() {
            assert!((x - y).abs() < 1e-9, "bond {i}: {x} -> {y}");
        }
    }

    /// Benzene's ring bonds are all non-rotors, so there is nothing to solve.
    /// Saying that plainly beats returning an unchanged structure and letting
    /// the user wonder whether the button worked.
    #[test]
    fn a_rigid_path_is_rejected() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = Molecule::parse("c1ccccc1", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        let carbons: Vec<u32> = mol
            .atoms()
            .filter(|x| x.atomic_number() == 6)
            .map(|x| x.index())
            .collect();

        let err = solve(
            &mol,
            &coordinates(&mol),
            carbons[0],
            carbons[3],
            Goal::AsFarAsPossible,
        )
        .unwrap_err();
        assert_eq!(err, Error::NoRotatableBond);
    }
}
