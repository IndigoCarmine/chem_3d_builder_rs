//! Geometry and graph primitives shared by the two 3D editing features:
//! interactive dihedral rotation (`dihedral.rs`) and the distance solver
//! (`ik.rs`).
//!
//! Coordinates live here as a plain `Vec<[f64; 3]>` rather than inside the
//! OpenBabel molecule. Every entry point of `openbabel_rs` takes a
//! process-global lock, so a solver that rotated through `Molecule::set_torsion`
//! would take that lock once per degree of freedom per sweep. Reading the
//! coordinates once, doing the arithmetic here, and writing back once keeps the
//! lock out of the inner loop.

use openbabel::Molecule;

/// Atom pairs closer than this fraction of the sum of their van der Waals radii
/// count as clashing. 0.75 is the usual "too close to be a real contact"
/// threshold — hydrogen bonds sit around 0.8, so a lower value would flag them.
const CLASH_FACTOR: f64 = 0.75;

/// How far out of plane an atom has to sit, in Å, for a structure to count as
/// three-dimensional. Well under a bond length, and well over the rounding in a
/// file that stores coordinates to three decimals.
const PLANAR_TOLERANCE: f64 = 0.01;

/// Subtract two points.
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

pub fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}

/// Unit vector, or `None` for a vector that has no direction — too short to
/// have one, which is what a "bond" between two atoms at the same position
/// would give, or not a number at all.
///
/// The NaN case is spelled out rather than left to `n < 1e-12`, which is *false*
/// for NaN: that spelling hands back `Some([NaN; 3])`, and the NaN then gets
/// rotated into the molecule and written back by `dihedral::commit`.
pub fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    if !n.is_finite() || n < 1e-12 {
        return None;
    }
    Some([a[0] / n, a[1] / n, a[2] / n])
}

/// Rotate `moving` about the axis through `origin` along `axis_unit` (which must
/// be normalized) by `radians`, in place. Rodrigues' rotation formula.
pub fn rotate_about_axis(
    coords: &mut [[f64; 3]],
    moving: &[u32],
    origin: [f64; 3],
    axis_unit: [f64; 3],
    radians: f64,
) {
    let (s, c) = radians.sin_cos();
    for &i in moving {
        let p = coords[i as usize];
        let v = sub(p, origin);
        let par = dot(v, axis_unit);
        // v = v_parallel + v_perp; only v_perp turns.
        let perp = [
            v[0] - par * axis_unit[0],
            v[1] - par * axis_unit[1],
            v[2] - par * axis_unit[2],
        ];
        let tang = cross(axis_unit, perp);
        coords[i as usize] = [
            origin[0] + par * axis_unit[0] + c * perp[0] + s * tang[0],
            origin[1] + par * axis_unit[1] + c * perp[1] + s * tang[1],
            origin[2] + par * axis_unit[2] + c * perp[2] + s * tang[2],
        ];
    }
}

/// Neighbour lists, built from the bond table in one pass.
///
/// `Atom::neighbors()` would answer the same question, but through the FFI once
/// per atom; the traversals here run over the whole molecule repeatedly.
pub fn adjacency(mol: &Molecule) -> Vec<Vec<u32>> {
    let n = mol.num_atoms() as usize;
    let mut adj = vec![Vec::new(); n];
    for b in mol.bonds() {
        let (i, j) = (b.begin_atom_index() as usize, b.end_atom_index() as usize);
        if i < n && j < n {
            adj[i].push(j as u32);
            adj[j].push(i as u32);
        }
    }
    adj
}

/// Shortest path from `a` to `b` inclusive of both, or `None` when they are in
/// different fragments. BFS, so the path is the fewest-bonds one.
pub fn shortest_path(adj: &[Vec<u32>], a: u32, b: u32) -> Option<Vec<u32>> {
    if a as usize >= adj.len() || b as usize >= adj.len() {
        return None;
    }
    if a == b {
        return Some(vec![a]);
    }

    let mut prev = vec![u32::MAX; adj.len()];
    let mut queue = std::collections::VecDeque::from([a]);
    prev[a as usize] = a;

    while let Some(cur) = queue.pop_front() {
        for &next in &adj[cur as usize] {
            if prev[next as usize] != u32::MAX {
                continue;
            }
            prev[next as usize] = cur;
            if next == b {
                let mut path = vec![b];
                let mut walk = b;
                while walk != a {
                    walk = prev[walk as usize];
                    path.push(walk);
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(next);
        }
    }
    None
}

/// One rotatable degree of freedom: turning `moving` about the `b`–`c` axis.
///
/// `moving` is `Molecule::find_children(b, c)` — the atoms reachable from `c`
/// without going back through `b`, excluding both endpoints. `c` itself sits on
/// the axis, so leaving it out changes nothing.
#[derive(Clone, Debug)]
pub struct Dof {
    pub b: u32,
    pub c: u32,
    pub moving: Vec<u32>,
}

impl Dof {
    /// The degree of freedom for the `b`–`c` bond, or `None` when the two are
    /// not separable — a ring bond leaves every atom reachable both ways, so
    /// `find_children` comes back with the whole rest of the molecule and
    /// "rotating" it would tear the ring apart. `Bond::is_rotor()` is the
    /// caller's guard for that; this is the belt-and-braces one.
    pub fn new(mol: &Molecule, b: u32, c: u32) -> Option<Self> {
        let moving = mol.find_children(b, c);
        if moving.contains(&b) {
            return None;
        }
        Some(Self { b, c, moving })
    }

    /// Apply `radians` of rotation to `coords`.
    pub fn apply(&self, coords: &mut [[f64; 3]], radians: f64) {
        let Some(axis) = unit(sub(coords[self.c as usize], coords[self.b as usize])) else {
            return;
        };
        rotate_about_axis(coords, &self.moving, coords[self.b as usize], axis, radians);
    }
}

/// A pair of atoms that ended up closer than their van der Waals radii allow,
/// with the distance in Å.
#[derive(Clone, Copy, Debug)]
pub struct Clash {
    pub a: u32,
    pub b: u32,
    pub distance: f64,
}

/// Non-bonded atom pairs that overlap, worst first.
///
/// Only 1-2 and 1-3 neighbours are excluded. Those two distances are fixed by
/// bond lengths and angles, which no rotation here changes, so reporting them
/// would be noise the user cannot act on. 1-4 pairs are deliberately kept:
/// their separation *is* the torsion angle, so they are the first thing a
/// dihedral rotation can drive into itself.
pub fn clashes(mol: &Molecule, coords: &[[f64; 3]]) -> Vec<Clash> {
    let n = coords.len();
    let radii: Vec<f64> = mol
        .atoms()
        .map(|a| openbabel::elements::vdw_radius(a.atomic_number()))
        .collect();
    if radii.len() != n {
        return Vec::new();
    }
    let atoms: Vec<_> = mol.atoms().collect();

    let mut out = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let limit = CLASH_FACTOR * (radii[i] + radii[j]);
            let d = distance(coords[i], coords[j]);
            if d >= limit {
                continue;
            }
            let (ai, aj) = (&atoms[i], &atoms[j]);
            if ai.is_connected(aj) || ai.is_one_three(aj) {
                continue;
            }
            out.push(Clash {
                a: i as u32,
                b: j as u32,
                distance: d,
            });
        }
    }
    out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
    out
}

/// Whether every atom lies in one plane, to within `PLANAR_TOLERANCE`.
///
/// A structure file can declare three dimensions and still hold a flat drawing:
/// the Mol2 and SDF readers set the dimension from the format, not from the
/// numbers, so a 2D sketch saved as Mol2 arrives with `has_3d()` true and every
/// z at zero. Nothing downstream survives that — the force field sees atoms on
/// top of each other, and a dihedral has no geometry to turn — so the caller
/// needs to know before adopting it.
///
/// The test is exact rather than a z-only check, because a flat structure can
/// lie in any plane, not just the xy one.
pub fn is_planar(coords: &[[f64; 3]]) -> bool {
    if coords.len() < 4 {
        return true; // three points always share a plane
    }
    let n = coords.len() as f64;
    let centre = coords.iter().fold([0.0; 3], |acc, p| {
        [acc[0] + p[0] / n, acc[1] + p[1] / n, acc[2] + p[2] / n]
    });

    // Two spanning directions: the atom farthest from the centre, then the one
    // with the largest component perpendicular to it. Their cross product is the
    // plane's normal if there is a plane at all.
    let offsets: Vec<[f64; 3]> = coords.iter().map(|p| sub(*p, centre)).collect();
    let Some(first) = offsets
        .iter()
        .max_by(|a, b| norm(**a).total_cmp(&norm(**b)))
        .copied()
        .and_then(unit)
    else {
        return true; // every atom at the same point
    };

    let mut normal = None;
    let mut best = 0.0;
    for v in &offsets {
        let n = cross(first, *v);
        let len = norm(n);
        if len > best {
            best = len;
            normal = unit(n);
        }
    }
    let Some(normal) = normal else {
        return true; // collinear
    };

    offsets
        .iter()
        .all(|v| dot(*v, normal).abs() < PLANAR_TOLERANCE)
}

/// Flatten to the `[x0, y0, z0, x1, …]` layout `Molecule::set_coordinates` wants.
pub fn flatten(coords: &[[f64; 3]]) -> Vec<f64> {
    coords.iter().flat_map(|p| p.iter().copied()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `unit()` is the gate every rotation goes through, so a NaN it lets past
    /// ends up written into the molecule by `dihedral::commit`. The natural
    /// `n < 1e-12` spelling passes NaN straight through, which is why this is
    /// pinned rather than left to the reader.
    #[test]
    fn a_direction_that_is_not_a_number_has_no_unit_vector() {
        assert_eq!(unit([f64::NAN, 0.0, 0.0]), None);
        assert_eq!(unit([0.0, f64::INFINITY, 0.0]), None);
        assert_eq!(unit([0.0, 0.0, 0.0]), None, "no direction at all");
        assert!(unit([3.0, 0.0, 4.0]).is_some());
    }

    /// A quarter turn about z sends +x to +y. Checking the rotation against a
    /// hand-computable case is what keeps the sign convention honest; every
    /// other test here compares a result against another result.
    #[test]
    fn quarter_turn_about_z() {
        let mut coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        rotate_about_axis(
            &mut coords,
            &[1],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2,
        );
        assert!(
            distance(coords[1], [0.0, 1.0, 0.0]) < 1e-12,
            "{:?}",
            coords[1]
        );
        assert_eq!(coords[0], [0.0, 0.0, 0.0], "an atom not in the set moved");
    }

    /// The case this exists for: a 2D sketch saved as Mol2 comes back with the
    /// dimension declared as 3 and every z at zero. Importing that as a
    /// conformation gives overlapping atoms and an energy in the hundreds of
    /// millions, so the importer tests the coordinates instead of the claim.
    #[test]
    fn a_flat_structure_is_recognised() {
        let flat = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [1.5, 1.5, 0.0],
            [0.0, 1.5, 0.0],
            [3.0, 0.7, 0.0],
        ];
        assert!(is_planar(&flat));

        // The same points in a tilted plane: still flat, and a z-only check
        // would miss it.
        let tilted: Vec<[f64; 3]> = flat.iter().map(|p| [p[0], p[1], p[0] * 0.5]).collect();
        assert!(is_planar(&tilted), "a tilted plane is still a plane");
    }

    /// A real conformation must not be mistaken for a drawing, or every import
    /// would be silently rebuilt and the user's geometry thrown away.
    #[test]
    fn a_real_conformation_is_not_planar() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = openbabel::Molecule::parse("CCCC", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        assert!(!is_planar(&crate::bridge::coordinates(&mol)));
    }

    /// `find_children` gives one side of a bond; the path between two atoms is
    /// what picks which bonds to turn. Both underpin the solver, so a chain of
    /// known length is worth pinning.
    #[test]
    fn shortest_path_follows_the_chain() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = openbabel::Molecule::parse("CCCC", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        let adj = adjacency(&mol);

        let carbons: Vec<u32> = mol
            .atoms()
            .filter(|a| a.atomic_number() == 6)
            .map(|a| a.index())
            .collect();
        let path = shortest_path(&adj, carbons[0], carbons[3]).expect("connected");
        assert_eq!(path.len(), 4, "butane's carbons are four bonds end to end");
        assert_eq!(path.first(), Some(&carbons[0]));
        assert_eq!(path.last(), Some(&carbons[3]));
    }
}
