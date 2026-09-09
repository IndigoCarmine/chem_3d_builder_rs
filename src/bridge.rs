//! Bridge layer between the three crates.
//!
//! The 2D editor (`chembuider_rs`) produces a flat structure with **no explicit
//! hydrogens** and z = 0. OpenBabel does the missing 2D → 3D work — hydrogens
//! from real valence rules and an initial geometry from fragment templates plus
//! a force-field cleanup (`gen3d`) — and stays the single source of truth from
//! then on: minimization, energies and every export format read the same
//! `openbabel::Molecule`. The 3D viewer (`moleucle_3dview_rs`) wants nm units
//! (Å × 0.1), so only that last hop converts.
//!
//! The structure reaches OpenBabel as an MDL V2000 mol block rather than through
//! `add_atom`/`add_bond`, because that is the only way to get its valence model.
//! An `OBAtom`'s implicit-hydrogen count is a stored property, not a derived one:
//! atoms built through the API have zero, `add_hydrogens` only makes an existing
//! count explicit, and there is no "assign typical valence" call. The MDL reader
//! is what applies the model, so this hands it a block and lets it do the
//! chemistry — including the charge-aware part the old valence table got wrong
//! (N → NH₃, but N⁺ → NH₄⁺).

use chembuider_rs::BondOrder;
use openbabel::Molecule;

fn order_num(o: &BondOrder) -> u32 {
    match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
    }
}

/// Serialize the editor's structure as an MDL V2000 mol block.
///
/// Every atom's valence field is left 0 — "no marking" — which is what makes
/// OpenBabel fill the default valence with implicit hydrogens. Writing the
/// atom's real valence there (as OpenBabel's own mol writer does) instead pins
/// it to the hydrogen-less skeleton we started from.
fn mol_block(src: &chembuider_rs::Molecule) -> String {
    let mut s = String::new();
    // Title / program / comment lines.
    s.push_str("\n  chem_3d_builder_rs\n\n");
    s.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        src.atoms.len(),
        src.bonds.len()
    ));

    for a in &src.atoms {
        s.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            a.pos[0] as f64, a.pos[1] as f64, 0.0, a.element,
        ));
    }

    // Bond block indices are 1-based positions in the atom block above.
    let index_of = |id: u32| src.atoms.iter().position(|a| a.id == id).map(|i| i + 1);
    for b in &src.bonds {
        let (Some(i), Some(j)) = (index_of(b.begin), index_of(b.end)) else {
            continue;
        };
        s.push_str(&format!(
            "{:>3}{:>3}{:>3}  0  0  0  0\n",
            i,
            j,
            order_num(&b.order)
        ));
    }

    // `M  CHG` properties rather than the atom block's legacy charge column,
    // eight per line.
    let charged: Vec<(usize, i8)> = src
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.charge != 0)
        .map(|(i, a)| (i + 1, a.charge))
        .collect();
    for chunk in charged.chunks(8) {
        s.push_str(&format!("M  CHG{:>3}", chunk.len()));
        for (idx, charge) in chunk {
            s.push_str(&format!("{idx:>4}{charge:>4}"));
        }
        s.push('\n');
    }

    s.push_str("M  END\n");
    s
}

/// A converted structure, and anything the caller should tell the user about it.
pub struct Converted {
    pub mol: Molecule,
    /// Set when the builder never managed a clean geometry and the structure
    /// handed back has atoms closer together than they should be. The molecule
    /// is still usable — this is a note, not a failure.
    pub warning: Option<String>,
}

/// Convert the editor's 2D structure into a 3D OpenBabel molecule: explicit
/// hydrogens added and an initial geometry generated.
///
/// `ff_id` is the force field the caller is about to show an energy for. A
/// less-than-clean geometry is only accepted if that particular field can make a
/// real number out of it: how close two atoms may get before `1/r` stops meaning
/// anything is the force field's business rather than something a distance
/// threshold here can decide, and the user can switch fields between presses.
pub fn to_ob_molecule(src: &chembuider_rs::Molecule, ff_id: &str) -> Result<Converted, String> {
    if src.atoms.is_empty() {
        return Err("構造が空です。まず左側で原子を描いてください。".to_string());
    }

    // Check symbols before serializing: `atomic_number` reports an unknown
    // element as 0 rather than an error, and the mol block would carry the bad
    // symbol into a vaguer parse failure.
    for a in &src.atoms {
        if openbabel::elements::atomic_number(&a.element) == 0 {
            return Err(format!("未知の元素記号です: {}", a.element));
        }
    }

    // Two atoms drawn on top of each other send OpenBabel's builder to NaN, and
    // it reports success anyway (see below). Catching it here is what lets us
    // name the actual mistake instead of describing the symptom.
    if let Some((a, b)) = coincident_atoms(src) {
        return Err(format!(
            "原子 {a} と {b} がほぼ同じ位置にあります。重なった原子を削除してください。"
        ));
    }

    let block = mol_block(src);

    // OpenBabel's 3D builder is not deterministic, and on some perfectly ordinary
    // structures it fails outright a fraction of the time — 5-benzylidenebarbituric
    // acid, for one, lands on overlapping atoms or NaN coordinates in roughly one
    // attempt in eight. The failures are independent, so simply building again
    // clears them; what is *not* acceptable is handing a broken geometry on,
    // which is what happened before `validate_geometry` existed (see it for why
    // OpenBabel's own success flag cannot be trusted).
    //
    // Each attempt re-parses, because a failed `generate_3d` leaves its wreckage
    // in the molecule it was given.
    const ATTEMPTS: usize = 5;
    let mut last_err = String::new();
    // The best structure seen so far that is usable even if not clean, kept so a
    // run of unlucky attempts still hands the user something to work with.
    let mut fallback: Option<(Molecule, String)> = None;

    for _ in 0..ATTEMPTS {
        let mut mol =
            Molecule::parse(&block, "mol").map_err(|e| format!("構造の解釈に失敗しました: {e}"))?;

        // The reader gave us implicit hydrogens; MM and the viewer need real atoms.
        mol.add_hydrogens();
        if !mol.generate_3d() {
            last_err = "3D 構造の生成に失敗しました。".to_string();
            continue;
        }

        let note = match inspect_geometry(&mol) {
            Geometry::Clean => None,
            // Close contacts are ugly, not fatal: the viewer draws the structure
            // and the force field pulls the atoms apart.
            Geometry::Cramped(msg) => Some(msg),
            Geometry::Unusable(msg) => {
                last_err = msg;
                continue;
            }
        };

        // The geometry checks catch the shapes we know about; the energy catches
        // the rest. Measured over 400 builds of a structure OBBuilder struggles
        // with, one press in two hundred came back geometrically *clean* and yet
        // carried an energy of 1e8 — a strain the distance tests have no way to
        // see. Since the energy is the number the user is about to be shown, it
        // is also the right thing to judge the attempt by.
        let Some(_) = usable_energy(&mol, ff_id) else {
            last_err = "3D 構造の生成に失敗しました（エネルギーが異常な値になりました）。\
                        もう一度お試しください。"
                .to_string();
            continue;
        };

        match note {
            None => return Ok(Converted { mol, warning: None }),
            Some(msg) => {
                last_err = msg.clone();
                fallback.get_or_insert((mol, msg));
            }
        }
    }

    // Every attempt was flawed. Prefer handing back a cramped structure over
    // refusing outright — the user asked for a 3D structure, and this one draws,
    // minimizes and exports; it just wants the force field run over it.
    if let Some((mol, msg)) = fallback {
        return Ok(Converted {
            mol,
            warning: Some(format!(
                "{msg}\n「▶ MM 最小化」で解消できることがあります。"
            )),
        });
    }

    Err(format!(
        "{last_err}\n（{ATTEMPTS} 回試行しました。もう一度お試しいただくか、2D 構造を描き直してください。）"
    ))
}

/// A structure this app builds is drawn by hand in the 2D editor, so it runs to
/// tens of atoms and lands in the tens or hundreds of kJ/mol. The failed builds
/// this guards against come back at 1e6 and above, so there are two orders of
/// magnitude of daylight and the exact ceiling is not delicate. Force fields
/// differ in unit (kcal/mol against kJ/mol) but only by a factor of four, which
/// at this scale changes nothing.
const ENERGY_CEILING: f64 = 1.0e5;

/// The molecule's energy under `ff_id`, if it is a number the app can show and
/// minimize from.
///
/// `Molecule::energy` returns `Some(NaN)` and `Some(inf)` as ordinary values —
/// the FFI marks the call successful before it has looked at the result — so
/// every caller has to check, and this is where that happens.
fn usable_energy(mol: &Molecule, ff_id: &str) -> Option<f64> {
    mol.energy(ff_id)
        .filter(|e| e.is_finite() && e.abs() <= ENERGY_CEILING)
}

/// What OpenBabel actually produced, as opposed to what it reported.
enum Geometry {
    Clean,
    /// Usable, but with atoms closer than any real contact. Draws, minimizes and
    /// exports; it just looks wrong until a force field has been run over it.
    Cramped(String),
    /// Nothing downstream can do anything with this.
    Unusable(String),
}

/// Sort a geometry OpenBabel reported as a success into the three cases above.
///
/// `generate_3d` returning true is not enough. When OBBuilder cannot place an
/// atom it fails in one of two ways, and calls both of them success:
///
/// * it writes **NaN coordinates** (printing "There exists NaN in calculated
///   coordinates" to stderr and nothing else), or
/// * it writes finite coordinates with **two atoms on the same point**, which
///   stays invisible until a force field divides by that zero distance and every
///   energy comes back NaN.
///
/// Either way the structure flows on to the viewer, which draws nothing, and to
/// the energy badge, which reads `E NaN` — with no error raised anywhere. That
/// silent failure is precisely what "the structure won't go 3D" looks like from
/// the outside, so both stay `Unusable`.
///
/// Atoms that are merely *too close* are a different matter. They come from the
/// same unlucky builder run, but the distance is still non-zero, so the energy
/// is finite, the viewer draws the structure and minimization pulls it apart.
/// Refusing those meant a structure the user could have fixed in one click was
/// thrown away instead, so they are `Cramped` and flow on with a warning.
/// Name an atom the way the rest of the app does — element symbol plus its
/// 1-based position, `C9` rather than `原子 9`.
///
/// A bare number is not something the user can find: nothing in the 2D editor
/// shows atom indices, and `add_hydrogens` runs before `generate_3d`, so a high
/// number is often a hydrogen they never drew. The symbol at least says which
/// kind of atom is meant, and marks the added hydrogens as hydrogens.
fn atom_name(atom: &openbabel::Atom<'_>, index: usize) -> String {
    format!(
        "{}{}",
        openbabel::elements::symbol(atom.atomic_number()),
        index + 1
    )
}

fn inspect_geometry(mol: &Molecule) -> Geometry {
    let coords = coordinates(mol);

    // A conformer read that comes back short leaves every loop below vacuously
    // satisfied, so the structure would be pronounced Clean and then panic on
    // the first index in `to_viewer_molecule`.
    if coords.len() != mol.num_atoms() as usize {
        return Geometry::Unusable(
            "3D 構造の生成に失敗しました（座標を読み出せませんでした）。".to_string(),
        );
    }

    if !coords.iter().flatten().all(|v| v.is_finite()) {
        return Geometry::Unusable(
            "3D 構造の生成に失敗しました（座標が NaN になりました）。\
             原子の重なりや無理な結合がないか確認してください。"
                .to_string(),
        );
    }

    // No molecule this app builds is anywhere near this wide. When OBBuilder
    // cannot place an atom it sometimes flings it off to a finite but absurd
    // coordinate rather than writing NaN — values up to 1e243 Å have been
    // observed — and the rest of the structure then collapses in on itself. That
    // collapse on its own reads as a merely tight contact, so without this bound
    // the runaway is misfiled as recoverable and handed on, and its energy is
    // `inf`.
    const MAX_COORD: f64 = 1.0e4;
    if let Some(v) = coords
        .iter()
        .flatten()
        .find(|v| v.abs() > MAX_COORD)
        .copied()
    {
        return Geometry::Unusable(format!(
            "3D 構造の生成に失敗しました（原子が {v:.3e} Å に飛びました）。\
             もう一度お試しください。"
        ));
    }

    // Below this a pair is effectively one point and every force field's `1/r`
    // terms stop meaning anything.
    const COINCIDENT: f64 = 1e-3;
    // Bonded atoms are legitimately this close — a C–H bond is about 1.09 Å — so
    // what actually signals a failed build is two atoms that are *not* bonded
    // sitting on top of each other. Nothing non-bonded comes near 1.2 Å in a
    // real structure. Measured over 8000 builds, the bad ones land between 0.5
    // and 1.1 Å, which the old flat 0.5 Å bound waved through as Clean while the
    // energy ran to 1e14.
    const MIN_NONBONDED: f64 = 1.2;

    let atoms: Vec<_> = mol.atoms().collect();
    let mut cramped: Option<String> = None;
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            let (p, q) = (coords[i], coords[j]);
            let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2);
            if d2 >= MIN_NONBONDED * MIN_NONBONDED {
                continue;
            }
            let d = d2.sqrt();
            let msg = format!(
                "{} と {} が近すぎます: {:.3} Å。",
                atom_name(&atoms[i], i),
                atom_name(&atoms[j], j),
                d
            );

            if d < COINCIDENT {
                return Geometry::Unusable(format!(
                    "3D 構造の生成に失敗しました（{msg}）。\
                     水素は自動で付加されるため、番号は描いた原子より大きくなります。\
                     2D 構造を描き直すか、原子の配置を広げてみてください。"
                ));
            }
            // Only reached for pairs already known to be close, so the cost of
            // asking OpenBabel about the topology is paid a handful of times.
            let bonded = atoms[i].is_connected(&atoms[j]) || atoms[i].is_one_three(&atoms[j]);
            if !bonded {
                cramped.get_or_insert(msg);
            }
        }
    }

    match cramped {
        Some(msg) => Geometry::Cramped(msg),
        None => Geometry::Clean,
    }
}

/// The first pair of atoms drawn close enough together to be the same point, as
/// 1-based indices for the message.
///
/// The editor has no fixed unit — a structure drawn by clicking is in the tens,
/// one built in a test is around 1.4 per bond — so the threshold is relative to
/// the shortest bond actually drawn rather than an absolute distance.
fn coincident_atoms(src: &chembuider_rs::Molecule) -> Option<(usize, usize)> {
    let dist = |p: [f32; 2], q: [f32; 2]| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();

    let index_of = |id: u32| src.atoms.iter().position(|a| a.id == id);
    let mut lengths: Vec<f32> = src
        .bonds
        .iter()
        .filter_map(|b| Some((index_of(b.begin)?, index_of(b.end)?)))
        .map(|(i, j)| dist(src.atoms[i].pos, src.atoms[j].pos))
        .collect();
    if lengths.is_empty() {
        return None; // nothing bonded yet: no scale to judge against
    }
    // The median, not the minimum: an overlapping pair is usually bonded to each
    // other, and a near-zero-length bond would otherwise set the scale by which
    // it is then judged normal.
    lengths.sort_by(|a, b| a.total_cmp(b));
    let typical = lengths[lengths.len() / 2];
    if typical <= 0.0 {
        return None;
    }
    // A tenth of a normal bond: far below any real one, far above the jitter
    // between two deliberate clicks.
    let limit = typical * 0.1;

    for i in 0..src.atoms.len() {
        for j in (i + 1)..src.atoms.len() {
            if dist(src.atoms[i].pos, src.atoms[j].pos) < limit {
                return Some((i + 1, j + 1));
            }
        }
    }
    None
}

/// Bulk coordinate read — one lock, versus a lock and three FFI calls per atom
/// through `Atom::coords()`.
pub(crate) fn coordinates(mol: &Molecule) -> Vec<[f64; 3]> {
    mol.conformer_coordinates(0).unwrap_or_else(|| {
        mol.atoms()
            .map(|a| {
                let (x, y, z) = a.coords();
                [x, y, z]
            })
            .collect()
    })
}

/// Build the viewer molecule from the current geometry. Positions are scaled
/// Å → nm (× 0.1) to match the viewer's loaders.
pub fn to_viewer_molecule(mol: &Molecule) -> moleucle_3dview_rs::Molecule {
    use lin_alg::f32::Vec3 as FVec3;
    use moleucle_3dview_rs::molecule::{Atom, Bond, Element};

    const ANGSTROM_TO_NM: f32 = 0.1;

    let coords = coordinates(mol);
    let atoms = mol
        .atoms()
        .enumerate()
        .map(|(i, a)| Atom {
            position: FVec3::new(
                coords[i][0] as f32 * ANGSTROM_TO_NM,
                coords[i][1] as f32 * ANGSTROM_TO_NM,
                coords[i][2] as f32 * ANGSTROM_TO_NM,
            ),
            element: Element::new(&openbabel::elements::symbol(a.atomic_number())),
            id: i,
            meta: None,
        })
        .collect();

    // OpenBabel keeps aromatic rings Kekulé in `order` (1/2) with aromaticity on
    // a separate flag, so these are already the 1–3 the viewer draws as parallel
    // lines. Order 5 is a file-format convention, not something `order` returns.
    let bonds = mol
        .bonds()
        .map(|b| Bond {
            atom_a: b.begin_atom_index() as usize,
            atom_b: b.end_atom_index() as usize,
            order: b.order() as u8,
        })
        .collect();

    moleucle_3dview_rs::Molecule::from_atoms_bonds(atoms, bonds)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(mol: &mut chembuider_rs::Molecule, el: &str, x: f32, y: f32, charge: i8) -> u32 {
        mol.add_atom(el.to_string(), [x, y], charge)
    }

    fn carbon(mol: &mut chembuider_rs::Molecule, x: f32, y: f32) -> u32 {
        atom(mol, "C", x, y, 0)
    }

    /// 5-benzylidenebarbituric acid: three carbonyls plus an exocyclic C=C on
    /// the ring carbon between two of them. This is the structure OpenBabel's
    /// builder fails on at random, so both the conversion test and the stress
    /// test above are built from it.
    fn benzylidene_barbituric_acid() -> chembuider_rs::Molecule {
        let mut m = chembuider_rs::Molecule::default();
        let c2 = carbon(&mut m, 0.0, 1.4);
        let n1 = atom(&mut m, "N", -1.21, 0.7, 0);
        let c6 = carbon(&mut m, -1.21, -0.7);
        let c5 = carbon(&mut m, 0.0, -1.4);
        let c4 = carbon(&mut m, 1.21, -0.7);
        let n3 = atom(&mut m, "N", 1.21, 0.7, 0);
        for (a, b) in [(c2, n1), (n1, c6), (c6, c5), (c5, c4), (c4, n3), (n3, c2)] {
            m.add_bond(a, b, BondOrder::Single);
        }
        for (c, x, y) in [(c2, 0.0, 2.8), (c6, -2.42, -1.4), (c4, 2.42, -1.4)] {
            let o = atom(&mut m, "O", x, y, 0);
            m.add_bond(c, o, BondOrder::Double);
        }
        let ch = carbon(&mut m, 0.0, -2.8);
        m.add_bond(c5, ch, BondOrder::Double);
        let ph: Vec<u32> = (0..6)
            .map(|i| {
                let a = std::f32::consts::FRAC_PI_3 * i as f32;
                carbon(&mut m, 1.21 + 1.4 * a.cos(), -3.5 - 1.4 * a.sin())
            })
            .collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            m.add_bond(ph[i], ph[(i + 1) % 6], order);
        }
        m.add_bond(ch, ph[2], BondOrder::Single);
        m
    }

    /// 5-benzylidenebarbituric acid: three carbonyls plus an exocyclic C=C on the
    /// ring carbon between two of them. Grounds that this shape converts at all,
    /// so a report of "it won't go 3D" points at the drawing, not the chemistry.
    #[test]
    fn benzylidene_barbituric_acid_converts() {
        let _ob = crate::test_support::ob_guard();
        let m = benzylidene_barbituric_acid();
        // Repeatedly, because this is exactly the structure OpenBabel's builder
        // fails on at random — about one attempt in eight lands on overlapping
        // atoms or NaN coordinates. A single pass here would go green while the
        // app still failed for the user roughly every eighth press.
        //
        // The energy assertion is the load-bearing one. A cramped geometry is
        // handed on rather than refused, so "it converted" no longer says much;
        // what has to hold is that whatever comes back is something the rest of
        // the app can use, and a NaN energy is exactly what it cannot.
        for run in 0..20 {
            let converted = to_ob_molecule(&m, "UFF")
                .unwrap_or_else(|e| panic!("run {run}: benzylidene barbituric acid failed: {e}"));
            assert_eq!(converted.mol.formula(), "C11H8N2O3");
            assert!(
                converted.mol.energy("UFF").is_some_and(|e| e.is_finite()),
                "run {run}: a converted structure must have a real energy, not NaN"
            );
            if let Some(warning) = &converted.warning {
                assert!(
                    warning.contains("近すぎます"),
                    "run {run}: a warning must say what is wrong, got {warning:?}"
                );
            }
        }
    }

    /// Does asking the same topology for its energy over and over poison
    /// OpenBabel's shared force-field object? `OBForceField::Setup` takes a fast
    /// path when the molecule looks topologically identical to the last one, and
    /// a setup that failed leaves the shared instance marked invalid — which is
    /// exactly the shape of this crate's retry loop.
    /// `cargo test -- --ignored --nocapture forcefield_singleton_poisoning`
    #[test]
    #[ignore = "measurement"]
    fn forcefield_singleton_poisoning() {
        let _ob = crate::test_support::ob_guard();
        let m = benzylidene_barbituric_acid();
        let block = mol_block(&m);
        let mut line = String::new();
        let (mut none, mut total) = (0, 0);
        for _ in 0..60 {
            let mut mol = Molecule::parse(&block, "mol").expect("parse");
            mol.add_hydrogens();
            if !mol.generate_3d() {
                line.push('g');
                continue;
            }
            total += 1;
            match mol.energy("UFF") {
                None => {
                    none += 1;
                    line.push('N');
                }
                Some(e) if !e.is_finite() => line.push('!'),
                Some(e) if e.abs() > 1.0e5 => line.push('#'),
                Some(_) => line.push('.'),
            }
        }
        println!(
            "total={total} none={none}
{line}"
        );
    }

    /// Prints the distribution of single-point energies OBBuilder produces for
    /// the awkward structure below, so `ENERGY_CEILING` is chosen from data
    /// rather than guessed. Ignored: it is a measurement, not an assertion.
    /// `cargo test -- --ignored --nocapture press_energy_distribution`
    #[test]
    #[ignore = "measurement; prints the distribution ENERGY_CEILING is drawn from"]
    fn press_energy_distribution() {
        let _ob = crate::test_support::ob_guard();
        let m = benzylidene_barbituric_acid();
        let block = mol_block(&m);

        let mut energies: Vec<f64> = Vec::new();
        let mut nonfinite = 0;
        for _ in 0..300 {
            let mut mol = Molecule::parse(&block, "mol").expect("parse");
            mol.add_hydrogens();
            if !mol.generate_3d() {
                continue;
            }
            if matches!(inspect_geometry(&mol), Geometry::Unusable(_)) {
                continue;
            }
            match mol.energy("UFF") {
                Some(e) if e.is_finite() => energies.push(e.abs()),
                _ => nonfinite += 1,
            }
        }
        energies.sort_by(f64::total_cmp);
        let at = |q: f64| energies[((energies.len() - 1) as f64 * q) as usize];
        println!(
            "n={} nonfinite={nonfinite} min={:.3e} p50={:.3e} p90={:.3e} p99={:.3e} max={:.3e}",
            energies.len(),
            energies[0],
            at(0.50),
            at(0.90),
            at(0.99),
            energies[energies.len() - 1],
        );
    }

    /// How often pressing "→ 3D 生成" gives back something with an absurd
    /// energy, measured rather than argued about.
    ///
    /// Ignored because it is a statistical check over 400 builds, not a unit
    /// test: run it with `cargo test -- --ignored --nocapture bad_press_rate`
    /// after touching `inspect_geometry`. The number it prints is the one the
    /// user experiences as "何回か押していると NaN になる"; before the coordinate
    /// bound and the non-bonded threshold it sat at roughly one press in
    /// seventy.
    #[test]
    #[ignore = "statistical; run explicitly after changing inspect_geometry"]
    fn bad_press_rate() {
        let _ob = crate::test_support::ob_guard();
        let m = benzylidene_barbituric_acid();

        const PRESSES: usize = 400;
        let (mut refused, mut warned, mut absurd, mut worst) = (0, 0, 0, 0.0f64);
        for _ in 0..PRESSES {
            match to_ob_molecule(&m, "UFF") {
                Err(_) => refused += 1,
                Ok(c) => {
                    if c.warning.is_some() {
                        warned += 1;
                    }
                    let e = c.mol.energy("UFF").unwrap_or(f64::NAN);
                    if !e.is_finite() || e.abs() > 1.0e6 {
                        absurd += 1;
                    }
                    if e.is_finite() && e.abs() > worst {
                        worst = e.abs();
                    }
                }
            }
        }
        println!(
            "presses {PRESSES}: refused {refused}, warned {warned}, \
             absurd energy {absurd}, worst |E| {worst:.3e}"
        );
        assert_eq!(absurd, 0, "a press must never hand back an absurd energy");
    }

    /// The failure that made pressing "→ 3D 生成" repeatedly produce a NaN
    /// energy: OBBuilder sometimes flings one atom to a finite but astronomical
    /// coordinate and lets the rest collapse together. Judged on separation
    /// alone that reads as a merely tight contact and gets handed on, and the
    /// force field then returns `inf`. The coordinate bound is what names it.
    #[test]
    fn a_runaway_coordinate_is_unusable() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = Molecule::parse("CCO", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");
        assert!(
            matches!(inspect_geometry(&mol), Geometry::Clean),
            "a freshly built ethanol should be clean"
        );

        let mut coords = coordinates(&mol);
        coords[0][0] = 1.0e208;
        assert!(
            mol.set_coordinates(&crate::geom3d::flatten(&coords)),
            "write"
        );

        match inspect_geometry(&mol) {
            Geometry::Unusable(msg) => assert!(msg.contains("飛びました"), "{msg}"),
            _ => panic!("an atom at 1e208 Å must not be handed on"),
        }
    }

    /// Bonded atoms are legitimately closer than the non-bonded threshold — a
    /// C–H bond is about 1.09 Å — so a check that did not consult the topology
    /// would call every ordinary structure cramped and warn on all of them.
    #[test]
    fn ordinary_bonds_are_not_mistaken_for_close_contacts() {
        let _ob = crate::test_support::ob_guard();
        for smiles in ["CCO", "c1ccccc1", "CC(=O)N"] {
            let mut mol = Molecule::parse(smiles, "smi").expect("parse");
            assert!(mol.generate_3d(), "gen3d {smiles}");
            assert!(
                matches!(inspect_geometry(&mol), Geometry::Clean),
                "{smiles} should be clean, its shortest bond is about 1.09 Å"
            );
        }
    }

    /// Two atoms on the same spot drive OpenBabel's builder to NaN coordinates
    /// while `generate_3d` still returns true, so the structure silently reaches
    /// the viewer as nothing at all. Both guards exist to stop that.
    #[test]
    fn overlapping_atoms_are_rejected_not_silently_nan() {
        let _ob = crate::test_support::ob_guard();
        let mut m = chembuider_rs::Molecule::default();
        let a = carbon(&mut m, 0.0, 0.0);
        let b = carbon(&mut m, 1.4, 0.0);
        m.add_bond(a, b, BondOrder::Single);
        // A third carbon all but on top of `b` — a double-click in the editor.
        let c = carbon(&mut m, 1.4 + 0.001, 0.0);
        m.add_bond(b, c, BondOrder::Single);

        let err = to_ob_molecule(&m, "UFF")
            .err()
            .expect("overlapping atoms must be reported");
        assert!(err.contains("同じ位置"), "unhelpful message: {err}");
    }

    #[test]
    fn methane_gets_four_hydrogens() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = chembuider_rs::Molecule::default();
        carbon(&mut mol, 0.0, 0.0);
        let ob = to_ob_molecule(&mol, "UFF").expect("build").mol;
        assert_eq!(ob.formula(), "CH4");
        assert_eq!(ob.num_atoms(), 5, "C + 4H");
    }

    #[test]
    fn ethene_double_bond_leaves_two_h_per_carbon() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, BondOrder::Double);
        let ob = to_ob_molecule(&mol, "UFF").expect("build").mol;
        assert_eq!(ob.formula(), "C2H4", "sp2 C contributes 2 H each");
    }

    #[test]
    fn generated_structure_is_3d() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, BondOrder::Single);
        let ob = to_ob_molecule(&mol, "UFF").expect("build").mol;
        assert_eq!(ob.dimension(), 3);
        assert!(
            coordinates(&ob).iter().any(|p| p[2].abs() > 1e-6),
            "expected non-zero z somewhere"
        );
    }

    /// The editor has always carried a formal charge; before OpenBabel it was
    /// dropped on the way to 3D.
    #[test]
    fn formal_charge_reaches_openbabel() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = chembuider_rs::Molecule::default();
        atom(&mut mol, "N", 0.0, 0.0, 1);
        let ob = to_ob_molecule(&mol, "UFF").expect("build").mol;
        assert_eq!(ob.total_charge(), 1);
        assert_eq!(ob.formula(), "H4N+", "ammonium takes four hydrogens");
    }

    #[test]
    fn unknown_element_is_rejected() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = chembuider_rs::Molecule::default();
        atom(&mut mol, "Xx", 0.0, 0.0, 0);
        assert!(to_ob_molecule(&mol, "UFF").is_err());
    }

    #[test]
    fn empty_structure_is_rejected() {
        let _ob = crate::test_support::ob_guard();
        let mol = chembuider_rs::Molecule::default();
        assert!(to_ob_molecule(&mol, "UFF").is_err());
    }

    /// Benzene drawn as an alternating ring: OpenBabel perceives the aromaticity
    /// the old exporter never did (which is what earns `C.ar` in MOL2), while
    /// `order` stays Kekulé so the viewer can draw it.
    #[test]
    fn benzene_is_perceived_aromatic_and_reaches_the_viewer() {
        let _ob = crate::test_support::ob_guard();
        let mut ring = chembuider_rs::Molecule::default();
        let ids: Vec<u32> = (0..6)
            .map(|i| {
                let angle = std::f32::consts::FRAC_PI_3 * i as f32;
                carbon(&mut ring, 1.4 * angle.cos(), 1.4 * angle.sin())
            })
            .collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            ring.add_bond(ids[i], ids[(i + 1) % 6], order);
        }

        let ob = to_ob_molecule(&ring, "UFF").expect("build").mol;
        assert_eq!(ob.formula(), "C6H6");
        assert_eq!(
            ob.bonds().filter(|b| b.is_aromatic()).count(),
            6,
            "the whole ring should come back aromatic"
        );

        let view = to_viewer_molecule(&ob);
        assert_eq!(view.atoms.len(), 12);
        assert!(
            view.bonds.iter().all(|b| (1..=3).contains(&b.order)),
            "the viewer draws `order` parallel lines, so it must stay Kekulé"
        );
        assert!(
            view.atoms.iter().any(|a| a.position.x.abs() > 1e-6),
            "expected real coordinates, not an all-zero geometry"
        );
        // Benzene spans ~3 Å, so in nm every atom sits well inside 1.0.
        assert!(
            view.atoms.iter().all(|a| a.position.x.abs() < 1.0),
            "positions should be nm, not Å"
        );
    }
}
