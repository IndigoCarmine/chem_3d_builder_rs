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
            a.pos[0] as f64,
            a.pos[1] as f64,
            0.0,
            a.element,
        ));
    }

    // Bond block indices are 1-based positions in the atom block above.
    let index_of = |id: u32| src.atoms.iter().position(|a| a.id == id).map(|i| i + 1);
    for b in &src.bonds {
        let (Some(i), Some(j)) = (index_of(b.begin), index_of(b.end)) else {
            continue;
        };
        s.push_str(&format!("{:>3}{:>3}{:>3}  0  0  0  0\n", i, j, order_num(&b.order)));
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

/// Convert the editor's 2D structure into a 3D OpenBabel molecule: explicit
/// hydrogens added and an initial geometry generated.
pub fn to_ob_molecule(src: &chembuider_rs::Molecule) -> Result<Molecule, String> {
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

    let mut mol = Molecule::parse(&mol_block(src), "mol")
        .map_err(|e| format!("構造の解釈に失敗しました: {e}"))?;

    // The reader gave us implicit hydrogens; MM and the viewer need real atoms.
    mol.add_hydrogens();
    if !mol.generate_3d() {
        return Err("3D 構造の生成に失敗しました。".to_string());
    }

    Ok(mol)
}

/// Bulk coordinate read — one lock, versus a lock and three FFI calls per atom
/// through `Atom::coords()`.
fn coordinates(mol: &Molecule) -> Vec<[f64; 3]> {
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

    let mut result = moleucle_3dview_rs::Molecule::default();
    result.atoms = atoms;
    result.bonds = bonds;
    result
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

    #[test]
    fn methane_gets_four_hydrogens() {
        let mut mol = chembuider_rs::Molecule::default();
        carbon(&mut mol, 0.0, 0.0);
        let ob = to_ob_molecule(&mol).expect("build");
        assert_eq!(ob.formula(), "CH4");
        assert_eq!(ob.num_atoms(), 5, "C + 4H");
    }

    #[test]
    fn ethene_double_bond_leaves_two_h_per_carbon() {
        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, BondOrder::Double);
        let ob = to_ob_molecule(&mol).expect("build");
        assert_eq!(ob.formula(), "C2H4", "sp2 C contributes 2 H each");
    }

    #[test]
    fn generated_structure_is_3d() {
        let mut mol = chembuider_rs::Molecule::default();
        let a = carbon(&mut mol, 0.0, 0.0);
        let b = carbon(&mut mol, 1.4, 0.0);
        mol.add_bond(a, b, BondOrder::Single);
        let ob = to_ob_molecule(&mol).expect("build");
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
        let mut mol = chembuider_rs::Molecule::default();
        atom(&mut mol, "N", 0.0, 0.0, 1);
        let ob = to_ob_molecule(&mol).expect("build");
        assert_eq!(ob.total_charge(), 1);
        assert_eq!(ob.formula(), "H4N+", "ammonium takes four hydrogens");
    }

    #[test]
    fn unknown_element_is_rejected() {
        let mut mol = chembuider_rs::Molecule::default();
        atom(&mut mol, "Xx", 0.0, 0.0, 0);
        assert!(to_ob_molecule(&mol).is_err());
    }

    #[test]
    fn empty_structure_is_rejected() {
        let mol = chembuider_rs::Molecule::default();
        assert!(to_ob_molecule(&mol).is_err());
    }

    /// Benzene drawn as an alternating ring: OpenBabel perceives the aromaticity
    /// the old exporter never did (which is what earns `C.ar` in MOL2), while
    /// `order` stays Kekulé so the viewer can draw it.
    #[test]
    fn benzene_is_perceived_aromatic_and_reaches_the_viewer() {
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

        let ob = to_ob_molecule(&ring).expect("build");
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
