//! Export the current 3D structure to common chemistry file formats.
//!
//! Both writers take the element list, explicit bonds and current Ångström
//! coordinates of a [`Mol3D`] and produce a single-model file:
//!   - [`to_mol2`]  : Tripos MOL2 — ATOM block (element, xyz, a coarse SYBYL
//!                    atom type) plus a BOND block carrying bond orders.
//!   - [`to_pdb`]   : PDB — column-exact `HETATM` records plus `CONECT` bond
//!                    records (one MODEL, `UNL` residue).
//!
//! Coordinates are written in Ångström (the units MM works in) with no unit
//! conversion. The SYBYL/PDB atom names are derived from a per-element counter
//! so they read as `C1`, `C2`, `H1`, … regardless of interleaving.

use std::collections::HashMap;

use crate::bridge::Mol3D;

/// Assigns per-element 1-based indices so atoms read as `C1`, `H1`, `C2`, …
struct ElementNamer {
    counts: HashMap<String, usize>,
}

impl ElementNamer {
    fn new() -> Self {
        Self {
            counts: HashMap::new(),
        }
    }

    /// Next 1-based index for `element`, advancing its counter.
    fn next_index(&mut self, element: &str) -> usize {
        let c = self.counts.entry(element.to_string()).or_insert(0);
        *c += 1;
        *c
    }

    /// A free-form atom name (`C1`, `H12`, …) for MOL2.
    fn name(&mut self, element: &str) -> String {
        format!("{}{}", element, self.next_index(element))
    }
}

/// Highest bond order incident on each atom — used to pick a coarse SYBYL type.
fn max_bond_orders(mol: &Mol3D) -> Vec<u8> {
    let mut max = vec![0u8; mol.elements.len()];
    for &(a, b, o) in &mol.bonds {
        if o > max[a] {
            max[a] = o;
        }
        if o > max[b] {
            max[b] = o;
        }
    }
    max
}

/// A coarse SYBYL atom type from the element and its highest incident bond
/// order. We do not detect aromaticity, so aromatic rings come out as `C.2`
/// / `N.2` rather than `C.ar`; adequate for most viewers.
fn sybyl_type(element: &str, max_order: u8) -> &'static str {
    match element {
        "C" => match max_order {
            3 => "C.1",
            2 => "C.2",
            _ => "C.3",
        },
        "N" => match max_order {
            3 => "N.1",
            2 => "N.2",
            _ => "N.3",
        },
        "O" => match max_order {
            2 => "O.2",
            _ => "O.3",
        },
        "S" => "S.3",
        "P" => "P.3",
        "H" => "H",
        "B" => "B",
        "F" => "F",
        "Si" => "Si",
        "Cl" => "Cl",
        "Br" => "Br",
        "I" => "I",
        // Fall back to a plain element for anything we don't special-case.
        _ => "Du",
    }
}

/// Serialize the structure as a Tripos MOL2 string.
pub fn to_mol2(mol: &Mol3D, positions: &[[f64; 3]]) -> String {
    let n_atoms = mol.elements.len();
    let n_bonds = mol.bonds.len();
    let max_order = max_bond_orders(mol);
    let mut names = ElementNamer::new();

    let mut s = String::new();
    s.push_str("@<TRIPOS>MOLECULE\n");
    s.push_str("MOL\n");
    s.push_str(&format!("{n_atoms} {n_bonds} 0 0 0\n"));
    s.push_str("SMALL\n");
    s.push_str("NO_CHARGES\n\n");

    s.push_str("@<TRIPOS>ATOM\n");
    for i in 0..n_atoms {
        let el = &mol.elements[i];
        let name = names.name(el);
        let [x, y, z] = positions[i];
        let ty = sybyl_type(el, max_order[i]);
        s.push_str(&format!(
            "{:>7} {:<8} {:>10.4} {:>10.4} {:>10.4} {:<5} {:>4} {:<8} {:>9.4}\n",
            i + 1,
            name,
            x,
            y,
            z,
            ty,
            1,
            "MOL",
            0.0,
        ));
    }

    s.push_str("@<TRIPOS>BOND\n");
    for (bi, &(a, b, o)) in mol.bonds.iter().enumerate() {
        s.push_str(&format!("{:>6} {:>5} {:>5} {:>4}\n", bi + 1, a + 1, b + 1, o));
    }

    s
}

/// Four-character PDB atom-name field (columns 13–16). Single-letter elements
/// keep column 13 blank (` C1 `); the field is truncated to four chars for
/// very high indices.
fn pdb_atom_name(element: &str, idx: usize) -> String {
    let raw = format!("{element}{idx}");
    if element.len() == 1 && raw.len() <= 3 {
        format!(" {raw:<3}")
    } else {
        let mut t = raw;
        t.truncate(4);
        format!("{t:<4}")
    }
}

/// Serialize the structure as a single-model PDB string.
pub fn to_pdb(mol: &Mol3D, positions: &[[f64; 3]]) -> String {
    let n = mol.elements.len();
    let mut names = ElementNamer::new();

    let mut s = String::new();
    s.push_str("HEADER    chem_3d_builder_rs export\n");

    for i in 0..n {
        let el = &mol.elements[i];
        let name = pdb_atom_name(el, names.next_index(el));
        let [x, y, z] = positions[i];
        // Column-exact HETATM record (PDB v3.3). The literal spaces reproduce
        // the blank columns (12, 17, 21, 27, 28–30, 67–76).
        s.push_str(&format!(
            "HETATM{serial:>5} {name} {res:<3} {chain}{resseq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{temp:>6.2}          {elem:>2}\n",
            serial = i + 1,
            name = name,
            res = "UNL",
            chain = "A",
            resseq = 1,
            x = x,
            y = y,
            z = z,
            occ = 1.0,
            temp = 0.0,
            elem = el.to_uppercase(),
        ));
    }

    // CONECT records: list every atom's bonded partners, up to four per line.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b, _) in &mol.bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    for (i, neighbors) in adj.iter().enumerate() {
        for chunk in neighbors.chunks(4) {
            let mut line = format!("CONECT{:>5}", i + 1);
            for &j in chunk {
                line.push_str(&format!("{:>5}", j + 1));
            }
            line.push('\n');
            s.push_str(&line);
        }
    }

    s.push_str("END\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Methane: C at origin with four hydrogens, four C–H single bonds.
    fn methane() -> (Mol3D, Vec<[f64; 3]>) {
        let mol = Mol3D {
            elements: vec![
                "C".into(),
                "H".into(),
                "H".into(),
                "H".into(),
                "H".into(),
            ],
            positions: vec![[0.0; 3]; 5],
            bonds: vec![(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)],
        };
        let pos = vec![
            [0.000, 0.000, 0.000],
            [0.629, 0.629, 0.629],
            [-0.629, -0.629, 0.629],
            [-0.629, 0.629, -0.629],
            [0.629, -0.629, -0.629],
        ];
        (mol, pos)
    }

    #[test]
    fn mol2_has_correct_counts_and_blocks() {
        let (mol, pos) = methane();
        let out = to_mol2(&mol, &pos);
        assert!(out.contains("@<TRIPOS>MOLECULE"));
        assert!(out.contains("@<TRIPOS>ATOM"));
        assert!(out.contains("@<TRIPOS>BOND"));
        // "5 4 0 0 0" counts line (5 atoms, 4 bonds).
        assert!(out.contains("5 4 0 0 0"), "counts line missing:\n{out}");
        // sp3 carbon → C.3
        assert!(out.contains("C.3"), "expected C.3 SYBYL type:\n{out}");
        // One ATOM line per atom, one BOND line per bond.
        assert_eq!(
            out.lines().filter(|l| l.contains("MOL ")).count(),
            5,
            "expected 5 atom lines"
        );
    }

    #[test]
    fn pdb_columns_and_conect() {
        let (mol, pos) = methane();
        let out = to_pdb(&mol, &pos);
        let hetatm: Vec<&str> = out.lines().filter(|l| l.starts_with("HETATM")).collect();
        assert_eq!(hetatm.len(), 5, "expected 5 HETATM records");
        // Column check on the carbon record: element right-justified in 77-78,
        // coordinates in the 8.3 fields.
        let carbon = hetatm[0];
        assert_eq!(&carbon[0..6], "HETATM");
        assert!(carbon.contains("UNL"), "residue name missing: {carbon}");
        assert!(carbon.trim_end().ends_with('C'), "element col: {carbon}");
        // Carbon bonds to 4 hydrogens → its CONECT lists all four.
        assert!(
            out.contains("CONECT    1    2    3    4    5"),
            "carbon CONECT missing:\n{out}"
        );
        assert!(out.trim_end().ends_with("END"));
    }
}
