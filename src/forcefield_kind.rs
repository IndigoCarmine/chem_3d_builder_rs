//! Selectable force fields, mapped to OpenBabel's force-field plugin ids.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FfKind {
    #[default]
    Uff,
    Mmff94,
    Mmff94s,
    Gaff,
    Ghemical,
}

impl FfKind {
    pub const ALL: [FfKind; 5] = [
        FfKind::Uff,
        FfKind::Mmff94,
        FfKind::Mmff94s,
        FfKind::Gaff,
        FfKind::Ghemical,
    ];

    pub fn label(&self) -> &'static str {
        self.ob_id()
    }

    /// The plugin id OpenBabel resolves this force field by. These five are the
    /// whole runnable line-up: OpenBabel 3.2.1 excludes `forcefieldmm2.cpp` from
    /// its build, so MM2 has no implementation to select.
    pub fn ob_id(&self) -> &'static str {
        match self {
            FfKind::Uff => "UFF",
            FfKind::Mmff94 => "MMFF94",
            FfKind::Mmff94s => "MMFF94s",
            FfKind::Gaff => "GAFF",
            FfKind::Ghemical => "Ghemical",
        }
    }

    /// The unit the force field reports energies in. Not uniform across the
    /// line-up, so an energy is only meaningful next to its own force field's
    /// unit. Hardcoded because reading it from OpenBabel takes the global lock,
    /// which the read-out cannot afford per frame — `energy_unit_matches_openbabel`
    /// keeps these honest.
    pub fn energy_unit(&self) -> &'static str {
        match self {
            FfKind::Mmff94 | FfKind::Mmff94s => "kcal/mol",
            FfKind::Uff | FfKind::Gaff | FfKind::Ghemical => "kJ/mol",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_id_resolves_and_unit_matches_openbabel() {
        let _ob = crate::test_support::ob_guard();
        for k in FfKind::ALL {
            assert_eq!(
                openbabel::forcefield_energy_unit(k.ob_id()).as_deref(),
                Some(k.energy_unit()),
                "{} reports a different unit than we label it with",
                k.label()
            );
        }
    }

    /// OpenBabel 3.2.1 ships no runnable MM2, which is why `FfKind` dropped it.
    #[test]
    fn mm2_is_not_available() {
        let _ob = crate::test_support::ob_guard();
        assert_eq!(openbabel::forcefield_energy_unit("MM2"), None);
    }
}
