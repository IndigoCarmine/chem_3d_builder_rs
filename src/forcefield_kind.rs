//! Selectable force fields, mapped to the `zunda_rs` implementations.

use zunda_rs::{
    ForceField, GaffForceField, GhemicalForceField, Mm2ForceField, Mmff94ForceField, UffForceField,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FfKind {
    Uff,
    Mmff94,
    Gaff,
    Ghemical,
    Mm2,
}

impl FfKind {
    pub const ALL: [FfKind; 5] = [
        FfKind::Uff,
        FfKind::Mmff94,
        FfKind::Gaff,
        FfKind::Ghemical,
        FfKind::Mm2,
    ];

    pub fn label(&self) -> &'static str {
        match self {
            FfKind::Uff => "UFF",
            FfKind::Mmff94 => "MMFF94",
            FfKind::Gaff => "GAFF",
            FfKind::Ghemical => "Ghemical",
            FfKind::Mm2 => "MM2",
        }
    }

    /// Construct a fresh, un-setup force field of this kind.
    pub fn make(&self) -> Box<dyn ForceField> {
        match self {
            FfKind::Uff => Box::new(UffForceField::new()),
            FfKind::Mmff94 => Box::new(Mmff94ForceField::new()),
            FfKind::Gaff => Box::new(GaffForceField::new()),
            FfKind::Ghemical => Box::new(GhemicalForceField::new()),
            FfKind::Mm2 => Box::new(Mm2ForceField::new()),
        }
    }
}

impl Default for FfKind {
    fn default() -> Self {
        FfKind::Uff
    }
}
