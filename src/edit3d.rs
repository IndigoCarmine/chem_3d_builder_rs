//! State and side panel for editing the 3D structure directly: the selection
//! made by clicking in the viewport, the dihedral rotation it arms, and the
//! distance solver that runs on a pair of atoms.
//!
//! This module only works out *what* the user asked for. Carrying it out needs
//! the molecule, the viewer and the energy badge at once, so the panel returns
//! an [`Action`] and `main.rs` applies it.

/// The most atoms any one consumer of `picked` needs: four, for a dihedral
/// restraint. The solver takes two of them, an angle restraint three.
const MAX_PICKED: usize = 4;

use crate::constraints::{Constraint, Kind};
use crate::dihedral::Selection;
use crate::ik::Goal;
use crate::theme::{self, PAL};
use crate::viewport3d::Hit;
use egui::RichText;
use moleucle_3dview_rs::{AtomGroup, AtomGroupState, AtomPairState};

/// How far one wheel notch turns a bond, in degrees. Modifiers scale it: the
/// default is for finding a rough conformation, Shift for settling on one,
/// Ctrl for sweeping through a full turn quickly.
pub const STEP_DEGREES: f64 = 5.0;
pub const STEP_FINE_DEGREES: f64 = 1.0;
pub const STEP_COARSE_DEGREES: f64 = 15.0;

/// How many geometry snapshots to keep. The 2D editor caps its undo stack at 50
/// whole-molecule clones; these are coordinate arrays for a structure that is
/// usually far smaller, but 20 is plenty for feeling out a conformation.
const UNDO_DEPTH: usize = 20;

const ANCHOR_COLOR: (f32, f32, f32) = (0.98, 0.66, 0.15);
const BOND_COLOR: (f32, f32, f32) = (0.20, 0.75, 0.95);
const IK_COLOR: (f32, f32, f32) = (0.85, 0.30, 0.80);
/// The half of the molecule the wheel would swing. Deliberately pale: it can
/// cover most of the structure, and it is context rather than a selection.
const MOVING_COLOR: (f32, f32, f32) = (0.45, 0.55, 0.70);

/// What `main.rs` should do once the panel and the viewport have been read.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Action {
    /// Turn the selected bond by this many degrees.
    Rotate(f64),
    /// Set the selected bond's torsion to this absolute angle, in degrees.
    SetAngle(f64),
    Solve(Goal),
    Undo,
}

/// Which goal the radio buttons are on. Kept apart from [`Goal`] so the typed
/// distance survives a trip through the other two options.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GoalKind {
    #[default]
    Distance,
    Far,
    Near,
}

/// What the panel needs to know about the molecule to draw itself.
///
/// Gathered by `main.rs`, because every read takes OpenBabel's global lock and
/// the panel should not be the thing deciding when that is safe.
#[derive(Clone, Debug, Default)]
pub struct Readout {
    /// Torsion of the armed rotation.
    pub torsion_degrees: Option<f64>,
    /// Why the selected bond cannot be turned, when it cannot.
    ///
    /// Clicking a ring or terminal bond is an easy mistake to make — nothing on
    /// screen distinguishes them — and without this the panel would answer by
    /// showing the same "pick an anchor and a bond" hint the user had just
    /// followed.
    pub rotation_error: Option<&'static str>,
    /// Current separation of the two solver atoms.
    pub pair_distance: Option<f64>,
    /// What a restraint of the panel's current kind would measure over the
    /// picked atoms, or `None` if too few are picked.
    pub pending_constraint: Option<f64>,
    /// What each existing restraint measures right now, in the same order as
    /// `Edit3d::constraints`, so the panel can show the drift from target.
    pub constraint_values: Vec<f64>,
}

#[derive(Default)]
pub struct Edit3d {
    pub selection: Selection,
    /// The two endpoints of `selection.bond`, resolved by `main.rs` so the
    /// panel and the overlays can name them without reaching into OpenBabel.
    bond_atoms: Option<(u32, u32)>,
    /// The atoms a turn of the wheel would move, likewise resolved upstream.
    moving: Vec<u32>,
    /// Atoms picked in the viewport, oldest first. Capped at the largest number
    /// any consumer needs — three, for an angle restraint — and a click on an
    /// atom already in the list removes it, so clicking toggles.
    ///
    /// The solver still works on a pair: it takes the newest two, which is the
    /// same pair the old two-deep list would have held after the same clicks.
    pub picked: Vec<u32>,
    /// Restraints the next minimization will run under. Kept here because the
    /// panel that edits them is this one; `main.rs` clones them into the worker.
    pub constraints: Vec<Constraint>,
    /// Which restraint the "追加" button would build.
    pub constraint_kind: Kind,
    /// Target box for the pending restraint.
    constraint_target: f64,
    /// What `constraint_target` was last seeded for: the restraint kind and the
    /// atoms picked at the time.
    ///
    /// Seeded once per selection rather than tracked every frame. A box rewritten
    /// from the molecule on every pass cannot be typed into at all — the digits
    /// are overwritten between keystrokes — and guarding that with the widget's
    /// focus does not hold, because a `DragValue` in text-entry mode does not
    /// report focus on the response the panel sees. It would also fight the user
    /// during a minimization, when the measured value is moving under them.
    constraint_seeded_for: Option<(Kind, Vec<u32>)>,
    pub goal_kind: GoalKind,
    pub target_distance: f64,
    /// Result of the last solve, shown until the next one.
    pub report: Option<String>,
    /// Geometry as it stood before each edit, newest last.
    undo: Vec<Vec<[f64; 3]>>,
    /// Absolute-angle box. Tracks the molecule except while the user is in the
    /// widget, where it would otherwise fight them for the value.
    angle_input: f64,
    angle_active: bool,
}

impl Edit3d {
    pub fn new() -> Self {
        Self {
            target_distance: 4.0,
            ..Default::default()
        }
    }

    /// Drop everything tied to one particular molecule. Called whenever the
    /// molecule is replaced rather than edited — the indices would otherwise
    /// point into a structure that no longer exists.
    pub fn reset(&mut self) {
        let (target, goal) = (self.target_distance, self.goal_kind);
        *self = Self::new();
        self.target_distance = target;
        self.goal_kind = goal;
    }

    /// Publish what `main.rs` resolved the current selection into: the bond's
    /// endpoints and the atoms that would swing. Both are cleared whenever the
    /// selection does not name a usable rotation.
    pub fn set_resolved(&mut self, bond_atoms: Option<(u32, u32)>, moving: Vec<u32>) {
        self.bond_atoms = bond_atoms;
        self.moving = moving;
    }

    /// A left click: an atom toggles in the solver pair, a bond becomes the
    /// rotation axis, empty space clears both.
    /// The newest two picked atoms — the pair the solver acts on.
    pub fn ik_pair(&self) -> Option<[u32; 2]> {
        match self.picked[..] {
            [.., a, b] => Some([a, b]),
            _ => None,
        }
    }

    /// The newest `arity` picked atoms, in pick order: what a restraint of
    /// `kind` would be built over. `None` until enough atoms are picked.
    pub fn constraint_atoms(&self, kind: Kind) -> Option<&[u32]> {
        let n = Constraint::arity(kind);
        (self.picked.len() >= n).then(|| &self.picked[self.picked.len() - n..])
    }

    pub fn on_primary(&mut self, hit: Hit) {
        match hit {
            Hit::Atom(i) => {
                let i = i as u32;
                match self.picked.iter().position(|&x| x == i) {
                    Some(pos) => {
                        self.picked.remove(pos);
                    }
                    None => {
                        self.picked.push(i);
                        if self.picked.len() > MAX_PICKED {
                            self.picked.remove(0);
                        }
                    }
                }
                self.report = None;
            }
            Hit::Bond(i) => self.selection.bond = Some(i as u32),
            Hit::Nothing => {
                self.selection.bond = None;
                self.bond_atoms = None;
                self.moving.clear();
                self.picked.clear();
            }
        }
    }

    /// A right click: an atom becomes the anchor, empty space clears it.
    pub fn on_secondary(&mut self, hit: Hit) {
        match hit {
            Hit::Atom(i) => self.selection.anchor = Some(i as u32),
            Hit::Bond(_) => {}
            Hit::Nothing => self.selection.anchor = None,
        }
    }

    /// Degrees to turn for `lines` of wheel travel — one line per notch — or
    /// `None` when nothing is armed or the wheel did not move.
    pub fn wheel_degrees(&self, lines: f32, modifiers: egui::Modifiers) -> Option<f64> {
        if lines.abs() <= f32::EPSILON || !self.selection.is_armed() {
            return None;
        }
        let step = if modifiers.shift {
            STEP_FINE_DEGREES
        } else if modifiers.command {
            STEP_COARSE_DEGREES
        } else {
            STEP_DEGREES
        };
        Some(step * lines as f64)
    }

    pub fn push_undo(&mut self, coords: Vec<[f64; 3]>) {
        if self.undo.len() == UNDO_DEPTH {
            self.undo.remove(0);
        }
        self.undo.push(coords);
    }

    pub fn pop_undo(&mut self) -> Option<Vec<[f64; 3]>> {
        self.undo.pop()
    }

    /// Forget the history but keep the selection — for when the geometry is
    /// replaced by something the snapshots no longer describe, such as the
    /// result of a minimization run.
    pub fn clear_undo(&mut self) {
        self.undo.clear();
    }

    pub fn goal(&self) -> Goal {
        match self.goal_kind {
            GoalKind::Distance => Goal::Distance(self.target_distance),
            GoalKind::Far => Goal::AsFarAsPossible,
            GoalKind::Near => Goal::AsCloseAsPossible,
        }
    }

    /// The highlight overlays for the current selection.
    ///
    /// One atom can hold two roles at once — the anchor is often an end of the
    /// solver pair. Overlapping groups stack spheres in the same place, so each
    /// atom is assigned to exactly one group, most specific first.
    pub fn overlays(&self) -> (AtomGroupState, AtomPairState) {
        let mut used: Vec<usize> = Vec::new();
        let mut groups: Vec<AtomGroup> = Vec::new();

        let mut add = |indices: &[u32], color: (f32, f32, f32)| {
            let mut fresh = Vec::new();
            for &i in indices {
                let i = i as usize;
                if !used.contains(&i) {
                    used.push(i);
                    fresh.push(i);
                }
            }
            if !fresh.is_empty() {
                groups.push(AtomGroup {
                    atom_indices: fresh,
                    color,
                });
            }
        };

        if let Some(anchor) = self.selection.anchor {
            add(&[anchor], ANCHOR_COLOR);
        }
        add(&self.picked, IK_COLOR);

        let pairs = match self.bond_atoms {
            Some((a, b)) if self.selection.bond.is_some() => {
                add(&[a, b], BOND_COLOR);
                vec![(a as usize, b as usize)]
            }
            _ => Vec::new(),
        };
        // Last, so anything already spoken for keeps its own colour.
        add(&self.moving, MOVING_COLOR);

        (
            AtomGroupState {
                groups,
                visible: true,
                opacity: 0.55,
            },
            AtomPairState::new(pairs),
        )
    }
}

/// Draw the panel. `label_of` turns an atom index into something like `C3`;
/// only `main.rs` can read the element out of the molecule.
pub fn panel(
    ui: &mut egui::Ui,
    edit: &mut Edit3d,
    readout: Readout,
    label_of: &dyn Fn(u32) -> String,
    enabled: bool,
) -> Option<Action> {
    let mut action = None;

    ui.add_enabled_ui(enabled, |ui| {
        ui.add_space(4.0);
        theme::group_label(ui, "3D 編集");
        ui.add_space(6.0);

        // ── selection ──
        row(ui, "アンカー", edit.selection.anchor.map(label_of));
        row(
            ui,
            "結合",
            edit.bond_atoms
                .filter(|_| edit.selection.bond.is_some())
                .map(|(a, b)| format!("{} – {}", label_of(a), label_of(b))),
        );
        row(
            ui,
            "対象原子",
            (!edit.picked.is_empty()).then(|| {
                edit.picked
                    .iter()
                    .map(|&i| label_of(i))
                    .collect::<Vec<_>>()
                    .join(" , ")
            }),
        );
        if ui.button("選択を解除").clicked() {
            edit.selection.clear();
            edit.bond_atoms = None;
            edit.moving.clear();
            edit.picked.clear();
            edit.report = None;
        }

        ui.add_space(10.0);
        ui.separator();
        ui.add_space(6.0);

        // ── dihedral ──
        ui.label(RichText::new("二面角").size(12.0).color(PAL.label));
        match readout.torsion_degrees {
            Some(current) => {
                // Follow the molecule unless the user is in the box; otherwise
                // every frame would snap their half-set value back.
                if !edit.angle_active {
                    edit.angle_input = current;
                }
                let dv = ui.add(
                    egui::DragValue::new(&mut edit.angle_input)
                        .speed(1.0)
                        .range(-180.0..=180.0)
                        .suffix("°"),
                );
                edit.angle_active = dv.has_focus() || dv.dragged();
                if dv.changed() {
                    action = Some(Action::SetAngle(edit.angle_input));
                }
                ui.label(
                    RichText::new(
                        "ビュー上でホイールを回すと、アンカー側を固定したまま反対側が回ります（Shift で細かく、Ctrl で粗く）。",
                    )
                    .size(11.0)
                    .color(PAL.label),
                );
            }
            None => {
                edit.angle_active = false;
                match readout.rotation_error {
                    Some(reason) => {
                        ui.label(RichText::new(reason).size(11.0).color(PAL.danger_text));
                    }
                    None => {
                        ui.label(
                            RichText::new(
                                "原子を右クリックしてアンカーに、結合を左クリックして回転軸にします。",
                            )
                            .size(11.0)
                            .color(PAL.label),
                        );
                    }
                }
            }
        }

        ui.add_space(10.0);
        ui.separator();
        ui.add_space(6.0);

        // ── inverse kinematics ──
        ui.label(
            RichText::new("逆運動学（2原子の距離）")
                .size(12.0)
                .color(PAL.label),
        );
        row(
            ui,
            "現在",
            readout.pair_distance.map(|d| format!("{d:.3} Å")),
        );

        ui.radio_value(&mut edit.goal_kind, GoalKind::Distance, "距離を指定");
        ui.horizontal(|ui| {
            ui.add_space(20.0);
            ui.add_enabled(
                edit.goal_kind == GoalKind::Distance,
                egui::DragValue::new(&mut edit.target_distance)
                    .speed(0.05)
                    .range(0.5..=100.0)
                    .suffix(" Å"),
            );
        });
        ui.radio_value(&mut edit.goal_kind, GoalKind::Far, "できるだけ遠く");
        ui.radio_value(&mut edit.goal_kind, GoalKind::Near, "できるだけ近く");

        ui.add_space(6.0);
        if ui
            .add_enabled(edit.ik_pair().is_some(), egui::Button::new("解く"))
            .on_disabled_hover_text("3D ビューで原子を2つ左クリックしてください。")
            .clicked()
        {
            action = Some(Action::Solve(edit.goal()));
        }

        if let Some(report) = &edit.report {
            ui.add_space(4.0);
            ui.label(RichText::new(report).size(11.0));
        }

        ui.add_space(10.0);
        ui.separator();
        ui.add_space(6.0);

        // ── MM constraints ──
        ui.label(
            RichText::new("MM 制約")
                .size(12.0)
                .color(PAL.label),
        );
        ui.horizontal(|ui| {
            for &kind in Kind::ALL {
                ui.radio_value(&mut edit.constraint_kind, kind, kind.label());
            }
        });

        let kind = edit.constraint_kind;
        let picked: Option<Vec<u32>> = edit.constraint_atoms(kind).map(|a| a.to_vec());
        row(
            ui,
            "対象",
            picked.as_ref().map(|a| {
                a.iter()
                    .map(|&i| label_of(i))
                    .collect::<Vec<_>>()
                    .join(" – ")
            }),
        );
        row(
            ui,
            "現在",
            readout
                .pending_constraint
                .map(|v| format!("{v:.3} {}", kind.unit())),
        );

        match (&picked, readout.pending_constraint) {
            (Some(atoms), Some(current)) => {
                // Seed on a change of selection, then leave it alone: what the
                // user typed is theirs until they pick different atoms.
                let key = (kind, atoms.clone());
                if edit.constraint_seeded_for.as_ref() != Some(&key) {
                    edit.constraint_target = current;
                    edit.constraint_seeded_for = Some(key);
                }
                ui.horizontal(|ui| {
                    ui.label(RichText::new("目標").size(11.0).color(PAL.label));
                    ui.add(
                        egui::DragValue::new(&mut edit.constraint_target)
                            .speed(if kind == Kind::Distance { 0.02 } else { 1.0 })
                            .range(kind.range())
                            .suffix(format!(" {}", kind.unit())),
                    );
                });

                ui.horizontal_wrapped(|ui| {
                    // Snapping back to the measured value is the one-click
                    // "hold this where it is", and it has to be a button now
                    // that the box no longer follows the molecule.
                    if ui
                        .small_button("現在値")
                        .on_hover_text(format!("{current:.3} {}", kind.unit()))
                        .clicked()
                    {
                        edit.constraint_target = current;
                    }
                    for (name, value) in kind.presets() {
                        if ui
                            .small_button(*name)
                            .on_hover_text(format!("{value} {}", kind.unit()))
                            .clicked()
                        {
                            edit.constraint_target = *value;
                        }
                    }
                });

                ui.add_space(4.0);
                if ui.button("＋ 制約を追加").clicked()
                    && let Some(c) = Constraint::new(kind, atoms, edit.constraint_target)
                {
                    edit.constraints.push(c);
                }
            }
            _ => {
                edit.constraint_seeded_for = None;
                ui.label(
                    RichText::new(format!(
                        "3D ビューで原子を{}つ左クリックしてください。",
                        Constraint::arity(kind)
                    ))
                    .size(11.0)
                    .color(PAL.label),
                );
            }
        }

        if !edit.constraints.is_empty() {
            ui.add_space(6.0);
            let mut remove = None;
            for (i, c) in edit.constraints.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    if ui
                        .small_button("×")
                        .on_hover_text("この制約を削除")
                        .clicked()
                    {
                        remove = Some(i);
                    }
                    let kind = c.kind();
                    let atoms = c
                        .atoms()
                        .iter()
                        .map(|&a| label_of(a))
                        .collect::<Vec<_>>()
                        .join("–");
                    ui.label(RichText::new(atoms).size(11.0));

                    // Editable in place: the target that was right when the
                    // restraint was added is often not the one wanted after
                    // seeing what the minimization did with it.
                    let mut target = c.target();
                    if ui
                        .add(
                            egui::DragValue::new(&mut target)
                                .speed(if kind == Kind::Distance { 0.02 } else { 1.0 })
                                .range(kind.range())
                                .suffix(format!(" {}", kind.unit())),
                        )
                        .changed()
                    {
                        c.set_target(target);
                    }

                    if let Some(v) = readout.constraint_values.get(i) {
                        ui.label(
                            RichText::new(format!("現在 {v:.2}"))
                                .size(11.0)
                                .color(PAL.label),
                        );
                    }
                });
            }
            if let Some(i) = remove {
                edit.constraints.remove(i);
            }
            match crate::constraints::conflict(&edit.constraints) {
                Some(message) => {
                    ui.label(RichText::new(message).size(11.0).color(PAL.danger_text));
                }
                None => {
                    ui.label(
                        RichText::new("「▶ MM 最小化」を実行すると、これらを保ったまま最適化します。")
                            .size(11.0)
                            .color(PAL.label),
                    );
                }
            }
        }

        ui.add_space(10.0);
        ui.separator();
        ui.add_space(6.0);

        if ui
            .add_enabled(!edit.undo.is_empty(), egui::Button::new("↶ 3D 編集を取り消す"))
            .clicked()
        {
            action = Some(Action::Undo);
        }
    });

    action
}

/// One `label   value` line, with an em dash where there is nothing selected.
fn row(ui: &mut egui::Ui, label: &str, value: Option<String>) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(12.0).color(PAL.label));
        ui.label(RichText::new(value.unwrap_or_else(|| "—".to_string())).size(12.0));
    });
}
