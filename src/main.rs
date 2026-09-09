//! chem_3d_builder_rs — draw a 2D structure on the left, turn it into a 3D
//! structure, and run molecular-mechanics minimization on the right.
//!
//! Glues together three crates:
//! - chembuider_rs      : 2D structure editor (egui widget)
//! - moleucle_3dview_rs : 3D molecule viewer (egui + wgpu)
//! - openbabel          : hydrogen addition, 3D generation, MM force fields,
//!   minimization and every import/export format
//!
//! See `bridge.rs` for the 2D → 3D conversion, `mm_session.rs` for how a
//! minimization run is kept off the UI thread, and `edit3d.rs` / `dihedral.rs` /
//! `ik.rs` for editing the 3D structure in place.

mod bridge;
mod constraints;
mod dihedral;
mod edit3d;
mod forcefield_kind;
mod geom3d;
mod ik;
mod mm_session;
#[cfg(test)]
mod test_support;
mod theme;
mod viewport3d;

use chembuider_rs::{ChemStructEditor, Config};
use edit3d::{Action, Edit3d, Readout};
use eframe::egui;
use egui::RichText;
use forcefield_kind::FfKind;
use mm_session::{MmState, PollOutcome};
use moleucle_3dview_rs::RenderStyle;
use theme::PAL;
use viewport3d::Viewport3d;

/// Export targets offered in the save dialog. OpenBabel picks the writer from
/// the extension, so this list and its inference must agree — `writes_every_offered_format`
/// pins that.
const FORMATS: &[(&str, &[&str])] = &[
    ("Tripos Mol2", &["mol2"]),
    ("PDB", &["pdb"]),
    ("MDL SDF", &["sdf", "mol"]),
    ("XYZ", &["xyz"]),
    ("SMILES", &["smi"]),
];

/// Import sources. Same deal in reverse — OpenBabel picks the reader from the
/// extension, and `reads_every_offered_format` pins that every one of these
/// resolves to a reader that is actually built into the library.
///
/// CML and CIF are deliberately absent: OpenBabel builds them against libxml2,
/// which `openbabel_rs` does not pull in, so offering them would put filters in
/// the dialog that fail on every file.
const READ_FORMATS: &[(&str, &[&str])] = &[
    ("MDL Molfile / SDF", &["mol", "sdf"]),
    ("Tripos Mol2", &["mol2"]),
    ("PDB", &["pdb"]),
    ("XYZ", &["xyz"]),
    ("SMILES", &["smi"]),
    ("InChI", &["inchi"]),
];

/// Image targets.
///
/// Unlike the structure formats these are not interchangeable: each is a
/// different renderer, and two of them write a PNG. `rfd` reports only the path
/// the user chose, never which filter was active, so the extension cannot tell
/// them apart — the choice has to be made before the dialog opens, which is why
/// the toolbar's 画像 button is a menu.
#[derive(Clone, Copy, PartialEq, Eq)]
enum ImageKind {
    /// The 3D viewport exactly as it is on screen.
    View3d,
    /// The 2D structure drawing, rasterised.
    Structure2dPng,
    /// The structure as a vector drawing, from OpenBabel's depiction writer.
    Structure2dSvg,
}

const IMAGE_FORMATS: &[(ImageKind, &str, &str)] = &[
    (ImageKind::View3d, "3D ビュー (PNG)", "png"),
    (ImageKind::Structure2dPng, "2D 構造式 (PNG)", "png"),
    (ImageKind::Structure2dSvg, "2D 構造式 (SVG)", "svg"),
];

impl ImageKind {
    fn label(self) -> &'static str {
        IMAGE_FORMATS
            .iter()
            .find(|(kind, ..)| *kind == self)
            .map(|(_, label, _)| *label)
            .unwrap_or("画像")
    }

    fn extension(self) -> &'static str {
        IMAGE_FORMATS
            .iter()
            .find(|(kind, ..)| *kind == self)
            .map(|(.., ext)| *ext)
            .unwrap_or("png")
    }
}

struct App {
    editor: ChemStructEditor,
    viewport: Viewport3d,
    render_state: Option<egui_wgpu::RenderState>,
    ff_kind: FfKind,
    mm: MmState,
    steps_per_frame: usize,
    /// Selection and controls for editing the 3D structure directly.
    edit: Edit3d,
    /// Last conversion error, surfaced in the energy badge (the refined design
    /// has no status bar — drawing tools and errors go through the widget's
    /// keyboard shortcuts / the badge instead).
    error: Option<String>,
    /// What `error` held the last time it was written to the log. The badge is
    /// the only place an error is shown, so nothing about a failure escapes the
    /// running process unless it is mirrored out — see `log_error_changes`.
    logged_error: Option<String>,
    /// The same, for the 3D viewport's draw error.
    logged_draw_error: Option<String>,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        Self::setup_fonts(&cc.egui_ctx);
        theme::install(&cc.egui_ctx);

        if cc.wgpu_render_state.is_none() {
            // Fatal for the 3D half of the app: the central panel can only show
            // an apology from here on. Worth a line in the log, because the
            // window still opens and the failure otherwise looks like a
            // rendering glitch rather than a missing backend.
            log::error!("no wgpu render state — the 3D viewport is unavailable");
        }

        let mut editor = ChemStructEditor::default();
        // Use built-in defaults so the app never touches the user's config files.
        editor.config = Config::embedded();
        Self {
            editor,
            viewport: Viewport3d::new(),
            render_state: cc.wgpu_render_state.clone(),
            ff_kind: FfKind::default(),
            mm: MmState::Empty,
            steps_per_frame: 8,
            edit: Edit3d::new(),
            error: None,
            logged_error: None,
            logged_draw_error: None,
        }
    }

    /// Mirror `error` to the log whenever it changes.
    ///
    /// Every failure path in the app ends in `self.error = Some(..)`, and the
    /// energy badge is the only thing that reads it — so a run that went wrong
    /// leaves no trace outside the window. Watching the field from one place
    /// beats logging at each of the two dozen assignments: it cannot be
    /// forgotten by a later one, and comparing against the last logged value
    /// keeps a message that persists across frames to a single line.
    ///
    /// `draw` is the 3D viewport's own failure, which never reaches `error` —
    /// it is painted straight into the central panel — but fails every frame
    /// once it starts, so it needs the same edge detection.
    fn log_error_changes(&mut self, draw: Option<String>) {
        if draw != self.logged_draw_error {
            if let Some(message) = &draw {
                log::error!("3D draw failed: {message}");
            }
            self.logged_draw_error = draw;
        }

        if self.error != self.logged_error {
            if let Some(message) = &self.error {
                log::warn!("{message}");
            }
            self.logged_error = self.error.clone();
        }
    }

    /// Load a system Japanese font into egui (default fonts lack CJK glyphs).
    /// Mirrors the approach in the sibling `Shiratama-rs` repo.
    fn setup_fonts(ctx: &egui::Context) {
        let mut fonts = egui::FontDefinitions::default();
        let mut loaded = false;

        #[cfg(target_os = "windows")]
        for path in [
            "C:\\Windows\\Fonts\\YuGothM.ttc",
            "C:\\Windows\\Fonts\\meiryo.ttc",
            "C:\\Windows\\Fonts\\msgothic.ttc",
        ] {
            if let Ok(bytes) = std::fs::read(path) {
                fonts.font_data.insert(
                    "japanese".to_owned(),
                    std::sync::Arc::new(egui::FontData::from_owned(bytes)),
                );
                loaded = true;
                break;
            }
        }

        #[cfg(target_os = "macos")]
        if let Ok(bytes) = std::fs::read("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc")
        {
            fonts.font_data.insert(
                "japanese".to_owned(),
                std::sync::Arc::new(egui::FontData::from_owned(bytes)),
            );
            loaded = true;
        }

        #[cfg(target_os = "linux")]
        for path in [
            "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
            "/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc",
        ] {
            if let Ok(bytes) = std::fs::read(path) {
                fonts.font_data.insert(
                    "japanese".to_owned(),
                    std::sync::Arc::new(egui::FontData::from_owned(bytes)),
                );
                loaded = true;
                break;
            }
        }

        if loaded {
            // Append as a fallback so Latin keeps the default font and CJK glyphs
            // resolve through the Japanese font.
            for family in [egui::FontFamily::Proportional, egui::FontFamily::Monospace] {
                fonts
                    .families
                    .entry(family)
                    .or_default()
                    .push("japanese".to_owned());
            }
        }

        ctx.set_fonts(fonts);
    }

    /// Hand a freshly built or loaded molecule to the viewer and the MM state.
    /// Every 3D selection is dropped: the indices in it belong to the structure
    /// being replaced.
    fn adopt(&mut self, mol: openbabel::Molecule) {
        self.viewport.set_molecule(bridge::to_viewer_molecule(&mol));
        self.viewport.focus_on_molecule_center();
        self.edit.reset();
        self.error = None;
        self.mm = MmState::ready(mol, self.ff_kind);
    }

    /// Convert the current 2D structure into a 3D structure ready to minimize.
    fn generate_3d(&mut self) {
        match bridge::to_ob_molecule(&self.editor.molecule, self.ff_kind.ob_id()) {
            Ok(converted) => {
                self.adopt(converted.mol);
                // After `adopt`, which clears the error: a cramped geometry is
                // still worth showing, and the badge is where it goes.
                self.error = converted.warning;
            }
            Err(e) => self.error = Some(e),
        }
    }

    /// Read a structure file into the 3D side. OpenBabel picks the reader from
    /// the extension, the same way the export path picks the writer.
    ///
    /// The 2D canvas is deliberately left alone: `bridge.rs` only converts 2D →
    /// OpenBabel, so there is nothing to draw the imported structure with. The
    /// button's hover text says so, because pressing "→ 3D 生成" afterwards
    /// throws the import away.
    fn import_structure(&mut self) {
        let mut dialog = rfd::FileDialog::new().set_title("構造を読み込み");
        for (name, exts) in READ_FORMATS {
            dialog = dialog.add_filter(*name, exts);
        }
        let Some(path) = dialog.pick_file() else {
            return; // user cancelled
        };
        let Some(path) = path.to_str() else {
            self.error = Some("パスに使用できない文字が含まれています。".to_string());
            return;
        };

        let mut mol = match openbabel::Molecule::read_file(path, None) {
            Ok(mol) => mol,
            Err(e) => {
                self.error = Some(format!("読み込みに失敗しました: {e}"));
                return;
            }
        };
        if mol.num_atoms() == 0 {
            self.error = Some("ファイルに原子がありません。".to_string());
            return;
        }
        // SMILES and InChI carry no geometry, and a 2D file has none either.
        // Everything downstream — the viewer, MM, the editing gestures — needs
        // real 3D coordinates, so build them here rather than failing later.
        //
        // `has_3d()` alone is not enough: it reports the dimension the *format*
        // declares, and Mol2 and SDF declare three whatever the numbers say. A
        // 2D sketch saved as Mol2 therefore arrives claiming to be 3D with every
        // z at zero, which loads as a pile of overlapping atoms with an energy in
        // the hundreds of millions. Checking the coordinates catches that.
        let flat = geom3d::is_planar(&bridge::coordinates(&mol));
        if !mol.has_3d() || flat {
            mol.add_hydrogens();
            if !mol.generate_3d() {
                self.error = Some("立体構造の生成に失敗しました。".to_string());
                return;
            }
        }
        self.adopt(mol);
        if flat {
            // Said plainly, because it means the file did not carry the
            // conformation the user may have expected to get back.
            self.error = Some("平面構造だったため、立体構造を生成しました。".to_string());
        }
    }

    /// Write the current 3D structure out. OpenBabel picks the writer from the
    /// extension the user chose in the native save dialog. A no-op with an error
    /// badge if there is no 3D structure to write yet. Shared by the toolbar's
    /// Export button and Ctrl+S.
    fn export_structure(&mut self) {
        // `mol` is None while a worker owns the molecule, but the button and the
        // shortcut are both gated on that, so this only fires before 3D generation.
        let Some(mol) = self.mm.mol() else {
            self.error = Some("3D構造がありません。先に「→ 3D 生成」してください。".to_string());
            return;
        };

        let mut dialog = rfd::FileDialog::new()
            .set_title("構造をエクスポート")
            .set_file_name("molecule.mol2");
        for (name, exts) in FORMATS {
            dialog = dialog.add_filter(*name, exts);
        }
        let Some(path) = dialog.save_file() else {
            return; // user cancelled
        };

        // `write_file` takes a &str, so a path we cannot render is a real error
        // rather than something to unwrap through.
        let Some(path) = path.to_str() else {
            self.error = Some("パスに使用できない文字が含まれています。".to_string());
            return;
        };

        match mol.write_file(path, None) {
            Ok(()) => self.error = None,
            Err(e) => self.error = Some(format!("保存に失敗しました: {e}")),
        }
    }

    /// Write a picture of the kind the user picked from the 画像 menu.
    fn export_image(&mut self, kind: ImageKind) {
        let ext = kind.extension();
        let Some(path) = rfd::FileDialog::new()
            .set_title(format!("{} を保存", kind.label()))
            .set_file_name(format!("molecule.{ext}"))
            .add_filter(kind.label(), &[ext])
            .save_file()
        else {
            return; // user cancelled
        };

        if let Err(e) = self.write_image(&path, kind) {
            self.error = Some(format!("画像の保存に失敗しました: {e}"));
        } else {
            self.error = None;
        }
    }

    fn write_image(&self, path: &std::path::Path, kind: ImageKind) -> Result<(), String> {
        match kind {
            ImageKind::View3d => {
                let render_state = self
                    .render_state
                    .as_ref()
                    .ok_or("wgpu レンダラが利用できません。")?;
                let (width, height, rgba) = self.viewport.read_rgba(render_state)?;
                write_png(path, width, height, &rgba)
            }
            ImageKind::Structure2dPng => {
                let png = chembuider_rs::molecule::image::molecule_to_png(&self.editor.molecule)
                    .ok_or("2D構造式を描画できませんでした。")?;
                std::fs::write(path, png).map_err(|e| e.to_string())
            }
            ImageKind::Structure2dSvg => {
                let mol = self
                    .mm
                    .mol()
                    .ok_or("3D構造がありません。先に「→ 3D 生成」してください。")?;
                let svg = mol.to_svg().ok_or("SVG を生成できませんでした。")?;
                std::fs::write(path, svg).map_err(|e| e.to_string())
            }
        }
    }

    /// Collect a finished minimization run, or show the next recorded frame.
    fn advance_mm(&mut self, ctx: &egui::Context) {
        match self.mm.poll() {
            PollOutcome::Frame(coords) => {
                let _ = self.viewport.update_positions_angstrom(&coords);
                // The snapshots describe a geometry the minimizer has now moved
                // past, so "undo the last rotation" would jump backwards
                // through a run the user did not ask to undo.
                self.edit.clear_undo();
                ctx.request_repaint();
            }
            // Without this the UI sleeps and never polls the worker again.
            PollOutcome::Waiting => ctx.request_repaint(),
            PollOutcome::Failed(e) => self.error = Some(e),
            PollOutcome::Idle => {}
        }
    }

    /// Push `coords` into OpenBabel and the viewer, and re-read the energy.
    fn commit_geometry(&mut self, coords: &[[f64; 3]]) -> Result<(), String> {
        let mol = self
            .mm
            .mol_mut()
            .ok_or("最小化の実行中は構造を編集できません。")?;
        dihedral::commit(mol, coords)?;
        let as_f32: Vec<[f32; 3]> = coords
            .iter()
            .map(|p| [p[0] as f32, p[1] as f32, p[2] as f32])
            .collect();
        self.viewport.update_positions_angstrom(&as_f32)?;
        self.mm.refresh_badge(self.ff_kind);
        Ok(())
    }

    /// Turn the armed bond. `absolute` sets the torsion outright rather than
    /// adding to it, which is what the panel's angle box does.
    fn rotate_selected_bond(&mut self, degrees: f64, absolute: bool) {
        let (Some(anchor), Some(bond)) = (self.edit.selection.anchor, self.edit.selection.bond)
        else {
            return;
        };
        let Some(mol) = self.mm.mol() else {
            self.error = Some("最小化の実行中は構造を編集できません。".to_string());
            return;
        };

        let rotation = match dihedral::Rotation::resolve(mol, bond, anchor) {
            Ok(r) => r,
            Err(reject) => {
                // Name the bond: the panel shows a selection, but the log did
                // not, which made a refusal impossible to tell from a misclick
                // after the fact.
                let named = match mol.bond(bond) {
                    Some(b) => format!(
                        "{} – {}: {}",
                        atom_label(mol, b.begin_atom_index()),
                        atom_label(mol, b.end_atom_index()),
                        reject.message()
                    ),
                    None => reject.message().to_string(),
                };
                self.error = Some(with_rotatable_hint(&named, mol));
                return;
            }
        };
        let delta = if absolute {
            let Some(current) = rotation.angle_degrees(mol) else {
                self.error = Some("この結合には角度を測る置換基がありません。".to_string());
                return;
            };
            degrees - current
        } else {
            degrees
        };

        let coords = bridge::coordinates(mol);
        let next = rotation.rotated(&coords, delta);
        // Only remember the previous geometry once the new one is actually in
        // place, so a failed edit does not leave an undo step that changes
        // nothing.
        match self.commit_geometry(&next) {
            Ok(()) => {
                self.edit.push_undo(coords);
                self.error = None;
            }
            Err(e) => self.error = Some(e),
        }
    }

    /// Solve for the requested separation of the two selected atoms.
    fn solve_ik(&mut self, goal: ik::Goal) {
        let Some([a, b]) = self.edit.ik_pair() else {
            return;
        };
        let Some(mol) = self.mm.mol() else {
            self.error = Some("最小化の実行中は構造を編集できません。".to_string());
            return;
        };

        let coords = bridge::coordinates(mol);
        let solution = match ik::solve(mol, &coords, a, b, goal) {
            Ok(s) => s,
            Err(e) => {
                // Same reason as the bond refusal: without the pair, neither the
                // user nor the log can tell a correct refusal from a misclick.
                let named = format!(
                    "{} – {}: {}",
                    atom_label(mol, a),
                    atom_label(mol, b),
                    e.message()
                );
                self.error = Some(with_rotatable_hint(&named, mol));
                self.edit.report = None;
                return;
            }
        };

        let clashes = geom3d::clashes(mol, &solution.coords);
        let mut report = format!(
            "{:.3} Å → {:.3} Å（{} 回の反復）",
            solution.start, solution.achieved, solution.sweeps
        );
        if solution.unreachable {
            report.push_str(&format!(
                "\n⚠ 指定の距離には届きません。この経路で到達できるのは {:.2}〜{:.2} Å です。",
                solution.reach.0, solution.reach.1
            ));
        }
        if let Some(worst) = clashes.first() {
            report.push_str(&format!(
                "\n⚠ 原子の重なり {} 件（最悪: {} – {} が {:.2} Å）",
                clashes.len(),
                atom_label(mol, worst.a),
                atom_label(mol, worst.b),
                worst.distance,
            ));
        }

        match self.commit_geometry(&solution.coords) {
            Ok(()) => {
                self.edit.push_undo(coords);
                self.edit.report = Some(report);
                self.error = None;
            }
            Err(e) => self.error = Some(e),
        }
    }

    fn undo_3d(&mut self) {
        let Some(coords) = self.edit.pop_undo() else {
            return;
        };
        if let Err(e) = self.commit_geometry(&coords) {
            self.error = Some(e);
        } else {
            self.error = None;
        }
    }

    /// Resolve the selection against the molecule once a frame: refresh the
    /// viewer's highlights, hand the panel the names it needs, and read back the
    /// numbers it displays.
    ///
    /// Doing all three in one pass matters because each is an OpenBabel call
    /// behind a process-global lock, and resolving the rotation is the expensive
    /// one — the panel wants the torsion and the highlights want the moving
    /// side, which is the same resolution twice.
    fn sync_selection(&mut self) -> Readout {
        let Some(mol) = self.mm.mol() else {
            self.edit.set_resolved(None, Vec::new());
            self.publish_overlays();
            return Readout::default();
        };

        let bond_atoms = self
            .edit
            .selection
            .bond
            .and_then(|b| mol.bond(b))
            .map(|b| (b.begin_atom_index(), b.end_atom_index()));

        let (rotation, rotation_error) =
            match (self.edit.selection.anchor, self.edit.selection.bond) {
                (Some(anchor), Some(bond)) => {
                    match dihedral::Rotation::resolve(mol, bond, anchor) {
                        Ok(r) => (Some(r), None),
                        Err(reject) => (None, Some(reject.message())),
                    }
                }
                _ => (None, None),
            };
        let torsion_degrees = rotation.as_ref().and_then(|r| r.angle_degrees(mol));
        let moving = rotation
            .as_ref()
            .map(|r| r.moving_atoms().to_vec())
            .unwrap_or_default();

        let pair_distance = self.edit.ik_pair().map(|[a, b]| mol.distance(a, b));

        // The restraint readouts: what the pending pick would measure, and where
        // every restraint already added stands right now.
        let kind = self.edit.constraint_kind;
        let pending_constraint = self
            .edit
            .constraint_atoms(kind)
            .and_then(|atoms| constraints::Constraint::from_current(kind, atoms, mol))
            .map(|c| c.target());
        let constraint_values = self
            .edit
            .constraints
            .iter()
            .map(|c| c.measure(mol))
            .collect();

        self.edit.set_resolved(bond_atoms, moving);
        self.publish_overlays();

        Readout {
            torsion_degrees,
            rotation_error,
            pair_distance,
            pending_constraint,
            constraint_values,
        }
    }

    fn publish_overlays(&mut self) {
        let (groups, pairs) = self.edit.overlays();
        self.viewport.set_state_by_type(groups);
        self.viewport.set_state_by_type(pairs);
    }

    /// The refined design's single compact toolbar:
    /// `編集 │ 力場 · 変換  ───  E badge │ 表示`.
    /// Drawing tools (element / bond-order / tool / delete) live in the 2D
    /// widget's own keyboard shortcuts, so they are intentionally absent here.
    fn toolbar(&mut self, ui: &mut egui::Ui) {
        // While a run is in flight a worker holds the molecule and OpenBabel's
        // global lock, so everything that would reach OpenBabel waits it out.
        let running = self.mm.is_running();
        // None in both Empty and Running, which is exactly what the MM and
        // export buttons want to be disabled on.
        let has_mol = self.mm.mol().is_some();

        ui.horizontal(|ui| {
            // ── 編集 ──
            if ui.button("↺ 元に戻す").clicked() {
                self.editor.undo();
            }
            let cleaning = self.editor.is_cleaning();
            if ui
                .selectable_label(cleaning, "整形")
                .on_hover_text("2Dレイアウトの自動整形（再クリックで停止）")
                .clicked()
            {
                self.editor.toggle_cleanup();
            }
            if ui
                .add_enabled_ui(!running, |ui| theme::danger_button(ui, "クリア"))
                .inner
                .clicked()
            {
                self.editor.molecule = chembuider_rs::Molecule::default();
                self.editor.selected_atoms.clear();
                self.mm = MmState::Empty;
                self.edit.reset();
                self.error = None;
            }

            theme::divider(ui);

            // ── 力場 ──
            theme::group_label(ui, "力場");
            ui.add_enabled_ui(!running, |ui| {
                egui::ComboBox::from_id_salt("ff_kind")
                    .selected_text(self.ff_kind.label())
                    .show_ui(ui, |ui| {
                        for k in FfKind::ALL {
                            ui.selectable_value(&mut self.ff_kind, k, k.label());
                        }
                    });
            });

            // ── 変換 ──
            if ui
                .add_enabled_ui(!running, |ui| theme::accent_button(ui, "→ 3D 生成"))
                .inner
                .on_hover_text("水素を付加し初期立体構造を生成（読み込んだ構造は破棄されます）")
                .clicked()
            {
                self.generate_3d();
            }

            // A run now continues until the geometry converges, so it can last
            // longer than the old fixed budget did — the button stays live and
            // stops it at the next segment boundary.
            let mm_button = if running {
                egui::Button::new(RichText::new("■ 停止").size(13.0).color(PAL.danger_text))
                    .fill(PAL.danger_bg)
                    .stroke(egui::Stroke::new(1.0, PAL.danger_border))
                    .corner_radius(egui::CornerRadius::same(9))
                    .min_size(theme::h(theme::CTRL_H))
            } else {
                egui::Button::new((
                    RichText::new("▶").size(10.0).color(PAL.mm_play),
                    RichText::new(" MM 最小化").size(13.0).color(PAL.text),
                ))
                .corner_radius(egui::CornerRadius::same(9))
                .min_size(theme::h(theme::CTRL_H))
            };
            // `has_mol` is false while running (the worker owns the molecule), so
            // the stop button has to be enabled on its own terms.
            if ui.add_enabled(has_mol || running, mm_button).clicked() {
                if running {
                    self.mm.cancel();
                } else {
                    let ctx = ui.ctx().clone();
                    let constraints = self.edit.constraints.clone();
                    // Refused rather than run: OpenBabel would drop the dihedral
                    // and hand back a structure that ignored what was asked for,
                    // with nothing on screen saying so.
                    match constraints::conflict(&constraints) {
                        Some(message) => self.error = Some(message.to_string()),
                        None => {
                            self.error = None;
                            self.mm
                                .start(self.ff_kind, self.steps_per_frame, constraints, &ctx);
                        }
                    }
                }
            }

            theme::divider(ui);

            // ── 入出力 ──
            let import_button = egui::Button::new(RichText::new("⭱ 読み込み").size(13.0))
                .corner_radius(egui::CornerRadius::same(9))
                .min_size(theme::h(theme::CTRL_H));
            if ui
                .add_enabled(!running, import_button)
                .on_hover_text(
                    "構造ファイルを 3D 側に読み込む (Ctrl+O)\n左の 2D キャンバスには反映されません",
                )
                .clicked()
            {
                self.import_structure();
            }

            let export_button = egui::Button::new(RichText::new("⭳ エクスポート").size(13.0))
                .corner_radius(egui::CornerRadius::same(9))
                .min_size(theme::h(theme::CTRL_H));
            if ui
                .add_enabled(has_mol, export_button)
                .on_hover_text("Mol2 / PDB / SDF / XYZ / SMILES で保存 (Ctrl+S)")
                .clicked()
            {
                self.export_structure();
            }

            // A menu rather than a button: two of the three targets are PNGs, and
            // the save dialog cannot report which filter was chosen, so the kind
            // has to be settled before it opens.
            let mut image_kind = None;
            ui.add_enabled_ui(has_mol, |ui| {
                // No dropdown arrow in the label: the CJK fallback fonts this app
                // loads have no glyph for one, so it renders as a tofu box.
                ui.menu_button(RichText::new("🖼 画像").size(13.0), |ui| {
                    for (kind, label, _) in IMAGE_FORMATS {
                        if ui.button(*label).clicked() {
                            image_kind = Some(*kind);
                            ui.close();
                        }
                    }
                })
                .response
                .on_hover_text("3D ビュー / 2D 構造式を画像で保存");
            });
            if let Some(kind) = image_kind {
                self.export_image(kind);
            }

            // ── right cluster (margin-auto): エネルギー … 表示 ──
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let mut style = self.viewport.render_style();
                egui::ComboBox::from_id_salt("view_style")
                    .selected_text(style_label(style))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut style, RenderStyle::BallStick, "Ball & Stick");
                        ui.selectable_value(&mut style, RenderStyle::BallOnly, "Ball Only");
                        ui.selectable_value(&mut style, RenderStyle::Circles, "Circles");
                        ui.selectable_value(&mut style, RenderStyle::Wireframe, "Wireframe");
                    });
                self.viewport.set_render_style(style);
                theme::group_label(ui, "表示");

                theme::divider(ui);

                theme::energy_badge(ui, self.mm.badge(), self.error.as_deref());
            });
        });
    }
}

/// `C3`-style label for an atom, for messages and the selection panel.
/// How many bonds in `mol` can be turned at all.
///
/// A refusal that only says "this bond cannot be turned" leaves the user
/// clicking bonds one by one to find out whether *any* of them can — on a fused
/// aromatic system the answer is often one, or none, and nothing on screen
/// distinguishes them. Saying how many exist turns a hunt into a decision.
fn rotatable_bond_count(mol: &openbabel::Molecule) -> usize {
    mol.bonds().filter(|b| b.is_rotor()).count()
}

/// Append the molecule-wide count to a refusal, so it says what *would* work.
fn with_rotatable_hint(message: &str, mol: &openbabel::Molecule) -> String {
    match rotatable_bond_count(mol) {
        0 => format!(
            "{message}
この分子に回転できる結合はありません。"
        ),
        n => format!(
            "{message}
この分子で回転できる結合は {n} 本です。"
        ),
    }
}

fn atom_label(mol: &openbabel::Molecule, index: u32) -> String {
    match mol.atom(index) {
        Some(a) => format!(
            "{}{}",
            openbabel::elements::symbol(a.atomic_number()),
            index + 1
        ),
        None => format!("#{}", index + 1),
    }
}

fn write_png(path: &std::path::Path, width: u32, height: u32, rgba: &[u8]) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().map_err(|e| e.to_string())?;
    writer.write_image_data(rgba).map_err(|e| e.to_string())
}

/// Human-friendly label for the 3D render-style combo.
fn style_label(style: RenderStyle) -> &'static str {
    match style {
        RenderStyle::BallStick => "Ball & Stick",
        RenderStyle::BallOnly => "Ball Only",
        RenderStyle::Circles => "Circles",
        RenderStyle::Wireframe => "Wireframe",
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.advance_mm(&ctx);

        // Ctrl+S / Ctrl+O (Cmd on macOS) export and import. Consumed so they do
        // not also reach the 2D editor widget — including while a minimization
        // holds the molecule, where they are inert like the buttons beside them.
        let running = self.mm.is_running();
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::S)) && !running {
            self.export_structure();
        }
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::O)) && !running {
            self.import_structure();
        }

        let toolbar_frame = egui::Frame::new()
            .fill(PAL.toolbar_bg)
            .inner_margin(egui::Margin::symmetric(16, 6));
        egui::Panel::top("toolbar")
            .frame(toolbar_frame)
            .show(ui, |ui| {
                self.toolbar(ui);
            });

        // Bare workspace: white 2D canvas on the left, dark 3D viewport on the
        // right (both widgets paint their own background). No panel chrome —
        // the refined design shows the canvases edge-to-edge.
        let bare = egui::Frame::new().fill(PAL.panel);
        egui::Panel::left("editor_2d")
            .resizable(true)
            .default_size(720.0)
            .frame(bare)
            .show(ui, |ui| {
                let _ = self.editor.ui(ui);
            });

        let readout = self.sync_selection();

        // The 3D editing controls only mean anything once there is a structure,
        // and the panel has to be laid out before the central viewport claims
        // the rest of the width.
        let mut action = None;
        let mut draw_error = None;
        if self.viewport.has_molecule() {
            let panel_frame = egui::Frame::new()
                .fill(PAL.toolbar_bg)
                .inner_margin(egui::Margin::symmetric(12, 8));
            let labels: Vec<String> = self
                .mm
                .mol()
                .map(|mol| (0..mol.num_atoms()).map(|i| atom_label(mol, i)).collect())
                .unwrap_or_default();
            let label_of = move |i: u32| {
                labels
                    .get(i as usize)
                    .cloned()
                    .unwrap_or_else(|| format!("#{}", i + 1))
            };
            let edit = &mut self.edit;
            egui::Panel::right("edit_3d")
                .resizable(true)
                .default_size(260.0)
                .frame(panel_frame)
                .show(ui, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        action = edit3d::panel(ui, edit, readout, &label_of, !running);
                    });
                });
        }

        egui::CentralPanel::default()
            .frame(bare)
            .show(ui, |ui| match &self.render_state {
                Some(rs) => {
                    // The wheel belongs to the camera unless a rotation is armed
                    // — that is the whole difference between zooming and turning
                    // a bond, and it has to be decided before the viewport reads
                    // the scroll.
                    let capture = self.edit.selection.is_armed() && !running;
                    match self.viewport.show(ui, rs, capture) {
                        Ok(response) => {
                            if let Some(hit) = response.primary {
                                self.edit.on_primary(hit);
                            }
                            if let Some(hit) = response.secondary {
                                self.edit.on_secondary(hit);
                            }
                            if action.is_none() {
                                let modifiers = ui.input(|i| i.modifiers);
                                if let Some(deg) =
                                    self.edit.wheel_degrees(response.scroll_lines, modifiers)
                                {
                                    action = Some(Action::Rotate(deg));
                                }
                            }
                        }
                        Err(e) => {
                            ui.colored_label(egui::Color32::RED, format!("3D描画エラー: {e}"));
                            draw_error = Some(e.to_string());
                        }
                    }
                }
                None => {
                    ui.heading("wgpu バックエンドが利用できません");
                    ui.label("eframe を wgpu レンダラで起動してください。");
                }
            });

        match action {
            Some(Action::Rotate(deg)) => self.rotate_selected_bond(deg, false),
            Some(Action::SetAngle(deg)) => self.rotate_selected_bond(deg, true),
            Some(Action::Solve(goal)) => self.solve_ik(goal),
            Some(Action::Undo) => self.undo_3d(),
            None => {}
        }

        // Last, so an error raised by the actions just above is carried out on
        // this frame rather than the next one.
        self.log_error_changes(draw_error);
    }

    fn on_exit(&mut self) {
        if let Some(rs) = &self.render_state {
            self.viewport.free_egui_texture(rs);
        }
    }
}

fn main() -> eframe::Result {
    // eframe, egui and wgpu all report through `log`, so without a logger every
    // adapter-selection failure and surface warning is silently discarded. Warn
    // is the useful default here — info floods with per-frame wgpu chatter — and
    // RUST_LOG still overrides it.
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    // Not strictly required — every OpenBabel entry point does this itself — but
    // it pins where the library's data directory is resolved, which is the part
    // that has to survive being packaged into an installer.
    openbabel::init();

    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1500.0, 900.0])
            .with_title("Chem 3D Builder — 2D → 3D → MM"),
        wgpu_options: Default::default(),
        ..Default::default()
    };

    eframe::run_native(
        "chem_3d_builder_rs",
        options,
        Box::new(|cc| Ok(Box::new(App::new(cc)))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(name);
        std::fs::create_dir_all(&dir).expect("temp dir");
        dir
    }

    /// `export_structure` has no format logic of its own — it hands the path to
    /// OpenBabel and lets the extension pick the writer. So what needs pinning is
    /// that every extension the dialog offers is one OpenBabel actually infers.
    #[test]
    fn writes_every_offered_format() {
        let _ob = crate::test_support::ob_guard();
        let mut mol = openbabel::Molecule::parse("CCO", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");

        let dir = temp_dir("chem_3d_builder_rs_export_test");

        for (name, exts) in FORMATS {
            for ext in *exts {
                let path = dir.join(format!("molecule.{ext}"));
                let p = path.to_str().expect("utf-8 temp path");
                mol.write_file(p, None)
                    .unwrap_or_else(|e| panic!("{name} (.{ext}): {e}"));
                let written = std::fs::read_to_string(&path).expect("read back");
                assert!(!written.trim().is_empty(), "{name} (.{ext}) wrote nothing");
            }
        }

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The mirror of the export test, and the reason it is worth having twice:
    /// readers and writers are separate OpenBabel plugins, and this build ships
    /// a different set of each. A format whose reader is missing would leave the
    /// import dialog offering a filter that fails on every file the user picks.
    ///
    /// Every format is tried before reporting, rather than panicking on the
    /// first: when the plugin set changes, one run should name everything that
    /// broke instead of one thing at a time.
    #[test]
    fn reads_every_offered_format() {
        let _ob = crate::test_support::ob_guard();
        let mut source = openbabel::Molecule::parse("CCO", "smi").expect("parse");
        assert!(source.generate_3d(), "gen3d");
        let heavy = source.num_heavy_atoms();

        let dir = temp_dir("chem_3d_builder_rs_import_test");
        let mut failures = Vec::new();

        for (name, exts) in READ_FORMATS {
            for ext in *exts {
                let path = dir.join(format!("molecule.{ext}"));
                let p = path.to_str().expect("utf-8 temp path");
                if let Err(e) = source.write_file(p, None) {
                    failures.push(format!("{name} (.{ext}) write: {e}"));
                    continue;
                }
                match openbabel::Molecule::read_file(p, None) {
                    // Hydrogens survive some formats and not others (SMILES and
                    // InChI drop them by design), so the heavy-atom count is the
                    // part that has to round-trip.
                    Ok(read) if read.num_heavy_atoms() == heavy => {}
                    Ok(read) => failures.push(format!(
                        "{name} (.{ext}): {} heavy atoms back, expected {heavy}",
                        read.num_heavy_atoms()
                    )),
                    Err(e) => failures.push(format!("{name} (.{ext}) read: {e}")),
                }
            }
        }

        let _ = std::fs::remove_dir_all(&dir);
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
