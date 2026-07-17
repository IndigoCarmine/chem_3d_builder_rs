//! chem_3d_builder_rs — draw a 2D structure on the left, turn it into a 3D
//! structure, and run molecular-mechanics minimization on the right.
//!
//! Glues together three crates:
//! - chembuider_rs      : 2D structure editor (egui widget)
//! - moleucle_3dview_rs : 3D molecule viewer (egui + wgpu)
//! - openbabel          : hydrogen addition, 3D generation, MM force fields,
//!   minimization and every export format
//!
//! See `bridge.rs` for the 2D → 3D conversion and `mm_session.rs` for how a
//! minimization run is kept off the UI thread.

mod bridge;
mod forcefield_kind;
mod mm_session;
mod theme;

use chembuider_rs::{ChemStructEditor, Config};
use eframe::egui;
use egui::RichText;
use forcefield_kind::FfKind;
use mm_session::{MmState, PollOutcome};
use moleucle_3dview_rs::{InteractiveMoleculeViewport, RenderStyle};
use theme::PAL;

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

struct App {
    editor: ChemStructEditor,
    viewport: InteractiveMoleculeViewport,
    render_state: Option<egui_wgpu::RenderState>,
    ff_kind: FfKind,
    mm: MmState,
    steps_per_frame: usize,
    /// Last conversion error, surfaced in the energy badge (the refined design
    /// has no status bar — drawing tools and errors go through the widget's
    /// keyboard shortcuts / the badge instead).
    error: Option<String>,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        Self::setup_fonts(&cc.egui_ctx);
        theme::install(&cc.egui_ctx);

        let mut editor = ChemStructEditor::default();
        // Use built-in defaults so the app never touches the user's config files.
        editor.config = Config::embedded();
        Self {
            editor,
            viewport: InteractiveMoleculeViewport::new(None),
            render_state: cc.wgpu_render_state.clone(),
            ff_kind: FfKind::default(),
            mm: MmState::Empty,
            steps_per_frame: 8,
            error: None,
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
        if let Ok(bytes) = std::fs::read("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc") {
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

    /// Convert the current 2D structure into a 3D structure ready to minimize.
    fn generate_3d(&mut self) {
        match bridge::to_ob_molecule(&self.editor.molecule) {
            Ok(mol) => {
                self.viewport.set_molecule(bridge::to_viewer_molecule(&mol));
                self.viewport.focus_on_molecule_center();
                self.error = None;
                self.mm = MmState::ready(mol, self.ff_kind);
            }
            Err(e) => self.error = Some(e),
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

    /// Collect a finished minimization run, or show the next recorded frame.
    fn advance_mm(&mut self, ctx: &egui::Context) {
        match self.mm.poll() {
            PollOutcome::Frame(coords) => {
                let _ = self.viewport.update_positions_angstrom(&coords);
                ctx.request_repaint();
            }
            // Without this the UI sleeps and never polls the worker again.
            PollOutcome::Waiting => ctx.request_repaint(),
            PollOutcome::Failed(e) => self.error = Some(e),
            PollOutcome::Idle => {}
        }
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
        let playing = self.mm.is_playing();

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
                .on_hover_text("水素を付加し初期立体構造を生成")
                .clicked()
            {
                self.generate_3d();
            }

            let mm_button = if playing {
                egui::Button::new(RichText::new("■ 停止").size(13.0).color(PAL.danger_text))
                    .fill(PAL.danger_bg)
                    .stroke(egui::Stroke::new(1.0, PAL.danger_border))
                    .corner_radius(egui::CornerRadius::same(9))
                    .min_size(theme::h(theme::CTRL_H))
            } else if running {
                // OpenBabel's minimizer has no cancel path, so this goes inert
                // rather than offering a stop that would not stop anything.
                egui::Button::new(RichText::new("⏳ 計算中…").size(13.0).color(PAL.text))
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
            if ui.add_enabled(has_mol, mm_button).clicked() {
                if playing {
                    self.mm.stop();
                } else {
                    let ctx = ui.ctx().clone();
                    self.mm.start(self.ff_kind, self.steps_per_frame, &ctx);
                }
            }

            // ── 出力 ──
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

        // Ctrl+S (Cmd+S on macOS) exports the current structure. Consume the key
        // so it doesn't also reach the 2D editor widget — including while a
        // minimization holds the molecule, where the shortcut is inert like the
        // button beside it.
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::S))
            && !self.mm.is_running()
        {
            self.export_structure();
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

        egui::CentralPanel::default().frame(bare).show(ui, |ui| {
            match &self.render_state {
                Some(rs) => {
                    if let Err(e) = self.viewport.show(ui, rs) {
                        ui.colored_label(egui::Color32::RED, format!("3D描画エラー: {e}"));
                    }
                }
                None => {
                    ui.heading("wgpu バックエンドが利用できません");
                    ui.label("eframe を wgpu レンダラで起動してください。");
                }
            }
        });
    }

    fn on_exit(&mut self) {
        if let Some(rs) = &self.render_state {
            self.viewport.free_egui_texture(rs);
        }
    }
}

fn main() -> eframe::Result {
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

    /// `export_structure` has no format logic of its own — it hands the path to
    /// OpenBabel and lets the extension pick the writer. So what needs pinning is
    /// that every extension the dialog offers is one OpenBabel actually infers.
    #[test]
    fn writes_every_offered_format() {
        let mut mol = openbabel::Molecule::parse("CCO", "smi").expect("parse");
        assert!(mol.generate_3d(), "gen3d");

        let dir = std::env::temp_dir().join("chem_3d_builder_rs_export_test");
        std::fs::create_dir_all(&dir).expect("temp dir");

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
}
