//! chem_3d_builder_rs — draw a 2D structure on the left, turn it into a 3D
//! structure, and run molecular-mechanics minimization on the right.
//!
//! Glues together three crates:
//!   - chembuider_rs        : 2D structure editor (egui widget)
//!   - moleucle_3dview_rs   : 3D molecule viewer (egui + wgpu)
//!   - zunda_rs             : MM force fields + minimizers
//!
//! See `bridge.rs` for the 2D → 3D conversion (hydrogen addition + initial 3D).

mod bridge;
mod export;
mod forcefield_kind;
mod mm_session;
mod theme;

use chembuider_rs::{ChemStructEditor, Config};
use eframe::egui;
use egui::RichText;
use forcefield_kind::FfKind;
use mm_session::MmSession;
use moleucle_3dview_rs::{InteractiveMoleculeViewport, RenderStyle};
use theme::PAL;

struct App {
    editor: ChemStructEditor,
    viewport: InteractiveMoleculeViewport,
    render_state: Option<egui_wgpu::RenderState>,
    ff_kind: FfKind,
    mm: Option<MmSession>,
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
            mm: None,
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

    /// Convert the current 2D structure into an initial 3D structure and set up MM.
    fn generate_3d(&mut self) {
        let mol3d = bridge::to_mol_3d(&self.editor.molecule);
        match MmSession::new(mol3d, self.ff_kind) {
            Ok(session) => {
                self.viewport.set_molecule(session.viewer_molecule());
                self.viewport.focus_on_molecule_center();
                self.error = None;
                self.mm = Some(session);
            }
            Err(e) => self.error = Some(e),
        }
    }

    /// Write the current 3D structure to a Mol2 or PDB file. The format is
    /// chosen from the extension the user picks in the native save dialog
    /// (defaults to Mol2). A no-op with an error badge if there is no 3D
    /// structure yet. Shared by the toolbar's Export button and Ctrl+S.
    fn export_structure(&mut self) {
        // Take owned copies so the immutable borrow of `self.mm` ends before we
        // touch `self.error`.
        let (mol, positions) = match self.mm.as_ref() {
            Some(mm) => (mm.mol().clone(), mm.positions_f64()),
            None => {
                self.error =
                    Some("3D構造がありません。先に「→ 3D 生成」してください。".to_string());
                return;
            }
        };

        let path = match rfd::FileDialog::new()
            .set_title("構造をエクスポート")
            .add_filter("Tripos Mol2", &["mol2"])
            .add_filter("PDB", &["pdb"])
            .set_file_name("molecule.mol2")
            .save_file()
        {
            Some(p) => p,
            None => return, // user cancelled
        };

        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let contents = match ext.as_str() {
            "pdb" => export::to_pdb(&mol, &positions),
            // Default to Mol2 for `.mol2` and anything unrecognized.
            _ => export::to_mol2(&mol, &positions),
        };

        match std::fs::write(&path, contents) {
            Ok(()) => self.error = None,
            Err(e) => self.error = Some(format!("保存に失敗しました: {e}")),
        }
    }

    /// Advance the live minimization by one frame, if running.
    fn advance_mm(&mut self, ctx: &egui::Context) {
        let spf = self.steps_per_frame;
        let mut new_coords: Option<Vec<[f32; 3]>> = None;
        let mut repaint = false;

        if let Some(mm) = self.mm.as_mut() {
            if mm.minimizing {
                mm.step(spf);
                new_coords = Some(mm.positions_angstrom());
                if mm.minimizing {
                    repaint = true;
                }
            }
        }

        if let Some(coords) = new_coords {
            let _ = self.viewport.update_positions_angstrom(&coords);
        }
        if repaint {
            ctx.request_repaint();
        }
    }

    /// The refined design's single compact toolbar:
    /// `編集 │ 力場 · 変換  ───  E badge │ 表示`.
    /// Drawing tools (element / bond-order / tool / delete) live in the 2D
    /// widget's own keyboard shortcuts, so they are intentionally absent here.
    fn toolbar(&mut self, ui: &mut egui::Ui) {
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
            if theme::danger_button(ui, "クリア").clicked() {
                self.editor.molecule = chembuider_rs::Molecule::default();
                self.editor.selected_atoms.clear();
                self.mm = None;
                self.error = None;
            }

            theme::divider(ui);

            // ── 力場 ──
            theme::group_label(ui, "力場");
            egui::ComboBox::from_id_salt("ff_kind")
                .selected_text(self.ff_kind.label())
                .show_ui(ui, |ui| {
                    for k in FfKind::ALL {
                        ui.selectable_value(&mut self.ff_kind, k, k.label());
                    }
                });

            // ── 変換 ──
            if theme::accent_button(ui, "→ 3D 生成")
                .on_hover_text("水素を付加し初期立体構造を生成")
                .clicked()
            {
                self.generate_3d();
            }

            let has_mm = self.mm.is_some();
            let minimizing = self.mm.as_ref().is_some_and(|m| m.minimizing);
            let mm_button = if minimizing {
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
            if ui.add_enabled(has_mm, mm_button).clicked() {
                if let Some(mm) = self.mm.as_mut() {
                    mm.minimizing = !mm.minimizing;
                }
            }

            // ── 出力 ──
            let export_button = egui::Button::new(RichText::new("⭳ エクスポート").size(13.0))
                .corner_radius(egui::CornerRadius::same(9))
                .min_size(theme::h(theme::CTRL_H));
            if ui
                .add_enabled(has_mm, export_button)
                .on_hover_text("Mol2 / PDB で保存 (Ctrl+S)")
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

                let energy = self.mm.as_ref().map(|m| (m.energy, m.total_steps));
                theme::energy_badge(ui, energy, self.error.as_deref());
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
        // so it doesn't also reach the 2D editor widget.
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::S)) {
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
