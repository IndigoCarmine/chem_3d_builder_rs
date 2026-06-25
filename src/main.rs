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
mod forcefield_kind;
mod mm_session;

use chembuider_rs::{BondOrder, ChemStructEditor, Config, Tool};
use eframe::egui;
use forcefield_kind::FfKind;
use mm_session::MmSession;
use moleucle_3dview_rs::{InteractiveMoleculeViewport, RenderStyle};

/// Elements offered as quick-pick buttons in the 2D toolbar.
const ELEMENTS: [&str; 9] = ["C", "N", "O", "H", "F", "Cl", "Br", "S", "P"];

struct App {
    editor: ChemStructEditor,
    viewport: InteractiveMoleculeViewport,
    render_state: Option<egui_wgpu::RenderState>,
    ff_kind: FfKind,
    mm: Option<MmSession>,
    steps_per_frame: usize,
    status: String,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
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
            status: "左で構造を描き、「→ 3D 生成」で立体化、「MM最小化」で構造緩和します。"
                .to_string(),
        }
    }

    /// Convert the current 2D structure into an initial 3D structure and set up MM.
    fn generate_3d(&mut self) {
        let mol3d = bridge::to_mol_3d(&self.editor.molecule);
        match MmSession::new(mol3d, self.ff_kind) {
            Ok(session) => {
                self.viewport.set_molecule(session.viewer_molecule());
                self.viewport.focus_on_molecule_center();
                self.status = format!(
                    "3D構造を生成しました（{} 原子, 水素付加済み, {}）。E = {:.3} kJ/mol",
                    session.atom_count(),
                    self.ff_kind.label(),
                    session.energy,
                );
                self.mm = Some(session);
            }
            Err(e) => self.status = e,
        }
    }

    /// Advance the live minimization by one frame, if running.
    fn advance_mm(&mut self, ctx: &egui::Context) {
        let spf = self.steps_per_frame;
        let mut new_mol = None;
        let mut finished_status = None;
        let mut repaint = false;

        if let Some(mm) = self.mm.as_mut() {
            if mm.minimizing {
                mm.step(spf);
                new_mol = Some(mm.viewer_molecule());
                if mm.minimizing {
                    repaint = true;
                } else {
                    finished_status = Some(format!(
                        "MM最小化終了: {} ステップ, E = {:.3} kJ/mol {}",
                        mm.total_steps,
                        mm.energy,
                        if mm.converged { "（収束）" } else { "（停止）" },
                    ));
                }
            }
        }

        if let Some(mol) = new_mol {
            self.viewport.set_molecule(mol);
        }
        if let Some(s) = finished_status {
            self.status = s;
        }
        if repaint {
            ctx.request_repaint();
        }
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            // ── 2D editing tools ──
            ui.label("ツール:");
            if ui
                .selectable_label(self.editor.tool == Tool::Bond, "✏ 結合")
                .clicked()
            {
                self.editor.tool = Tool::Bond;
            }
            if ui
                .selectable_label(self.editor.tool == Tool::Select, "⬚ 選択")
                .clicked()
            {
                self.editor.tool = Tool::Select;
            }
            if ui
                .selectable_label(self.editor.tool == Tool::Eraser, "✖ 消去")
                .clicked()
            {
                self.editor.tool = Tool::Eraser;
            }

            ui.separator();
            ui.label("元素:");
            for el in ELEMENTS {
                if ui
                    .selectable_label(self.editor.current_element == el, el)
                    .clicked()
                {
                    self.editor.current_element = el.to_string();
                }
            }

            ui.separator();
            ui.label("次数:");
            for (label, order) in [
                ("単", BondOrder::Single),
                ("二重", BondOrder::Double),
                ("三重", BondOrder::Triple),
            ] {
                if ui
                    .selectable_label(self.editor.current_bond_order == order, label)
                    .clicked()
                {
                    self.editor.current_bond_order = order;
                }
            }

            ui.separator();
            if ui.button("↶ 元に戻す").clicked() {
                self.editor.undo();
            }
            let cleaning = self.editor.is_cleaning();
            if ui
                .selectable_label(cleaning, "✨ 整形")
                .on_hover_text("2Dレイアウトの自動整形（再クリックで停止）")
                .clicked()
            {
                self.editor.toggle_cleanup();
            }
            if ui.button("🗑 クリア").clicked() {
                self.editor.molecule = chembuider_rs::Molecule::default();
                self.editor.selected_atoms.clear();
                self.mm = None;
                self.status = "キャンバスをクリアしました。".to_string();
            }
        });

        ui.separator();

        ui.horizontal_wrapped(|ui| {
            // ── 3D / MM controls ──
            egui::ComboBox::from_label("力場")
                .selected_text(self.ff_kind.label())
                .show_ui(ui, |ui| {
                    for k in FfKind::ALL {
                        ui.selectable_value(&mut self.ff_kind, k, k.label());
                    }
                });

            if ui
                .button("→ 3D 生成")
                .on_hover_text("水素を付加し初期立体構造を生成")
                .clicked()
            {
                self.generate_3d();
            }

            let has_mm = self.mm.is_some();
            let minimizing = self.mm.as_ref().is_some_and(|m| m.minimizing);
            let mm_label = if minimizing { "⛔ 停止" } else { "▶ MM最小化" };
            if ui
                .add_enabled(has_mm, egui::Button::new(mm_label))
                .clicked()
            {
                if let Some(mm) = self.mm.as_mut() {
                    mm.minimizing = !mm.minimizing;
                }
            }

            if let Some(mm) = self.mm.as_ref() {
                ui.separator();
                ui.label(format!(
                    "E = {:.3} kJ/mol | {} steps{}",
                    mm.energy,
                    mm.total_steps,
                    if mm.converged { " | 収束" } else { "" },
                ));
            }

            ui.separator();
            let mut style = self.viewport.render_style();
            egui::ComboBox::from_label("表示")
                .selected_text(format!("{style:?}"))
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut style, RenderStyle::BallStick, "BallStick");
                    ui.selectable_value(&mut style, RenderStyle::BallOnly, "BallOnly");
                    ui.selectable_value(&mut style, RenderStyle::Circles, "Circles");
                    ui.selectable_value(&mut style, RenderStyle::Wireframe, "Wireframe");
                });
            self.viewport.set_render_style(style);
        });
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.advance_mm(&ctx);

        egui::Panel::top("toolbar").show_inside(ui, |ui| {
            self.toolbar(ui);
        });

        egui::Panel::bottom("status").show_inside(ui, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!(
                    "原子: {}  結合: {}  |  {}",
                    self.editor.molecule.atoms.len(),
                    self.editor.molecule.bonds.len(),
                    self.status,
                ));
            });
        });

        egui::Panel::left("editor_2d")
            .resizable(true)
            .default_size(560.0)
            .show_inside(ui, |ui| {
                ui.heading("2D 構造エディタ");
                let _ = self.editor.ui(ui);
            });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            ui.heading("3D ビュー + MM");
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
