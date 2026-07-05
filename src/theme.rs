//! Visual theme for the app — the light "studio ribbon" look from the refined
//! `Chem 3D Builder.dc.html` design: warm neutrals, one compact toolbar, a blue
//! primary action and a green energy read-out.
//!
//! The design specifies every colour in OKLCH, so we convert OKLCH → sRGB at
//! startup (`oklch`) instead of eyeballing hex values, then apply them through
//! egui's `Style`/`Visuals`. Accent-coloured elements (primary button, energy
//! badge) are drawn by the small helpers here.

use egui::{
    Button, Color32, CornerRadius, FontId, Frame, Margin, Response, RichText, Stroke, TextStyle,
    Ui, Vec2, vec2,
};
use std::sync::LazyLock;

/// Uniform height (px) for every toolbar control — buttons, combo boxes and the
/// energy badge — so the single toolbar row lines up cleanly.
pub const CTRL_H: f32 = 32.0;

/// Convert an OKLCH colour (as written in the design's CSS) to sRGB.
/// `l` in 0..1, `c` chroma, `h` hue in degrees. Follows Björn Ottosson's
/// OKLab → linear-sRGB matrices, then applies the sRGB transfer function.
fn oklch(l: f32, c: f32, h_deg: f32) -> Color32 {
    let h = h_deg.to_radians();
    let a = c * h.cos();
    let b = c * h.sin();

    // OKLab → LMS (cube-rooted)
    let l_ = l + 0.396_337_78 * a + 0.215_803_76 * b;
    let m_ = l - 0.105_561_35 * a - 0.063_854_17 * b;
    let s_ = l - 0.089_484_18 * a - 1.291_485_5 * b;
    let (l3, m3, s3) = (l_ * l_ * l_, m_ * m_ * m_, s_ * s_ * s_);

    // LMS → linear sRGB
    let r = 4.076_741_7 * l3 - 3.307_711_6 * m3 + 0.230_969_94 * s3;
    let g = -1.268_438 * l3 + 2.609_757_4 * m3 - 0.341_319_4 * s3;
    let bl = -0.004_196_086 * l3 - 0.703_418_6 * m3 + 1.707_614_7 * s3;

    Color32::from_rgb(to_srgb(r), to_srgb(g), to_srgb(bl))
}

fn to_srgb(c: f32) -> u8 {
    let c = c.clamp(0.0, 1.0);
    let s = if c <= 0.003_130_8 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    };
    (s.clamp(0.0, 1.0) * 255.0).round() as u8
}

/// The palette, computed once from the design's OKLCH values.
pub struct Palette {
    pub panel: Color32,
    pub toolbar_bg: Color32,
    pub surface: Color32,
    pub hover: Color32,
    pub border: Color32,
    pub border_soft: Color32,
    pub border_hover: Color32,

    pub text: Color32,
    pub text_strong: Color32,
    pub label: Color32,

    pub accent: Color32,
    pub accent_border: Color32,
    pub sel_bg: Color32,

    pub danger_bg: Color32,
    pub danger_border: Color32,
    pub danger_text: Color32,

    pub energy_bg: Color32,
    pub energy_border: Color32,
    pub energy_e: Color32,
    pub energy_text: Color32,
    pub energy_unit: Color32,
    pub energy_steps: Color32,
    pub energy_div: Color32,
    pub mm_play: Color32,
}

pub static PAL: LazyLock<Palette> = LazyLock::new(|| Palette {
    panel: oklch(0.995, 0.002, 100.0),
    toolbar_bg: oklch(0.985, 0.003, 100.0),
    surface: Color32::WHITE,
    hover: oklch(0.965, 0.004, 100.0),
    border: oklch(0.88, 0.005, 100.0),
    border_soft: oklch(0.91, 0.005, 100.0),
    border_hover: oklch(0.80, 0.02, 250.0),

    text: oklch(0.34, 0.01, 260.0),
    text_strong: oklch(0.30, 0.01, 260.0),
    label: oklch(0.60, 0.012, 260.0),

    accent: oklch(0.59, 0.135, 251.0),
    accent_border: oklch(0.53, 0.14, 250.0),
    sel_bg: oklch(0.94, 0.03, 250.0),

    danger_bg: oklch(0.98, 0.012, 25.0),
    danger_border: oklch(0.90, 0.02, 25.0),
    danger_text: oklch(0.50, 0.14, 25.0),

    energy_bg: oklch(0.965, 0.006, 145.0),
    energy_border: oklch(0.88, 0.02, 150.0),
    energy_e: oklch(0.50, 0.05, 150.0),
    energy_text: oklch(0.35, 0.06, 150.0),
    energy_unit: oklch(0.55, 0.03, 150.0),
    energy_steps: oklch(0.52, 0.02, 150.0),
    energy_div: oklch(0.85, 0.02, 150.0),
    mm_play: oklch(0.55, 0.16, 145.0),
});

/// Install the global egui style: light neutrals, rounded white controls,
/// blue selection, and the design's type scale. Applied to every theme slot
/// so the app always renders light regardless of the OS preference.
pub fn install(ctx: &egui::Context) {
    let p = &*PAL;
    let cr = CornerRadius::same(8);

    ctx.all_styles_mut(move |style| {
        let mut v = egui::Visuals::light();

        v.panel_fill = p.panel;
        v.window_fill = p.panel;
        v.window_stroke = Stroke::new(1.0, p.border);
        v.faint_bg_color = p.toolbar_bg;
        v.extreme_bg_color = p.surface;

        v.selection.bg_fill = p.sel_bg;
        v.selection.stroke = Stroke::new(1.0, p.accent_border);

        // Plain labels / separators.
        v.widgets.noninteractive.fg_stroke = Stroke::new(1.0, p.text);
        v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, p.border_soft);
        v.widgets.noninteractive.corner_radius = cr;

        // Idle controls: white face, soft border, rounded.
        v.widgets.inactive.bg_fill = p.surface;
        v.widgets.inactive.weak_bg_fill = p.surface;
        v.widgets.inactive.bg_stroke = Stroke::new(1.0, p.border);
        v.widgets.inactive.fg_stroke = Stroke::new(1.0, p.text);
        v.widgets.inactive.corner_radius = cr;

        // Hover.
        v.widgets.hovered.bg_fill = p.hover;
        v.widgets.hovered.weak_bg_fill = p.hover;
        v.widgets.hovered.bg_stroke = Stroke::new(1.0, p.border_hover);
        v.widgets.hovered.fg_stroke = Stroke::new(1.0, p.text_strong);
        v.widgets.hovered.corner_radius = cr;

        // Pressed / selected.
        v.widgets.active.bg_fill = p.sel_bg;
        v.widgets.active.weak_bg_fill = p.sel_bg;
        v.widgets.active.bg_stroke = Stroke::new(1.0, p.accent_border);
        v.widgets.active.fg_stroke = Stroke::new(1.0, p.text_strong);
        v.widgets.active.corner_radius = cr;

        v.widgets.open.corner_radius = cr;

        style.visuals = v;

        // One uniform control height for the whole toolbar so buttons, combos
        // and the energy badge all line up. `interact_size.y` is the floor that
        // buttons/combos/selectable-labels snap to; keep `button_padding.y`
        // small enough that it never exceeds this (otherwise heights diverge).
        style.spacing.item_spacing = vec2(7.0, 7.0);
        style.spacing.button_padding = vec2(12.0, 5.0);
        style.spacing.interact_size.y = CTRL_H;
        style.spacing.menu_margin = Margin::same(6);

        style.text_styles = [
            (TextStyle::Heading, FontId::proportional(16.0)),
            (TextStyle::Body, FontId::proportional(13.5)),
            (TextStyle::Button, FontId::proportional(13.0)),
            (TextStyle::Small, FontId::proportional(10.5)),
            (TextStyle::Monospace, FontId::monospace(12.5)),
        ]
        .into();
    });
}

/// Small section label, e.g. `力場` / `表示`.
pub fn group_label(ui: &mut Ui, text: &str) {
    ui.label(RichText::new(text).size(10.5).strong().color(PAL.label));
}

/// The primary blue call-to-action ("→ 3D 生成").
pub fn accent_button(ui: &mut Ui, label: &str) -> Response {
    let p = &*PAL;
    ui.add(
        Button::new(RichText::new(label).size(13.0).strong().color(Color32::WHITE))
            .fill(p.accent)
            .stroke(Stroke::new(1.0, p.accent_border))
            .corner_radius(CornerRadius::same(9))
            .min_size(vec2(0.0, CTRL_H)),
    )
}

/// A red-tinted destructive button ("クリア").
pub fn danger_button(ui: &mut Ui, label: &str) -> Response {
    let p = &*PAL;
    ui.add(
        Button::new(RichText::new(label).size(13.0).color(p.danger_text))
            .fill(p.danger_bg)
            .stroke(Stroke::new(1.0, p.danger_border))
            .corner_radius(CornerRadius::same(8))
            .min_size(vec2(0.0, CTRL_H)),
    )
}

/// The right-aligned read-out badge. Priority: `error` (red) → `energy` (green
/// `E value kJ/mol │ steps`) → muted `E 未計算`.
///
/// NOTE: must be placed inside a right-to-left toolbar cluster. `ui.horizontal`
/// then inherits that direction (which keeps the frame content-sized rather
/// than stretched), so the pieces are added in reverse to read left-to-right.
pub fn energy_badge(ui: &mut Ui, energy: Option<(f64, usize)>, error: Option<&str>) {
    let p = &*PAL;
    let (bg, border) = if error.is_some() {
        (p.danger_bg, p.danger_border)
    } else {
        (p.energy_bg, p.energy_border)
    };
    Frame::new()
        .fill(bg)
        .stroke(Stroke::new(1.0, border))
        .corner_radius(CornerRadius::same(9))
        // Vertical margin 0: the inner row already snaps to `CTRL_H`
        // (interact_size), so the badge ends up exactly control-height and
        // lines up with the buttons/combos. Text is centred within that row.
        .inner_margin(Margin::symmetric(13, 0))
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 8.0;
                if let Some(msg) = error {
                    // reversed → "⚠ {msg}"
                    ui.label(RichText::new(msg).size(11.5).color(p.danger_text));
                    ui.label(RichText::new("⚠").size(12.0).color(p.danger_text));
                } else if let Some((e, steps)) = energy {
                    // reversed → "E {value} kJ/mol │ {steps} steps"
                    ui.label(
                        RichText::new(format!("{steps} steps"))
                            .monospace()
                            .size(11.5)
                            .color(p.energy_steps),
                    );
                    ui.label(RichText::new("│").size(12.0).color(p.energy_div));
                    ui.label(
                        RichText::new("kJ/mol")
                            .monospace()
                            .size(11.0)
                            .color(p.energy_unit),
                    );
                    ui.label(
                        RichText::new(format!("{e:.3}"))
                            .monospace()
                            .size(13.5)
                            .strong()
                            .color(p.energy_text),
                    );
                    ui.label(RichText::new("E").monospace().size(11.0).color(p.energy_e));
                } else {
                    ui.label(
                        RichText::new("未計算")
                            .monospace()
                            .size(11.5)
                            .color(p.energy_steps),
                    );
                    ui.label(RichText::new("E").monospace().size(11.0).color(p.energy_e));
                }
            });
        });
}

/// A thin vertical divider matching the toolbar group separators.
pub fn divider(ui: &mut Ui) -> Response {
    ui.add(egui::Separator::default().vertical().spacing(14.0))
}

/// Convenience for a height-only min-size `Vec2`.
pub fn h(height: f32) -> Vec2 {
    vec2(0.0, height)
}
