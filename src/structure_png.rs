//! Rasterise the 2D structure editor's molecule to a PNG on every platform.
//!
//! chembuider-rs has its own `molecule::image::molecule_to_png`, but that module
//! is `#[cfg(windows)]` — it lives next to the ChemDraw clipboard / EMF plumbing
//! that genuinely needs the `windows`/`gdi32` stack, and its `resvg` dependency
//! is declared only under `[target.'cfg(windows)'.dependencies]`. The PNG path
//! itself is pure resvg and works anywhere, so we rebuild it here from the
//! crate's *public* scene API (`chembuider_rs::widget::scene`) and depend on
//! resvg ourselves. This keeps "2D 構造式 (PNG)" available on Linux and macOS
//! without patching the upstream crate.
//!
//! The geometry is produced exactly as upstream does — same `SCALE_FACTOR`, same
//! `build_labels`/`build_bonds`/`build_ring_fills`, same egui font stack for
//! label layout — so the exported picture matches what the on-screen editor and
//! the Windows exporter draw.

use std::sync::{Arc, OnceLock};

use chembuider_rs::config::StyleConfig;
use chembuider_rs::molecule::Molecule;
use chembuider_rs::widget::SCALE_FACTOR;
use chembuider_rs::widget::scene::{self, MoleculeScene, RenderTransform};
use eframe::egui;
use egui::Vec2;
use resvg::usvg::fontdb;

/// Extra pixels of empty space around the drawn structure.
const MARGIN_PX: f32 = 12.0;

/// System fonts plus egui's own bundled proportional font, loaded once. The
/// returned family name is egui's font (so the SVG labels rasterise with the
/// same glyphs the widget shows); it falls back to `"sans-serif"` if the bundled
/// font can't be extracted.
fn egui_font() -> (Arc<fontdb::Database>, String) {
    static DB: OnceLock<(Arc<fontdb::Database>, String)> = OnceLock::new();
    DB.get_or_init(|| {
        let mut db = fontdb::Database::new();
        db.load_system_fonts();
        let mut family = String::from("sans-serif");
        let defs = egui::FontDefinitions::default();
        if let Some(first) = defs
            .families
            .get(&egui::FontFamily::Proportional)
            .and_then(|names| names.first())
            && let Some(data) = defs.font_data.get(first)
        {
            let before: Vec<fontdb::ID> = db.faces().map(|f| f.id).collect();
            db.load_font_data(data.font.to_vec());
            if let Some(face) = db.faces().find(|f| !before.contains(&f.id))
                && let Some((fam, _)) = face.families.first()
            {
                family = fam.clone();
            }
        }
        (Arc::new(db), family)
    })
    .clone()
}

/// Build the shared scene at `SCALE_FACTOR` px/unit (offset 0), plus the pixel
/// bounding rect of everything in it (ring fills, bonds, labels, bare atoms).
fn build_scene(
    mol: &Molecule,
    style: &StyleConfig,
    layout: impl FnMut(egui::text::LayoutJob) -> Arc<egui::Galley>,
) -> (MoleculeScene, egui::Rect) {
    let xf = RenderTransform {
        coord_scale: SCALE_FACTOR,
        offset: Vec2::ZERO,
        visual_scale: 1.0,
    };
    let labels = scene::build_labels(mol, style, xf, layout);
    let bonds = scene::build_bonds(mol, style, xf, &labels);
    let ring_fills = scene::build_ring_fills(mol, xf);

    let mut min = egui::pos2(f32::INFINITY, f32::INFINITY);
    let mut max = egui::pos2(f32::NEG_INFINITY, f32::NEG_INFINITY);
    let mut grow = |p: egui::Pos2| {
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
    };
    for a in &mol.atoms {
        grow(xf.apply(a.pos));
    }
    for poly in &ring_fills {
        poly.points.iter().for_each(|p| grow(*p));
    }
    for b in &bonds {
        for s in &b.strokes {
            grow(s.a);
            grow(s.b);
        }
        for f in &b.fills {
            f.iter().for_each(|p| grow(*p));
        }
    }
    for l in labels.values() {
        grow(l.rect.min);
        grow(l.rect.max);
    }
    let bounds = egui::Rect::from_min_max(min, max);
    (
        MoleculeScene {
            ring_fills,
            bonds,
            labels,
            label_bg_pad: style.label_bg_pad,
        },
        bounds,
    )
}

fn xml_escape(c: char) -> String {
    match c {
        '&' => "&amp;".into(),
        '<' => "&lt;".into(),
        '>' => "&gt;".into(),
        '"' => "&quot;".into(),
        '\'' => "&apos;".into(),
        other => other.to_string(),
    }
}

/// Render the molecule to an SVG string from the shared scene: ring-fill
/// polygons, then bond strokes / wedge fills, then white label masks, then each
/// label glyph placed at egui's exact layout position.
fn molecule_to_svg(mol: &Molecule) -> String {
    let style = StyleConfig::default();
    let family = egui_font().1;

    // Standalone egui font stack: the same layout egui does on-screen, so label
    // geometry matches.
    let mut fonts = egui::epaint::text::Fonts::new(
        egui::epaint::text::TextOptions::default(),
        egui::FontDefinitions::default(),
    );
    let mut view = fonts.with_pixels_per_point(1.0);
    let (scene, bounds) = build_scene(mol, &style, |job| view.layout_job(job));

    let w = bounds.width() + 2.0 * MARGIN_PX;
    let h = bounds.height() + 2.0 * MARGIN_PX;
    let dx = MARGIN_PX - bounds.min.x;
    let dy = MARGIN_PX - bounds.min.y;

    let mut s = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" viewBox="0 0 {w:.0} {h:.0}"><g transform="translate({dx:.3},{dy:.3})">"#
    );

    // Ring fills (backmost). Translucent RGBA → fill + fill-opacity.
    for poly in &scene.ring_fills {
        let pts = poly
            .points
            .iter()
            .map(|p| format!("{:.2},{:.2}", p.x, p.y))
            .collect::<Vec<_>>()
            .join(" ");
        let [r, g, b, a] = poly.rgba;
        s.push_str(&format!(
            r#"<polygon points="{pts}" fill="rgb({r},{g},{b})" fill-opacity="{:.3}"/>"#,
            a as f32 / 255.0
        ));
    }

    // Bonds: stroked lines + filled (black) stereo wedges.
    for bond in &scene.bonds {
        for st in &bond.strokes {
            s.push_str(&format!(
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="black" stroke-width="{:.3}" stroke-linecap="round"/>"#,
                st.a.x, st.a.y, st.b.x, st.b.y, st.width
            ));
        }
        for fill in &bond.fills {
            let pts = fill
                .iter()
                .map(|p| format!("{:.2},{:.2}", p.x, p.y))
                .collect::<Vec<_>>()
                .join(" ");
            s.push_str(&format!(r#"<polygon points="{pts}" fill="black"/>"#));
        }
    }

    // White masks behind labels (hide any crossing bonds), then the label glyphs.
    for label in scene.labels.values() {
        let r = label.rect.expand(scene.label_bg_pad);
        s.push_str(&format!(
            r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="white"/>"#,
            r.min.x,
            r.min.y,
            r.width(),
            r.height()
        ));
    }
    for label in scene.labels.values() {
        let base = label.rect.min;
        for g in label.galley.glyph_placements() {
            let ch = xml_escape(g.chr);
            let style_attr = if g.italic {
                r#" font-style="italic""#
            } else {
                ""
            };
            let x = base.x + g.x;
            let y = base.y + g.baseline_y;
            let mut emit = |gx: f32| {
                s.push_str(&format!(
                    r#"<text x="{gx:.3}" y="{y:.3}" font-size="{:.3}" font-family="{family}"{style_attr} fill="black" xml:space="preserve">{ch}</text>"#,
                    g.size
                ));
            };
            emit(x);
            if let Some(dxb) = g.bold_offset {
                emit(x + dxb); // faux-bold double strike, as egui does on-screen
            }
        }
    }

    s.push_str("</g></svg>");
    s
}

/// PNG bytes of the structure with a transparent background, or `None` when the
/// molecule is empty or the render fails.
pub fn molecule_to_png(mol: &Molecule) -> Option<Vec<u8>> {
    if mol.atoms.is_empty() {
        return None;
    }
    let svg = molecule_to_svg(mol);
    let opt = resvg::usvg::Options {
        fontdb: egui_font().0,
        ..Default::default()
    };
    let tree = resvg::usvg::Tree::from_str(&svg, &opt).ok()?;

    let size = tree.size();
    let w = size.width().ceil() as u32;
    let h = size.height().ceil() as u32;
    if w == 0 || h == 0 || w > 4000 || h > 4000 {
        return None;
    }
    let mut pixmap = resvg::tiny_skia::Pixmap::new(w, h)?;
    resvg::render(
        &tree,
        resvg::tiny_skia::Transform::identity(),
        &mut pixmap.as_mut(),
    );
    pixmap.encode_png().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chembuider_rs::molecule::BondOrder;

    /// The whole reason this module exists: the PNG path must run on the current
    /// (non-Windows) host, not just on Windows via `chembuider_rs::molecule::image`.
    #[test]
    fn png_has_a_valid_signature() {
        let mut mol = Molecule::default();
        let c = mol.add_atom("C".to_string(), [0.0, 0.0], 0);
        let o = mol.add_atom("O".to_string(), [1.5, 0.0], 0);
        mol.add_bond(c, o, BondOrder::Double);

        let png = molecule_to_png(&mol).expect("png rendered");
        assert_eq!(
            &png[0..8],
            &[0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A],
            "PNG signature"
        );
    }

    #[test]
    fn empty_molecule_renders_nothing() {
        assert!(molecule_to_png(&Molecule::default()).is_none());
    }
}
