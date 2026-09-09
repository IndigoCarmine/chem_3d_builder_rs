//! The 3D viewport, assembled here instead of using
//! `moleucle_3dview_rs::InteractiveMoleculeViewport`.
//!
//! The upstream widget picks atoms only: its `handle_interaction` calls
//! `MoleculeViewer::pick`, matches `ViewerEvent::AtomClicked` and drops the
//! `BondClicked` the same call already computed. It also binds the wheel
//! straight to `camera.dolly()` and the secondary button to panning, so a host
//! cannot see a bond click, a right-click or a scroll. The 3D editing gestures
//! need all three.
//!
//! So the render half of `InteractiveMoleculeViewport::show` (0.11.0) is
//! reproduced below against the crate's public parts — `MoleculeViewer`,
//! `OrbitalCamera`, `OffscreenRenderer`, `RenderFrameState` — and only the
//! input half differs. Keep the gestures identical to upstream's when updating
//! the dependency; the three deliberate differences are marked `DIFFERS`.

use eframe::egui::{self, PointerButton, Sense};
use lin_alg::f32::Vec3;
use moleucle_3dview_rs::{
    AtomGroupRender, AtomPairRender, Molecule, MoleculeViewer, OffscreenRenderer, RenderFrameState,
    RenderStyle, SharedRenderStates,
    camera::{Camera, OrbitalCamera},
    new_shared_states, set_state_by_type,
    viewer::ViewerEvent,
};

/// What a ray through the pointer hit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hit {
    Atom(usize),
    Bond(usize),
    Nothing,
}

/// What the user did to the 3D view this frame.
#[derive(Clone, Copy, Debug, Default)]
pub struct Frame3dResponse {
    pub primary: Option<Hit>,
    /// A right-click, not a right-drag. egui does not report a click after a
    /// drag, so this never fires while the user is panning.
    pub secondary: Option<Hit>,
    pub hovered: Option<Hit>,
    /// Wheel travel in lines — one notch on a mouse wheel — and non-zero only
    /// when `show` was called with `capture_scroll`, in which case the camera
    /// did not consume it.
    pub scroll_lines: f32,
}

/// Everything a hover pick depends on. Parking the pointer over the view would
/// otherwise re-run the whole ray test every frame; upstream caches on the same
/// key and the reason carries over unchanged.
#[derive(Clone, Copy, PartialEq, Eq)]
struct HoverKey {
    pointer: [u32; 2],
    view_proj: [u32; 16],
    geometry_revision: u64,
}

pub struct Viewport3d {
    viewer: MoleculeViewer,
    camera: OrbitalCamera,
    offscreen: OffscreenRenderer,
    shared_states: SharedRenderStates,
    last_hover: Option<(HoverKey, Hit)>,
    /// Size of the last rendered frame, for `read_rgba` to label its output.
    last_size: (u32, u32),
}

impl Viewport3d {
    pub fn new() -> Self {
        let mut offscreen = OffscreenRenderer::new();
        // Atom highlights (anchor / bond ends / the IK pair) and the bond axis.
        offscreen.add_additional_render(Box::new(AtomGroupRender::new()));
        offscreen.add_additional_render(Box::new(AtomPairRender::new()));
        // 0.11.0 defaults to `Circles`; this app has always opened on ball and
        // stick, and bonds have to be visible to be clicked.
        offscreen.set_render_style(RenderStyle::BallStick);

        Self {
            viewer: MoleculeViewer::new(),
            camera: OrbitalCamera::default(),
            offscreen,
            shared_states: new_shared_states(),
            last_hover: None,
            last_size: (1, 1),
        }
    }

    pub fn set_molecule(&mut self, molecule: Molecule) {
        self.viewer.set_molecule(molecule);
        self.last_hover = None;
    }

    pub fn has_molecule(&self) -> bool {
        self.viewer.molecule.is_some()
    }

    /// Update atom positions in place from Ångström coordinates.
    pub fn update_positions_angstrom(&mut self, coords: &[[f32; 3]]) -> Result<(), String> {
        self.last_hover = None;
        self.viewer.update_positions_angstrom(coords)
    }

    pub fn focus_on_molecule_center(&mut self) {
        if let Some(molecule) = self.viewer.molecule.as_ref() {
            self.camera.center = molecule.center();
            self.camera.radius = molecule.radius() * 2.0;
        }
    }

    /// Publish an overlay's state to the renderer.
    pub fn set_state_by_type<T: 'static + Send + Sync>(&mut self, state: T) {
        set_state_by_type(&self.shared_states, state);
    }

    pub fn render_style(&self) -> RenderStyle {
        self.offscreen.render_style()
    }

    pub fn set_render_style(&mut self, render_style: RenderStyle) {
        self.offscreen.set_render_style(render_style);
    }

    pub fn free_egui_texture(&mut self, render_state: &egui_wgpu::RenderState) {
        self.offscreen.free_egui_texture(render_state);
    }

    /// The last rendered frame as `(width, height, rgba)`. Blocks on the GPU
    /// copy, so this is for an explicit export, not per frame.
    pub fn read_rgba(
        &self,
        render_state: &egui_wgpu::RenderState,
    ) -> Result<(u32, u32, Vec<u8>), String> {
        let rgba = self.offscreen.read_rgba(render_state)?;
        Ok((self.last_size.0, self.last_size.1, rgba))
    }

    pub fn show(
        &mut self,
        ui: &mut egui::Ui,
        render_state: &egui_wgpu::RenderState,
        capture_scroll: bool,
    ) -> Result<Frame3dResponse, String> {
        let available = ui.available_size_before_wrap();
        let width = available.x.max(1.0) as u32;
        let height = available.y.max(1.0) as u32;
        self.camera.set_aspect(width as f32 / height as f32);
        self.last_size = (width, height);

        if let Some(molecule) = self.viewer.molecule.as_ref() {
            let distance = (self.camera.position() - molecule.center()).magnitude();
            self.offscreen.submit_lod_distance(distance);
        }

        self.offscreen
            .ensure_resources(render_state, width, height)?;

        let view_proj = self.camera.view_projection().data;
        let cam_rot = self.camera.camera_rotation();
        let camera_right = cam_rot.rotate_vec(Vec3::new(1.0, 0.0, 0.0));
        let camera_up = cam_rot.rotate_vec(Vec3::new(0.0, 1.0, 0.0));
        let camera_forward = cam_rot.rotate_vec(Vec3::new(0.0, 0.0, 1.0));

        let frame = RenderFrameState::new(
            self.viewer.molecule.as_ref(),
            view_proj,
            Some(self.camera.position()),
            self.camera.fov_y(),
            camera_right,
            camera_up,
            camera_forward,
            self.viewer.color_fn,
            Some(&self.shared_states),
            self.offscreen.render_style(),
            self.offscreen.mesh_resolution(),
            self.offscreen.is_low_mode(),
            self.viewer.molecule_opacity,
        )
        .with_geometry_revision(self.viewer.revision())
        .with_visible_atoms(self.viewer.visible_atoms())
        .with_periodic_images(self.viewer.periodic_images())
        .with_atom_attrs(
            self.viewer.atom_radii.as_deref(),
            self.viewer.atom_colors.as_deref(),
        );

        self.offscreen
            .render_frame_with_state(render_state, &frame)?;

        let texture_id = self
            .offscreen
            .texture_id()
            .ok_or_else(|| "No texture id registered".to_string())?;

        let response = ui.add(
            egui::Image::from_texture(egui::load::SizedTexture::new(
                texture_id,
                egui::vec2(width as f32, height as f32),
            ))
            .sense(Sense::click_and_drag()),
        );

        let ctx = ui.ctx().clone();
        Ok(self.handle_interaction(&ctx, &response, capture_scroll))
    }

    /// Ray-pick at a pointer position local to `rect`.
    fn pick_at(&self, pointer: egui::Pos2, rect: egui::Rect) -> Hit {
        let local = pointer - rect.min;
        let (origin, dir) = self.camera.ray_from_screen(
            local.x,
            local.y,
            rect.width().max(1.0),
            rect.height().max(1.0),
        );
        match self.viewer.pick(origin, dir) {
            // DIFFERS (1/3): upstream matches `AtomClicked` only and throws the
            // bond away.
            Some(ViewerEvent::AtomClicked(i)) => Hit::Atom(i),
            Some(ViewerEvent::BondClicked(i)) => Hit::Bond(i),
            _ => Hit::Nothing,
        }
    }

    fn hover_key(&self, pointer: egui::Pos2) -> HoverKey {
        let view_proj = self.camera.view_projection().data;
        let mut view_proj_bits = [0u32; 16];
        for (dst, src) in view_proj_bits.iter_mut().zip(view_proj.iter()) {
            *dst = src.to_bits();
        }
        HoverKey {
            pointer: [pointer.x.to_bits(), pointer.y.to_bits()],
            view_proj: view_proj_bits,
            geometry_revision: self.viewer.revision(),
        }
    }

    fn handle_interaction(
        &mut self,
        ctx: &egui::Context,
        response: &egui::Response,
        capture_scroll: bool,
    ) -> Frame3dResponse {
        let mut out = Frame3dResponse::default();

        // Orbiting is not pointing: the hover key includes the camera, so a drag
        // misses the cache every frame, exactly when the frame budget is
        // tightest.
        if response.hovered()
            && !response.dragged()
            && let Some(pointer) = response.hover_pos()
        {
            let key = self.hover_key(pointer);
            let hit = match self.last_hover {
                Some((cached, hit)) if cached == key => hit,
                _ => {
                    let hit = self.pick_at(pointer, response.rect);
                    self.last_hover = Some((key, hit));
                    hit
                }
            };
            out.hovered = Some(hit);
        }

        // DIFFERS (2/3): the host gets first refusal on the wheel, so a dihedral
        // rotation can claim it while the camera keeps it the rest of the time.
        //
        // The two want different readings of the same wheel. Camera motion wants
        // `smooth_scroll_delta`, which eases one notch out over several frames --
        // that is what upstream uses and why dollying feels continuous. A host
        // stepping something notch by notch wants the raw events, because the
        // smoothed value only adds up to a whole notch if every one of those
        // frames gets painted, and the tail of the easing arrives after the app
        // has gone back to sleep. Measured, that lost about a fifth of every
        // turn of the wheel.
        if response.hovered() {
            if capture_scroll {
                out.scroll_lines = wheel_lines(ctx);
            } else {
                let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
                if scroll.abs() > f32::EPSILON {
                    self.camera.dolly(scroll * 0.02);
                }
            }
        }

        // Primary drag: orbit (Shift+drag or middle/right drag: pan)
        if response.dragged_by(PointerButton::Primary) {
            let delta = response.drag_delta();
            if ctx.input(|i| i.modifiers.shift) {
                self.camera
                    .pan(lin_alg::f32::Vec2::new(delta.x * 0.01, delta.y * 0.01));
            } else {
                self.camera.orbit(delta.x * 0.005, delta.y * 0.005);
            }
        }

        if response.dragged_by(PointerButton::Secondary)
            || response.dragged_by(PointerButton::Middle)
        {
            let delta = response.drag_delta();
            self.camera
                .pan(lin_alg::f32::Vec2::new(delta.x * 0.01, delta.y * 0.01));
        }

        if response.clicked_by(PointerButton::Primary)
            && let Some(pointer) = response.interact_pointer_pos()
        {
            out.primary = Some(self.pick_at(pointer, response.rect));
        }

        // DIFFERS (3/3): upstream has no secondary click at all. A right *drag*
        // pans, and egui does not report a click after a drag, so the two do not
        // collide.
        if response.clicked_by(PointerButton::Secondary)
            && let Some(pointer) = response.interact_pointer_pos()
        {
            out.secondary = Some(self.pick_at(pointer, response.rect));
        }

        out
    }
}

impl Default for Viewport3d {
    fn default() -> Self {
        Self::new()
    }
}

/// This frame's unsmoothed wheel travel, in lines.
///
/// A mouse reports whole lines, one per notch, which is the unit a stepped
/// gesture wants. Trackpads report points instead, so those are divided back
/// through egui's own points-per-line to keep one "notch" meaning the same
/// amount of turn on both.
fn wheel_lines(ctx: &egui::Context) -> f32 {
    let points_per_line = ctx
        .options(|o| o.input_options.line_scroll_speed)
        .max(f32::EPSILON);
    ctx.input(|i| {
        i.raw
            .events
            .iter()
            .filter_map(|event| match event {
                egui::Event::MouseWheel { unit, delta, .. } => Some(match unit {
                    egui::MouseWheelUnit::Line => delta.y,
                    egui::MouseWheelUnit::Point => delta.y / points_per_line,
                    // A page is a viewport-full. Nothing sends this for a
                    // dihedral, but treating it as one line would make a
                    // page-scroll do almost nothing at all.
                    egui::MouseWheelUnit::Page => delta.y * 20.0,
                }),
                _ => None,
            })
            .sum()
    })
}
