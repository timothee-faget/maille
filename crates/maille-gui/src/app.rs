//! Fenêtre de l'interface : la vue 3D plein écran (le viewer, dessiné dans une texture montrée comme
//! une image) et, par-dessus, des panneaux flottants arrondis et semi-transparents. Ne décide de
//! rien : l'état est dans `state`.

use std::path::{Path, PathBuf};

use eframe::egui::{
    self, Align2, Color32, CornerRadius, Frame, Margin, Pos2, Rect, RichText, Sense, Stroke, Vec2,
};
use eframe::egui_wgpu::wgpu;
use egui_lucide::Lucide;
use maille_viewer::{Cut, ResultsView};

use crate::results;
use crate::state::{Action, Item, Space, State, Status, Tool};
use crate::triad;

/// Fond des panneaux : sombre et translucide, la vue reste lisible derrière.
const PANEL: Color32 = Color32::from_rgba_premultiplied(26, 28, 34, 224);
const GAP: f32 = 12.0;
const TOOLBAR_WIDTH: f32 = 36.0;
const ICON_SIZE: f32 = 18.0;
/// Marge du trièdre d'orientation au coin bas droit de la vue, et sa taille.
const TRIAD_MARGIN: f32 = 56.0;
const TRIAD_LENGTH: f32 = 34.0;
const INSPECTOR_WIDTH: f32 = 250.0;
/// Facteur de zoom par pixel de défilement, comme `maille view`.
const SCROLL_ZOOM: f64 = 1.005;

/// Résultats ouverts : la vue, sa texture côté egui, et le réglage de la coupe qu'on édite.
struct Open {
    view: ResultsView,
    texture: egui::TextureId,
    /// Taille de la texture, en pixels.
    size: [u32; 2],
    name: String,
    cut_axis: Option<usize>,
    cut_fraction: f32,
    cut_flipped: bool,
}

pub struct App {
    state: State,
    open: Option<Open>,
    /// Chemin tapé dans l'inspecteur « Ouvrir ».
    path: String,
    /// À ouvrir dès la première image (argument de la ligne de commande).
    pending: Option<PathBuf>,
    error: Option<String>,
}

impl App {
    pub fn new(open: Option<PathBuf>) -> Self {
        Self {
            state: State::default(),
            open: None,
            path: String::new(),
            pending: open,
            error: None,
        }
    }

    fn load(&mut self, path: &Path, frame: &eframe::Frame) {
        let Some(gpu) = frame.wgpu_render_state() else {
            self.error = Some("l'interface a besoin de wgpu".to_string());
            return;
        };
        let loaded = results::open(path, None)
            .and_then(|r| results::to_nodal(&r))
            .and_then(|nodal| {
                ResultsView::new(gpu.device.clone(), gpu.queue.clone(), &nodal, [1, 1], 1.0)
            });
        match loaded {
            Ok(view) => {
                let texture = gpu.renderer.write().register_native_texture(
                    &gpu.device,
                    &view.texture().create_view(&Default::default()),
                    wgpu::FilterMode::Linear,
                );
                if let Some(old) = self.open.take() {
                    gpu.renderer.write().free_texture(&old.texture);
                }
                self.open = Some(Open {
                    view,
                    texture,
                    size: [1, 1],
                    name: path
                        .file_name()
                        .map(|n| n.to_string_lossy().into_owned())
                        .unwrap_or_default(),
                    cut_axis: None,
                    cut_fraction: 0.5,
                    cut_flipped: false,
                });
                self.state.progress.results_open = true;
                self.state.go(Space::Results);
                self.state.set_tool(Tool::Navigate);
                self.error = None;
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    fn act(&mut self, action: Action) {
        if let Some(open) = self.open.as_mut() {
            match action {
                Action::Fit => open.view.fit(),
                Action::View(axis) => open.view.look_along(axis),
            }
        }
    }

    /// Dessine la vue dans la zone libre et lit la souris. Sans résultats : une invitation.
    fn viewport(&mut self, ui: &mut egui::Ui, frame: &eframe::Frame) {
        let rect = ui.max_rect();
        let ppp = ui.ctx().pixels_per_point();
        let (Some(open), Some(gpu)) = (self.open.as_mut(), frame.wgpu_render_state()) else {
            ui.painter().rect_filled(rect, 0.0, Color32::from_gray(242));
            ui.painter().text(
                rect.center(),
                Align2::CENTER_CENTER,
                "Dépose un fichier .pvd ou .vtu ici,\nou ouvre-le avec l'outil « Ouvrir » (Résultats)",
                egui::FontId::proportional(18.0),
                Color32::from_gray(110),
            );
            return;
        };
        let pixels = [
            (rect.width() * ppp).round().max(1.0) as u32,
            (rect.height() * ppp).round().max(1.0) as u32,
        ];
        if pixels != open.size {
            open.view.resize(pixels, f64::from(ppp));
            open.size = pixels;
            gpu.renderer.write().update_egui_texture_from_wgpu_texture(
                &gpu.device,
                &open.view.texture().create_view(&Default::default()),
                wgpu::FilterMode::Linear,
                open.texture,
            );
        }
        open.view.draw();
        ui.painter().image(
            open.texture,
            rect,
            Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
            Color32::WHITE,
        );

        draw_triad(ui, rect, open.view.orientation());
        let response = ui.allocate_rect(rect, Sense::click_and_drag());
        let delta = response.drag_delta();
        if delta != Vec2::ZERO {
            let pan = response.dragged_by(egui::PointerButton::Secondary)
                || ui.input(|i| i.modifiers.shift);
            open.view
                .drag([f64::from(delta.x * ppp), f64::from(delta.y * ppp)], pan);
        }
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if let (true, Some(at)) = (scroll != 0.0, response.hover_pos()) {
                let cursor = [
                    f64::from((at.x - rect.min.x) * ppp),
                    f64::from((at.y - rect.min.y) * ppp),
                ];
                open.view
                    .zoom_at(SCROLL_ZOOM.powf(f64::from(scroll * ppp)), cursor);
            }
        }
    }

    fn spaces(&mut self, ctx: &egui::Context) {
        panel("spaces", Align2::CENTER_TOP, [0.0, GAP]).show(ctx, |ui| {
            floating().show(ui, |ui| {
                ui.horizontal(|ui| {
                    for space in Space::ALL {
                        status_dot(ui, self.state.space_status(space));
                        let selected = self.state.space == space;
                        if ui.selectable_label(selected, space.label()).clicked() {
                            self.state.go(space);
                        }
                        ui.add_space(6.0);
                    }
                });
            });
        });
    }

    fn toolbar(&mut self, ctx: &egui::Context) {
        let groups = crate::state::toolbar(self.state.space);
        let mut clicked = None;
        panel("toolbar", Align2::RIGHT_CENTER, [-GAP, 0.0]).show(ctx, |ui| {
            floating().show(ui, |ui| {
                ui.set_width(TOOLBAR_WIDTH);
                for (g, group) in groups.iter().enumerate() {
                    if g > 0 {
                        ui.separator();
                    }
                    for &item in &group.items {
                        let (label, selected) = match item {
                            Item::Tool(t) => (t.label(), self.state.tool() == t),
                            Item::Action(a) => (a.label(), false),
                        };
                        let button = match item {
                            Item::Action(Action::View(axis)) => egui::Button::new(
                                RichText::new(["X", "Y", "Z"][axis])
                                    .color(AXIS_COLORS[axis])
                                    .strong(),
                            ),
                            _ => egui::Button::image(
                                icon(item)
                                    .size(ICON_SIZE)
                                    .color(Color32::from_gray(228))
                                    .stroke_width(1.75)
                                    .image(),
                            ),
                        }
                        .selected(selected)
                        .min_size(Vec2::splat(TOOLBAR_WIDTH));
                        if ui.add(button).on_hover_text(label).clicked() {
                            clicked = Some(item);
                        }
                    }
                }
            });
        });
        match clicked {
            Some(Item::Tool(t)) => self.state.pick(t),
            Some(Item::Action(a)) => self.act(a),
            None => {}
        }
    }

    fn inspector(&mut self, ctx: &egui::Context, frame: &eframe::Frame) {
        if !self.state.inspector {
            return;
        }
        let tool = self.state.tool();
        let mut to_open = None;
        let mut hide = false;
        panel(
            "inspector",
            Align2::RIGHT_CENTER,
            [-(GAP * 2.0 + TOOLBAR_WIDTH + 22.0), 0.0],
        )
        .show(ctx, |ui| {
            floating().show(ui, |ui| {
                ui.set_width(INSPECTOR_WIDTH);
                ui.horizontal(|ui| {
                    ui.heading(tool.label());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let close = egui::Button::image(
                            Lucide::PanelRightClose
                                .size(ICON_SIZE)
                                .color(Color32::from_gray(200))
                                .image(),
                        )
                        .frame(false);
                        let tip = "Masquer (un clic sur l'outil le remontre)";
                        if ui.add(close).on_hover_text(tip).clicked() {
                            hide = true;
                        }
                    });
                });
                ui.separator();
                match tool {
                    Tool::Navigate => {
                        ui.label(
                            "Glisser : tourner (3D) ou déplacer (2D).\n\
                             Maj + glisser ou clic droit : déplacer.\n\
                             Molette : zoomer. F : cadrer. X, Y, Z : vues.",
                        );
                    }
                    Tool::Open => {
                        ui.label("Fichier .pvd ou .vtu :");
                        ui.text_edit_singleline(&mut self.path);
                        if ui.button("Ouvrir").clicked() && !self.path.trim().is_empty() {
                            to_open = Some(PathBuf::from(self.path.trim()));
                        }
                        ui.weak("Ou dépose-le sur la fenêtre.");
                    }
                    Tool::Field | Tool::Deformation | Tool::Cut => match self.open.as_mut() {
                        Some(open) => results_tool(ui, tool, open),
                        None => {
                            ui.label("Ouvre d'abord un calcul (outil « Ouvrir »).");
                        }
                    },
                    _ => {
                        ui.label(
                            "Cet outil n'est pas encore disponible : il arrive dans une \
                             prochaine étape de M13.",
                        );
                    }
                }
                if let Some(e) = &self.error {
                    ui.add_space(6.0);
                    ui.colored_label(Color32::from_rgb(240, 110, 100), e);
                }
            });
        });
        if hide {
            self.state.inspector = false;
        }
        if let Some(path) = to_open {
            self.load(&path, frame);
        }
    }

    /// Carte du parcours : la liste de contrôle de l'espace, désactivable.
    fn card(&mut self, ctx: &egui::Context) {
        if !self.state.cards {
            panel("card-off", Align2::LEFT_BOTTOM, [GAP, -GAP]).show(ctx, |ui| {
                if ui.button("? Parcours").clicked() {
                    self.state.cards = true;
                }
            });
            return;
        }
        let space = self.state.space;
        let steps = self.state.checklist(space);
        let mut clicked = None;
        let mut close = false;
        panel("card", Align2::LEFT_BOTTOM, [GAP, -GAP]).show(ctx, |ui| {
            floating().show(ui, |ui| {
                ui.set_width(220.0);
                ui.horizontal(|ui| {
                    ui.label(RichText::new(space.label()).strong());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .small_button("×")
                            .on_hover_text("Masquer les cartes")
                            .clicked()
                        {
                            close = true;
                        }
                    });
                });
                ui.separator();
                for step in &steps {
                    ui.horizontal(|ui| {
                        status_dot(ui, step.status);
                        let text = if step.optional {
                            format!("{} (facultatif)", step.label)
                        } else {
                            step.label.to_string()
                        };
                        if ui.link(text).clicked() {
                            clicked = Some(step.tool);
                        }
                    });
                }
            });
        });
        if let Some(tool) = clicked {
            self.state.set_tool(tool);
        }
        if close {
            self.state.cards = false;
        }
    }
}

/// Couleurs des axes x, y, z (rouge, vert, bleu, comme partout).
const AXIS_COLORS: [Color32; 3] = [
    Color32::from_rgb(226, 80, 80),
    Color32::from_rgb(90, 190, 100),
    Color32::from_rgb(80, 130, 235),
];

/// Icône d'un élément de la barre d'outils.
fn icon(item: Item) -> Lucide {
    match item {
        Item::Tool(t) => match t {
            Tool::Navigate => Lucide::MousePointer2,
            Tool::Import => Lucide::FileInput,
            Tool::MeshSize => Lucide::Grid3x3,
            Tool::FaceSets => Lucide::Square,
            Tool::VertexSets => Lucide::Dot,
            Tool::ManageSets => Lucide::Layers,
            Tool::Quality => Lucide::Gauge,
            Tool::Material => Lucide::Package,
            Tool::Support => Lucide::Anchor,
            Tool::Load => Lucide::ArrowDownToLine,
            Tool::Thermal => Lucide::Flame,
            Tool::Contact => Lucide::CircleDot,
            Tool::StepAmplitude => Lucide::ChartLine,
            Tool::Run => Lucide::Play,
            Tool::Open => Lucide::FolderOpen,
            Tool::Field => Lucide::Palette,
            Tool::Deformation => Lucide::Spline,
            Tool::Cut => Lucide::Scissors,
            Tool::Probe => Lucide::Crosshair,
            Tool::Increments => Lucide::Clock,
            Tool::Curve => Lucide::ChartSpline,
            Tool::Capture => Lucide::Camera,
        },
        Item::Action(Action::Fit) => Lucide::Maximize,
        // les vues X, Y, Z portent leur lettre, pas une icône
        Item::Action(Action::View(_)) => Lucide::Axis3d,
    }
}

/// Trièdre d'orientation au coin bas droit de la vue : x rouge, y vert, z bleu, les axes qui
/// s'éloignent en plus pâle. Il suit la caméra : on voit toujours d'où l'on regarde.
fn draw_triad(ui: &egui::Ui, rect: Rect, (right, up, eye): ([f64; 3], [f64; 3], [f64; 3])) {
    let painter = ui.painter();
    let center = rect.right_bottom() - Vec2::splat(TRIAD_MARGIN);
    painter.circle_filled(center, TRIAD_LENGTH + 14.0, Color32::from_black_alpha(18));
    for axis in triad::axes(right, up, eye) {
        let index = ["X", "Y", "Z"]
            .iter()
            .position(|n| *n == axis.name)
            .unwrap_or(0);
        let color = if axis.toward < -0.05 {
            AXIS_COLORS[index].gamma_multiply(0.45)
        } else {
            AXIS_COLORS[index]
        };
        let tip = center + Vec2::new(axis.screen[0] as f32, axis.screen[1] as f32) * TRIAD_LENGTH;
        painter.line_segment([center, tip], Stroke::new(2.5, color));
        // un axe qui pointe droit vers nous n'a pas de direction : la lettre reste sur le centre
        let at = if (tip - center).length() < 4.0 {
            tip
        } else {
            tip + (tip - center).normalized() * 9.0
        };
        painter.text(
            at,
            Align2::CENTER_CENTER,
            axis.name,
            egui::FontId::proportional(13.0),
            color,
        );
    }
}

/// Un panneau flottant, collé à `anchor` de la fenêtre.
fn panel(id: &str, anchor: Align2, offset: [f32; 2]) -> egui::Area {
    egui::Area::new(egui::Id::unique(id))
        .anchor(anchor, offset)
        .order(egui::Order::Foreground)
}

/// Cadre commun des panneaux flottants.
fn floating() -> Frame {
    Frame::new()
        .fill(PANEL)
        .corner_radius(CornerRadius::same(14))
        .inner_margin(Margin::same(10))
        .stroke(Stroke::new(1.0, Color32::from_white_alpha(24)))
}

fn status_dot(ui: &mut egui::Ui, status: Status) {
    let color = match status {
        Status::Done => Color32::from_rgb(90, 200, 120),
        Status::Warning => Color32::from_rgb(240, 180, 70),
        Status::Missing => Color32::from_rgb(220, 90, 90),
    };
    let (rect, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
    ui.painter().circle_filled(rect.center(), 4.5, color);
}

/// Inspecteur des outils Champ, Déformée et Coupe : ils règlent la vue des résultats ouverts.
fn results_tool(ui: &mut egui::Ui, tool: Tool, open: &mut Open) {
    match tool {
        Tool::Field => {
            let current = open.view.quantity();
            for (k, label) in open.view.quantities().into_iter().enumerate() {
                if ui.selectable_label(current == k, label).clicked() {
                    open.view.set_quantity(k);
                }
            }
            ui.separator();
            if ui.button("Changer de colormap").clicked() {
                open.view.next_colormap();
            }
        }
        Tool::Deformation => {
            if open.view.has_displacements() {
                let mut deformed = open.view.deformed();
                if ui.checkbox(&mut deformed, "Montrer la déformée").changed() {
                    open.view.set_deformed(deformed);
                }
            } else {
                ui.label("Ce calcul n'a pas de déplacements (thermique).");
            }
        }
        Tool::Cut => {
            if !open.view.is_solid() {
                ui.label("Le plan de coupe n'existe qu'en 3D.");
                return;
            }
            let mut changed = false;
            ui.horizontal(|ui| {
                changed |= ui.radio_value(&mut open.cut_axis, None, "Aucune").changed();
                for (axis, name) in ["X", "Y", "Z"].into_iter().enumerate() {
                    changed |= ui
                        .radio_value(&mut open.cut_axis, Some(axis), name)
                        .changed();
                }
            });
            if open.cut_axis.is_some() {
                changed |= ui
                    .add(egui::Slider::new(&mut open.cut_fraction, 0.0..=1.0).text("position"))
                    .changed();
                changed |= ui
                    .checkbox(&mut open.cut_flipped, "Garder l'autre côté")
                    .changed();
            }
            if changed {
                open.view.set_cut(open.cut_axis.map(|axis| Cut {
                    axis,
                    fraction: f64::from(open.cut_fraction),
                    flipped: open.cut_flipped,
                }));
            }
        }
        _ => {}
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        if let Some(path) = self.pending.take() {
            self.load(&path, frame);
        }
        let dropped = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .map(|f| f.path().to_path_buf())
                .next()
        });
        if let Some(path) = dropped {
            self.load(&path, frame);
        }
        if !ctx.egui_wants_keyboard_input() {
            let (fit, axis) = ctx.input(|i| {
                (
                    i.key_pressed(egui::Key::F),
                    [egui::Key::X, egui::Key::Y, egui::Key::Z]
                        .iter()
                        .position(|&k| i.key_pressed(k)),
                )
            });
            if fit {
                self.act(Action::Fit);
            }
            if let Some(axis) = axis {
                self.act(Action::View(axis));
            }
        }

        egui::CentralPanel::default()
            .frame(Frame::NONE)
            .show(ui, |ui| self.viewport(ui, frame));
        self.spaces(&ctx);
        self.card(&ctx);
        self.toolbar(&ctx);
        self.inspector(&ctx, frame);
        if let Some(open) = &self.open {
            ctx.send_viewport_cmd(egui::ViewportCommand::Title(format!(
                "Maille — {}",
                open.name
            )));
        }
    }
}
