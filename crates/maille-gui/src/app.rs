//! Fenêtre de l'interface : la vue 3D plein écran (le viewer, dessiné dans une texture montrée comme
//! une image) et, par-dessus, des panneaux flottants arrondis et semi-transparents. Ne décide de
//! rien : l'état est dans `state`.

use std::error::Error;
use std::path::{Path, PathBuf};

use eframe::egui::{
    self, Align2, Color32, CornerRadius, Frame, Margin, Pos2, Rect, RichText, Sense, Stroke, Vec2,
};
use eframe::egui_wgpu::{wgpu, RenderState};
use egui_lucide::Lucide;
use maille_io::meshing::Source;
use maille_viewer::{Cut, NodalResults, ResultsView, Theme};

use crate::compute::{ComputeState, Job, Running};
use crate::files::{self, Written};
use crate::library;
use crate::mesh_pick::{self, Rect as PickRect};
use crate::mesh_space::{self, Built, Cad, Imported};
use crate::model_space::{self, Diagnostic, Draft, Severity};
use crate::model_tools::{self, Forms, MeshInfo};
use crate::quality::POOR;
use crate::results;
use crate::results_ui::{self, Playback};
use crate::series::Series;
use crate::sets::{self, PickFace};
use crate::settings::{self, ThemeChoice};
use crate::state::{Action, Item, Space, State, Status, Tool};
use crate::triad;

/// Fond des panneaux : translucide, la vue reste lisible derrière.
const PANEL_DARK: Color32 = Color32::from_rgba_premultiplied(26, 28, 34, 224);
const PANEL_LIGHT: Color32 = Color32::from_rgba_premultiplied(238, 239, 242, 230);
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
pub(crate) struct Open {
    pub(crate) view: ResultsView,
    texture: egui::TextureId,
    /// Taille de la texture, en pixels.
    size: [u32; 2],
    name: String,
    cut_axis: Option<usize>,
    cut_fraction: f32,
    cut_flipped: bool,
    /// Les Increments du calcul, la lecture et la sonde (espace Résultats seulement).
    pub(crate) play: Option<Playback>,
}

/// Réglages de l'espace Maillage, tels qu'on les édite.
struct MeshUi {
    imported: Option<Imported>,
    built: Option<Built>,
    /// Taille h choisie, et maillage du second ordre.
    size: f64,
    quadratic: bool,
    /// Dernier message de l'espace (maillage écrit…).
    note: Option<String>,
    /// Faces de bord désignables à la souris, et entités choisies (outil « Faces → Set »).
    pick_faces: Vec<PickFace>,
    selection: Vec<String>,
    /// Nom du Set à créer, et Set examiné dans « Sets » (avec le nouveau nom qu'on lui tape).
    set_name: String,
    inspected: Option<String>,
    rename: String,
    /// Entités surlignées dans la vue en ce moment.
    shown: Vec<String>,
    /// Sélection faite sur le Mesh lui-même par l'outil `tool`, et son numéro de version (le
    /// surlignage ne se refait que s'il change).
    region: Option<(Tool, Region)>,
    region_version: u64,
    /// Outil « Cadre → Set » : ce que le cadre prend, s'il s'ajoute à la sélection, coin de départ.
    box_mode: BoxMode,
    box_add: bool,
    box_start: Option<[f64; 2]>,
    /// Outil « Faces → Set » sur un `.msh` : angle (degrés) au-delà duquel la région s'arrête.
    angle: f64,
}

/// Ce qu'une sélection sur le Mesh désigne.
#[derive(Debug, Clone, PartialEq)]
enum Region {
    Nodes(Vec<usize>),
    Elements(Vec<usize>),
    Faces(maille_io::msh::Set),
}

impl Region {
    fn count(&self) -> String {
        match self {
            Region::Nodes(n) => format!("{} Nodes", n.len()),
            Region::Elements(e) => format!("{} Elements", e.len()),
            Region::Faces(s) => format!("{} faces, {} Nodes", s.faces.len(), s.nodes.len()),
        }
    }

    /// Le Set du `.msh` qu'elle donnera.
    fn to_set(&self) -> maille_io::msh::Set {
        match self {
            Region::Nodes(n) => mesh_pick::node_set(n),
            Region::Elements(e) => mesh_pick::element_set(e),
            Region::Faces(s) => s.clone(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum BoxMode {
    Nodes,
    Elements,
}

/// Ce dont dépendent les diagnostics du modèle : si rien n'a changé, on ne les refait pas.
#[derive(Debug, Clone, PartialEq)]
struct ModelKey {
    draft: Draft,
    sets: Vec<String>,
    user: Vec<String>,
    nodes: usize,
    cells: usize,
}

/// Réglages de l'espace Modèle, tels qu'on les édite.
struct ModelUi {
    draft: Draft,
    library: model_tools::LibraryUi,
    forms: Forms,
    diags: Vec<Diagnostic>,
    key: Option<ModelKey>,
    note: Option<String>,
    /// Surfaces rigides dessinées en ce moment dans la vue.
    shown: Vec<maille_io::model::ContactDef>,
    /// Le calcul en cours, son état, et le modèle (clé) du dernier calcul dont les résultats sont
    /// ouverts : s'il change ensuite, les résultats sont périmés.
    job: Option<Job>,
    compute: ComputeState,
    launched: Option<ModelKey>,
    fresh: Option<ModelKey>,
}

/// Ce que l'inspecteur demande à l'application, une fois les panneaux dessinés.
enum Todo {
    Open(PathBuf),
    Import(PathBuf),
    Remesh,
    Save,
    CreateSet,
    CreateMeshSet,
    RemainingSet,
    DeleteSet(String),
    RenameSet(String, String),
    SaveModel,
    Run,
}

pub struct App {
    state: State,
    /// La vue des résultats (espace Résultats).
    results: Option<Open>,
    /// La vue du maillage (espaces Maillage et Modèle).
    mesh_view: Option<Open>,
    mesh: MeshUi,
    model: ModelUi,
    /// Le thème voulu, et le fichier où on le garde.
    theme: ThemeChoice,
    settings_path: Option<PathBuf>,
    /// Chemins tapés dans les inspecteurs « Ouvrir » et « Importer ».
    path: String,
    mesh_path: String,
    /// À ouvrir dès la première image (argument de la ligne de commande).
    pending: Option<PathBuf>,
    error: Option<String>,
    /// Ce que l'interface a écrit pendant la session, pour ne pas écraser le travail de l'utilisateur.
    written: Written,
    /// Des fichiers à écrire qui existent et ont été modifiés en dehors de l'interface : on attend
    /// l'accord de l'utilisateur.
    confirm: Option<Confirm>,
}

/// Écriture en attente d'accord, voir `App::write_or_ask`.
struct Confirm {
    files: Vec<files::File>,
    conflicts: Vec<PathBuf>,
    then: AfterWrite,
}

/// Ce qui suit l'écriture des fichiers.
#[derive(Clone, Copy)]
enum AfterWrite {
    /// Le maillage seul (« Écrire le .msh »).
    Mesh,
    /// Le modèle (et son maillage) ; le dernier fichier est le TOML.
    Model,
    /// Le modèle, puis le calcul sur ce TOML.
    Run,
}

impl App {
    pub fn new(open: Option<PathBuf>) -> Self {
        let settings_path = settings::default_path();
        Self {
            theme: settings::load(settings_path.as_deref()),
            settings_path,
            state: State::default(),
            results: None,
            mesh_view: None,
            mesh: MeshUi {
                imported: None,
                built: None,
                size: 1.0,
                quadratic: false,
                note: None,
                pick_faces: Vec::new(),
                selection: Vec::new(),
                set_name: String::new(),
                inspected: None,
                rename: String::new(),
                shown: Vec::new(),
                region: None,
                region_version: 0,
                box_mode: BoxMode::Nodes,
                box_add: false,
                box_start: None,
                angle: 30.0,
            },
            model: ModelUi {
                draft: Draft::default(),
                library: model_tools::LibraryUi::new(library::default_user_path()),
                forms: Forms::default(),
                diags: Vec::new(),
                key: None,
                note: None,
                shown: Vec::new(),
                job: None,
                compute: ComputeState::Idle,
                launched: None,
                fresh: None,
            },
            path: String::new(),
            mesh_path: String::new(),
            pending: open,
            error: None,
            written: Written::default(),
            confirm: None,
        }
    }

    /// Un fichier déposé ou ouvert : des résultats (`.pvd`, `.vtu`), sinon un maillage ou une géométrie.
    fn route(&mut self, path: &Path, frame: &eframe::Frame) {
        let is_results = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| matches!(e.to_ascii_lowercase().as_str(), "pvd" | "vtu"));
        if is_results {
            self.load(path, frame);
        } else {
            self.import(path, frame);
        }
    }

    fn load(&mut self, path: &Path, frame: &eframe::Frame) {
        self.model.fresh = None;
        let Some(gpu) = frame.wgpu_render_state() else {
            self.error = Some("l'interface a besoin de wgpu".to_string());
            return;
        };
        let loaded = Series::open(path).and_then(|mut series| {
            let last = series.len() - 1;
            let nodal = results::to_nodal(series.get(last)?)?;
            let mut open = open_view(gpu, &nodal, path)?;
            open.play = Some(Playback::new(series));
            Ok(open)
        });
        match loaded {
            Ok(open) => {
                if let Some(old) = self.results.replace(open) {
                    free_view(gpu, &old);
                }
                self.state.progress.results_open = true;
                self.state.progress.probed = false;
                self.state.go(Space::Results);
                self.state.set_tool(Tool::Navigate);
                self.error = None;
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Importe un dessin, une surface, une pièce ou un maillage. Une géométrie est maillée tout de suite
    /// à la taille proposée, pour qu'il y ait quelque chose à voir ; on règle h ensuite.
    fn import(&mut self, path: &Path, frame: &eframe::Frame) {
        match mesh_space::import(path) {
            Ok(Imported::Cad(cad)) => {
                self.mesh.size = cad.suggested;
                self.mesh.quadratic = false;
                self.mesh.note = None;
                self.mesh.built = None;
                self.mesh.imported = Some(Imported::Cad(cad));
                self.state.progress.geometry = true;
                self.state.progress.meshed = false;
                self.state.progress.sets = 0;
                self.state.go(Space::Mesh);
                self.state.set_tool(Tool::MeshSize);
                self.error = None;
                self.remesh(frame);
            }
            Ok(Imported::Mesh(file, built)) => {
                self.mesh.imported = Some(Imported::Mesh(file.clone(), built.clone()));
                self.mesh.note = None;
                self.state.progress.geometry = true;
                self.state.go(Space::Mesh);
                self.state.set_tool(Tool::Quality);
                self.error = None;
                self.show_built(*built, &file, frame);
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Remaille la géométrie importée avec la taille et l'ordre choisis.
    fn remesh(&mut self, frame: &eframe::Frame) {
        let Some(Imported::Cad(cad)) = &self.mesh.imported else {
            return;
        };
        let user = self
            .mesh
            .built
            .as_ref()
            .map(|b| b.user.clone())
            .unwrap_or_default();
        let lost = self
            .mesh
            .built
            .as_ref()
            .map(|b| b.mesh_sets.iter().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        match mesh_space::mesh(cad, self.mesh.size, self.mesh.quadratic, &user) {
            Ok(built) => {
                let path = cad.path.clone();
                self.mesh.note = (!lost.is_empty()).then(|| {
                    format!(
                        "Le maillage a changé : les Sets tirés sur le Mesh sont perdus ({}).",
                        lost.join(", ")
                    )
                });
                self.set_region(None);
                self.error = None;
                self.show_built(built, &path, frame);
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    fn show_built(&mut self, built: Built, path: &Path, frame: &eframe::Frame) {
        let Some(gpu) = frame.wgpu_render_state() else {
            self.error = Some("l'interface a besoin de wgpu".to_string());
            return;
        };
        match open_view(gpu, &built.nodal(), path) {
            Ok(open) => {
                if let Some(old) = self.mesh_view.replace(open) {
                    free_view(gpu, &old);
                }
                self.state.progress.meshed = true;
                self.state.progress.sets = built.sets().len();
                self.mesh.pick_faces = sets::pickable(&built);
                self.mesh.selection.clear();
                self.mesh.shown.clear();
                self.mesh.built = Some(built);
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// La vue à montrer : les résultats dans l'espace Résultats, le maillage dans les deux autres.
    fn view_mut(&mut self) -> Option<&mut Open> {
        match self.state.space {
            Space::Results => self.results.as_mut(),
            Space::Mesh | Space::Model => self.mesh_view.as_mut(),
        }
    }

    fn act(&mut self, action: Action) {
        if let Some(open) = self.view_mut() {
            match action {
                Action::Fit => open.view.fit(),
                Action::View(axis) => open.view.look_along(axis),
            }
        }
    }

    /// Dessine la vue dans la zone libre et lit la souris. Sans résultats : une invitation.
    /// Rend le pixel et le Maj d'un clic dans la vue, pour désigner une entité.
    fn viewport(&mut self, ui: &mut egui::Ui, frame: &eframe::Frame) -> Option<([f64; 2], bool)> {
        let rect = ui.max_rect();
        let ppp = ui.ctx().pixels_per_point();
        // « Cadre → Set » : le glisser trace un cadre au lieu de tourner la vue
        let box_tool = self.state.space == Space::Mesh && self.state.tool() == Tool::VertexSets;
        let box_start = self.mesh.box_start;
        let invitation = match self.state.space {
            Space::Results => {
                "Dépose un fichier .pvd ou .vtu ici,\nou ouvre-le avec l'outil « Ouvrir »"
            }
            _ => {
                "Dépose un dessin .dxf, une surface .stl, une pièce .step\nou un maillage .msh ici, ou utilise l'outil « Importer »"
            }
        };
        let dark = ui.visuals().dark_mode;
        let (Some(open), Some(gpu)) = (self.view_mut(), frame.wgpu_render_state()) else {
            let fill = if dark {
                Color32::from_rgb(37, 39, 43)
            } else {
                Color32::from_gray(204)
            };
            ui.painter().rect_filled(rect, 0.0, fill);
            ui.painter().text(
                rect.center(),
                Align2::CENTER_CENTER,
                invitation,
                egui::FontId::proportional(18.0),
                ui.visuals().weak_text_color(),
            );
            return None;
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
        open.view
            .set_theme(if dark { Theme::Dark } else { Theme::Light });
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
        let boxing = box_tool && response.dragged_by(egui::PointerButton::Primary);
        if delta != Vec2::ZERO && !boxing {
            let pan = response.dragged_by(egui::PointerButton::Secondary)
                || ui.input(|i| i.modifiers.shift);
            open.view
                .drag([f64::from(delta.x * ppp), f64::from(delta.y * ppp)], pan);
        }
        let pixel = |at: Pos2| {
            [
                f64::from((at.x - rect.min.x) * ppp),
                f64::from((at.y - rect.min.y) * ppp),
            ]
        };
        if response.hovered() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if let (true, Some(at)) = (scroll != 0.0, response.hover_pos()) {
                open.view
                    .zoom_at(SCROLL_ZOOM.powf(f64::from(scroll * ppp)), pixel(at));
            }
        }
        let click = if response.clicked() {
            let shift = ui.input(|i| i.modifiers.shift);
            response.interact_pointer_pos().map(|at| (pixel(at), shift))
        } else {
            None
        };
        // le cadre : départ au début du glisser, rendu à la fin
        let mut start = box_start;
        let mut finished = None;
        if box_tool {
            let primary = egui::PointerButton::Primary;
            let here = response.interact_pointer_pos();
            if response.drag_started_by(primary) {
                start = here.map(pixel);
            }
            if let (Some(s), Some(h)) = (start, here) {
                let corner = Pos2::new(
                    rect.min.x + (s[0] / f64::from(ppp)) as f32,
                    rect.min.y + (s[1] / f64::from(ppp)) as f32,
                );
                let band = Rect::from_two_pos(corner, h);
                ui.painter().rect_filled(
                    band,
                    0.0,
                    Color32::from_rgba_unmultiplied(255, 140, 0, 40),
                );
                ui.painter().rect_stroke(
                    band,
                    0.0,
                    egui::Stroke::new(1.5, Color32::from_rgb(255, 140, 0)),
                    egui::StrokeKind::Inside,
                );
                if response.drag_stopped_by(primary) {
                    finished = Some(PickRect::from_corners(s, pixel(h)));
                }
            }
            if response.drag_stopped_by(primary) {
                start = None;
            }
        }
        self.mesh.box_start = start;
        if let Some(band) = finished {
            self.finish_box(band);
        }
        click
    }

    fn spaces(&mut self, ctx: &egui::Context) {
        panel("spaces", Align2::CENTER_TOP, [0.0, GAP]).show(ctx, |ui| {
            floating(ui).show(ui, |ui| {
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
        let mut theme_clicked = false;
        panel("toolbar", Align2::RIGHT_CENTER, [-GAP, 0.0]).show(ctx, |ui| {
            floating(ui).show(ui, |ui| {
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
                                    .color(ui.visuals().text_color())
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
                ui.separator();
                let lucide = match self.theme {
                    ThemeChoice::System => Lucide::Monitor,
                    ThemeChoice::Light => Lucide::Sun,
                    ThemeChoice::Dark => Lucide::Moon,
                };
                let button = egui::Button::image(
                    lucide
                        .size(ICON_SIZE)
                        .color(ui.visuals().text_color())
                        .stroke_width(1.75)
                        .image(),
                )
                .min_size(Vec2::splat(TOOLBAR_WIDTH));
                let tip = format!("Thème : {} (clic pour changer)", self.theme.label());
                theme_clicked = ui.add(button).on_hover_text(tip).clicked();
            });
        });
        if theme_clicked {
            self.theme = self.theme.next();
            if let Some(path) = &self.settings_path {
                if let Err(e) = settings::save(path, self.theme) {
                    self.error = Some(format!("{} : {e}", path.display()));
                }
            }
        }
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
        let mut todo = None;
        let mut goto = None;
        let mut hide = false;
        panel(
            "inspector",
            Align2::RIGHT_CENTER,
            [-(GAP * 2.0 + TOOLBAR_WIDTH + 22.0), 0.0],
        )
        .show(ctx, |ui| {
            floating(ui).show(ui, |ui| {
                ui.set_width(INSPECTOR_WIDTH);
                ui.horizontal(|ui| {
                    ui.heading(tool.label());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let close = egui::Button::image(
                            Lucide::PanelRightClose
                                .size(ICON_SIZE)
                                .color(ui.visuals().weak_text_color())
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
                            todo = Some(Todo::Open(PathBuf::from(self.path.trim())));
                        }
                        ui.weak("Ou dépose-le sur la fenêtre.");
                    }
                    Tool::Import => {
                        ui.label("Dessin .dxf, surface .stl, pièce .step ou maillage .msh :");
                        ui.text_edit_singleline(&mut self.mesh_path);
                        if ui.button("Importer").clicked() && !self.mesh_path.trim().is_empty() {
                            todo = Some(Todo::Import(PathBuf::from(self.mesh_path.trim())));
                        }
                        ui.weak("Ou dépose-le sur la fenêtre.");
                        if let Some(built) = &self.mesh.built {
                            ui.separator();
                            ui.label(mesh_space::summary(built));
                        }
                    }
                    Tool::MeshSize => {
                        if let Some(t) = mesh_size_tool(ui, &mut self.mesh) {
                            todo = Some(t);
                        }
                    }
                    Tool::Quality => quality_tool(ui, &self.mesh),
                    Tool::FaceSets => {
                        if let Some(t) = face_sets_tool(ui, &mut self.mesh) {
                            todo = Some(t);
                        }
                    }
                    Tool::ManageSets => {
                        if let Some(t) = manage_sets_tool(ui, &mut self.mesh) {
                            todo = Some(t);
                        }
                    }
                    Tool::VertexSets => {
                        if let Some(t) = region_tool(ui, &mut self.mesh, Tool::VertexSets) {
                            todo = Some(t);
                        }
                    }
                    Tool::Material
                    | Tool::Support
                    | Tool::Load
                    | Tool::Thermal
                    | Tool::Contact
                    | Tool::StepAmplitude
                    | Tool::Run => match &self.mesh.built {
                        None => {
                            ui.label(
                                "Importe et maille d'abord une géométrie, ou importe un .msh.",
                            );
                        }
                        Some(built) => {
                            let sets = built.sets();
                            let element_sets = built.element_sets();
                            let info = MeshInfo {
                                sets: &sets,
                                element_sets: &element_sets,
                                kind: built.file.kind,
                            };
                            let m = &mut self.model;
                            match tool {
                                Tool::Material => {
                                    let asked = model_tools::material_tool(
                                        ui,
                                        &mut m.draft,
                                        &info,
                                        &mut m.library,
                                    );
                                    if let Some(model_tools::MaterialAction::RemainingSet) = asked {
                                        todo = Some(Todo::RemainingSet);
                                    }
                                }
                                Tool::Support => model_tools::support_tool(
                                    ui,
                                    &mut m.draft,
                                    &info,
                                    &mut m.forms.bc,
                                ),
                                Tool::Load => {
                                    model_tools::load_tool(ui, &mut m.draft, &info, &mut m.forms)
                                }
                                Tool::Thermal => {
                                    model_tools::thermal_tool(ui, &mut m.draft, &info, &mut m.forms)
                                }
                                Tool::Contact => model_tools::contact_tool(
                                    ui,
                                    &mut m.draft,
                                    &info,
                                    &mut m.forms.contact,
                                ),
                                Tool::StepAmplitude => model_tools::step_tool(ui, &mut m.draft),
                                _ => {
                                    let calculation = model_space::run_label(&m.draft, &built.file);
                                    let request = model_tools::run_tool(
                                        ui,
                                        calculation,
                                        &m.diags,
                                        m.note.as_deref(),
                                        &m.compute,
                                    );
                                    if request.run {
                                        todo = Some(Todo::Run);
                                    } else if request.save {
                                        todo = Some(Todo::SaveModel);
                                    }
                                    goto = request.goto;
                                }
                            }
                            if tool != Tool::Run {
                                ui.separator();
                                goto = model_tools::diagnostics_ui(ui, &m.diags, Some(tool));
                            }
                        }
                    },
                    Tool::Probe | Tool::Increments | Tool::Curve | Tool::Capture => {
                        match self.results.as_mut() {
                            Some(open) => {
                                let shown = match tool {
                                    Tool::Probe => results_ui::probe_tool(ui, open),
                                    Tool::Increments => results_ui::increments_tool(ui, open),
                                    Tool::Curve => results_ui::curve_tool(ui, open),
                                    _ => results_ui::capture_tool(ui, open),
                                };
                                if let Err(e) = shown {
                                    self.error = Some(e.to_string());
                                }
                            }
                            None => {
                                ui.label("Ouvre d'abord un calcul (outil « Ouvrir »).");
                            }
                        }
                    }
                    Tool::Field | Tool::Deformation | Tool::Cut => match self.results.as_mut() {
                        Some(open) => results_tool(ui, tool, open),
                        None => {
                            ui.label("Ouvre d'abord un calcul (outil « Ouvrir »).");
                        }
                    },
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
        if let Some(target) = goto {
            if target == Tool::MeshSize {
                self.state.go(Space::Mesh);
            }
            self.state.set_tool(target);
        }
        match todo {
            Some(Todo::SaveModel) => self.save_model(AfterWrite::Model, ctx),
            Some(Todo::Run) => self.start_run(ctx),
            Some(Todo::Open(path)) => self.load(&path, frame),
            Some(Todo::Import(path)) => self.import(&path, frame),
            Some(Todo::Remesh) => self.remesh(frame),
            Some(Todo::Save) => self.save(ctx),
            Some(Todo::CreateSet) => self.create_set(),
            Some(Todo::CreateMeshSet) => self.create_mesh_set(),
            Some(Todo::RemainingSet) => self.create_remaining_set(),
            Some(Todo::DeleteSet(name)) => self.edit_set(&name, None),
            Some(Todo::RenameSet(name, new_name)) => self.edit_set(&name, Some(&new_name)),
            None => {}
        }
    }

    /// Clic dans la vue avec l'outil « Faces → Set » : désigne l'entité sous le curseur (Maj : l'ajoute
    /// ou la retire de la sélection ; un clic dans le vide vide la sélection).
    fn pick_at(&mut self, cursor: [f64; 2], add: bool) {
        let on_geometry = matches!(self.mesh.imported, Some(Imported::Cad(_)));
        if self.state.space == Space::Mesh
            && self.state.tool() == Tool::FaceSets
            && matches!(self.mesh.imported, Some(Imported::Mesh(..)))
        {
            self.pick_skin_at(cursor, add);
            return;
        }
        if !(on_geometry && self.state.space == Space::Mesh && self.state.tool() == Tool::FaceSets)
        {
            return;
        }
        let (Some(open), Some(built)) = (&self.mesh_view, &self.mesh.built) else {
            return;
        };
        let hit = sets::pick(
            built,
            &self.mesh.pick_faces,
            &|p| open.view.project(p),
            cursor,
        )
        .map(str::to_string);
        let selection = &mut self.mesh.selection;
        match hit {
            Some(name) if add => match selection.iter().position(|n| *n == name) {
                Some(k) => {
                    selection.remove(k);
                }
                None => selection.push(name),
            },
            Some(name) => *selection = vec![name],
            None if !add => selection.clear(),
            None => {}
        }
    }

    /// Clic avec « Faces → Set » sur un `.msh` : la face de la peau sous le curseur, étendue aux voisines
    /// tant que l'angle reste sous le réglage. Maj : ajoute la région à la sélection.
    fn pick_skin_at(&mut self, cursor: [f64; 2], add: bool) {
        let (Some(open), Some(built)) = (&self.mesh_view, &self.mesh.built) else {
            return;
        };
        let skin = mesh_pick::skin(&built.file);
        let hit = mesh_pick::pick_skin(&built.file, &skin, &|p| open.view.project(p), cursor);
        let tool = Tool::FaceSets;
        let region = match hit {
            Some(i) => {
                let mut set = mesh_pick::angle_region(&built.file, &skin, i, self.mesh.angle);
                if let (true, Some((t, Region::Faces(old)))) = (add, &self.mesh.region) {
                    if *t == tool {
                        set.faces.extend(&old.faces);
                        set.nodes.extend(&old.nodes);
                        set.faces.sort_unstable();
                        set.faces.dedup();
                        set.nodes.sort_unstable();
                        set.nodes.dedup();
                    }
                }
                Some((tool, Region::Faces(set)))
            }
            None if add => return,
            None => None,
        };
        self.set_region(region);
    }

    /// Fin du glisser de « Cadre → Set » : les Nodes ou les Elements du cadre.
    fn finish_box(&mut self, band: PickRect) {
        let (Some(open), Some(built)) = (&self.mesh_view, &self.mesh.built) else {
            return;
        };
        let project = |p| open.view.project(p);
        let tool = Tool::VertexSets;
        let region = match self.mesh.box_mode {
            BoxMode::Nodes => {
                let mut ids = mesh_pick::nodes_in_rect(&built.file, &project, band);
                if let (true, Some((t, Region::Nodes(old)))) =
                    (self.mesh.box_add, &self.mesh.region)
                {
                    if *t == tool {
                        ids.extend(old);
                    }
                }
                ids.sort_unstable();
                ids.dedup();
                Region::Nodes(ids)
            }
            BoxMode::Elements => {
                let mut ids = mesh_pick::elements_in_rect(&built.file, &project, band);
                if let (true, Some((t, Region::Elements(old)))) =
                    (self.mesh.box_add, &self.mesh.region)
                {
                    if *t == tool {
                        ids.extend(old);
                    }
                }
                ids.sort_unstable();
                ids.dedup();
                Region::Elements(ids)
            }
        };
        self.set_region(Some((tool, region)));
    }

    fn set_region(&mut self, region: Option<(Tool, Region)>) {
        self.mesh.region = region;
        self.mesh.region_version += 1;
    }

    /// Crée le Set des Elements que les autres Materials ne couvrent pas, et le donne au premier.
    fn create_remaining_set(&mut self) {
        let Some(built) = self.mesh.built.as_mut() else {
            return;
        };
        let elements = self.model.draft.remaining_elements(&built.file);
        let name = (1..)
            .map(|k| {
                if k == 1 {
                    "reste".to_string()
                } else {
                    format!("reste{k}")
                }
            })
            .find(|n| !built.file.sets.contains_key(n))
            .unwrap_or_default();
        match built.add_mesh_set(&name, mesh_pick::element_set(&elements)) {
            Ok(()) => {
                self.model.draft.material_set = name.clone();
                self.mesh.note = Some(format!(
                    "Set « {name} » créé ({} Elements).",
                    elements.len()
                ));
                self.state.progress.sets = built.sets().len();
                self.error = None;
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Crée un Set sur le Mesh avec la sélection et le nom tapé.
    fn create_mesh_set(&mut self) {
        let (Some(built), Some((_, region))) = (self.mesh.built.as_mut(), &self.mesh.region) else {
            return;
        };
        match built.add_mesh_set(&self.mesh.set_name, region.to_set()) {
            Ok(()) => {
                self.mesh.note = Some(format!(
                    "Set « {} » créé ({}). Il est écrit dans le .msh.",
                    self.mesh.set_name.trim(),
                    region.count()
                ));
                self.mesh.set_name.clear();
                self.state.progress.sets = built.sets().len();
                self.error = None;
                self.set_region(None);
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Met la surbrillance de la vue du maillage à jour : la sélection pour « Faces → Set », le Set
    /// examiné pour « Sets », rien sinon.
    fn sync_highlight(&mut self) {
        // l'outil Contact dessine les surfaces rigides ; en le quittant, on rend le surlignage des Sets
        let rigid = match (self.state.space, self.state.tool()) {
            (Space::Model, Tool::Contact) => self.model.draft.contacts.clone(),
            _ => Vec::new(),
        };
        if rigid != self.model.shown {
            if let (Some(open), Some(built)) = (self.mesh_view.as_mut(), &self.mesh.built) {
                let triangles = if rigid.is_empty() {
                    sets::highlight(built, &self.mesh.shown)
                } else {
                    rigid
                        .iter()
                        .flat_map(|c| model_space::rigid_triangles(c, &built.file))
                        .collect()
                };
                open.view.set_highlight(triangles);
            }
            self.model.shown = rigid;
        }
        if !self.model.shown.is_empty() {
            return;
        }
        let wanted: Vec<String> = match (self.state.space, self.state.tool(), &self.mesh.built) {
            // une sélection sur le Mesh : une clé qui change avec elle (« # » : pas des entités)
            (Space::Mesh, tool @ (Tool::FaceSets | Tool::VertexSets), _)
                if self.mesh.region.as_ref().is_some_and(|(t, _)| *t == tool) =>
            {
                vec![format!("#region:{}", self.mesh.region_version)]
            }
            (Space::Mesh, Tool::FaceSets, _) => self.mesh.selection.clone(),
            (Space::Mesh, Tool::ManageSets, Some(built)) => match self.mesh.inspected.as_deref() {
                Some(name) if built.mesh_sets.contains(name) => vec![format!("#set:{name}")],
                Some(name) => built.entities_of(name),
                None => Vec::new(),
            },
            _ => Vec::new(),
        };
        if wanted == self.mesh.shown {
            return;
        }
        if let (Some(open), Some(built)) = (self.mesh_view.as_mut(), &self.mesh.built) {
            let triangles = match wanted.first().map(String::as_str) {
                Some(key) if key.starts_with('#') => {
                    let skin = mesh_pick::skin(&built.file);
                    match (&self.mesh.region, key.strip_prefix("#set:")) {
                        (_, Some(name)) => built
                            .file
                            .sets
                            .get(name)
                            .map(|set| mesh_pick::highlight_set(&built.file, &skin, set))
                            .unwrap_or_default(),
                        (Some((_, region)), None) => {
                            mesh_pick::highlight_set(&built.file, &skin, &region.to_set())
                        }
                        _ => Vec::new(),
                    }
                }
                _ => sets::highlight(built, &wanted),
            };
            open.view.set_highlight(triangles);
        }
        self.mesh.shown = wanted;
    }

    /// Crée un Set avec la sélection et le nom tapé.
    fn create_set(&mut self) {
        let Some(built) = self.mesh.built.as_mut() else {
            return;
        };
        match built.add_set(&self.mesh.set_name, &self.mesh.selection) {
            Ok(()) => {
                self.mesh.note = Some(format!("Set « {} » créé.", self.mesh.set_name.trim()));
                self.mesh.set_name.clear();
                self.mesh.selection.clear();
                self.state.progress.sets = built.sets().len();
                self.error = None;
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Supprime ou renomme un Set créé ici.
    fn edit_set(&mut self, name: &str, new_name: Option<&str>) {
        let Some(built) = self.mesh.built.as_mut() else {
            return;
        };
        match new_name {
            None => {
                built.remove_set(name);
                if self.mesh.inspected.as_deref() == Some(name) {
                    self.mesh.inspected = None;
                }
                self.error = None;
            }
            Some(new_name) => match built.rename_set(name, new_name) {
                Ok(()) => {
                    self.mesh.inspected = Some(new_name.trim().to_string());
                    self.error = None;
                }
                Err(e) => self.error = Some(e.to_string()),
            },
        }
        self.state.progress.sets = built.sets().len();
    }

    /// Refait les diagnostics du modèle si le brouillon ou le maillage ont changé, et met à jour la
    /// progression. Rend vrai si quelque chose a changé.
    fn refresh_model(&mut self) -> bool {
        let draft = &self.model.draft;
        let progress = &mut self.state.progress;
        progress.material = draft.material.is_some();
        progress.bcs = !draft.bcs.is_empty();
        progress.loaded = draft.has_loading();
        let Some(built) = &self.mesh.built else {
            progress.errors = 1;
            self.model.diags.clear();
            self.model.key = None;
            return false;
        };
        let key = ModelKey {
            draft: draft.clone(),
            sets: built.file.sets.keys().cloned().collect(),
            user: built.user.keys().cloned().collect(),
            nodes: built.file.nodes.len(),
            cells: built.file.cells.len(),
        };
        if self.model.key.as_ref() == Some(&key) {
            return false;
        }
        let user: Vec<&str> = built.user.keys().map(String::as_str).collect();
        self.model.diags = model_space::diagnose(draft, &built.file, &user);
        progress.errors = self
            .model
            .diags
            .iter()
            .filter(|d| d.severity == Severity::Error)
            .count();
        self.model.key = Some(key);
        true
    }

    /// Écrit le modèle (et le maillage, pour une géométrie) à côté de la source : `<source>_modele.toml`.
    fn save_model(&mut self, then: AfterWrite, ctx: &egui::Context) {
        match self.model_files() {
            Ok(files) => self.write_or_ask(files, then, ctx),
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    /// Les fichiers du modèle, à côté de la source : le maillage s'il faut l'écrire, puis
    /// `<source>_modele.toml` en dernier.
    fn model_files(&self) -> Result<Vec<files::File>, Box<dyn Error>> {
        let built = self.mesh.built.as_ref().ok_or("rien n'est maillé")?;
        let mut out = Vec::new();
        let (source, msh_path, meshing) = match &self.mesh.imported {
            Some(Imported::Cad(cad)) => {
                let msh = mesh_space::mesh_file(cad, built)?;
                let msh_path = msh.0.clone();
                out.push(msh);
                let name = cad
                    .path
                    .file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_default();
                let recipe =
                    mesh_space::recipe(&name, self.mesh.size, self.mesh.quadratic, &built.user);
                (cad.path.clone(), msh_path, Some(recipe))
            }
            Some(Imported::Mesh(path, _)) if built.mesh_sets.is_empty() => {
                (path.clone(), path.clone(), None)
            }
            // les Sets tirés sur le Mesh : un `.msh` à côté, l'original reste intact
            Some(Imported::Mesh(path, _)) => {
                let stem = path
                    .file_stem()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default();
                let with_sets = path.with_file_name(format!("{stem}_sets.msh"));
                out.push((
                    with_sets.clone(),
                    maille_io::msh::write(&built.file)?.into_bytes(),
                ));
                (path.clone(), with_sets, None)
            }
            None => return Err("rien n'est importé".into()),
        };
        let mesh_name = msh_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let text = self.model.draft.to_toml(&mesh_name, &built.file, meshing)?;
        let stem = source
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        out.push((
            source.with_file_name(format!("{stem}_modele.toml")),
            text.into_bytes(),
        ));
        Ok(out)
    }

    /// Écrit les fichiers, sauf si l'un d'eux a été modifié en dehors de l'interface : on demande
    /// alors avant de l'écraser (voir `confirm_dialog`).
    fn write_or_ask(&mut self, files: Vec<files::File>, then: AfterWrite, ctx: &egui::Context) {
        let conflicts = self.written.conflicts(&files);
        if conflicts.is_empty() {
            self.write(files, then, ctx);
        } else {
            self.confirm = Some(Confirm {
                files,
                conflicts,
                then,
            });
        }
    }

    fn write(&mut self, files: Vec<files::File>, then: AfterWrite, ctx: &egui::Context) {
        if let Err(e) = self.written.write(&files) {
            self.error = Some(e);
            return;
        }
        self.error = None;
        let Some((last, _)) = files.last() else {
            return;
        };
        let note = Some(format!("Écrit : {}", last.display()));
        match then {
            AfterWrite::Mesh => self.mesh.note = note,
            AfterWrite::Model => self.model.note = note,
            AfterWrite::Run => {
                self.model.note = note;
                let ctx = ctx.clone();
                self.model.job = Some(Job::start(last, move || ctx.request_repaint()));
                self.model.compute = ComputeState::Running(Running::default());
                self.model.launched = self.model.key.clone();
            }
        }
    }

    /// Demande avant d'écraser des fichiers modifiés en dehors de l'interface.
    fn confirm_dialog(&mut self, ctx: &egui::Context) {
        let Some(confirm) = &self.confirm else {
            return;
        };
        let mut choice = None;
        let modal = egui::Modal::new(egui::Id::unique("confirm-overwrite")).show(ctx, |ui| {
            ui.set_max_width(420.0);
            ui.heading("Écraser des fichiers ?");
            ui.label("Ces fichiers existent déjà et ont été modifiés en dehors de Maille :");
            for path in &confirm.conflicts {
                ui.monospace(path.display().to_string());
            }
            ui.label("Les écraser fait perdre ces modifications.");
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Écraser").clicked() {
                    choice = Some(true);
                }
                if ui.button("Annuler").clicked() {
                    choice = Some(false);
                }
            });
        });
        if modal.should_close() && choice.is_none() {
            choice = Some(false);
        }
        match choice {
            Some(true) => {
                if let Some(Confirm { files, then, .. }) = self.confirm.take() {
                    self.write(files, then, ctx);
                }
            }
            Some(false) => {
                if let Some(Confirm { then, .. }) = self.confirm.take() {
                    let note =
                        Some("Rien n'a été écrit : les fichiers existants sont gardés.".into());
                    match then {
                        AfterWrite::Mesh => self.mesh.note = note,
                        AfterWrite::Model | AfterWrite::Run => self.model.note = note,
                    }
                }
            }
            None => {}
        }
    }

    /// Enregistre le modèle puis lance le calcul dans un thread.
    fn start_run(&mut self, ctx: &egui::Context) {
        if self.model.job.is_some() || self.confirm.is_some() {
            return;
        }
        self.save_model(AfterWrite::Run, ctx);
    }

    /// Prend les messages du calcul ; à la fin, ouvre les résultats dans l'espace Résultats.
    fn poll_job(&mut self, frame: &eframe::Frame) {
        let Some(job) = &mut self.model.job else {
            return;
        };
        for message in job.poll() {
            self.model.compute.apply(message);
        }
        if !job.is_finished() {
            return;
        }
        self.model.job = None;
        if let ComputeState::Finished(pvd) = self.model.compute.clone() {
            self.load(&pvd, frame);
            if self.error.is_none() {
                self.model.fresh = self.model.launched.clone();
            }
        }
    }

    /// Panneau flottant du calcul en cours, dans tous les espaces.
    fn compute_panel(&mut self, ctx: &egui::Context) {
        let ComputeState::Running(running) = &self.model.compute else {
            return;
        };
        let mut cancel = false;
        panel("compute", Align2::CENTER_BOTTOM, [0.0, -GAP]).show(ctx, |ui| {
            floating(ui).show(ui, |ui| {
                ui.set_width(260.0);
                cancel = model_tools::progress_ui(ui, running);
            });
        });
        if cancel {
            if let Some(job) = &self.model.job {
                job.cancel();
            }
            if let ComputeState::Running(r) = &mut self.model.compute {
                r.stopping = true;
            }
        }
    }

    /// Bandeau de l'espace Résultats quand le modèle a changé depuis le calcul dont on voit les résultats.
    fn stale_banner(&mut self, ctx: &egui::Context) {
        let stale = self.model.fresh.is_some() && self.model.fresh != self.model.key;
        if !(stale && self.state.space == Space::Results) {
            return;
        }
        let mut relaunch = false;
        panel("stale", Align2::CENTER_TOP, [0.0, 60.0]).show(ctx, |ui| {
            floating(ui).show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.colored_label(
                        Color32::from_rgb(230, 180, 80),
                        "Résultats périmés : le modèle a changé depuis ce calcul.",
                    );
                    relaunch = ui.button("Relancer le calcul").clicked();
                });
            });
        });
        if relaunch {
            self.state.go(Space::Model);
            self.state.set_tool(Tool::Run);
        }
    }

    /// Écrit le maillage à côté de la source, comme `maille mesh`.
    fn save(&mut self, ctx: &egui::Context) {
        if let (Some(Imported::Cad(cad)), Some(built)) = (&self.mesh.imported, &self.mesh.built) {
            match mesh_space::mesh_file(cad, built) {
                Ok(file) => self.write_or_ask(vec![file], AfterWrite::Mesh, ctx),
                Err(e) => self.error = Some(e.to_string()),
            }
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
            floating(ui).show(ui, |ui| {
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

/// Montre `nodal` : une vue du viewer sur le device d'eframe, et sa texture côté egui.
fn open_view(gpu: &RenderState, nodal: &NodalResults, path: &Path) -> Result<Open, Box<dyn Error>> {
    let view = ResultsView::new(gpu.device.clone(), gpu.queue.clone(), nodal, [1, 1], 1.0)?;
    let texture = gpu.renderer.write().register_native_texture(
        &gpu.device,
        &view.texture().create_view(&Default::default()),
        wgpu::FilterMode::Linear,
    );
    Ok(Open {
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
        play: None,
    })
}

/// Libère la texture d'une vue remplacée.
fn free_view(gpu: &RenderState, old: &Open) {
    gpu.renderer.write().free_texture(&old.texture);
}

/// Inspecteur de la taille h : curseur, ordre des Elements, bouton Mailler, écriture du `.msh`.
fn mesh_size_tool(ui: &mut egui::Ui, mesh: &mut MeshUi) -> Option<Todo> {
    let Some(imported) = &mesh.imported else {
        ui.label("Importe d'abord un dessin, une surface ou une pièce (outil « Importer »).");
        return None;
    };
    let Imported::Cad(cad) = imported else {
        ui.label("Un .msh est déjà maillé. Importe le dessin, la surface ou la pièce pour changer la taille.");
        return None;
    };
    let mut todo = None;
    ui.label(format!(
        "{} ({})",
        cad.path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        format_name(cad)
    ));
    ui.label("Taille h des Elements, dans les unités du fichier :");
    let (lo, hi) = (cad.suggested / 8.0, cad.suggested * 4.0);
    ui.add(egui::Slider::new(&mut mesh.size, lo..=hi).logarithmic(true));
    ui.checkbox(
        &mut mesh.quadratic,
        "Elements du second ordre (Tri6, Tet10)",
    );
    if ui.button("Mailler").clicked() {
        todo = Some(Todo::Remesh);
    }
    if let Some(built) = &mesh.built {
        ui.separator();
        ui.label(mesh_space::summary(built));
        if let Some(unit) = &built.unit {
            ui.label(format!(
                "Unité du fichier STEP : {unit} (rien n'est converti)."
            ));
        }
        if ui.button("Enregistrer le .msh").clicked() {
            todo = Some(Todo::Save);
        }
    }
    if let Some(note) = &mesh.note {
        ui.weak(note);
    }
    todo
}

fn format_name(cad: &Cad) -> &'static str {
    match cad.source {
        Source::Dxf => "dessin DXF",
        Source::Stl => "surface STL",
        Source::Step => "pièce STEP",
    }
}

/// Vrai si une géométrie est importée : les Sets se créent dessus, pas sur un `.msh`.
fn needs_geometry(ui: &mut egui::Ui, mesh: &MeshUi) -> bool {
    match (&mesh.imported, &mesh.built) {
        (Some(Imported::Cad(_)), Some(_)) => false,
        (Some(Imported::Mesh(..)), _) => {
            ui.label("Un .msh n'a pas de géométrie : désigne ses faces avec « Faces → Set » ou un cadre.");
            true
        }
        _ => {
            ui.label("Importe et maille d'abord une géométrie (outil « Importer »).");
            true
        }
    }
}

/// Inspecteur « Faces → Set » : la sélection faite à la souris, et le nom du Set à créer.
fn face_sets_tool(ui: &mut egui::Ui, mesh: &mut MeshUi) -> Option<Todo> {
    if matches!(mesh.imported, Some(Imported::Mesh(..))) {
        return region_tool(ui, mesh, Tool::FaceSets);
    }
    if needs_geometry(ui, mesh) {
        return None;
    }
    let mut todo = None;
    ui.label("Clique une face (3D) ou un bord (2D). Maj + clic : en ajouter ou en retirer.");
    ui.separator();
    if mesh.selection.is_empty() {
        ui.weak("Rien de sélectionné.");
    } else {
        ui.label(format!("Sélection : {}", mesh.selection.join(", ")));
        if ui.button("Effacer la sélection").clicked() {
            mesh.selection.clear();
        }
    }
    ui.separator();
    ui.label("Nom du Set :");
    ui.text_edit_singleline(&mut mesh.set_name);
    let ready = !mesh.selection.is_empty() && !mesh.set_name.trim().is_empty();
    if ui
        .add_enabled(ready, egui::Button::new("Créer le Set"))
        .clicked()
    {
        todo = Some(Todo::CreateSet);
    }
    if let Some(note) = &mesh.note {
        ui.weak(note);
    }
    todo
}

/// Inspecteurs des Sets tirés sur le Mesh : « Faces → Set » sur un `.msh` (région par angle) et
/// « Cadre → Set » (Nodes ou Elements dans un cadre). `tool` dit lequel.
fn region_tool(ui: &mut egui::Ui, mesh: &mut MeshUi, tool: Tool) -> Option<Todo> {
    if mesh.built.is_none() {
        ui.label("Importe et maille d'abord une géométrie, ou importe un .msh.");
        return None;
    }
    if tool == Tool::FaceSets {
        ui.label("Clique une face de la surface : la région s'étend aux faces voisines.");
        ui.label("Angle maximal entre voisines (degrés) :");
        ui.add(egui::Slider::new(&mut mesh.angle, 1.0..=90.0));
        ui.weak("Un plan se prend d'un clic ; un coin plus vif que l'angle arrête la région. Maj + clic : ajouter.");
    } else {
        ui.label("Glisse un cadre dans la vue. Oriente la vue d'abord avec l'outil flèche.");
        ui.horizontal(|ui| {
            ui.radio_value(&mut mesh.box_mode, BoxMode::Nodes, "Nodes");
            ui.radio_value(&mut mesh.box_mode, BoxMode::Elements, "Elements");
        });
        ui.checkbox(&mut mesh.box_add, "Ajouter à la sélection");
        ui.weak(
            "En 3D, le cadre prend aussi ce qui est caché derrière. Un Set d'Elements sert à \
             donner un Material à une partie de la pièce.",
        );
    }
    ui.separator();
    let mut todo = None;
    let current = mesh.region.as_ref().filter(|(t, _)| *t == tool);
    match current {
        None => {
            ui.weak("Rien de sélectionné.");
        }
        Some((_, region)) => {
            ui.label(format!("Sélection : {}", region.count()));
        }
    }
    let selected = current.is_some();
    if ui
        .add_enabled(selected, egui::Button::new("Effacer la sélection"))
        .clicked()
    {
        mesh.region = None;
        mesh.region_version += 1;
    }
    ui.separator();
    ui.label("Nom du Set :");
    ui.text_edit_singleline(&mut mesh.set_name);
    let ready = selected && !mesh.set_name.trim().is_empty();
    if ui
        .add_enabled(ready, egui::Button::new("Créer le Set"))
        .clicked()
    {
        todo = Some(Todo::CreateMeshSet);
    }
    ui.weak("Ce Set est écrit dans le .msh ; un remaillage le perd.");
    if let Some(note) = &mesh.note {
        ui.weak(note);
    }
    todo
}

/// Inspecteur « Sets » : les entités de la géométrie et les Sets créés ; cliquer en surligne un.
fn manage_sets_tool(ui: &mut egui::Ui, mesh: &mut MeshUi) -> Option<Todo> {
    let Some(built) = &mesh.built else {
        ui.label("Importe et maille d'abord une géométrie, ou importe un .msh.");
        return None;
    };
    let mut todo = None;
    let pick = |ui: &mut egui::Ui, name: &str, inspected: &mut Option<String>| {
        let selected = inspected.as_deref() == Some(name);
        if ui.selectable_label(selected, name).clicked() {
            *inspected = (!selected).then(|| name.to_string());
        }
    };
    ui.label("Entités de la géométrie (celles du fichier, pour un .msh) :");
    for name in &built.entities {
        pick(ui, name, &mut mesh.inspected);
    }
    ui.separator();
    ui.label("Sets créés :");
    if built.user.is_empty() && built.mesh_sets.is_empty() {
        ui.weak("Aucun : utilise l'outil « Faces → Set » ou « Cadre → Set ».");
    }
    for (name, members) in &built.user {
        pick(ui, name, &mut mesh.inspected);
        ui.weak(format!("    {}", members.join(", ")));
    }
    for name in &built.mesh_sets {
        pick(ui, name, &mut mesh.inspected);
        let set = &built.file.sets[name];
        ui.weak(format!(
            "    sur le maillage : {}",
            if !set.elements.is_empty() {
                format!("{} Elements", set.elements.len())
            } else if !set.faces.is_empty() {
                format!("{} faces", set.faces.len())
            } else {
                format!("{} Nodes", set.nodes.len())
            }
        ));
    }
    if let Some(name) = mesh
        .inspected
        .clone()
        .filter(|n| built.user.contains_key(n) || built.mesh_sets.contains(n))
    {
        ui.separator();
        ui.label(format!("« {name} » :"));
        ui.text_edit_singleline(&mut mesh.rename);
        ui.horizontal(|ui| {
            if ui.button("Renommer").clicked() {
                todo = Some(Todo::RenameSet(name.clone(), mesh.rename.clone()));
            }
            if ui.button("Supprimer").clicked() {
                todo = Some(Todo::DeleteSet(name.clone()));
            }
        });
    }
    todo
}

/// Inspecteur de la qualité : ce que le viewer colorie, en chiffres.
fn quality_tool(ui: &mut egui::Ui, mesh: &MeshUi) {
    let Some(built) = &mesh.built else {
        ui.label("Importe et maille d'abord une géométrie, ou importe un .msh.");
        return;
    };
    let q = &built.quality;
    ui.label(mesh_space::summary(built));
    ui.separator();
    ui.label(format!("Qualité minimale : {:.2}", q.min));
    ui.label(format!("Qualité moyenne : {:.2}", q.mean));
    if q.poor == 0 {
        ui.colored_label(
            Color32::from_rgb(90, 200, 120),
            format!("Aucun Element sous {POOR}."),
        );
    } else {
        ui.colored_label(
            Color32::from_rgb(240, 180, 70),
            format!(
                "{} Element(s) sous {POOR} : le calcul y sera moins précis, réduire la taille h.",
                q.poor
            ),
        );
    }
    ui.weak("1 : triangle équilatéral ou tétraèdre régulier. La vue colore le défaut de forme (1 − qualité) : rouge pour un Element écrasé, au pire Element de chaque Node.");
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
    painter.circle_filled(
        center,
        TRIAD_LENGTH + 14.0,
        ui.visuals().text_color().gamma_multiply(0.08),
    );
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
fn floating(ui: &egui::Ui) -> Frame {
    let dark = ui.visuals().dark_mode;
    Frame::new()
        .fill(if dark { PANEL_DARK } else { PANEL_LIGHT })
        .corner_radius(CornerRadius::same(14))
        .inner_margin(Margin::same(10))
        .stroke(Stroke::new(
            1.0,
            if dark {
                Color32::from_white_alpha(24)
            } else {
                Color32::from_black_alpha(36)
            },
        ))
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
                    open.refresh_marker();
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
            self.route(&path, frame);
        }
        let dropped = ctx.input(|i| {
            i.raw
                .dropped_files
                .iter()
                .map(|f| f.path().to_path_buf())
                .next()
        });
        if let Some(path) = dropped {
            self.route(&path, frame);
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

        let preference = match self.theme {
            ThemeChoice::System => egui::ThemePreference::System,
            ThemeChoice::Light => egui::ThemePreference::Light,
            ThemeChoice::Dark => egui::ThemePreference::Dark,
        };
        if ctx.options(|o| o.theme_preference) != preference {
            ctx.set_theme(preference);
        }
        self.poll_job(frame);
        self.refresh_model();
        let mut clicked = None;
        egui::CentralPanel::default()
            .frame(Frame::NONE)
            .show(ui, |ui| clicked = self.viewport(ui, frame));
        if let Some((cursor, add)) = clicked {
            self.pick_at(cursor, add);
            if self.state.space == Space::Results && self.state.tool() == Tool::Probe {
                if let Some(open) = self.results.as_mut() {
                    open.probe_at(cursor);
                    self.state.progress.probed =
                        open.play.as_ref().is_some_and(|p| p.probe.is_some());
                }
            }
        }
        if let Some(open) = self.results.as_mut() {
            if let Err(e) = open.tick(&ctx) {
                self.error = Some(e.to_string());
                if let Some(play) = open.play.as_mut() {
                    play.playing = false;
                }
            }
        }
        self.sync_highlight();
        self.spaces(&ctx);
        self.card(&ctx);
        self.toolbar(&ctx);
        self.inspector(&ctx, frame);
        self.compute_panel(&ctx);
        self.stale_banner(&ctx);
        self.confirm_dialog(&ctx);
        if self.refresh_model() {
            ctx.request_repaint();
        }
        let shown = match self.state.space {
            Space::Results => &self.results,
            _ => &self.mesh_view,
        };
        if let Some(open) = shown {
            ctx.send_viewport_cmd(egui::ViewportCommand::Title(format!(
                "Maille — {}",
                open.name
            )));
        }
    }
}
