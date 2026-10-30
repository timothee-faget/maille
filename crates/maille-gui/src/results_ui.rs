//! Les outils de l'espace Résultats qui parcourent un calcul : Sets (montrés sur la vue), Increments
//! (curseur et animation), Sonde (valeurs en un Node), Courbe XY (dans le temps, forces de contact) et
//! Capture PNG. Les calculs sont dans `series` (pur, testé) ; ici, seulement l'affichage et la vue.

use std::collections::BTreeSet;
use std::error::Error;
use std::fs;
use std::path::Path;
use std::time::{Duration, Instant};

use eframe::egui::{self, Color32, Pos2, Sense, Stroke, Vec2};
use maille_io::msh::MshFile;

use crate::app::Open;
use crate::mesh_pick::{self, SkinFace};
use crate::results;
use crate::series::{self, Series};

/// Couleurs des Sets, attribuées dans l'ordre des noms (la même couleur à chaque ouverture).
const SET_COLORS: [[u8; 3]; 8] = [
    [230, 25, 75],
    [60, 180, 75],
    [255, 225, 25],
    [0, 130, 200],
    [245, 130, 48],
    [145, 30, 180],
    [70, 240, 240],
    [240, 50, 230],
];

/// Rayon, en pixels, dans lequel un clic désigne un Node.
pub(crate) const PICK_RADIUS: f64 = 16.0;
const PLOT_HEIGHT: f32 = 150.0;

/// Une courbe tracée : ce qu'elle montre, pour quel Node, et ses points (temps, valeur).
pub(crate) struct Curve {
    pub label: String,
    pub node: Option<usize>,
    pub points: Vec<[f64; 2]>,
}

/// Ce qui change d'un Increment à l'autre dans la vue d'un calcul.
pub(crate) struct Playback {
    pub series: Series,
    pub index: usize,
    pub playing: bool,
    pub last_step: Instant,
    pub fps: f32,
    /// La même échelle de couleurs pour tous les Increments.
    pub fixed_range: bool,
    /// Rang du champ pour lequel l'échelle a été calculée.
    pub locked_quantity: Option<usize>,
    pub probe: Option<usize>,
    pub curve: Option<Curve>,
    pub note: Option<String>,
    /// Le `.msh` du calcul (outil « Sets ») : lu, ou son chemin (trouvé ou à taper) et, s'il n'a pas
    /// pu être lu, pourquoi.
    pub sets: Option<MeshSets>,
    pub mesh_path: String,
    pub mesh_note: Option<String>,
}

/// Les Sets du `.msh` du calcul et ceux qu'on montre. Les Nodes du fichier suivent la vue (déformée,
/// Increment affiché) : le surlignage se calcule dessus.
pub(crate) struct MeshSets {
    file: MshFile,
    skin: Vec<SkinFace>,
    shown: BTreeSet<String>,
}

impl Playback {
    pub(crate) fn new(series: Series) -> Self {
        let index = series.len() - 1;
        let fixed_range = series.len() > 1;
        // le `.msh` ne se lit qu'à la première ouverture de l'outil « Sets »
        let found = series.mesh_path();
        let mesh_note = found.is_none().then(|| {
            "Le modèle de ces résultats n'est pas à côté : tape le chemin du .msh.".to_string()
        });
        Self {
            series,
            index,
            playing: false,
            last_step: Instant::now(),
            fps: 6.0,
            fixed_range,
            locked_quantity: None,
            probe: None,
            curve: None,
            note: None,
            sets: None,
            mesh_path: found.map(|p| p.display().to_string()).unwrap_or_default(),
            mesh_note,
        }
    }
}

impl MeshSets {
    fn new(file: MshFile) -> Self {
        Self {
            skin: mesh_pick::skin(&file),
            file,
            shown: BTreeSet::new(),
        }
    }
}

impl Open {
    /// Montre l'Increment `k` (depuis 0).
    pub(crate) fn goto(&mut self, k: usize) -> Result<(), Box<dyn Error>> {
        let Some(play) = self.play.as_mut() else {
            return Ok(());
        };
        let nodal = results::to_nodal(play.series.get(k)?)?;
        play.index = k;
        self.view.set_results(&nodal)?;
        self.refresh_marker();
        self.refresh_sets();
        Ok(())
    }

    /// Les Sets cochés, chacun dans sa couleur, sur la pièce telle qu'elle est affichée.
    pub(crate) fn refresh_sets(&mut self) {
        let Some(sets) = self.play.as_mut().and_then(|p| p.sets.as_mut()) else {
            return;
        };
        if sets.shown.is_empty() {
            self.view.set_sets(Vec::new());
            return;
        }
        for (node, p) in sets.file.nodes.iter_mut().zip(self.view.node_positions()) {
            (node.x, node.y, node.z) = (p[0], p[1], p[2]);
        }
        let shown = sets
            .file
            .sets
            .iter()
            .zip(SET_COLORS.iter().cycle())
            .filter(|((name, _), _)| sets.shown.contains(*name))
            .map(|((_, set), &color)| {
                (color, mesh_pick::highlight_set(&sets.file, &sets.skin, set))
            })
            .collect();
        self.view.set_sets(shown);
    }

    /// Lit le `.msh` tapé dans l'outil « Sets ».
    fn load_mesh(&mut self) {
        let Some(play) = self.play.as_mut() else {
            return;
        };
        let path = Path::new(play.mesh_path.trim()).to_path_buf();
        match play.series.read_mesh(&path) {
            Ok(file) => {
                play.sets = Some(MeshSets::new(file));
                play.mesh_note = None;
            }
            Err(e) => play.mesh_note = Some(e.to_string()),
        }
        self.refresh_sets();
    }

    /// Fixe (ou libère) l'intervalle de la légende : celui du champ affiché sur tous les Increments.
    fn lock_range(&mut self) -> Result<(), Box<dyn Error>> {
        let Some(play) = self.play.as_mut() else {
            return Ok(());
        };
        if !play.fixed_range || play.series.len() < 2 {
            self.view.set_range_lock(None);
            play.locked_quantity = None;
            return Ok(());
        }
        let quantity = self.view.quantity();
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for k in 0..play.series.len() {
            let nodal = results::to_nodal(play.series.get(k)?)?;
            let (a, b) = maille_viewer::quantity_range(&nodal, quantity)?;
            lo = lo.min(a);
            hi = hi.max(b);
        }
        self.view.set_range_lock(Some((lo, hi)));
        play.locked_quantity = Some(quantity);
        Ok(())
    }

    /// Un petit repère sur le Node de la sonde (trois triangles, un par plan).
    pub(crate) fn refresh_marker(&mut self) {
        let Some(node) = self.play.as_ref().and_then(|p| p.probe) else {
            self.view.set_highlight(Vec::new());
            return;
        };
        let positions = self.view.node_positions();
        let Some(&p) = positions.get(node) else {
            return;
        };
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for q in &positions {
            for k in 0..3 {
                lo[k] = lo[k].min(q[k]);
                hi[k] = hi[k].max(q[k]);
            }
        }
        let size = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max) * 0.012;
        let at = |a: usize, b: usize, (da, db): (f64, f64)| {
            let mut v = p;
            v[a] += da * size;
            v[b] += db * size;
            v
        };
        let tri = |a: usize, b: usize| {
            [
                at(a, b, (-1.0, -1.0)),
                at(a, b, (1.0, -1.0)),
                at(a, b, (0.0, 1.5)),
            ]
        };
        self.view
            .set_highlight(vec![tri(0, 1), tri(0, 2), tri(1, 2)]);
    }

    /// Une image de plus de l'animation, et l'échelle de couleurs si elle a besoin d'être refaite.
    pub(crate) fn tick(&mut self, ctx: &egui::Context) -> Result<(), Box<dyn Error>> {
        let Some(play) = self.play.as_ref() else {
            return Ok(());
        };
        if play.fixed_range && play.locked_quantity != Some(self.view.quantity()) {
            self.lock_range()?;
        }
        let Some(play) = self.play.as_mut() else {
            return Ok(());
        };
        if play.playing {
            let step = Duration::from_secs_f32(1.0 / play.fps.max(1.0));
            if play.last_step.elapsed() >= step {
                let next = (play.index + 1) % play.series.len();
                play.last_step = Instant::now();
                self.goto(next)?;
            }
            ctx.request_repaint_after(Duration::from_millis(16));
        }
        Ok(())
    }

    /// Désigne le Node sous le curseur (pixels de la vue) comme sonde.
    pub(crate) fn probe_at(&mut self, cursor: [f64; 2]) {
        if self.play.is_none() {
            return;
        }
        let projected: Vec<[f64; 3]> = self
            .view
            .node_positions()
            .into_iter()
            .map(|p| self.view.project(p))
            .collect();
        let node = series::nearest_node(&projected, cursor, PICK_RADIUS);
        if let Some(play) = self.play.as_mut() {
            play.probe = node;
            play.curve = None;
        }
        self.refresh_marker();
    }
}

/// Inspecteur « Sets » : une case par Set du `.msh` du calcul, dans sa couleur.
pub(crate) fn sets_tool(ui: &mut egui::Ui, open: &mut Open) {
    if open
        .play
        .as_ref()
        .is_some_and(|p| p.sets.is_none() && p.mesh_note.is_none())
    {
        open.load_mesh();
    }
    let Some(play) = open.play.as_mut() else {
        return;
    };
    let mut changed = false;
    match play.sets.as_mut() {
        Some(sets) if sets.file.sets.is_empty() => {
            ui.label("Le maillage de ce calcul n'a pas de Set.");
        }
        Some(sets) => {
            ui.horizontal(|ui| {
                if ui.button("Tous").clicked() {
                    sets.shown = sets.file.sets.keys().cloned().collect();
                    changed = true;
                }
                if ui.button("Aucun").clicked() {
                    sets.shown.clear();
                    changed = true;
                }
            });
            for ((name, set), &[r, g, b]) in sets.file.sets.iter().zip(SET_COLORS.iter().cycle()) {
                ui.horizontal(|ui| {
                    let (rect, _) = ui.allocate_exact_size(Vec2::splat(12.0), Sense::hover());
                    ui.painter()
                        .rect_filled(rect, 2.0, Color32::from_rgb(r, g, b));
                    let mut on = sets.shown.contains(name);
                    if ui.checkbox(&mut on, name).changed() {
                        if on {
                            sets.shown.insert(name.clone());
                        } else {
                            sets.shown.remove(name);
                        }
                        changed = true;
                    }
                    ui.weak(mesh_pick::describe(set));
                });
            }
        }
        None => {
            if let Some(note) = &play.mesh_note {
                ui.label(note);
            }
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut play.mesh_path);
                if ui.button("Lire").clicked() {
                    changed = true;
                }
            });
            if changed {
                open.load_mesh();
            }
            return;
        }
    }
    if changed {
        open.refresh_sets();
    }
}

/// Inspecteur « Increments » : le curseur, la lecture et l'échelle de couleurs.
pub(crate) fn increments_tool(ui: &mut egui::Ui, open: &mut Open) -> Result<(), Box<dyn Error>> {
    let Some(play) = open.play.as_mut() else {
        return Ok(());
    };
    let n = play.series.len();
    if n < 2 {
        ui.label(format!(
            "Ce calcul n'a qu'un incrément (t = {}).",
            play.series.time(0)
        ));
        return Ok(());
    }
    let mut k = play.index;
    ui.add(egui::Slider::new(&mut k, 0..=n - 1).text("incrément"));
    ui.label(format!("t = {}", play.series.time(k)));
    ui.horizontal(|ui| {
        if ui.button("<<").clicked() {
            k = 0;
        }
        let label = if play.playing { "Pause" } else { "Lecture" };
        if ui.button(label).clicked() {
            play.playing = !play.playing;
            play.last_step = Instant::now();
        }
        if ui.button(">>").clicked() {
            k = n - 1;
        }
    });
    ui.add(egui::Slider::new(&mut play.fps, 1.0..=30.0).text("images/s"));
    let before = play.fixed_range;
    ui.checkbox(&mut play.fixed_range, "Même échelle de couleurs")
        .on_hover_text("La légende couvre tous les incréments, pour comparer les images");
    let changed = play.fixed_range != before;
    let (index, playing) = (play.index, play.playing);
    if changed {
        open.lock_range()?;
    }
    if k != index {
        if let Some(play) = open.play.as_mut() {
            play.playing = playing && k != n - 1;
        }
        open.goto(k)?;
    }
    Ok(())
}

/// Inspecteur « Sonde » : les valeurs de l'Increment affiché au Node cliqué.
pub(crate) fn probe_tool(ui: &mut egui::Ui, open: &mut Open) -> Result<(), Box<dyn Error>> {
    let Some(play) = open.play.as_mut() else {
        return Ok(());
    };
    ui.label("Clique un Node dans la vue.");
    let Some(node) = play.probe else {
        return Ok(());
    };
    let index = play.index;
    let time = play.series.time(index);
    let r = play.series.get(index)?;
    let p = &r.nodes[node];
    ui.separator();
    ui.label(format!("Node {node} : ({}, {}, {})", p.x, p.y, p.z));
    ui.weak(format!("incrément {index}, t = {time}"));
    egui::Grid::new("probe").num_columns(2).show(ui, |ui| {
        for (label, value) in series::readings(r, node) {
            ui.label(label);
            ui.monospace(format!("{value:.6e}"));
            ui.end_row();
        }
    });
    Ok(())
}

/// Inspecteur « Courbe XY » : une grandeur du Node de la sonde au fil du temps, ou la force sur une
/// surface rigide.
pub(crate) fn curve_tool(ui: &mut egui::Ui, open: &mut Open) -> Result<(), Box<dyn Error>> {
    let Some(play) = open.play.as_mut() else {
        return Ok(());
    };
    let mut choices: Vec<(String, Option<usize>)> = Vec::new();
    if let Some(node) = play.probe {
        let labels = series::readings(play.series.get(play.index)?, node);
        choices.extend(labels.into_iter().map(|(l, _)| (l, Some(node))));
    }
    let contacts = match play.series.contacts_csv() {
        Some(path) => series::contact_curves(
            &fs::read_to_string(&path)
                .map_err(|e| format!("impossible de lire {} : {e}", path.display()))?,
        )?,
        None => Vec::new(),
    };
    choices.extend(contacts.iter().map(|(l, _)| (l.clone(), None)));
    if choices.is_empty() {
        ui.label("Choisis d'abord un Node avec l'outil « Sonde ».");
        return Ok(());
    }
    if play.probe.is_some() {
        ui.weak(format!("Node de la sonde : {}", play.probe.unwrap_or(0)));
    }
    let current = play.curve.as_ref().map(|c| (c.label.clone(), c.node));
    let mut wanted = current.clone();
    egui::ComboBox::from_id_salt("curve")
        .selected_text(wanted.as_ref().map_or("Choisir…", |(l, _)| l.as_str()))
        .width(200.0)
        .show_ui(ui, |ui| {
            for (label, node) in &choices {
                let selected = wanted.as_ref() == Some(&(label.clone(), *node));
                if ui.selectable_label(selected, label).clicked() {
                    wanted = Some((label.clone(), *node));
                }
            }
        });
    if wanted != current {
        if let Some((label, node)) = wanted {
            let points = match node {
                Some(node) => play.series.curve(node, &label)?,
                None => contacts
                    .iter()
                    .find(|(l, _)| *l == label)
                    .map(|(_, p)| p.clone())
                    .unwrap_or_default(),
            };
            play.curve = Some(Curve {
                label,
                node,
                points,
            });
        }
    }
    if let Some(curve) = &play.curve {
        ui.label(&curve.label);
        plot(ui, &curve.points, Some(play.series.time(play.index)));
    }
    Ok(())
}

/// Trace une courbe (temps, valeur) avec ses bornes ; `marker` : le temps de l'Increment affiché.
fn plot(ui: &mut egui::Ui, points: &[[f64; 2]], marker: Option<f64>) {
    let width = ui.available_width();
    let (response, painter) = ui.allocate_painter(Vec2::new(width, PLOT_HEIGHT), Sense::hover());
    let rect = response.rect;
    let dark = ui.visuals().dark_mode;
    painter.rect_filled(
        rect,
        3.0,
        Color32::from_black_alpha(if dark { 90 } else { 24 }),
    );
    if points.len() < 2 {
        return;
    }
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for p in points {
        x0 = x0.min(p[0]);
        x1 = x1.max(p[0]);
        y0 = y0.min(p[1]);
        y1 = y1.max(p[1]);
    }
    if (x1 - x0).abs() < f64::EPSILON {
        x1 = x0 + 1.0;
    }
    if (y1 - y0).abs() < f64::EPSILON {
        y0 -= 0.5;
        y1 += 0.5;
    }
    let inner = rect.shrink2(Vec2::new(8.0, 14.0));
    let at = |p: [f64; 2]| {
        Pos2::new(
            inner.min.x + ((p[0] - x0) / (x1 - x0)) as f32 * inner.width(),
            inner.max.y - ((p[1] - y0) / (y1 - y0)) as f32 * inner.height(),
        )
    };
    if let Some(t) = marker.filter(|t| (x0..=x1).contains(t)) {
        let x = at([t, y0]).x;
        painter.line_segment(
            [Pos2::new(x, inner.min.y), Pos2::new(x, inner.max.y)],
            Stroke::new(1.0, Color32::from_gray(150)),
        );
    }
    painter.add(egui::Shape::line(
        points.iter().map(|p| at(*p)).collect(),
        Stroke::new(
            1.6,
            if dark {
                Color32::from_rgb(120, 190, 255)
            } else {
                Color32::from_rgb(20, 90, 190)
            },
        ),
    ));
    let font = egui::FontId::monospace(10.0);
    let weak = ui.visuals().weak_text_color();
    let corner = |p: Pos2, anchor, text: String| {
        painter.text(p, anchor, text, font.clone(), weak);
    };
    corner(
        rect.left_top() + Vec2::new(4.0, 1.0),
        egui::Align2::LEFT_TOP,
        format!("{y1:.3e}"),
    );
    corner(
        rect.left_bottom() + Vec2::new(4.0, -1.0),
        egui::Align2::LEFT_BOTTOM,
        format!("{y0:.3e}"),
    );
    corner(
        rect.right_bottom() + Vec2::new(-4.0, -1.0),
        egui::Align2::RIGHT_BOTTOM,
        format!("t = {x1}"),
    );
}

/// Inspecteur « Capture PNG » : l'image de la vue, à côté du `.pvd`.
pub(crate) fn capture_tool(ui: &mut egui::Ui, open: &mut Open) -> Result<(), Box<dyn Error>> {
    ui.label("Enregistre l'image de la vue, légende comprise.");
    if ui.button("Capturer").clicked() {
        let png = open.view.png()?;
        if let Some(play) = open.play.as_mut() {
            let stem = play
                .series
                .path()
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_default();
            let path = play
                .series
                .path()
                .with_file_name(format!("{stem}_{:04}.png", play.index));
            fs::write(&path, png)
                .map_err(|e| format!("impossible d'écrire {} : {e}", path.display()))?;
            play.note = Some(format!("Écrit : {}", path.display()));
        }
    }
    if let Some(note) = open.play.as_ref().and_then(|p| p.note.as_deref()) {
        ui.weak(note);
    }
    Ok(())
}
