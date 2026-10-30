//! Les inspecteurs de l'espace Modèle : des formulaires qui éditent le `Draft` de `model_space`. Ils ne
//! décident de rien d'autre que de ce qu'on affiche ; les règles sont dans `model_space`.

use std::path::PathBuf;

use eframe::egui::{self, Color32};
use maille_io::materials::{Library, LibraryMaterial};
use maille_io::model::{
    ContactDef, ContactKind, ConvectionDef, DofDef, EdgeLoadDef, HeatSourceDef, PlaneStateDef,
    PressureDef, ValueDef,
};
use maille_io::msh::CellKind;

use crate::compute::ComputeState;
use crate::library;
use crate::model_space::{Diagnostic, Draft, MaterialForm, Region, Severity, UnitSystem};
use crate::state::Tool;

/// Un formulaire « Set, DOF, valeur » pour ajouter une BC ou un Load.
#[derive(Debug, Clone)]
pub struct ValueForm {
    pub set: String,
    pub dof: DofDef,
    pub value: f64,
}

impl Default for ValueForm {
    fn default() -> Self {
        Self {
            set: String::new(),
            dof: DofDef::Ux,
            value: 0.0,
        }
    }
}

/// Les valeurs tapées dans les inspecteurs du Modèle (pas encore ajoutées au brouillon).
#[derive(Debug, Clone, Default)]
pub struct Forms {
    pub bc: ValueForm,
    pub load: ValueForm,
    pub pressure_set: String,
    pub pressure: f64,
    pub edge_set: String,
    /// Force linéique (fx, fy, fz) ; fz ne sert qu'en 3D.
    pub edge_load: [f64; 3],
    /// Température imposée : `set` et `value` de la BC en `t`.
    pub temperature: ValueForm,
    pub heat_set: String,
    pub heat: f64,
    pub convection_set: String,
    pub film: f64,
    pub sink: f64,
    pub contact: ContactForm,
}

/// Le formulaire d'un Contact : tous les paramètres des trois surfaces, on n'en garde que ceux du `kind`.
#[derive(Debug, Clone)]
pub struct ContactForm {
    pub set: String,
    pub kind: ContactKind,
    pub point: [f64; 3],
    pub normal: [f64; 3],
    pub axis: [f64; 3],
    pub center: [f64; 3],
    pub radius: f64,
    pub motion: [f64; 3],
}

impl Default for ContactForm {
    fn default() -> Self {
        Self {
            set: String::new(),
            kind: ContactKind::Plane,
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            center: [0.0; 3],
            radius: 1.0,
            motion: [0.0; 3],
        }
    }
}

impl ContactForm {
    /// Le Contact du formulaire pour un maillage de dimension `dim` (2 ou 3).
    fn to_def(&self, dim: usize) -> ContactDef {
        let v = |a: &[f64; 3]| Some(a[..dim].to_vec());
        let (plane, cylinder, sphere) = (
            self.kind == ContactKind::Plane,
            self.kind == ContactKind::Cylinder,
            self.kind == ContactKind::Sphere,
        );
        ContactDef {
            set: self.set.clone(),
            kind: self.kind,
            point: (plane || cylinder).then(|| self.point[..dim].to_vec()),
            normal: if plane { v(&self.normal) } else { None },
            axis: if cylinder && dim == 3 {
                v(&self.axis)
            } else {
                None
            },
            center: if sphere { v(&self.center) } else { None },
            radius: (cylinder || sphere).then_some(self.radius),
            motion: self
                .motion
                .iter()
                .any(|m| *m != 0.0)
                .then(|| self.motion[..dim].to_vec()),
        }
    }
}

/// Ce que le maillage dit au formulaire.
pub struct MeshInfo<'a> {
    /// Sets de bord du maillage.
    pub sets: &'a [&'a str],
    /// Sets d'Elements (le domaine compris) : ce qu'un Material peut recevoir.
    pub element_sets: &'a [&'a str],
    /// Sets qu'une force linéique peut viser : des bords en 2D, des arêtes en 3D.
    pub edge_sets: &'a [&'a str],
    pub kind: CellKind,
}

impl MeshInfo<'_> {
    fn dim(&self) -> usize {
        if self.planar() {
            2
        } else {
            3
        }
    }

    fn planar(&self) -> bool {
        matches!(self.kind, CellKind::Tri3 | CellKind::Tri6)
    }

    /// Les DOF de déplacement que le maillage accepte.
    fn dofs(&self) -> &'static [DofDef] {
        if self.planar() {
            &[DofDef::Ux, DofDef::Uy]
        } else {
            &[DofDef::Ux, DofDef::Uy, DofDef::Uz]
        }
    }
}

fn dof_label(dof: DofDef) -> &'static str {
    match dof {
        DofDef::Ux => "ux",
        DofDef::Uy => "uy",
        DofDef::Uz => "uz",
        DofDef::T => "t",
    }
}

fn with_unit(text: &str, unit: &str) -> String {
    if unit.is_empty() {
        text.to_string()
    } else {
        format!("{text} [{unit}]")
    }
}

/// Liste déroulante des Sets ; garde le choix précédent s'il existe encore.
fn set_picker(ui: &mut egui::Ui, id: &str, sets: &[&str], chosen: &mut String) {
    if !sets.contains(&chosen.as_str()) {
        *chosen = sets.first().map(|s| s.to_string()).unwrap_or_default();
    }
    egui::ComboBox::from_id_salt(id)
        .selected_text(chosen.as_str())
        .show_ui(ui, |ui| {
            for &set in sets {
                ui.selectable_value(chosen, set.to_string(), set);
            }
        });
}

fn dof_picker(ui: &mut egui::Ui, id: &str, dofs: &[DofDef], chosen: &mut DofDef) {
    if !dofs.contains(chosen) {
        *chosen = dofs[0];
    }
    egui::ComboBox::from_id_salt(id)
        .selected_text(dof_label(*chosen))
        .show_ui(ui, |ui| {
            for &dof in dofs {
                ui.selectable_value(chosen, dof, dof_label(dof));
            }
        });
}

/// La bibliothèque de Materials affichée dans l'inspecteur « Material ».
pub struct LibraryUi {
    pub library: Library,
    /// Le fichier de l'utilisateur, s'il a un dossier de configuration.
    pub path: Option<PathBuf>,
    /// Dernier message (fichier illisible, enregistrement fait ou refusé).
    pub note: Option<String>,
    /// Le Material choisi dans la liste.
    pub chosen: String,
    /// La limite d'élasticité du dernier Material pris, dans les unités du brouillon : proposée quand
    /// l'utilisateur coche la plasticité.
    pub yield_hint: Option<f64>,
}

impl LibraryUi {
    pub fn new(path: Option<PathBuf>) -> Self {
        let (library, note) = library::load(path.as_deref());
        Self {
            library,
            path,
            note,
            chosen: String::new(),
            yield_hint: None,
        }
    }
}

/// Inspecteur « Material » : unités, Material élastique, état plan et épaisseur.
pub fn material_tool(
    ui: &mut egui::Ui,
    draft: &mut Draft,
    mesh: &MeshInfo,
    lib: &mut LibraryUi,
) -> Option<MaterialAction> {
    ui.label("Unités (cohérentes, rien n'est converti) :");
    let mut system = match &draft.units {
        UnitSystem::Si => 0,
        UnitSystem::MmNMpa => 1,
        UnitSystem::Custom(_) => 2,
    };
    let before = system;
    ui.horizontal(|ui| {
        ui.radio_value(&mut system, 0, "N, m, Pa");
        ui.radio_value(&mut system, 1, "N, mm, MPa");
        ui.radio_value(&mut system, 2, "libre");
    });
    if system != before {
        draft.units = match system {
            0 => UnitSystem::Si,
            1 => UnitSystem::MmNMpa,
            _ => UnitSystem::Custom(String::new()),
        };
    }
    if let UnitSystem::Custom(text) = &mut draft.units {
        ui.text_edit_singleline(text);
        ui.weak("Aucune unité n'est affichée.");
    }
    ui.separator();
    let units = draft.units.clone();
    let library_units = library::library_units(&units);
    let mut action = None;
    if draft.material.is_none() {
        if ui.button("Définir un Material").clicked() {
            draft.material = Some(draft.steel());
        }
        ui.weak("Un acier pour commencer : règle ensuite E et ν.");
        ui.label("Ou prends-en un dans la bibliothèque :");
        if let Some(form) = library_picker(ui, "library-first", library_units, lib) {
            draft.material = Some(form);
        }
    } else {
        let large = draft.large;
        let several = !draft.regions.is_empty();
        let mut remove = None;
        // le premier Material, puis les autres
        egui::CollapsingHeader::new(format!(
            "{} · {}",
            draft.material.as_ref().map_or("", |m| &m.name),
            draft.material_set
        ))
        .id_salt("material-0")
        .default_open(true)
        .show(ui, |ui| {
            if several || draft.material_set != maille_io::msh::DOMAIN_SET {
                ui.label("Set d'Elements :");
                set_picker(
                    ui,
                    "material-0-set",
                    mesh.element_sets,
                    &mut draft.material_set,
                );
                if several && ui.button("Elements restants").clicked() {
                    action = Some(MaterialAction::RemainingSet);
                }
                ui.weak("« Elements restants » crée le Set de ce qui n'a pas d'autre Material.");
            }
            if let Some(m) = &mut draft.material {
                material_editor(ui, m, &units, large, library_units, lib, "library-0");
                ui.separator();
                if !several && ui.button("Retirer le Material").clicked() {
                    draft.material = None;
                }
            }
        });
        for (i, region) in draft.regions.iter_mut().enumerate() {
            egui::CollapsingHeader::new(format!("{} · {}", region.material.name, region.set))
                .id_salt(format!("material-{}", i + 1))
                .default_open(true)
                .show(ui, |ui| {
                    ui.label("Set d'Elements :");
                    set_picker(
                        ui,
                        &format!("material-{}-set", i + 1),
                        mesh.element_sets,
                        &mut region.set,
                    );
                    material_editor(
                        ui,
                        &mut region.material,
                        &units,
                        large,
                        library_units,
                        lib,
                        &format!("library-{}", i + 1),
                    );
                    ui.separator();
                    if ui.button("Retirer ce Material").clicked() {
                        remove = Some(i);
                    }
                });
        }
        if let Some(i) = remove {
            draft.regions.remove(i);
        }
        ui.separator();
        if ui.button("Ajouter un Material sur un autre Set").clicked() {
            let mut material = draft.steel();
            material.name = format!("materiau{}", draft.regions.len() + 2);
            let used: Vec<&str> = draft.sections().iter().map(|(set, _)| *set).collect();
            let set = mesh
                .element_sets
                .iter()
                .find(|s| !used.contains(s))
                .map(|s| s.to_string())
                .unwrap_or_default();
            draft.regions.push(Region { set, material });
        }
        ui.weak(
            "Chaque Element a un seul Material : les Sets d'Elements se créent avec « Cadre → Set » \
             (mode Elements).",
        );
        ui.separator();
        ui.checkbox(&mut draft.large, "Grandes déformations");
    }
    if mesh.planar() {
        ui.separator();
        ui.label("Triangles :");
        ui.horizontal(|ui| {
            ui.radio_value(
                &mut draft.plane_state,
                PlaneStateDef::Stress,
                "contraintes planes",
            );
            ui.radio_value(
                &mut draft.plane_state,
                PlaneStateDef::Strain,
                "déformations planes",
            );
        });
        ui.label(with_unit("Épaisseur", units.length()));
        let speed = (draft.thickness.abs() * 0.01).max(1e-9);
        ui.add(
            egui::DragValue::new(&mut draft.thickness)
                .speed(speed)
                .range(0.0..=f64::MAX),
        );
    }
    action
}

/// Ce que l'inspecteur « Material » demande à l'application, qui détient le Mesh.
pub enum MaterialAction {
    /// Créer le Set des Elements que les autres Materials ne couvrent pas, et le donner au premier.
    RemainingSet,
}

/// Le formulaire d'un Material : nom, E, ν, propriétés facultatives, bibliothèque, enregistrement.
fn material_editor(
    ui: &mut egui::Ui,
    m: &mut MaterialForm,
    units: &UnitSystem,
    large: bool,
    library_units: Option<maille_io::materials::UnitSystem>,
    lib: &mut LibraryUi,
    id: &str,
) {
    ui.label("Nom :");
    ui.text_edit_singleline(&mut m.name);
    ui.label(with_unit("Module d'Young E", units.stress()));
    let speed = (m.young.abs() * 0.01).max(1e-9);
    ui.add(
        egui::DragValue::new(&mut m.young)
            .speed(speed)
            .range(0.0..=f64::MAX),
    );
    ui.label("Coefficient de Poisson ν :");
    ui.add(egui::DragValue::new(&mut m.poisson).speed(0.005));
    let k = if *units == UnitSystem::Si { 1e6 } else { 1.0 };
    ui.separator();
    ui.label("Propriétés facultatives :");
    optional(ui, "conductivité k", &mut m.conductivity, 1.0);
    optional(ui, "masse volumique ρ", &mut m.density, 1.0);
    optional(ui, "chaleur massique c", &mut m.specific_heat, 1.0);
    optional(ui, "dilatation α [1/K]", &mut m.expansion, 1.2e-5);
    optional(
        ui,
        "limite élastique σy (plasticité)",
        &mut m.yield_stress,
        lib.yield_hint.unwrap_or(250.0 * k),
    );
    if m.yield_stress.is_some() {
        let mut h = m.hardening.unwrap_or(0.0);
        ui.horizontal(|ui| {
            ui.label(with_unit("écrouissage H", units.stress()));
            ui.add(
                egui::DragValue::new(&mut h)
                    .speed(1.0)
                    .range(0.0..=f64::MAX),
            );
        });
        m.hardening = Some(h);
    }
    if large {
        ui.checkbox(&mut m.neo_hookean, "loi néo-hookéenne");
    }
    ui.separator();
    ui.label("Remplacer par la bibliothèque :");
    if let Some(form) = library_picker(ui, id, library_units, lib) {
        *m = form;
    }
    if let (Some(path), Some(lu)) = (&lib.path, library_units) {
        if ui.button("Enregistrer dans ma bibliothèque").clicked() {
            let saved = LibraryMaterial::from_def(&m.to_def(large), lu, "");
            lib.note = Some(match library::save(path, &m.name, saved) {
                Ok((library, note)) => {
                    lib.library = library;
                    note.unwrap_or_else(|| format!("« {} » est dans {}.", m.name, path.display()))
                }
                Err(e) => format!("Enregistrement impossible : {e}"),
            });
        }
    }
}

/// Choisir un Material de la bibliothèque : rend son formulaire, converti dans les unités du
/// brouillon, quand on le prend.
fn library_picker(
    ui: &mut egui::Ui,
    id: &str,
    units: Option<maille_io::materials::UnitSystem>,
    lib: &mut LibraryUi,
) -> Option<MaterialForm> {
    let Some(units) = units else {
        ui.weak("Pour la bibliothèque, écris les unités « longueur, force, temps », par exemple « mm, N, s ».");
        return None;
    };
    egui::ComboBox::from_id_salt(id)
        .selected_text(if lib.chosen.is_empty() {
            "choisir…"
        } else {
            &lib.chosen
        })
        .show_ui(ui, |ui| {
            for name in lib.library.names() {
                ui.selectable_value(&mut lib.chosen, name.to_string(), name);
            }
        });
    let mut taken = None;
    if let Some(material) = lib.library.get(&lib.chosen) {
        if !material.description.is_empty() {
            ui.weak(&material.description);
        }
        if ui.button("Prendre ce Material").clicked() {
            let (form, yield_stress) = library::material_form(&lib.chosen, material, units);
            lib.yield_hint = yield_stress;
            taken = Some(form);
        }
        ui.weak("Valeurs usuelles : vérifie celles de ta nuance.");
    }
    if let Some(note) = &lib.note {
        ui.weak(note);
    }
    taken
}

/// Inspecteur « Encastrement » : les BC (déplacement imposé sur un Set ; 0 = bloqué).
pub fn support_tool(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo, form: &mut ValueForm) {
    ui.label("Un déplacement imposé sur chaque Node d'un Set ; 0 bloque la direction.");
    ui.separator();
    if draft.bcs.is_empty() {
        ui.weak("Aucune BC.");
    }
    let mut remove = None;
    for (i, bc) in draft.bcs.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} = {}", bc.set, dof_label(bc.dof), bc.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        draft.bcs.remove(i);
    }
    ui.separator();
    if mesh.sets.is_empty() {
        ui.weak("Maille d'abord, et crée des Sets (espace Maillage).");
        return;
    }
    ui.label("Ajouter :");
    ui.horizontal(|ui| {
        set_picker(ui, "bc-set", mesh.sets, &mut form.set);
        dof_picker(ui, "bc-dof", mesh.dofs(), &mut form.dof);
    });
    ui.horizontal(|ui| {
        ui.label(with_unit("valeur", draft.units.length()));
        ui.add(egui::DragValue::new(&mut form.value).speed(0.01));
    });
    ui.horizontal(|ui| {
        if ui.button("Ajouter la BC").clicked() && !form.set.is_empty() {
            draft.bcs.push(ValueDef {
                set: form.set.clone(),
                dof: form.dof,
                value: form.value,
            });
        }
        if ui.button("Encastrer ce Set").clicked() && !form.set.is_empty() {
            for &dof in mesh.dofs() {
                draft.bcs.push(ValueDef {
                    set: form.set.clone(),
                    dof,
                    value: 0.0,
                });
            }
        }
    });
}

/// Inspecteur « Load, Pressure » : efforts par Node et pressions sur un bord.
pub fn load_tool(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo, forms: &mut Forms) {
    if mesh.sets.is_empty() {
        ui.weak("Maille d'abord, et crée des Sets (espace Maillage).");
        return;
    }
    let units = draft.units.clone();
    ui.label(with_unit(
        "Pressure sur un bord (positive : elle comprime)",
        units.stress(),
    ));
    let mut remove = None;
    for (i, p) in draft.pressures.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {}", p.set, p.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove.take() {
        draft.pressures.remove(i);
    }
    ui.horizontal(|ui| {
        set_picker(ui, "pressure-set", mesh.sets, &mut forms.pressure_set);
        ui.add(egui::DragValue::new(&mut forms.pressure).speed(0.1));
    });
    if ui.button("Ajouter la Pressure").clicked() && !forms.pressure_set.is_empty() {
        draft.pressures.push(PressureDef {
            set: forms.pressure_set.clone(),
            value: forms.pressure,
        });
    }
    ui.separator();
    edge_load_form(ui, draft, mesh, forms);
    ui.separator();
    ui.label(with_unit("Load : effort par Node d'un Set", units.force()));
    ui.weak("C'est une force par Node, pas la force totale.");
    for (i, l) in draft.loads.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {} = {}", l.set, dof_label(l.dof), l.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        draft.loads.remove(i);
    }
    let form = &mut forms.load;
    ui.horizontal(|ui| {
        set_picker(ui, "load-set", mesh.sets, &mut form.set);
        dof_picker(ui, "load-dof", mesh.dofs(), &mut form.dof);
        ui.add(egui::DragValue::new(&mut form.value).speed(0.1));
    });
    if ui.button("Ajouter le Load").clicked() && !form.set.is_empty() {
        draft.loads.push(ValueDef {
            set: form.set.clone(),
            dof: form.dof,
            value: form.value,
        });
    }
}

/// Force linéique sur un bord (2D) ou des arêtes (3D) : composantes fixes, par unité de longueur.
fn edge_load_form(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo, forms: &mut Forms) {
    let unit = match (draft.units.force(), draft.units.length()) {
        ("", _) | (_, "") => String::new(),
        (force, length) => format!("{force}/{length}"),
    };
    let dim = mesh.dim();
    ui.label(with_unit("Force linéique sur un bord ou des arêtes", &unit));
    ui.weak(if dim == 2 {
        "Force par longueur de bord, sur toute l'épaisseur : σ · épaisseur pour une contrainte σ."
    } else {
        "Force par longueur d'arête, composantes fixes. Les Sets d'arêtes se tirent avec « Arêtes → Set »."
    });
    let mut remove = None;
    for (i, l) in draft.edge_loads.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {:?}", l.set, l.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        draft.edge_loads.remove(i);
    }
    if mesh.edge_sets.is_empty() {
        ui.weak(if dim == 2 {
            "Aucun Set de bord."
        } else {
            "Aucun Set d'arêtes : trace-en un avec l'outil « Arêtes → Set » (espace Maillage)."
        });
        return;
    }
    ui.horizontal(|ui| {
        set_picker(ui, "edge-load-set", mesh.edge_sets, &mut forms.edge_set);
        for (axis, v) in ["x", "y", "z"]
            .iter()
            .zip(forms.edge_load.iter_mut())
            .take(dim)
        {
            ui.label(*axis);
            ui.add(egui::DragValue::new(v).speed(0.1));
        }
    });
    if ui.button("Ajouter la force linéique").clicked() && !forms.edge_set.is_empty() {
        draft.edge_loads.push(EdgeLoadDef {
            set: forms.edge_set.clone(),
            value: forms.edge_load[..dim].to_vec(),
        });
    }
}

/// Ce que l'inspecteur « Calculer » demande à l'application.
#[derive(Debug, Default)]
pub struct RunRequest {
    pub save: bool,
    pub run: bool,
    pub goto: Option<Tool>,
}

/// Inspecteur « Calculer » : le calcul déduit, les diagnostics, le lancement du calcul (en tâche de
/// fond, progression dans le panneau flottant) et l'écriture du modèle.
pub fn run_tool(
    ui: &mut egui::Ui,
    calculation: Option<String>,
    diags: &[Diagnostic],
    note: Option<&str>,
    compute: &ComputeState,
) -> RunRequest {
    let mut request = RunRequest::default();
    match calculation {
        Some(text) => ui.label(format!("Calcul déduit : {text}")),
        None => ui.weak("Calcul déduit : pas encore (il manque des BC ou un chargement)."),
    };
    ui.separator();
    request.goto = diagnostics_ui(ui, diags, None);
    let blocked = diags.iter().any(|d| d.severity == Severity::Error);
    if diags.is_empty() {
        ui.label("✅ Le modèle est prêt.");
    }
    ui.separator();
    let idle = !compute.is_running();
    request.run = ui
        .add_enabled(!blocked && idle, egui::Button::new("▶ Calculer"))
        .on_hover_text("Enregistre le modèle puis lance le calcul de `maille solve`")
        .clicked();
    match compute {
        ComputeState::Finished(pvd) => {
            ui.weak(format!("Résultats : {}", pvd.display()));
        }
        ComputeState::Cancelled => {
            ui.weak("Calcul annulé.");
        }
        ComputeState::Failed(why) => {
            ui.colored_label(Color32::from_rgb(240, 110, 100), why);
        }
        _ => {}
    }
    request.save = ui
        .add_enabled(!blocked, egui::Button::new("Enregistrer le modèle (.toml)"))
        .clicked();
    ui.weak("Équivalent en ligne de commande : maille solve modele.toml.");
    if let Some(note) = note {
        ui.weak(note);
    }
    request
}

/// Le contenu du panneau « Calcul en cours » : étape, avancement, itération, coupes. Rend vrai si on
/// clique sur « Annuler ».
pub fn progress_ui(ui: &mut egui::Ui, running: &crate::compute::Running) -> bool {
    ui.label(egui::RichText::new("Calcul en cours").strong());
    if !running.phase.is_empty() {
        ui.weak(&running.phase);
    }
    if running.total > 0 {
        ui.add(
            egui::ProgressBar::new(running.done as f32 / running.total as f32)
                .text(format!("incrément {} / {}", running.done, running.total)),
        );
    } else {
        ui.add(egui::Spinner::new());
    }
    if let Some((iteration, residual)) = running.iteration {
        ui.weak(format!(
            "Newton : itération {iteration}, résidu {residual:.1e}"
        ));
    }
    if running.cuts > 0 {
        ui.weak(format!("incrément coupé {} fois", running.cuts));
    }
    ui.add_enabled(
        !running.stopping,
        egui::Button::new(if running.stopping {
            "Arrêt en cours…"
        } else {
            "Annuler"
        }),
    )
    .clicked()
}

/// Les diagnostics, une ligne chacun ; `only` ne garde que ceux d'un outil. Rend l'outil du
/// diagnostic cliqué.
pub fn diagnostics_ui(ui: &mut egui::Ui, diags: &[Diagnostic], only: Option<Tool>) -> Option<Tool> {
    let mut go = None;
    for d in diags.iter().filter(|d| only.is_none_or(|t| d.tool == t)) {
        let (mark, color) = match d.severity {
            Severity::Error => ("❌", Color32::from_rgb(240, 110, 100)),
            Severity::Warning => ("⚠", Color32::from_rgb(230, 180, 80)),
        };
        let response = ui.add(
            egui::Label::new(egui::RichText::new(format!("{mark} {}", d.message)).color(color))
                .sense(egui::Sense::click()),
        );
        if response
            .on_hover_text("Aller à l'outil qui corrige")
            .clicked()
        {
            go = Some(d.tool);
        }
    }
    go
}

/// Une propriété facultative : une case, et sa valeur quand elle est cochée.
fn optional(ui: &mut egui::Ui, label: &str, value: &mut Option<f64>, default: f64) {
    ui.horizontal(|ui| {
        let mut on = value.is_some();
        if ui.checkbox(&mut on, label).changed() {
            *value = on.then_some(default);
        }
        if let Some(v) = value {
            let speed = (v.abs() * 0.01).max(1e-12);
            ui.add(egui::DragValue::new(v).speed(speed));
        }
    });
}

/// Trois composantes (ou deux en 2D) d'un vecteur.
fn vector(ui: &mut egui::Ui, label: &str, v: &mut [f64; 3], dim: usize) {
    ui.horizontal(|ui| {
        ui.label(label);
        for c in v.iter_mut().take(dim) {
            ui.add(egui::DragValue::new(c).speed(0.05));
        }
    });
}

/// Inspecteur « Thermique » : températures imposées, sources, convections et couplage.
pub fn thermal_tool(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo, forms: &mut Forms) {
    if mesh.sets.is_empty() {
        ui.weak("Maille d'abord, et crée des Sets (espace Maillage).");
        return;
    }
    ui.label("Température imposée sur un Set :");
    let mut remove = None;
    for (i, b) in draft
        .bcs
        .iter()
        .enumerate()
        .filter(|(_, b)| b.dof == DofDef::T)
    {
        ui.horizontal(|ui| {
            ui.label(format!("{} · T = {}", b.set, b.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove.take() {
        draft.bcs.remove(i);
    }
    ui.horizontal(|ui| {
        set_picker(ui, "t-set", mesh.sets, &mut forms.temperature.set);
        ui.add(egui::DragValue::new(&mut forms.temperature.value).speed(1.0));
    });
    if ui.button("Ajouter la température").clicked() && !forms.temperature.set.is_empty() {
        draft.bcs.push(ValueDef {
            set: forms.temperature.set.clone(),
            dof: DofDef::T,
            value: forms.temperature.value,
        });
    }
    ui.separator();
    ui.label("HeatSource (puissance volumique) :");
    for (i, h) in draft.heat_sources.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {}", h.set, h.value));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove.take() {
        draft.heat_sources.remove(i);
    }
    ui.horizontal(|ui| {
        let mut sets = vec![maille_io::msh::DOMAIN_SET];
        sets.extend_from_slice(mesh.sets);
        set_picker(ui, "heat-set", &sets, &mut forms.heat_set);
        ui.add(egui::DragValue::new(&mut forms.heat).speed(0.1));
    });
    if ui.button("Ajouter la HeatSource").clicked() && !forms.heat_set.is_empty() {
        draft.heat_sources.push(HeatSourceDef {
            set: forms.heat_set.clone(),
            value: forms.heat,
        });
    }
    ui.separator();
    ui.label("Convection h (T − T∞) sur un bord :");
    for (i, c) in draft.convections.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · h = {}, T∞ = {}", c.set, c.film, c.sink));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        draft.convections.remove(i);
    }
    set_picker(ui, "conv-set", mesh.sets, &mut forms.convection_set);
    ui.horizontal(|ui| {
        ui.label("h");
        ui.add(egui::DragValue::new(&mut forms.film).speed(0.1));
        ui.label("T∞");
        ui.add(egui::DragValue::new(&mut forms.sink).speed(1.0));
    });
    if ui.button("Ajouter la Convection").clicked() && !forms.convection_set.is_empty() {
        draft.convections.push(ConvectionDef {
            set: forms.convection_set.clone(),
            film: forms.film,
            sink: forms.sink,
        });
    }
    ui.separator();
    let mut coupled = draft.reference_temperature.is_some();
    if ui
        .checkbox(&mut coupled, "Dilater la pièce par les températures")
        .on_hover_text(
            "Couplage thermo-mécanique : le thermique, puis le mécanique avec α (T − T_ref)",
        )
        .changed()
    {
        draft.reference_temperature = coupled.then_some(20.0);
    }
    if let Some(reference) = &mut draft.reference_temperature {
        ui.horizontal(|ui| {
            ui.label("T_ref (sans contrainte)");
            ui.add(egui::DragValue::new(reference).speed(1.0));
        });
    }
}

/// Inspecteur « Contact » : Nodes d'un bord contre une surface rigide, dessinée dans la vue.
pub fn contact_tool(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo, form: &mut ContactForm) {
    if mesh.sets.is_empty() {
        ui.weak("Maille d'abord, et crée des Sets (espace Maillage).");
        return;
    }
    ui.label("Sans frottement, contre une surface rigide (pénalité). Elle est dessinée en orange.");
    ui.separator();
    let mut remove = None;
    for (i, c) in draft.contacts.iter().enumerate() {
        ui.horizontal(|ui| {
            ui.label(format!("{} · {:?}", c.set, c.kind));
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        draft.contacts.remove(i);
    }
    ui.separator();
    ui.label("Ajouter :");
    let dim = mesh.dim();
    set_picker(ui, "contact-set", mesh.sets, &mut form.set);
    let three_d = dim == 3;
    ui.horizontal(|ui| {
        ui.radio_value(&mut form.kind, ContactKind::Plane, "plan");
        ui.radio_value(
            &mut form.kind,
            ContactKind::Cylinder,
            if three_d { "cylindre" } else { "cercle" },
        );
        if three_d {
            ui.radio_value(&mut form.kind, ContactKind::Sphere, "sphère");
        }
    });
    if !three_d && form.kind == ContactKind::Sphere {
        form.kind = ContactKind::Cylinder;
    }
    match form.kind {
        ContactKind::Plane => {
            vector(ui, "point", &mut form.point, dim);
            vector(ui, "normale (vers la pièce)", &mut form.normal, dim);
        }
        ContactKind::Cylinder => {
            vector(
                ui,
                if three_d { "point de l'axe" } else { "centre" },
                &mut form.point,
                dim,
            );
            if three_d {
                vector(ui, "axe", &mut form.axis, dim);
            }
        }
        ContactKind::Sphere => vector(ui, "centre", &mut form.center, dim),
    }
    if form.kind != ContactKind::Plane {
        ui.horizontal(|ui| {
            ui.label("rayon");
            ui.add(
                egui::DragValue::new(&mut form.radius)
                    .speed(0.05)
                    .range(0.0..=f64::MAX),
            );
        });
    }
    vector(ui, "déplacement à λ = 1", &mut form.motion, dim);
    if ui.button("Ajouter le Contact").clicked() && !form.set.is_empty() {
        draft.contacts.push(form.to_def(dim));
    }
}

/// Inspecteur « Step, Amplitude » : le quasi-statique (avec son amplitude) et le transitoire thermique.
pub fn step_tool(ui: &mut egui::Ui, draft: &mut Draft) {
    let mut on = draft.quasi_step.is_some();
    if ui.checkbox(&mut on, "Step quasi-statique").changed() {
        draft.quasi_step = on.then(Draft::default_quasi_step);
    }
    if let Some(step) = &mut draft.quasi_step {
        ui.horizontal(|ui| {
            ui.label("durée");
            ui.add(
                egui::DragValue::new(&mut step.duration)
                    .speed(0.05)
                    .range(0.0..=f64::MAX),
            );
            ui.label("incréments");
            ui.add(egui::DragValue::new(&mut step.increments).range(1..=100_000));
        });
        let mut custom = step.amplitude.is_some();
        if ui
            .checkbox(&mut custom, "Amplitude λ(t) personnalisée")
            .on_hover_text("Sans elle : une rampe de 0 à 1 sur la durée")
            .changed()
        {
            step.amplitude = custom.then(|| vec![[0.0, 0.0], [step.duration, 1.0]]);
        }
        if let Some(points) = &mut step.amplitude {
            let mut remove = None;
            for (i, p) in points.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label("t");
                    ui.add(egui::DragValue::new(&mut p[0]).speed(0.05));
                    ui.label("λ");
                    ui.add(egui::DragValue::new(&mut p[1]).speed(0.05));
                    if ui.small_button("×").clicked() {
                        remove = Some(i);
                    }
                });
            }
            if let Some(i) = remove {
                points.remove(i);
            }
            if ui.button("Ajouter un point").clicked() {
                let last = points.last().map_or(0.0, |p| p[0]);
                points.push([last + 1.0, 1.0]);
            }
        }
    }
    ui.separator();
    let mut thermal = draft.thermal_step.is_some();
    if ui.checkbox(&mut thermal, "Transitoire thermique").changed() {
        draft.thermal_step = thermal.then(Draft::default_thermal_step);
    }
    if let Some(t) = &mut draft.thermal_step {
        ui.horizontal(|ui| {
            ui.label("durée");
            ui.add(
                egui::DragValue::new(&mut t.duration)
                    .speed(1.0)
                    .range(0.0..=f64::MAX),
            );
            ui.label("incréments");
            ui.add(egui::DragValue::new(&mut t.increments).range(1..=100_000));
        });
        ui.horizontal(|ui| {
            ui.label("T initiale");
            ui.add(egui::DragValue::new(&mut t.initial_temperature).speed(1.0));
        });
        let mut theta = t.theta.unwrap_or(1.0);
        ui.horizontal(|ui| {
            ui.label("θ (1 : Euler implicite, 0,5 : Crank-Nicolson)");
            ui.add(
                egui::DragValue::new(&mut theta)
                    .speed(0.05)
                    .range(0.0..=1.0),
            );
        });
        t.theta = ((theta - 1.0).abs() > f64::EPSILON).then_some(theta);
    }
}
