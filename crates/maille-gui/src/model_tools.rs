//! Les inspecteurs de l'espace Modèle : des formulaires qui éditent le `Draft` de `model_space`. Ils ne
//! décident de rien d'autre que de ce qu'on affiche ; les règles sont dans `model_space`.

use eframe::egui::{self, Color32};
use maille_io::model::Physics;
use maille_io::model::{DofDef, PlaneStateDef, PressureDef, ValueDef};
use maille_io::msh::CellKind;

use crate::model_space::{physics_label, Diagnostic, Draft, Severity, UnitSystem};
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
}

/// Ce que le maillage dit au formulaire.
pub struct MeshInfo<'a> {
    /// Sets de bord du maillage.
    pub sets: &'a [&'a str],
    pub kind: CellKind,
}

impl MeshInfo<'_> {
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

/// Inspecteur « Material » : unités, Material élastique, état plan et épaisseur.
pub fn material_tool(ui: &mut egui::Ui, draft: &mut Draft, mesh: &MeshInfo) {
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
    match &mut draft.material {
        None => {
            if ui.button("Définir un Material").clicked() {
                draft.material = Some(draft.steel());
            }
            ui.weak("Un acier pour commencer : règle ensuite E et ν.");
        }
        Some(m) => {
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
            if ui.button("Retirer le Material").clicked() {
                draft.material = None;
            }
        }
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

/// Inspecteur « Calculer » : le calcul déduit, les diagnostics et l'écriture du modèle. Rend vrai si
/// on demande d'enregistrer.
pub fn run_tool(
    ui: &mut egui::Ui,
    physics: Option<Physics>,
    diags: &[Diagnostic],
    note: Option<&str>,
) -> (bool, Option<Tool>) {
    match physics {
        Some(p) => ui.label(format!("Calcul déduit : {}", physics_label(p))),
        None => ui.weak("Calcul déduit : pas encore (il manque des BC ou un chargement)."),
    };
    ui.separator();
    let go = diagnostics_ui(ui, diags, None);
    let blocked = diags.iter().any(|d| d.severity == Severity::Error);
    if diags.is_empty() {
        ui.label("✅ Le modèle est prêt.");
    }
    ui.separator();
    let save = ui
        .add_enabled(!blocked, egui::Button::new("Enregistrer le modèle (.toml)"))
        .clicked();
    ui.weak("Puis : maille solve modele.toml. Le calcul depuis l'interface arrive ensuite.");
    if let Some(note) = note {
        ui.weak(note);
    }
    (save, go)
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
