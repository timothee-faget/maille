//! L'espace Modèle sans fenêtre (partie statique) : le brouillon du modèle que l'utilisateur remplit
//! (Material, BC, Load, Pressure), ses diagnostics en langage clair avant le calcul, et son écriture
//! en TOML par `model::write`. `app` ne fait que l'afficher.

use maille_io::model::{
    self, ContactDef, ContactKind, ConvectionDef, DofDef, EdgeLoadDef, HeatSourceDef,
    HyperelasticityDef, KinematicsDef, MeshingDef, ModelFile, Physics, PlaneStateDef, PressureDef,
    SectionDef, StepDef, StepKind, TemperaturesDef, ThermalStepDef, ValueDef,
};
use maille_io::msh::{CellKind, MshFile, DOMAIN_SET};
use maille_io::IoError;

use crate::state::Tool;

/// Système d'unités cohérent choisi au départ : il s'écrit dans la clé `units` et s'affiche à côté de
/// chaque champ. Le cœur ne convertit rien.
#[derive(Debug, Clone, PartialEq)]
pub enum UnitSystem {
    /// N, m, Pa.
    Si,
    /// N, mm, MPa.
    MmNMpa,
    /// Texte libre : aucune unité affichée.
    Custom(String),
}

impl UnitSystem {
    /// Le texte de la clé `units`.
    pub fn text(&self) -> String {
        match self {
            UnitSystem::Si => "N, m, Pa".to_string(),
            UnitSystem::MmNMpa => "N, mm, MPa".to_string(),
            UnitSystem::Custom(text) => text.clone(),
        }
    }

    /// Retrouve le système d'un texte de `units` ; un texte inconnu reste libre.
    pub fn from_text(text: &str) -> Self {
        match text.trim() {
            "N, m, Pa" => UnitSystem::Si,
            "N, mm, MPa" => UnitSystem::MmNMpa,
            other => UnitSystem::Custom(other.to_string()),
        }
    }

    pub fn length(&self) -> &'static str {
        match self {
            UnitSystem::Si => "m",
            UnitSystem::MmNMpa => "mm",
            UnitSystem::Custom(_) => "",
        }
    }

    pub fn force(&self) -> &'static str {
        match self {
            UnitSystem::Si | UnitSystem::MmNMpa => "N",
            UnitSystem::Custom(_) => "",
        }
    }

    /// Contrainte, module d'Young, pression.
    pub fn stress(&self) -> &'static str {
        match self {
            UnitSystem::Si => "Pa",
            UnitSystem::MmNMpa => "MPa",
            UnitSystem::Custom(_) => "",
        }
    }
}

/// Le Material élastique du modèle (un seul pour l'instant : il couvre tout le domaine).
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialForm {
    pub name: String,
    pub young: f64,
    pub poisson: f64,
    /// Thermique : conductivité k, masse volumique ρ et chaleur massique c (transitoire).
    pub conductivity: Option<f64>,
    pub density: Option<f64>,
    pub specific_heat: Option<f64>,
    /// Coefficient de dilatation α (couplage thermo-mécanique).
    pub expansion: Option<f64>,
    /// Plasticité de von Mises : limite σy0 et module d'écrouissage H.
    pub yield_stress: Option<f64>,
    pub hardening: Option<f64>,
    /// Loi néo-hookéenne (grandes déformations) à la place de Saint-Venant–Kirchhoff.
    pub neo_hookean: bool,
}

impl MaterialForm {
    /// Le Material du fichier de modèle ; `large` : le modèle est en grandes déformations.
    pub fn to_def(&self, large: bool) -> model::MaterialDef {
        model::MaterialDef {
            young: Some(self.young),
            poisson: Some(self.poisson),
            conductivity: self.conductivity,
            density: self.density,
            specific_heat: self.specific_heat,
            expansion: self.expansion,
            yield_stress: self.yield_stress,
            hardening: self.yield_stress.and(self.hardening),
            hyperelasticity: (large && self.neo_hookean).then_some(HyperelasticityDef::NeoHookean),
        }
    }
}

/// Un Material de plus, sur un Set d'Elements.
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub set: String,
    pub material: MaterialForm,
}

/// Ce que l'utilisateur a rempli dans l'espace Modèle.
#[derive(Debug, Clone, PartialEq)]
pub struct Draft {
    pub units: UnitSystem,
    /// Le Material, sur le Set d'Elements `material_set` (tout le domaine tant qu'il est seul).
    pub material: Option<MaterialForm>,
    pub material_set: String,
    /// D'autres Materials sur d'autres Sets d'Elements : ensemble, ils couvrent chaque Element une
    /// seule fois. L'épaisseur et l'état plan sont les mêmes pour tous.
    pub regions: Vec<Region>,
    /// Triangles seulement.
    pub plane_state: PlaneStateDef,
    pub thickness: f64,
    pub bcs: Vec<ValueDef>,
    /// Efforts par Node (`Load`).
    pub loads: Vec<ValueDef>,
    pub pressures: Vec<PressureDef>,
    /// Forces linéiques (N/mm) sur des bords (2D) ou des arêtes (3D).
    pub edge_loads: Vec<EdgeLoadDef>,
    pub heat_sources: Vec<HeatSourceDef>,
    pub convections: Vec<ConvectionDef>,
    pub contacts: Vec<ContactDef>,
    /// Grandes déformations (Lagrangien total).
    pub large: bool,
    /// T_ref du couplage thermo-mécanique : les températures de la partie thermique dilatent le modèle.
    pub reference_temperature: Option<f64>,
    /// Transitoire thermique (durée, incréments, θ, température initiale).
    pub thermal_step: Option<ThermalStepDef>,
    /// Calcul quasi-statique (Newton par incrément) avec son amplitude.
    pub quasi_step: Option<StepDef>,
}

impl Default for Draft {
    fn default() -> Self {
        Self {
            units: UnitSystem::MmNMpa,
            material: None,
            material_set: DOMAIN_SET.to_string(),
            regions: Vec::new(),
            plane_state: PlaneStateDef::Stress,
            thickness: 1.0,
            bcs: Vec::new(),
            loads: Vec::new(),
            pressures: Vec::new(),
            edge_loads: Vec::new(),
            heat_sources: Vec::new(),
            convections: Vec::new(),
            contacts: Vec::new(),
            large: false,
            reference_temperature: None,
            thermal_step: None,
            quasi_step: None,
        }
    }
}

impl Draft {
    /// Tous les Materials du modèle avec leur Set d'Elements : le premier, puis les autres.
    pub fn sections(&self) -> Vec<(&str, &MaterialForm)> {
        self.material
            .iter()
            .map(|m| (self.material_set.as_str(), m))
            .chain(self.regions.iter().map(|r| (r.set.as_str(), &r.material)))
            .collect()
    }

    /// Les Elements que les Sets des autres Materials (`regions`) ne couvrent pas : ce qui reste au
    /// premier Material quand le modèle en a plusieurs.
    pub fn remaining_elements(&self, mesh: &MshFile) -> Vec<usize> {
        let mut taken = vec![false; mesh.cells.len()];
        for r in &self.regions {
            for &e in mesh.sets.get(&r.set).map_or(&[][..], |s| &s.elements[..]) {
                taken[e] = true;
            }
        }
        (0..taken.len()).filter(|&e| !taken[e]).collect()
    }

    /// Un acier, dans le système d'unités choisi.
    pub fn steel(&self) -> MaterialForm {
        MaterialForm {
            name: "acier".to_string(),
            young: if self.units == UnitSystem::Si {
                210e9
            } else {
                210e3
            },
            poisson: 0.3,
            conductivity: None,
            density: None,
            specific_heat: None,
            expansion: None,
            yield_stress: None,
            hardening: None,
            neo_hookean: false,
        }
    }

    /// Une Step quasi-statique par défaut : rampe de 0 à 1 en 10 incréments.
    pub fn default_quasi_step() -> StepDef {
        StepDef {
            kind: StepKind::QuasiStatic,
            duration: 1.0,
            increments: 10,
            theta: None,
            initial_temperature: None,
            amplitude: None,
        }
    }

    /// Un transitoire thermique par défaut : 100 incréments d'Euler implicite à 20 degrés.
    pub fn default_thermal_step() -> ThermalStepDef {
        ThermalStepDef {
            duration: 100.0,
            increments: 20,
            theta: None,
            initial_temperature: 20.0,
        }
    }

    /// Le modèle a une partie mécanique (déplacements, efforts, pression, contact, dilatation…).
    pub fn has_mechanical(&self) -> bool {
        self.bcs.iter().any(|b| b.dof != DofDef::T)
            || self.loads.iter().any(|l| l.dof != DofDef::T)
            || !self.pressures.is_empty()
            || !self.edge_loads.is_empty()
            || !self.contacts.is_empty()
            || self.large
            || self.reference_temperature.is_some()
            || self.quasi_step.is_some()
    }

    /// Le modèle a une partie thermique.
    pub fn has_thermal(&self) -> bool {
        self.bcs.iter().any(|b| b.dof == DofDef::T)
            || self.loads.iter().any(|l| l.dof == DofDef::T)
            || !self.heat_sources.is_empty()
            || !self.convections.is_empty()
            || self.thermal_step.is_some()
    }

    /// Plasticité, grandes déformations ou contact : le calcul demande une Step quasi-statique.
    pub fn is_nonlinear(&self) -> bool {
        self.large
            || !self.contacts.is_empty()
            || self
                .sections()
                .iter()
                .any(|(_, m)| m.yield_stress.is_some())
    }

    /// Le fichier de modèle correspondant. `meshing` : la recette du maillage si on la connaît.
    pub fn to_file(
        &self,
        mesh_name: &str,
        mesh: &MshFile,
        meshing: Option<MeshingDef>,
    ) -> ModelFile {
        let planar = is_planar(mesh);
        let mut materials = std::collections::BTreeMap::new();
        let mut sections = Vec::new();
        for (set, m) in self.sections() {
            materials.insert(m.name.clone(), m.to_def(self.large));
            sections.push(SectionDef {
                set: set.to_string(),
                material: m.name.clone(),
                thickness: (planar && (self.thickness - 1.0).abs() > f64::EPSILON)
                    .then_some(self.thickness),
                plane_state: planar.then_some(self.plane_state),
            });
        }
        // sans partie mécanique, le transitoire thermique est le `[step]` du fichier
        let thermal_only = !self.has_mechanical();
        let (step, thermal_step) = match (&self.thermal_step, thermal_only) {
            (Some(t), true) => (
                Some(StepDef {
                    kind: StepKind::Transient,
                    duration: t.duration,
                    increments: t.increments,
                    theta: t.theta,
                    initial_temperature: Some(t.initial_temperature),
                    amplitude: None,
                }),
                None,
            ),
            (t, _) => (self.quasi_step.clone(), t.clone()),
        };
        ModelFile {
            units: self.units.text(),
            mesh: mesh_name.to_string(),
            analysis: None,
            kinematics: self.large.then_some(KinematicsDef::Large),
            materials,
            sections,
            bcs: self.bcs.clone(),
            loads: self.loads.clone(),
            heat_sources: self.heat_sources.clone(),
            convections: self.convections.clone(),
            pressures: self.pressures.clone(),
            edge_loads: self.edge_loads.clone(),
            contacts: self.contacts.clone(),
            temperatures: self.reference_temperature.map(|reference| TemperaturesDef {
                from: None,
                increment: None,
                reference,
            }),
            step,
            meshing,
            thermal_step,
        }
    }

    /// Le TOML du modèle.
    pub fn to_toml(
        &self,
        mesh_name: &str,
        mesh: &MshFile,
        meshing: Option<MeshingDef>,
    ) -> Result<String, IoError> {
        model::write(&self.to_file(mesh_name, mesh, meshing))
    }

    /// Le calcul que le contenu demande, s'il est assez rempli pour qu'on le déduise.
    pub fn physics(&self, mesh: &MshFile) -> Option<Physics> {
        model::physics(&self.to_file("", mesh, None)).ok()
    }

    /// Il y a de quoi charger la pièce : un effort, une pression, une source, une convection, un
    /// contact, ou une valeur imposée non nulle.
    pub fn has_loading(&self) -> bool {
        !self.loads.is_empty()
            || !self.pressures.is_empty()
            || !self.edge_loads.is_empty()
            || !self.heat_sources.is_empty()
            || !self.convections.is_empty()
            || !self.contacts.is_empty()
            || self.bcs.iter().any(|b| b.value != 0.0)
    }

    /// Noms des Sets que le modèle utilise.
    fn used_sets(&self) -> Vec<&str> {
        self.bcs
            .iter()
            .chain(&self.loads)
            .map(|v| v.set.as_str())
            .chain(self.pressures.iter().map(|p| p.set.as_str()))
            .chain(self.edge_loads.iter().map(|l| l.set.as_str()))
            .chain(self.heat_sources.iter().map(|h| h.set.as_str()))
            .chain(self.convections.iter().map(|c| c.set.as_str()))
            .chain(self.contacts.iter().map(|c| c.set.as_str()))
            .chain(self.sections().into_iter().map(|(set, _)| set))
            .collect()
    }
}

type Tri = [[f64; 3]; 3];

fn add(a: [f64; 3], b: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    (n > 1e-12).then(|| [a[0] / n, a[1] / n, a[2] / n])
}

/// Deux vecteurs unitaires orthogonaux entre eux et à `n` (unitaire).
fn plane_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let helper = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = unit(cross(n, helper)).unwrap_or([1.0, 0.0, 0.0]);
    (u, cross(n, u))
}

/// Un quadrilatère `a b c d` en deux triangles.
fn quad(out: &mut Vec<Tri>, a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) {
    out.push([a, b, c]);
    out.push([a, c, d]);
}

fn vec3(v: &Option<Vec<f64>>) -> Option<[f64; 3]> {
    let v = v.as_ref()?;
    Some([*v.first()?, *v.get(1)?, v.get(2).copied().unwrap_or(0.0)])
}

/// La surface rigide d'un Contact, en triangles à surligner dans la vue (position à λ = 0). En 2D, une
/// bande mince du côté opposé au corps ; en 3D, un carré (plan), un tube court (cylindre) ou une
/// sphère. Vide si les paramètres manquent.
pub fn rigid_triangles(contact: &ContactDef, mesh: &MshFile) -> Vec<Tri> {
    let mut lo = [f64::MAX; 3];
    let mut hi = [f64::MIN; 3];
    for p in &mesh.nodes {
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[k] = lo[k].min(v);
            hi[k] = hi[k].max(v);
        }
    }
    if mesh.nodes.is_empty() {
        return Vec::new();
    }
    let size = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-12);
    let planar = is_planar(mesh);
    let band = size * 0.01;
    let segments = 48;
    let mut out = Vec::new();
    match contact.kind {
        ContactKind::Plane => {
            let (Some(p), Some(n)) = (vec3(&contact.point), vec3(&contact.normal).and_then(unit))
            else {
                return out;
            };
            let l = size * 1.5;
            if planar {
                let t = [-n[1], n[0], 0.0];
                let (a, b) = (add(p, t, -l), add(p, t, l));
                quad(&mut out, a, b, add(b, n, -band), add(a, n, -band));
            } else {
                let (u, v) = plane_basis(n);
                let corner = |su: f64, sv: f64| add(add(p, u, su * l), v, sv * l);
                quad(
                    &mut out,
                    corner(-1.0, -1.0),
                    corner(1.0, -1.0),
                    corner(1.0, 1.0),
                    corner(-1.0, 1.0),
                );
            }
        }
        ContactKind::Cylinder => {
            let (Some(p), Some(r)) = (vec3(&contact.point), contact.radius) else {
                return out;
            };
            let angle = |k: usize| std::f64::consts::TAU * k as f64 / segments as f64;
            if planar {
                let ring = |k: usize, radius: f64| {
                    [
                        p[0] + radius * angle(k).cos(),
                        p[1] + radius * angle(k).sin(),
                        0.0,
                    ]
                };
                for k in 0..segments {
                    quad(
                        &mut out,
                        ring(k, r - band),
                        ring(k, r),
                        ring(k + 1, r),
                        ring(k + 1, r - band),
                    );
                }
            } else {
                let Some(axis) = vec3(&contact.axis).and_then(unit) else {
                    return out;
                };
                let (u, v) = plane_basis(axis);
                let l = size;
                let at = |k: usize, along: f64| {
                    add(
                        add(add(p, u, r * angle(k).cos()), v, r * angle(k).sin()),
                        axis,
                        along,
                    )
                };
                for k in 0..segments {
                    quad(&mut out, at(k, -l), at(k + 1, -l), at(k + 1, l), at(k, l));
                }
            }
        }
        ContactKind::Sphere => {
            let (Some(c), Some(r)) = (vec3(&contact.center), contact.radius) else {
                return out;
            };
            let rings = 16;
            let at = |i: usize, k: usize| {
                let theta = std::f64::consts::PI * i as f64 / rings as f64;
                let phi = std::f64::consts::TAU * k as f64 / segments as f64;
                [
                    c[0] + r * theta.sin() * phi.cos(),
                    c[1] + r * theta.sin() * phi.sin(),
                    c[2] + r * theta.cos(),
                ]
            };
            for i in 0..rings {
                for k in 0..segments {
                    quad(
                        &mut out,
                        at(i, k),
                        at(i + 1, k),
                        at(i + 1, k + 1),
                        at(i, k + 1),
                    );
                }
            }
        }
    }
    out
}

/// La phrase qui dit quel calcul est déduit, avec sa nature quand elle n'est pas linéaire.
pub fn run_label(draft: &Draft, mesh: &MshFile) -> Option<String> {
    let physics = draft.physics(mesh)?;
    let mut text = physics_label(physics).to_string();
    if draft.has_mechanical() && draft.quasi_step.is_some() {
        text.push_str(", quasi-statique");
    }
    if draft.thermal_step.is_some() {
        text.push_str(", transitoire");
    }
    if draft.large {
        text.push_str(", grandes déformations");
    }
    Some(text)
}

/// Phrase qui décrit le calcul déduit.
pub fn physics_label(physics: Physics) -> &'static str {
    match physics {
        Physics::Static => "statique",
        Physics::Thermal => "thermique",
        Physics::Coupled => "thermo-mécanique",
    }
}

fn is_planar(mesh: &MshFile) -> bool {
    matches!(mesh.kind, CellKind::Tri3 | CellKind::Tri6)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Bloque le calcul.
    Error,
    /// Le calcul se lance, mais le résultat risque d'être décevant.
    Warning,
}

/// Un problème dit en clair : ce qui ne va pas et quoi faire. Un clic mène à l'outil qui le corrige.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub message: String,
    pub tool: Tool,
}

fn error(tool: Tool, message: String) -> Diagnostic {
    Diagnostic {
        severity: Severity::Error,
        message,
        tool,
    }
}

fn warning(tool: Tool, message: String) -> Diagnostic {
    Diagnostic {
        severity: Severity::Warning,
        message,
        tool,
    }
}

/// Ce qu'une charge demande d'un Set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Needs {
    Any,
    Faces,
    Edges,
}

/// Ce qui ne va pas dans le brouillon, avant le calcul. `user_sets` : les Sets créés par l'utilisateur,
/// dont on signale ceux que rien n'utilise.
pub fn diagnose(draft: &Draft, mesh: &MshFile, user_sets: &[&str]) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    if draft.material.is_none() {
        out.push(error(
            Tool::Material,
            "Aucun Material : définis-en un (module d'Young et coefficient de Poisson).".into(),
        ));
    }
    section_diagnostics(draft, mesh, &mut out);
    for (_, m) in draft.sections() {
        if m.young <= 0.0 {
            out.push(error(
                Tool::Material,
                format!("Le module d'Young de « {} » doit être positif.", m.name),
            ));
        }
        if m.poisson <= -1.0 || m.poisson >= 0.5 {
            out.push(error(
                Tool::Material,
                format!(
                    "Le coefficient de Poisson de « {} » doit être entre -1 et 0,5 (exclus).",
                    m.name
                ),
            ));
        } else if m.poisson > 0.45 && matches!(mesh.kind, CellKind::Tri3 | CellKind::Tet4) {
            out.push(warning(
                Tool::MeshSize,
                format!(
                    "ν = {} avec des éléments linéaires : le résultat risque d'être trop raide, coche quadratique.",
                    m.poisson
                ),
            ));
        }
    }
    if is_planar(mesh) && draft.thickness <= 0.0 {
        out.push(error(
            Tool::Material,
            "L'épaisseur doit être positive.".into(),
        ));
    }

    let mut check_set = |tool: Tool, what: &str, set: &str, needs: Needs| {
        match mesh
        .sets
        .get(set)
    {
        None => out.push(error(
            tool,
            format!("{what} renvoie au Set « {set} », qui n'existe pas dans le maillage (remaillage ?) : choisis un autre Set."),
        )),
        Some(s) if needs == Needs::Faces && s.faces.is_empty() => out.push(error(
            tool,
            format!("{what} demande un bord, et le Set « {set} » n'a pas de faces : choisis un Set de bord."),
        )),
        Some(s) if needs == Needs::Edges && s.edges.is_empty() => out.push(error(
            tool,
            format!("{what} demande des arêtes, et le Set « {set} » n'en a pas : trace-le avec l'outil « Arêtes → Set »."),
        )),
        Some(_) => {}
    }
    };
    for bc in &draft.bcs {
        check_set(
            if bc.dof == DofDef::T {
                Tool::Thermal
            } else {
                Tool::Support
            },
            &format!("La BC « {} »", dof_name(bc.dof)),
            &bc.set,
            Needs::Any,
        );
    }
    for load in &draft.loads {
        check_set(
            Tool::Load,
            &format!("Le Load « {} »", dof_name(load.dof)),
            &load.set,
            Needs::Any,
        );
    }
    for p in &draft.pressures {
        check_set(Tool::Load, "La Pressure", &p.set, Needs::Faces);
    }
    for l in &draft.edge_loads {
        let needs = if mesh.kind.dim() == 2 {
            Needs::Faces
        } else {
            Needs::Edges
        };
        check_set(Tool::Load, "La force linéique", &l.set, needs);
    }
    for h in &draft.heat_sources {
        check_set(Tool::Thermal, "La HeatSource", &h.set, Needs::Any);
    }
    for c in &draft.convections {
        check_set(Tool::Thermal, "La Convection", &c.set, Needs::Faces);
    }
    for c in &draft.contacts {
        check_set(Tool::Contact, "Le Contact", &c.set, Needs::Any);
    }

    for name in user_sets {
        if !draft.used_sets().contains(name) {
            out.push(warning(
                Tool::Support,
                format!("Le Set « {name} » n'est utilisé par aucune BC ni aucun chargement."),
            ));
        }
    }

    let t_bcs = draft.bcs.iter().filter(|b| b.dof == DofDef::T).count();
    let mechanical_bcs: Vec<ValueDef> = draft
        .bcs
        .iter()
        .filter(|b| b.dof != DofDef::T)
        .cloned()
        .collect();
    if draft.bcs.is_empty() && draft.contacts.is_empty() && draft.convections.is_empty() {
        out.push(error(
            Tool::Support,
            "Aucune BC : la pièce flotte. Ajoute un encastrement ou un déplacement imposé.".into(),
        ));
    } else if draft.has_mechanical()
        && draft.contacts.is_empty()
        && out.iter().all(|d| d.tool != Tool::Support)
    {
        // un contact retient la pièce contre la surface rigide : on ne sait pas dire à l'avance
        let free = free_rigid_modes(mesh, &mechanical_bcs);
        if !free.is_empty() {
            out.push(error(
                Tool::Support,
                format!(
                    "La pièce peut encore bouger d'un bloc ({}) : ajoute une BC qui bloque ce mouvement.",
                    free.join(", ")
                ),
            ));
        }
    }

    if draft.has_thermal() {
        let materials: Vec<&MaterialForm> = draft.sections().into_iter().map(|(_, m)| m).collect();
        if materials
            .iter()
            .any(|m| m.conductivity.is_none_or(|k| k <= 0.0))
        {
            out.push(error(
                Tool::Material,
                "Le calcul thermique demande une conductivité k > 0 dans le Material.".into(),
            ));
        }
        if t_bcs == 0 && draft.convections.is_empty() {
            out.push(error(
                Tool::Thermal,
                "Aucune BC thermique : fixe une température sur un Set, ou ajoute une Convection."
                    .into(),
            ));
        } else if t_bcs == 0 {
            out.push(warning(
                Tool::Thermal,
                "Il y a une Convection mais aucune BC thermique : la température tendra vers celle du milieu."
                    .into(),
            ));
        }
        if draft.thermal_step.is_some()
            && materials.iter().any(|m| {
                m.density.is_none_or(|r| r <= 0.0) || m.specific_heat.is_none_or(|c| c <= 0.0)
            })
        {
            out.push(error(
                Tool::Material,
                "Le transitoire thermique demande une masse volumique ρ et une chaleur massique c dans le Material."
                    .into(),
            ));
        }
        if draft.has_mechanical() {
            if draft.reference_temperature.is_none() {
                out.push(error(
                    Tool::Thermal,
                    "Le modèle a une partie thermique et une partie mécanique : donne la température de référence pour les coupler."
                        .into(),
                ));
            }
            if materials.iter().any(|m| m.expansion.is_none()) {
                out.push(error(
                    Tool::Material,
                    "Le couplage thermo-mécanique demande un coefficient de dilatation α dans le Material."
                        .into(),
                ));
            }
        }
    }
    if draft.is_nonlinear() && draft.quasi_step.is_none() {
        out.push(error(
            Tool::StepAmplitude,
            "Contact, plasticité et grandes déformations demandent une Step quasi-statique : ajoute-la."
                .into(),
        ));
    }
    if let Some(step) = &draft.quasi_step {
        if step.duration <= 0.0 || step.increments == 0 {
            out.push(error(
                Tool::StepAmplitude,
                "La Step doit avoir une durée positive et au moins un incrément.".into(),
            ));
        }
    }
    if let Some(t) = &draft.thermal_step {
        if t.duration <= 0.0 || t.increments == 0 {
            out.push(error(
                Tool::StepAmplitude,
                "Le transitoire thermique doit avoir une durée positive et au moins un incrément."
                    .into(),
            ));
        }
    }
    if (!draft.bcs.is_empty() || !draft.convections.is_empty()) && !draft.has_loading() {
        out.push(warning(
            Tool::Load,
            "Aucun chargement : sans effort ni déplacement imposé, la pièce ne se déforme pas."
                .into(),
        ));
    }

    // les règles du cœur, quand rien d'autre n'a été trouvé
    if out.iter().all(|d| d.severity != Severity::Error) {
        if let Err(e) = check_core(&draft.to_file("", mesh, None), mesh) {
            out.push(error(Tool::Run, e.to_string()));
        }
    }
    out
}

/// Les Materials et leurs Sets d'Elements : chaque Set existe et porte des Elements, chaque Element
/// a un Material et un seul, deux Materials de même nom ont les mêmes valeurs.
fn section_diagnostics(draft: &Draft, mesh: &MshFile, out: &mut Vec<Diagnostic>) {
    let sections = draft.sections();
    if sections.is_empty() {
        return;
    }
    if !draft.regions.is_empty() && draft.material_set == DOMAIN_SET {
        out.push(error(
            Tool::Material,
            "Le premier Material couvre tout le domaine et chevauche les autres : donne-lui le Set des Elements restants.".into(),
        ));
    }
    let mut count = vec![0usize; mesh.cells.len()];
    for (set, m) in &sections {
        match mesh.sets.get(*set) {
            None => out.push(error(
                Tool::Material,
                format!("Le Material « {} » renvoie au Set « {set} », qui n'existe pas dans le maillage : choisis un autre Set.", m.name),
            )),
            Some(s) if s.elements.is_empty() => out.push(error(
                Tool::Material,
                format!("Le Material « {} » demande un Set d'Elements, et « {set} » n'en a pas : crée-le avec « Cadre → Set » en mode Elements.", m.name),
            )),
            Some(s) => s.elements.iter().for_each(|&e| count[e] += 1),
        }
    }
    let twice = count.iter().filter(|&&c| c > 1).count();
    if twice > 0 && draft.material_set != DOMAIN_SET {
        out.push(error(
            Tool::Material,
            format!("{twice} Elements ont plusieurs Materials : les Sets d'Elements des Materials ne doivent pas se chevaucher."),
        ));
    }
    let none = count.iter().filter(|&&c| c == 0).count();
    if none > 0 {
        out.push(error(
            Tool::Material,
            format!("{none} Elements n'ont aucun Material : ajoute-en un sur un Set qui les couvre (ou « Elements restants »)."),
        ));
    }
    for (i, (_, a)) in sections.iter().enumerate() {
        if sections[..i]
            .iter()
            .any(|(_, b)| b.name == a.name && b != a)
        {
            out.push(error(
                Tool::Material,
                format!("Deux Materials s'appellent « {} » mais n'ont pas les mêmes valeurs : change un nom.", a.name),
            ));
            break;
        }
    }
}

/// Construit le modèle avec les règles du cœur ; un calcul couplé est coupé, et la partie mécanique
/// reçoit des températures de remplissage (seule la construction nous intéresse).
fn check_core(file: &ModelFile, mesh: &MshFile) -> Result<(), IoError> {
    if model::physics(file)? == Physics::Coupled {
        let (thermal, mechanical) = model::split(file)?;
        model::build(&thermal, mesh, None)?;
        let filler = vec![0.0; mesh.nodes.len()];
        model::build(&mechanical, mesh, Some(&filler))?;
    } else {
        model::build(file, mesh, None)?;
    }
    Ok(())
}

fn dof_name(dof: DofDef) -> &'static str {
    match dof {
        DofDef::Ux => "ux",
        DofDef::Uy => "uy",
        DofDef::Uz => "uz",
        DofDef::T => "t",
    }
}

/// Les mouvements de corps rigide que les BC ne bloquent pas, en clair ; vide si la pièce est tenue.
///
/// On projette les modes rigides (3 en 2D, 6 en 3D) sur les DOF bloqués : un mouvement rigide reste
/// libre si une combinaison des modes s'annule sur tous les DOF bloqués, soit un rang de la matrice
/// des modes inférieur à leur nombre.
pub fn free_rigid_modes(mesh: &MshFile, bcs: &[ValueDef]) -> Vec<String> {
    let planar = is_planar(mesh);
    let names: &[&str] = if planar {
        &[
            "translation en X",
            "translation en Y",
            "rotation dans le plan",
        ]
    } else {
        &[
            "translation en X",
            "translation en Y",
            "translation en Z",
            "rotation autour de X",
            "rotation autour de Y",
            "rotation autour de Z",
        ]
    };
    let n = names.len();
    let nodes = &mesh.nodes;
    if nodes.is_empty() {
        return Vec::new();
    }
    let count = nodes.len() as f64;
    let c = [
        nodes.iter().map(|p| p.x).sum::<f64>() / count,
        nodes.iter().map(|p| p.y).sum::<f64>() / count,
        nodes.iter().map(|p| p.z).sum::<f64>() / count,
    ];
    let extent = nodes
        .iter()
        .flat_map(|p| [(p.x - c[0]).abs(), (p.y - c[1]).abs(), (p.z - c[2]).abs()])
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);

    // valeur de chaque mode sur le DOF `axis` du Node `i`
    let mode = |m: usize, i: usize, axis: usize| -> f64 {
        let p = &nodes[i];
        let r = [
            (p.x - c[0]) / extent,
            (p.y - c[1]) / extent,
            (p.z - c[2]) / extent,
        ];
        if planar {
            match m {
                0 => f64::from(axis == 0),
                1 => f64::from(axis == 1),
                _ => [-r[1], r[0]][axis],
            }
        } else {
            match m {
                0..=2 => f64::from(axis == m),
                // ω × r pour une rotation unitaire autour de l'axe m - 3
                _ => {
                    let mut w = [0.0; 3];
                    w[m - 3] = 1.0;
                    [
                        w[1] * r[2] - w[2] * r[1],
                        w[2] * r[0] - w[0] * r[2],
                        w[0] * r[1] - w[1] * r[0],
                    ][axis]
                }
            }
        }
    };

    // Gram des modes sur les DOF bloqués
    let mut gram = vec![vec![0.0; n]; n];
    let dim = if planar { 2 } else { 3 };
    for bc in bcs {
        let axis = match bc.dof {
            DofDef::Ux => 0,
            DofDef::Uy => 1,
            DofDef::Uz => 2,
            DofDef::T => continue,
        };
        if axis >= dim {
            continue;
        }
        let Some(set) = mesh.sets.get(&bc.set) else {
            continue;
        };
        for &i in &set.nodes {
            let v: Vec<f64> = (0..n).map(|m| mode(m, i, axis)).collect();
            for a in 0..n {
                for b in 0..n {
                    gram[a][b] += v[a] * v[b];
                }
            }
        }
    }

    let scale = (0..n).map(|a| gram[a][a]).fold(0.0, f64::max);
    let tol = 1e-9 * scale.max(1e-300);
    let mut free: Vec<String> = (0..n)
        .filter(|&a| gram[a][a] <= tol)
        .map(|a| names[a].to_string())
        .collect();
    if free.is_empty() && rank(gram, tol) < n {
        free.push(
            "un mouvement combiné, par exemple une rotation autour d'un point d'appui".into(),
        );
    }
    free
}

/// Rang d'une matrice symétrique positive, par élimination avec pivot sur la diagonale.
fn rank(mut a: Vec<Vec<f64>>, tol: f64) -> usize {
    let n = a.len();
    let mut rank = 0;
    for _ in 0..n {
        let Some((p, &pivot)) = (0..n)
            .map(|i| (i, &a[i][i]))
            .max_by(|x, y| x.1.total_cmp(y.1))
        else {
            break;
        };
        if pivot <= tol {
            break;
        }
        rank += 1;
        let row = a[p].clone();
        for (i, line) in a.iter_mut().enumerate() {
            if i == p {
                continue;
            }
            let f = line[p] / pivot;
            for (value, r) in line.iter_mut().zip(&row) {
                *value -= f * r;
            }
        }
        // le pivot ne sert plus
        a[p].iter_mut().for_each(|v| *v = 0.0);
        a.iter_mut().for_each(|line| line[p] = 0.0);
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_core::mesh::Node;
    use maille_io::msh::Set;

    /// Plaque 2 × 1 en deux triangles : Sets `gauche` (nodes 0, 3), `droite` (1, 2), `bas` (0, 1).
    fn plate() -> MshFile {
        let node = |x, y| Node { x, y, z: 0.0 };
        let mut sets = std::collections::BTreeMap::new();
        let nodes_set = |nodes: Vec<usize>, faces: Vec<(usize, usize)>| Set {
            nodes,
            faces,
            ..Set::default()
        };
        sets.insert("gauche".to_string(), nodes_set(vec![0, 3], vec![(1, 2)]));
        sets.insert("droite".to_string(), nodes_set(vec![1, 2], vec![(0, 1)]));
        sets.insert("bas".to_string(), nodes_set(vec![0, 1], vec![(0, 0)]));
        sets.insert(
            DOMAIN_SET.to_string(),
            Set {
                nodes: vec![0, 1, 2, 3],
                faces: Vec::new(),
                elements: vec![0, 1],
                ..Set::default()
            },
        );
        MshFile {
            nodes: vec![
                node(0.0, 0.0),
                node(2.0, 0.0),
                node(2.0, 1.0),
                node(0.0, 1.0),
            ],
            kind: CellKind::Tri3,
            cells: vec![vec![0, 1, 2], vec![0, 2, 3]],
            sets,
        }
    }

    fn bc(set: &str, dof: DofDef) -> ValueDef {
        ValueDef {
            set: set.to_string(),
            dof,
            value: 0.0,
        }
    }

    fn filled() -> Draft {
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.bcs = vec![bc("gauche", DofDef::Ux), bc("gauche", DofDef::Uy)];
        d.pressures = vec![PressureDef {
            set: "droite".into(),
            value: -10.0,
        }];
        d
    }

    fn messages(d: &Draft, mesh: &MshFile, user: &[&str]) -> Vec<String> {
        diagnose(d, mesh, user)
            .into_iter()
            .map(|d| d.message)
            .collect()
    }

    #[test]
    fn a_clamped_plate_with_a_pressure_has_no_diagnostic() {
        assert_eq!(messages(&filled(), &plate(), &[]), Vec::<String>::new());
    }

    #[test]
    fn a_missing_material_is_reported() {
        let mut d = filled();
        d.material = None;
        let m = messages(&d, &plate(), &[]);
        assert!(m.iter().any(|s| s.contains("Aucun Material")), "{m:?}");
    }

    #[test]
    fn a_part_that_can_still_move_as_a_block_is_reported() {
        // blocage de ux seulement sur le bord gauche : la plaque glisse en Y et tourne
        let mut d = filled();
        d.bcs = vec![bc("gauche", DofDef::Ux)];
        let m = messages(&d, &plate(), &[]);
        assert!(
            m.iter()
                .any(|s| s.contains("bouger d'un bloc") && s.contains("translation en Y")),
            "{m:?}"
        );
        // deux DOF sur un seul Node d'un bord de deux Nodes : la rotation reste libre
        d.bcs = vec![bc("gauche", DofDef::Ux), bc("bas", DofDef::Uy)];
        assert!(
            messages(&d, &plate(), &[]).is_empty(),
            "ux gauche + uy bas tient la plaque"
        );
    }

    #[test]
    fn a_single_support_point_leaves_a_combined_rotation() {
        // un seul Node encastré : les deux translations sont bloquées, la rotation reste libre
        let mut mesh = plate();
        mesh.sets.insert(
            "coin".to_string(),
            Set {
                nodes: vec![0],
                faces: Vec::new(),
                ..Set::default()
            },
        );
        let free = free_rigid_modes(&mesh, &[bc("coin", DofDef::Ux), bc("coin", DofDef::Uy)]);
        assert_eq!(free.len(), 1);
        assert!(free[0].contains("mouvement combiné"), "{free:?}");
    }

    #[test]
    fn supports_that_cross_block_everything() {
        // ux sur le bord x = 0 et uy en un point : tout est bloqué ; uy seul sur deux Nodes d'un bord
        // horizontal laisse la translation en X
        let mut mesh = plate();
        mesh.sets.insert(
            "coin".to_string(),
            Set {
                nodes: vec![0],
                faces: Vec::new(),
                ..Set::default()
            },
        );
        let free = free_rigid_modes(&mesh, &[bc("gauche", DofDef::Ux), bc("coin", DofDef::Uy)]);
        assert!(
            free.is_empty(),
            "ux sur x = 0 et uy en un point bloquent tout : {free:?}"
        );
        let free = free_rigid_modes(&mesh, &[bc("bas", DofDef::Uy), bc("coin", DofDef::Uy)]);
        assert!(
            free.iter().any(|s| s.contains("translation en X")),
            "{free:?}"
        );
    }

    #[test]
    fn a_three_dimensional_part_needs_six_blocked_motions() {
        let node = |x, y, z| Node { x, y, z };
        let mut sets = std::collections::BTreeMap::new();
        sets.insert(
            "pied".to_string(),
            Set {
                nodes: vec![0, 1, 2],
                faces: vec![(0, 3)],
                ..Set::default()
            },
        );
        let mesh = MshFile {
            nodes: vec![
                node(0.0, 0.0, 0.0),
                node(1.0, 0.0, 0.0),
                node(0.0, 1.0, 0.0),
                node(0.0, 0.0, 1.0),
            ],
            kind: CellKind::Tet4,
            cells: vec![vec![0, 1, 2, 3]],
            sets,
        };
        let all = [
            bc("pied", DofDef::Ux),
            bc("pied", DofDef::Uy),
            bc("pied", DofDef::Uz),
        ];
        assert!(free_rigid_modes(&mesh, &all).is_empty());
        assert!(
            !free_rigid_modes(&mesh, &all[..2]).is_empty(),
            "sans uz, la pièce bouge encore en Z"
        );
    }

    #[test]
    fn a_set_that_vanished_is_reported_by_name() {
        let mut d = filled();
        d.bcs.push(bc("disparu", DofDef::Ux));
        let m = messages(&d, &plate(), &[]);
        assert!(
            m.iter()
                .any(|s| s.contains("disparu") && s.contains("n'existe pas")),
            "{m:?}"
        );
    }

    #[test]
    fn a_pressure_needs_a_boundary_set() {
        let mut d = filled();
        d.pressures = vec![PressureDef {
            set: DOMAIN_SET.into(),
            value: 1.0,
        }];
        let m = messages(&d, &plate(), &[]);
        assert!(m.iter().any(|s| s.contains("n'a pas de faces")), "{m:?}");
    }

    #[test]
    fn an_unused_user_set_is_a_warning() {
        let diags = diagnose(&filled(), &plate(), &["bas", "gauche"]);
        let unused: Vec<_> = diags
            .iter()
            .filter(|d| d.message.contains("« bas »"))
            .collect();
        assert_eq!(unused.len(), 1);
        assert_eq!(unused[0].severity, Severity::Warning);
        assert!(!diags.iter().any(|d| d.message.contains("« gauche »")));
    }

    #[test]
    fn a_nearly_incompressible_material_with_linear_elements_is_a_warning() {
        let mut d = filled();
        d.material.as_mut().unwrap().poisson = 0.49;
        let diags = diagnose(&d, &plate(), &[]);
        assert!(diags
            .iter()
            .any(|d| d.severity == Severity::Warning && d.message.contains("quadratique")));
        assert!(diags.iter().all(|d| d.severity != Severity::Error));
    }

    #[test]
    fn an_impossible_material_is_an_error() {
        let mut d = filled();
        d.material.as_mut().unwrap().young = -1.0;
        d.material.as_mut().unwrap().poisson = 0.5;
        let diags = diagnose(&d, &plate(), &[]);
        assert_eq!(
            diags
                .iter()
                .filter(|d| d.severity == Severity::Error)
                .count(),
            2
        );
    }

    #[test]
    fn no_load_at_all_is_a_warning_but_an_imposed_displacement_is_a_load() {
        let mut d = filled();
        d.pressures.clear();
        assert!(messages(&d, &plate(), &[])
            .iter()
            .any(|s| s.contains("Aucun chargement")));
        d.bcs.push(ValueDef {
            set: "droite".into(),
            dof: DofDef::Ux,
            value: 0.1,
        });
        assert!(messages(&d, &plate(), &[]).is_empty());
    }

    #[test]
    fn the_calculation_is_deduced_from_the_content() {
        let mesh = plate();
        assert_eq!(Draft::default().physics(&mesh), None);
        assert_eq!(filled().physics(&mesh), Some(Physics::Static));
    }

    #[test]
    fn the_draft_round_trips_through_the_toml() {
        let mesh = plate();
        let d = filled();
        let text = d.to_toml("plaque.msh", &mesh, None).unwrap();
        let file = model::parse(&text).unwrap();
        assert_eq!(file, d.to_file("plaque.msh", &mesh, None));
        assert_eq!(file.units, "N, mm, MPa");
        assert_eq!(file.sections[0].plane_state, Some(PlaneStateDef::Stress));
        // et le modèle se construit
        model::build(&file, &mesh, None).unwrap();
    }

    /// Le parcours de l'interface sur `plaque.dxf` : maillage, Material, BC, calcul, comparé à la
    /// solution exacte de `examples/traction.toml` (ux = 0,001 x, uy = −ν 0,001 y).
    #[test]
    fn the_traction_plate_built_in_the_model_space_gives_the_exact_solution() {
        use crate::mesh_space::{import, mesh, Imported};
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/plaque.dxf");
        let Imported::Cad(cad) = import(&path).unwrap() else {
            panic!("plaque.dxf est un dessin");
        };
        let built = mesh(&cad, 0.5, false, &Default::default()).unwrap();
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.thickness = 2.0;
        d.bcs = vec![
            bc("gauche", DofDef::Ux),
            ValueDef {
                set: "droite".into(),
                dof: DofDef::Ux,
                value: 0.004,
            },
            bc("bas", DofDef::Uy),
        ];
        assert_eq!(messages(&d, &built.file, &[]), Vec::<String>::new());
        let text = d.to_toml("plaque.msh", &built.file, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &built.file, None).unwrap();
        let field = loaded.model.solve().unwrap();
        for (i, p) in built.file.nodes.iter().enumerate() {
            let (ux, uy) = (field.displacements[2 * i], field.displacements[2 * i + 1]);
            assert!((ux - 0.001 * p.x).abs() < 1e-9, "ux du Node {i}");
            assert!((uy + 0.3 * 0.001 * p.y).abs() < 1e-9, "uy du Node {i}");
        }
    }

    fn example(name: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(name)
    }

    fn example_mesh(name: &str) -> MshFile {
        maille_io::msh::read(&std::fs::read_to_string(example(name)).unwrap()).unwrap()
    }

    fn value(set: &str, dof: DofDef, value: f64) -> ValueDef {
        ValueDef {
            set: set.into(),
            dof,
            value,
        }
    }

    /// Le tube thermo-mécanique de `examples/tube_thermomeca.toml`, rempli comme l'interface le fait.
    fn tube() -> Draft {
        let mut d = Draft::default();
        d.units = UnitSystem::Custom("N, mm, MPa, K, W".into());
        d.material = Some(MaterialForm {
            expansion: Some(1.2e-5),
            conductivity: Some(0.05),
            ..d.steel()
        });
        d.plane_state = PlaneStateDef::Strain;
        d.bcs = vec![
            value("interieur", DofDef::T, 120.0),
            value("exterieur", DofDef::T, 20.0),
            value("sym_x", DofDef::Uy, 0.0),
            value("sym_y", DofDef::Ux, 0.0),
        ];
        d.reference_temperature = Some(20.0);
        d
    }

    #[test]
    fn the_coupled_tube_gives_the_example_file() {
        let mesh = example_mesh("tube.msh");
        let mut expected =
            model::parse(&std::fs::read_to_string(example("tube_thermomeca.toml")).unwrap())
                .unwrap();
        expected.temperatures.as_mut().unwrap().increment = None;
        let d = tube();
        assert_eq!(d.to_file("tube.msh", &mesh, None), expected);
        assert_eq!(d.physics(&mesh), Some(Physics::Coupled));
        assert_eq!(messages(&d, &mesh, &[]), Vec::<String>::new());
        assert_eq!(run_label(&d, &mesh).unwrap(), "thermo-mécanique");
    }

    #[test]
    fn the_hertz_contact_gives_the_example_file() {
        let mesh = example_mesh("hertz.msh");
        let mut expected =
            model::parse(&std::fs::read_to_string(example("hertz.toml")).unwrap()).unwrap();
        expected.analysis = None;
        expected.sections[0].thickness = None;
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.plane_state = PlaneStateDef::Strain;
        d.bcs = vec![
            value("axe", DofDef::Ux, 0.0),
            value("droite", DofDef::Ux, 0.0),
            value("bas", DofDef::Uy, 0.0),
        ];
        d.contacts = vec![ContactDef {
            set: "dessus".into(),
            kind: ContactKind::Cylinder,
            point: Some(vec![0.0, 10.0]),
            normal: None,
            axis: None,
            center: None,
            radius: Some(10.0),
            motion: Some(vec![0.0, -0.08]),
        }];
        // sans Step, un contact est refusé en clair
        let m = messages(&d, &mesh, &[]);
        assert!(m.iter().any(|s| s.contains("Step quasi-statique")), "{m:?}");
        d.quasi_step = Some(StepDef {
            increments: 8,
            ..Draft::default_quasi_step()
        });
        assert_eq!(d.to_file("hertz.msh", &mesh, None), expected);
        assert_eq!(messages(&d, &mesh, &[]), Vec::<String>::new());
        assert_eq!(run_label(&d, &mesh).unwrap(), "statique, quasi-statique");
        // la surface rigide se dessine : un anneau de 48 quadrilatères autour du centre (0, 10)
        let ring = rigid_triangles(&d.contacts[0], &mesh);
        assert_eq!(ring.len(), 96);
        for t in &ring {
            for p in t {
                let r = ((p[0]).powi(2) + (p[1] - 10.0).powi(2)).sqrt();
                assert!((9.0..=10.0 + 1e-9).contains(&r), "r = {r}");
            }
        }
    }

    #[test]
    fn a_rigid_plane_is_drawn_behind_its_surface() {
        let mesh = plate();
        let contact = ContactDef {
            set: "gauche".into(),
            kind: ContactKind::Plane,
            point: Some(vec![0.0, 1.5]),
            normal: Some(vec![0.0, -1.0]),
            axis: None,
            center: None,
            radius: None,
            motion: None,
        };
        let tris = rigid_triangles(&contact, &mesh);
        assert_eq!(tris.len(), 2);
        // le corps est du côté de la normale (y < 1,5) ; la bande est au-dessus de y = 1,5
        assert!(tris.iter().flatten().all(|p| p[1] >= 1.5 - 1e-12));
        assert!(rigid_triangles(
            &ContactDef {
                point: None,
                ..contact
            },
            &mesh
        )
        .is_empty());
    }

    #[test]
    fn thermal_diagnostics_name_what_is_missing() {
        let mesh = plate();
        let mut d = filled();
        d.heat_sources = vec![HeatSourceDef {
            set: DOMAIN_SET.into(),
            value: 1.0,
        }];
        let m = messages(&d, &mesh, &[]);
        assert!(m.iter().any(|s| s.contains("conductivité")), "{m:?}");
        assert!(
            m.iter().any(|s| s.contains("température de référence")),
            "{m:?}"
        );
        assert!(m.iter().any(|s| s.contains("dilatation")), "{m:?}");

        // thermique seul : une température imposée, un transitoire sans ρ ni c
        let mut t = Draft::default();
        t.material = Some(MaterialForm {
            conductivity: Some(1.0),
            ..t.steel()
        });
        t.bcs = vec![value("gauche", DofDef::T, 100.0)];
        t.thermal_step = Some(Draft::default_thermal_step());
        let m = messages(&t, &mesh, &[]);
        assert!(m.iter().any(|s| s.contains("masse volumique")), "{m:?}");
        let material = t.material.as_mut().unwrap();
        material.density = Some(7.8e-9);
        material.specific_heat = Some(4.6e8);
        assert_eq!(messages(&t, &mesh, &[]), Vec::<String>::new());
        // sans partie mécanique, le transitoire devient le `[step]` du fichier
        let file = t.to_file("p.msh", &mesh, None);
        assert_eq!(file.step.unwrap().kind, StepKind::Transient);
        assert!(file.thermal_step.is_none());
        assert_eq!(run_label(&t, &mesh).unwrap(), "thermique, transitoire");
    }

    #[test]
    fn a_convection_without_a_thermal_bc_is_a_warning() {
        let mut t = Draft::default();
        t.material = Some(MaterialForm {
            conductivity: Some(1.0),
            ..t.steel()
        });
        t.convections = vec![ConvectionDef {
            set: "droite".into(),
            film: 10.0,
            sink: 20.0,
        }];
        let diags = diagnose(&t, &plate(), &[]);
        assert!(diags.iter().any(|d| d.severity == Severity::Warning
            && d.message.contains("Convection mais aucune BC thermique")));
        assert!(
            diags.iter().all(|d| d.severity != Severity::Error),
            "{diags:?}"
        );
    }

    #[test]
    fn plasticity_and_large_strains_ask_for_a_quasi_static_step() {
        let mesh = plate();
        let mut d = filled();
        d.material.as_mut().unwrap().yield_stress = Some(250.0);
        d.plane_state = PlaneStateDef::Strain;
        let m = messages(&d, &mesh, &[]);
        assert!(m.iter().any(|s| s.contains("Step quasi-statique")), "{m:?}");
        d.quasi_step = Some(Draft::default_quasi_step());
        assert_eq!(messages(&d, &mesh, &[]), Vec::<String>::new());
        d.large = true;
        assert!(d.to_file("p.msh", &mesh, None).kinematics.is_some());
    }

    #[test]
    fn units_follow_the_chosen_system() {
        assert_eq!(UnitSystem::from_text("N, m, Pa"), UnitSystem::Si);
        assert_eq!(UnitSystem::from_text("N, mm, MPa").stress(), "MPa");
        let custom = UnitSystem::from_text("lbf, in, psi");
        assert_eq!(custom.stress(), "");
        assert_eq!(custom.text(), "lbf, in, psi");
        let d = Draft {
            units: UnitSystem::Si,
            ..Draft::default()
        };
        assert_eq!(d.steel().young, 210e9);
    }

    /// `plate()` avec deux Sets d'Elements : « bas » côté triangle 0, « haut » côté triangle 1.
    fn two_region_plate() -> MshFile {
        let mut mesh = plate();
        for (name, element) in [("moitie_a", 0), ("moitie_b", 1)] {
            mesh.sets.insert(
                name.to_string(),
                Set {
                    elements: vec![element],
                    ..Set::default()
                },
            );
        }
        mesh
    }

    fn two_materials(d: &mut Draft) {
        d.material_set = "moitie_a".into();
        let mut other = d.steel();
        other.name = "alu".into();
        other.young = 70e3;
        d.regions = vec![Region {
            set: "moitie_b".into(),
            material: other,
        }];
    }

    #[test]
    fn two_materials_on_two_element_sets_have_no_diagnostic_and_two_sections() {
        let mut d = filled();
        two_materials(&mut d);
        let mesh = two_region_plate();
        assert_eq!(messages(&d, &mesh, &[]), Vec::<String>::new());
        let file = d.to_file("plate.msh", &mesh, None);
        assert_eq!(file.materials.len(), 2);
        let sets: Vec<(&str, &str)> = file
            .sections
            .iter()
            .map(|s| (s.set.as_str(), s.material.as_str()))
            .collect();
        assert_eq!(sets, [("moitie_a", "acier"), ("moitie_b", "alu")]);
        // un Material n'est pas un Set inutile
        assert!(!messages(&d, &mesh, &["moitie_a"])
            .iter()
            .any(|m| m.contains("n'est utilisé")));
    }

    #[test]
    fn badly_placed_materials_are_said_in_plain_words() {
        let mesh = two_region_plate();
        // le premier garde tout le domaine et chevauche l'autre
        let mut d = filled();
        two_materials(&mut d);
        d.material_set = DOMAIN_SET.into();
        let m = messages(&d, &mesh, &[]);
        assert!(
            m.iter().any(|s| s.contains("couvre tout le domaine")),
            "{m:?}"
        );
        // deux Materials sur le même Set d'Elements
        let mut d = filled();
        two_materials(&mut d);
        d.regions[0].set = "moitie_a".into();
        let m = messages(&d, &mesh, &[]);
        assert!(
            m.iter().any(|s| s.contains("1 Elements ont plusieurs")),
            "{m:?}"
        );
        assert!(
            m.iter().any(|s| s.contains("1 Elements n'ont aucun")),
            "{m:?}"
        );
        // un Set qui n'existe pas, un Set de bord
        let mut d = filled();
        two_materials(&mut d);
        d.regions[0].set = "fantome".into();
        let m = messages(&d, &mesh, &[]);
        assert!(
            m.iter()
                .any(|s| s.contains("« fantome », qui n'existe pas")),
            "{m:?}"
        );
        d.regions[0].set = "gauche".into();
        let m = messages(&d, &mesh, &[]);
        assert!(
            m.iter().any(|s| s.contains("demande un Set d'Elements")),
            "{m:?}"
        );
        // deux Materials de même nom, valeurs différentes
        let mut d = filled();
        two_materials(&mut d);
        d.regions[0].material.name = "acier".into();
        let m = messages(&d, &mesh, &[]);
        assert!(
            m.iter().any(|s| s.contains("s'appellent « acier »")),
            "{m:?}"
        );
    }

    #[test]
    fn the_remaining_elements_are_those_the_other_materials_leave() {
        let mesh = two_region_plate();
        let mut d = filled();
        assert_eq!(d.remaining_elements(&mesh), vec![0, 1]);
        two_materials(&mut d);
        assert_eq!(d.remaining_elements(&mesh), vec![0]);
    }

    #[test]
    fn plasticity_in_any_material_asks_for_a_quasi_static_step() {
        let mut d = filled();
        two_materials(&mut d);
        assert!(!d.is_nonlinear());
        d.regions[0].material.yield_stress = Some(250.0);
        assert!(d.is_nonlinear());
    }

    /// Deux Materials en série dans la plaque 4 × 2 (ν = 0 pour que les moitiés ne se gênent pas) :
    /// σ = u / (L1/E1 + L2/E2), et la réaction vaut σ fois la section.
    #[test]
    fn two_materials_in_series_follow_the_series_formula() {
        let mut mesh = example_mesh("plaque.msh");
        let left: Vec<usize> = (0..mesh.cells.len())
            .filter(|&e| {
                mesh.cells[e][..3]
                    .iter()
                    .map(|&n| mesh.nodes[n].x)
                    .sum::<f64>()
                    / 3.0
                    <= 2.0
            })
            .collect();
        let right: Vec<usize> = (0..mesh.cells.len())
            .filter(|e| !left.contains(e))
            .collect();
        for (name, elements) in [("gauche_mat", left), ("droite_mat", right)] {
            mesh.sets.insert(
                name.to_string(),
                Set {
                    elements,
                    ..Set::default()
                },
            );
        }
        let mut d = Draft {
            thickness: 2.0,
            ..Draft::default()
        };
        let mut first = d.steel();
        first.poisson = 0.0;
        d.material = Some(first);
        d.material_set = "gauche_mat".into();
        let mut second = d.steel();
        (second.name, second.young, second.poisson) = ("alu".into(), 70e3, 0.0);
        d.regions = vec![Region {
            set: "droite_mat".into(),
            material: second,
        }];
        d.bcs = vec![
            bc("gauche", DofDef::Ux),
            value("droite", DofDef::Ux, 0.004),
            bc("bas", DofDef::Uy),
        ];
        assert_eq!(messages(&d, &mesh, &[]), Vec::<String>::new());
        let text = d.to_toml("plaque.msh", &mesh, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &mesh, None).unwrap();
        let field = loaded.model.solve().unwrap();

        let stress = 0.004 / (2.0 / 210e3 + 2.0 / 70e3);
        let force: f64 = mesh.sets["droite"]
            .nodes
            .iter()
            .map(|&n| field.reactions[2 * n])
            .sum();
        let expected = stress * 2.0 * 2.0;
        println!("force {force} attendue {expected}");
        assert!(
            (force - expected).abs() < 0.02 * expected,
            "{force} au lieu de {expected}"
        );
    }
}
