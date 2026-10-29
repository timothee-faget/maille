//! Forme du fichier de modèle TOML : des données, aucune règle. La construction du `Model` est
//! dans `build`. Une clé inconnue ou mal écrite est une erreur (`deny_unknown_fields`).
//!
//! ```toml
//! units = "N, mm, MPa"            # obligatoire, texte libre, recopié dans les résultats
//! mesh  = "tube.msh"
//! analysis = "static"             # ou "thermal"
//! kinematics = "small"            # ou "large" : grandes déformations (statique, avec [step])
//!
//! [materials.acier]
//! young = 210e3
//! poisson = 0.3
//!
//! [[sections]]                    # Material + épaisseur + état plan sur un Set d'Elements
//! set = "domaine"
//! material = "acier"
//! plane_state = "strain"          # "stress" ou "strain" (triangles, analyse statique)
//!
//! [[pressures]]                   # pression sur un bord : positive, elle comprime
//! set = "interieur"
//! value = 100.0
//!
//! [[contacts]]                    # contact sans frottement contre une surface rigide (statique)
//! set = "dessus"                  # bord du corps déformable (Set de Nodes)
//! kind = "cylinder"               # "plane", "cylinder" ou "sphere"
//! point = [0.0, 10.0]             # plan et cylindre ; `center` pour la sphère
//! radius = 10.0                   # cylindre et sphère ; `normal` pour le plan, `axis` en 3D
//! motion = [0.0, -0.08]           # déplacement de la surface à λ = 1 (nul par défaut)
//!
//! [[bcs]]                         # Dirichlet sur un Set de Nodes
//! set = "sym_y"
//! dof = "ux"                      # ux, uy, uz (statique), t (thermique)
//! value = 0.0
//! ```

use std::collections::BTreeMap;

use serde::Deserialize;

/// Fichier de modèle.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModelFile {
    /// Unités cohérentes de l'utilisateur : le cœur ne convertit rien.
    pub units: String,
    /// Nom du fichier de maillage `.msh`, ouvert par le CLI.
    pub mesh: String,
    pub analysis: AnalysisKind,
    /// Cinématique du calcul statique, petites déformations par défaut.
    pub kinematics: Option<KinematicsDef>,
    #[serde(default)]
    pub materials: BTreeMap<String, MaterialDef>,
    #[serde(default)]
    pub sections: Vec<SectionDef>,
    #[serde(default)]
    pub bcs: Vec<ValueDef>,
    #[serde(default)]
    pub loads: Vec<ValueDef>,
    #[serde(default)]
    pub heat_sources: Vec<HeatSourceDef>,
    #[serde(default)]
    pub convections: Vec<ConvectionDef>,
    #[serde(default)]
    pub pressures: Vec<PressureDef>,
    /// Contacts contre des surfaces rigides : demandent `[step] kind = "quasi_static"`.
    #[serde(default)]
    pub contacts: Vec<ContactDef>,
    /// Couplage thermo-mécanique : les températures d'un calcul thermique dilatent le modèle
    /// statique. Le CLI lit le fichier `from` et les passe à `build`.
    pub temperatures: Option<TemperaturesDef>,
    /// Sans `[step]` : calcul linéaire.
    pub step: Option<StepDef>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AnalysisKind {
    Static,
    Thermal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum KinematicsDef {
    Small,
    Large,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HyperelasticityDef {
    NeoHookean,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialDef {
    pub young: Option<f64>,
    pub poisson: Option<f64>,
    pub conductivity: Option<f64>,
    pub density: Option<f64>,
    pub specific_heat: Option<f64>,
    pub expansion: Option<f64>,
    /// Active la plasticité de von Mises.
    pub yield_stress: Option<f64>,
    /// Module d'écrouissage H, 0 par défaut (parfaitement plastique).
    pub hardening: Option<f64>,
    /// Loi hyperélastique, grandes déformations seulement ; Saint-Venant–Kirchhoff sans elle.
    pub hyperelasticity: Option<HyperelasticityDef>,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SectionDef {
    pub set: String,
    pub material: String,
    /// Triangles seulement, 1 par défaut.
    pub thickness: Option<f64>,
    /// Triangles seulement ; obligatoire en analyse statique.
    pub plane_state: Option<PlaneStateDef>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PlaneStateDef {
    Stress,
    Strain,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DofDef {
    Ux,
    Uy,
    Uz,
    T,
}

/// BC ou Load : une valeur sur un DOF de chaque Node d'un Set (force par Node, pas force totale).
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ValueDef {
    pub set: String,
    pub dof: DofDef,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeatSourceDef {
    pub set: String,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConvectionDef {
    pub set: String,
    pub film: f64,
    pub sink: f64,
}

/// Pression sur les faces d'un bord (statique) : positive, elle comprime ; négative, elle tire.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PressureDef {
    pub set: String,
    pub value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ContactKind {
    Plane,
    Cylinder,
    Sphere,
}

/// Contact des Nodes d'un bord contre une surface rigide. Les vecteurs ont 2 composantes (maillage 2D)
/// ou 3 ; seules les clés de `kind` sont admises :
/// - `plane` : `point` et `normal` (vers le corps) ;
/// - `cylinder` : `point`, `radius` et, en 3D, `axis` ; en 2D, c'est un cercle ;
/// - `sphere` : `center` et `radius` (3D).
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContactDef {
    pub set: String,
    pub kind: ContactKind,
    pub point: Option<Vec<f64>>,
    pub normal: Option<Vec<f64>>,
    pub axis: Option<Vec<f64>>,
    pub center: Option<Vec<f64>>,
    pub radius: Option<f64>,
    /// Déplacement de la surface rigide quand λ(t) = 1, nul par défaut.
    pub motion: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StepKind {
    Transient,
    QuasiStatic,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StepDef {
    pub kind: StepKind,
    pub duration: f64,
    pub increments: usize,
    /// Transitoire seulement, 1 par défaut (Euler implicite).
    pub theta: Option<f64>,
    /// Transitoire seulement, obligatoire : température uniforme au départ.
    pub initial_temperature: Option<f64>,
    /// Quasi-statique seulement : points (t, λ), rampe de 0 à 1 sur la durée par défaut.
    pub amplitude: Option<Vec<[f64; 2]>>,
}

/// `[temperatures]` : champ de températures lu dans un fichier de résultats thermiques (`.pvd`).
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemperaturesDef {
    /// Fichier `.pvd` du calcul thermique, sur le même maillage.
    pub from: String,
    /// `"last"` (par défaut) ou l'indice de l'incrément, depuis 0.
    pub increment: Option<IncrementDef>,
    /// T_ref : température sans contrainte, dans la même unité que les températures du fichier.
    pub reference: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(untagged)]
pub enum IncrementDef {
    Index(usize),
    Name(String),
}
