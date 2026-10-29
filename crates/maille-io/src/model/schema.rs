//! Forme du fichier de modèle TOML : des données, aucune règle. La construction du `Model` est
//! dans `build`. Une clé inconnue ou mal écrite est une erreur (`deny_unknown_fields`).
//!
//! ```toml
//! units = "N, mm, MPa"            # obligatoire, texte libre, recopié dans les résultats
//! mesh  = "tube.msh"
//! analysis = "static"             # ou "thermal" ; facultative : sans elle, déduite du contenu
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
//! [meshing]                       # facultative : de quoi refaire le .msh (maille mesh modele.toml)
//! source = "piece.step"
//! size = 2.0
//! sets = { encastrement = ["face_3", "face_5"] }   # un Set = plusieurs entités de la géométrie
//!
//! [[bcs]]                         # Dirichlet sur un Set de Nodes
//! set = "sym_y"
//! dof = "ux"                      # ux, uy, uz (statique), t (thermique)
//! value = 0.0
//! ```

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// Fichier de modèle.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ModelFile {
    /// Unités cohérentes de l'utilisateur : le cœur ne convertit rien.
    pub units: String,
    /// Nom du fichier de maillage `.msh`, ouvert par le CLI.
    pub mesh: String,
    /// Facultative : sans elle, le calcul se déduit du contenu (voir `physics`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub analysis: Option<AnalysisKind>,
    /// Cinématique du calcul statique, petites déformations par défaut.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kinematics: Option<KinematicsDef>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub materials: BTreeMap<String, MaterialDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sections: Vec<SectionDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bcs: Vec<ValueDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub loads: Vec<ValueDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub heat_sources: Vec<HeatSourceDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub convections: Vec<ConvectionDef>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub pressures: Vec<PressureDef>,
    /// Contacts contre des surfaces rigides : demandent `[step] kind = "quasi_static"`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contacts: Vec<ContactDef>,
    /// Couplage thermo-mécanique : les températures d'un calcul thermique dilatent le modèle
    /// statique. Le CLI lit le fichier `from` et les passe à `build`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperatures: Option<TemperaturesDef>,
    /// Sans `[step]` : calcul linéaire.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub step: Option<StepDef>,
    /// Recette du maillage : source CAO, taille h, Sets sur la géométrie. Sans elle, `mesh` est un
    /// `.msh` fait ailleurs.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub meshing: Option<MeshingDef>,
    /// Calcul couplé seulement : le thermique transitoire qui précède le mécanique ; sans lui, le
    /// thermique est stationnaire.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub thermal_step: Option<ThermalStepDef>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum AnalysisKind {
    Static,
    Thermal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum KinematicsDef {
    Small,
    Large,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum HyperelasticityDef {
    NeoHookean,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialDef {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub young: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub poisson: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expansion: Option<f64>,
    /// Active la plasticité de von Mises.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub yield_stress: Option<f64>,
    /// Module d'écrouissage H, 0 par défaut (parfaitement plastique).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hardening: Option<f64>,
    /// Loi hyperélastique, grandes déformations seulement ; Saint-Venant–Kirchhoff sans elle.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hyperelasticity: Option<HyperelasticityDef>,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SectionDef {
    pub set: String,
    pub material: String,
    /// Triangles seulement, 1 par défaut.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub thickness: Option<f64>,
    /// Triangles seulement ; obligatoire en analyse statique.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plane_state: Option<PlaneStateDef>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum PlaneStateDef {
    Stress,
    Strain,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum DofDef {
    Ux,
    Uy,
    Uz,
    T,
}

/// BC ou Load : une valeur sur un DOF de chaque Node d'un Set (force par Node, pas force totale).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ValueDef {
    pub set: String,
    pub dof: DofDef,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct HeatSourceDef {
    pub set: String,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ConvectionDef {
    pub set: String,
    pub film: f64,
    pub sink: f64,
}

/// Pression sur les faces d'un bord (statique) : positive, elle comprime ; négative, elle tire.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PressureDef {
    pub set: String,
    pub value: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
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
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ContactDef {
    pub set: String,
    pub kind: ContactKind,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub point: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub normal: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub axis: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub center: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// Déplacement de la surface rigide quand λ(t) = 1, nul par défaut.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motion: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum StepKind {
    Transient,
    QuasiStatic,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StepDef {
    pub kind: StepKind,
    pub duration: f64,
    pub increments: usize,
    /// Transitoire seulement, 1 par défaut (Euler implicite).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theta: Option<f64>,
    /// Transitoire seulement, obligatoire : température uniforme au départ.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub initial_temperature: Option<f64>,
    /// Quasi-statique seulement : points (t, λ), rampe de 0 à 1 sur la durée par défaut.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub amplitude: Option<Vec<[f64; 2]>>,
}

/// `[meshing]` : de quoi refaire le maillage. Un Set peut regrouper plusieurs entités de la géométrie
/// (calque DXF, `solid` ou `surface_N` du STL, face du STEP) : il survit au remaillage, puisqu'il les
/// désigne par leur nom.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct MeshingDef {
    /// Fichier CAO (`.dxf`, `.stl`, `.step`), relatif au dossier du modèle.
    pub source: String,
    /// Taille h des Elements, dans les unités du fichier ; facultative pour un STL seulement.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub size: Option<f64>,
    /// Elements du second ordre (Tri6, Tet10), faux par défaut.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quadratic: Option<bool>,
    /// Angle vif d'un STL, en degrés (30 par défaut).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub feature_angle: Option<f64>,
    /// Nom du Set → noms des entités de la géométrie qu'il regroupe.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub sets: BTreeMap<String, Vec<String>>,
}

/// `[thermal_step]` : le transitoire thermique d'un calcul couplé (mêmes clés que `[step]` transitoire).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ThermalStepDef {
    pub duration: f64,
    pub increments: usize,
    /// 1 par défaut (Euler implicite).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theta: Option<f64>,
    /// Obligatoire : température uniforme au départ.
    pub initial_temperature: f64,
}

/// `[temperatures]` : champ de températures lu dans un fichier de résultats thermiques (`.pvd`).
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct TemperaturesDef {
    /// Fichier `.pvd` du calcul thermique, sur le même maillage. Sans `from`, le modèle a sa propre
    /// partie thermique (calcul couplé) et `increment` choisit parmi ses Increments.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub from: Option<String>,
    /// `"last"` (par défaut) ou l'indice de l'incrément, depuis 0.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub increment: Option<IncrementDef>,
    /// T_ref : température sans contrainte, dans la même unité que les températures du fichier.
    pub reference: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum IncrementDef {
    Index(usize),
    Name(String),
}
