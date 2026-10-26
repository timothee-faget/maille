//! Forme du fichier de modèle TOML : des données, aucune règle. La construction du `Model` est
//! dans `build`. Une clé inconnue ou mal écrite est une erreur (`deny_unknown_fields`).
//!
//! ```toml
//! units = "N, mm, MPa"            # obligatoire, texte libre, recopié dans les résultats
//! mesh  = "tube.msh"
//! analysis = "static"             # ou "thermal"
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
    /// Couplage thermo-mécanique : pas encore pris en charge, refusé par `build`.
    pub temperatures: Option<toml::Value>,
    /// Sans `[step]` : calcul linéaire.
    pub step: Option<StepDef>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AnalysisKind {
    Static,
    Thermal,
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
