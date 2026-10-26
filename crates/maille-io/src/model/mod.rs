//! Fichier de modèle TOML.
//!
//! Deux temps : `parse` lit le TOML (`ModelFile`, dont le nom du `.msh`), puis `build` construit le
//! `Model` du cœur à partir d'un `MshFile` déjà lu. Lire et écrire les fichiers est l'affaire du
//! CLI, qui reçoit aussi le `Run` : quel calcul lancer.

mod build;
mod schema;

use maille_core::model::{Amplitude, Model, Step};

use crate::IoError;

pub use build::build;
pub use schema::*;

/// Lit un fichier de modèle. Les erreurs de syntaxe donnent la ligne et la colonne ; une clé
/// inconnue est refusée.
pub fn parse(text: &str) -> Result<ModelFile, IoError> {
    toml::from_str(text).map_err(|e| IoError::Toml(e.to_string()))
}

/// Modèle construit, prêt à résoudre.
#[derive(Debug, Clone, PartialEq)]
pub struct Loaded {
    /// Unités du fichier, recopiées telles quelles dans les résultats.
    pub units: String,
    pub model: Model,
    pub run: Run,
}

/// Calcul demandé par le fichier.
#[derive(Debug, Clone, PartialEq)]
pub enum Run {
    /// Sans `[step]` : `Model::solve`.
    Linear,
    /// `Model::solve_step` depuis une température uniforme (un Node : une valeur).
    Transient { step: Step, initial: Vec<f64> },
    /// `Model::solve_quasi_static`.
    QuasiStatic { step: Step, amplitude: Amplitude },
}
