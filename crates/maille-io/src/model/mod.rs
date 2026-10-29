//! Fichier de modèle TOML.
//!
//! Deux temps : `parse` lit le TOML (`ModelFile`, dont le nom du `.msh`), puis `build` construit le
//! `Model` du cœur à partir d'un `MshFile` déjà lu. Lire et écrire les fichiers est l'affaire du
//! CLI, qui reçoit aussi le `Run` : quel calcul lancer. Pour `[temperatures]`, le CLI lit le `.pvd`
//! et le `.vtu` (`TemperaturesDef::pick`, `vtk::read_point_data`) et passe les températures à `build`.

mod build;
mod physics;
mod schema;

use maille_core::mesh::Mesh;
use maille_core::model::{Amplitude, Model, Step};

use crate::vtk::Entry;
use crate::IoError;

pub use build::build;
pub use physics::{physics, split, Physics};
pub use schema::*;

/// Lit un fichier de modèle. Les erreurs de syntaxe donnent la ligne et la colonne ; une clé
/// inconnue est refusée.
pub fn parse(text: &str) -> Result<ModelFile, IoError> {
    toml::from_str(text).map_err(|e| IoError::Toml(e.to_string()))
}

/// Écrit un fichier de modèle : `parse(&write(&m)?)` redonne `m`. Les commentaires d'un fichier écrit
/// à la main ne survivent pas (un `ModelFile` n'en porte pas) ; l'ordre des clés est celui du schéma.
pub fn write(model: &ModelFile) -> Result<String, IoError> {
    toml::to_string(model).map_err(|e| IoError::Toml(e.to_string()))
}

/// Modèle construit, prêt à résoudre.
#[derive(Debug, Clone, PartialEq)]
pub struct Loaded {
    /// Unités du fichier, recopiées telles quelles dans les résultats.
    pub units: String,
    /// Le Mesh du modèle, nécessaire aux contraintes lissées des résultats (`vtk::result_arrays`).
    pub mesh: Mesh,
    pub model: Model,
    pub run: Run,
    /// `kinematics = "large"` : les résultats portent la marque des grandes déformations.
    pub large_strain: bool,
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

impl TemperaturesDef {
    /// Incrément du calcul thermique choisi par `increment` parmi ceux du `.pvd`.
    pub fn pick<'a>(&self, entries: &'a [Entry]) -> Result<&'a Entry, IoError> {
        let source = self.from.as_deref().unwrap_or("le calcul thermique");
        Ok(&entries[self.index(entries.len(), source)?])
    }

    /// Rang de l'Increment choisi parmi `count` Increments d'un calcul thermique, `source` nommant
    /// le fichier ou le calcul pour le message d'erreur.
    pub fn index(&self, count: usize, source: &str) -> Result<usize, IoError> {
        let last = count.saturating_sub(1);
        let index = match &self.increment {
            None => last,
            Some(IncrementDef::Name(n)) if n == "last" => last,
            Some(IncrementDef::Name(n)) => {
                return Err(IoError::InvalidModel {
                    problem: format!(
                        "[temperatures] increment = \"{n}\" : écrire \"last\" ou un indice entier, depuis 0"
                    ),
                })
            }
            Some(IncrementDef::Index(i)) => *i,
        };
        if index >= count {
            return Err(IoError::InvalidModel {
                problem: format!(
                    "[temperatures] increment = {index} : {source} contient {count} incrément(s), numérotés depuis 0"
                ),
            });
        }
        Ok(index)
    }
}
