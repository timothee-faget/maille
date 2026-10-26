//! Fichiers de Maille : maillage Gmsh `.msh` 4.1 ASCII (`msh`), modèle TOML (`model`).
//!
//! Ce crate travaille sur des `&str` et des `String` : lire et écrire les fichiers sur le disque
//! est l'affaire du CLI. Le cœur (`maille-core`) reste sans I/O.

mod error;
pub mod model;
pub mod msh;

pub use error::IoError;
