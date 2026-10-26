//! Fichiers de Maille : dessin DXF 2D (`dxf`), maillage Gmsh `.msh` 4.1 ASCII (`msh`), modèle TOML (`model`),
//! résultats VTK `.vtu` et `.pvd` (`vtk`).
//!
//! Ce crate travaille sur des `&str` et des `String` : lire et écrire les fichiers sur le disque
//! est l'affaire du CLI. Le cœur (`maille-core`) reste sans I/O.

pub mod dxf;
mod error;
pub mod model;
pub mod msh;
pub mod vtk;

pub use error::IoError;
