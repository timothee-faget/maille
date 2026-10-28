//! Fichiers de Maille : dessin DXF 2D (`dxf`), maillage Gmsh `.msh` 4.1 ASCII (`msh`), modèle TOML (`model`),
//! résultats VTK `.vtu` et `.pvd` (`vtk`), surface STL 3D (`stl`), pièce STEP (`step`, en cours : M11).
//!
//! Ce crate travaille sur des `&str` et des `String` : lire et écrire les fichiers sur le disque
//! est l'affaire du CLI. Le cœur (`maille-core`) reste sans I/O.

pub mod dxf;
mod error;
pub mod model;
pub mod msh;
#[cfg_attr(not(test), allow(dead_code))] // lu par la topologie B-rep (M11, PR 2)
mod step;
pub mod stl;
pub mod vtk;

pub use error::IoError;
