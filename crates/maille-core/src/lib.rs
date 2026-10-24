//! Cœur de calcul de Maille : maillage, éléments, assemblage, solveur.
//!
//! Pur : aucune I/O. Unités : le cœur ne convertit rien, l'utilisateur fournit
//! un système d'unités cohérent (ex. N, mm, MPa).

pub mod assembly;
pub mod element;
pub mod error;
pub mod material;
pub mod mesh;
pub mod model;
pub mod quadrature;
pub mod results;
pub mod shape;
pub mod solve;

pub use error::MailleError;
