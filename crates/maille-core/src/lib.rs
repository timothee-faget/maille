//! Cœur de calcul de Maille : maillage, éléments, assemblage, solveur.
//!
//! Pur : aucune I/O. Unités : le cœur ne convertit rien, l'utilisateur fournit
//! un système d'unités cohérent. Exemples :
//! - mécanique : N, mm, MPa ;
//! - thermique : W, mm, K (ou °C), d'où k en W/(mm·K), source Q en W/mm³, coefficient
//!   d'échange h en W/(mm²·K) ;
//! - thermique transitoire : en plus, le temps en s, donc ρ en kg/mm³ et c en J/(kg·K) (acier :
//!   ρ = 7,85·10⁻⁶ kg/mm³, c = 460 J/(kg·K)).
//!
//! - couplage thermo-mécanique : α en 1/K, et T dans la même unité que T_ref (K ou °C).
//!
//! Un même Material sert aux deux analyses (E, ν, α et k) : N et mm donnent des énergies en mJ,
//! donc des puissances en mW si le temps est en s. Mélanger MPa et W/(mm·K) reste cohérent : en
//! couplage, le champ T sort du calcul thermique et sert de donnée au calcul statique, et seul α
//! (en 1/K) relie les deux systèmes d'unités, chacun cohérent de son côté.

pub mod assembly;
pub mod element;
pub mod error;
pub mod material;
pub mod mesh;
pub mod mesher;
pub mod model;
pub mod plasticity;
pub mod quadrature;
pub mod results;
pub mod shape;
pub mod solve;

pub use error::MailleError;
