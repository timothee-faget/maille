//! Mailleur 2D : triangulation de Delaunay (Bowyer-Watson), bords respectés, raffinement de Ruppert.

// utilisés par le mailleur à partir de l'étape 3 (bords respectés)
#[allow(dead_code)]
mod delaunay;
#[allow(dead_code)]
mod geometry;
#[allow(dead_code)]
mod predicates;

pub use geometry::{Curve, CurveId, CurveKind, Geometry, MIN_CORNER};

/// Plafond de Nodes produits par le mailleur : au-delà, `MeshingFailed` plutôt qu'un calcul sans fin.
pub const MAX_NODES: usize = 1_000_000;
