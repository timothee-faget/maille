//! Pièce STEP (ISO 10303) vers `Brep` : le fichier texte Part 21 (`part21`), puis la topologie
//! B-rep, les surfaces des faces et les courbes des arêtes (`brep`). La pièce se maille ensuite
//! avec `maille_core::mesher::tetrahedralize_brep`.
//!
//! Les unités du fichier sont reprises telles quelles, comme pour le STL : le cœur ne convertit
//! rien. `read_brep` donne l'unité de longueur lue, pour que l'utilisateur choisisse la taille h et
//! son modèle en conséquence.

mod brep;
mod part21;

use crate::IoError;
pub use maille_core::mesher::{BSpline, Brep, BrepEdge, BrepFace, EdgeCurve, FaceSurface};

/// Lit un fichier STEP texte (AP203, AP214, AP242) : un seul solide (MANIFOLD_SOLID_BREP), sa
/// coque fermée, ses faces, leurs surfaces (PLANE, CYLINDRICAL_SURFACE, CONICAL_SURFACE,
/// SPHERICAL_SURFACE, TOROIDAL_SURFACE) et leurs
/// contours, les arêtes et leurs courbes (LINE, CIRCLE, ELLIPSE, B_SPLINE_CURVE_WITH_KNOTS,
/// rationnelle ou non). Renvoie la pièce et l'unité de longueur du fichier (`mm`, `m`, `inch`…,
/// ou `inconnue`) ; les angles sont convertis en radians.
///
/// Refus avec message : plusieurs solides (assemblage), aucun solide, cavités internes, coque
/// ouverte ou mal orientée, surface ou courbe non prise en charge, sommet hors de la courbe de son
/// arête.
pub fn read_brep(text: &str) -> Result<(Brep, String), IoError> {
    brep::build(&part21::parse(text)?)
}

fn error(message: String) -> IoError {
    IoError::Step { message }
}
