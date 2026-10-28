//! Pièce STEP (ISO 10303) vers `Brep` : le fichier texte Part 21 (`part21`), puis la topologie
//! B-rep et les courbes des arêtes (`brep`). À venir (M11) : les surfaces des faces, maillées face
//! par face.
//!
//! Les unités du fichier sont reprises telles quelles, comme pour le STL : le cœur ne convertit
//! rien. `read_brep` donne l'unité de longueur lue, pour que l'utilisateur choisisse la taille h et
//! son modèle en conséquence.

mod brep;
mod part21;

use crate::IoError;
pub use maille_core::mesher::{BSpline, Brep, BrepEdge, BrepFace, EdgeCurve};

/// Lit un fichier STEP texte (AP203, AP214, AP242) : un seul solide (MANIFOLD_SOLID_BREP), sa
/// coque fermée, ses faces et leurs contours, les arêtes et leurs courbes (LINE, CIRCLE, ELLIPSE,
/// B_SPLINE_CURVE_WITH_KNOTS, rationnelle ou non). Renvoie la pièce et l'unité de longueur du
/// fichier (`mm`, `m`, `inch`…, ou `inconnue`).
///
/// Refus avec message : plusieurs solides (assemblage), aucun solide, cavités internes, coque
/// ouverte ou mal orientée, courbe non prise en charge, sommet hors de la courbe de son arête.
pub fn read_brep(text: &str) -> Result<(Brep, String), IoError> {
    brep::build(&part21::parse(text)?)
}

fn error(message: String) -> IoError {
    IoError::Step { message }
}
