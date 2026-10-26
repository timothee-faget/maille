//! Dessin DXF 2D (ASCII) vers `Geometry` du mailleur. Entités lues : LINE, ARC, CIRCLE et
//! LWPOLYLINE (arrondis compris). Le calque d'une entité donne le nom du bord, pour y poser des BC.
//!
//! Les unités du dessin sont reprises telles quelles : le cœur ne convertit rien.

mod geometry;
mod reader;

use maille_core::mesher::Geometry;

use crate::IoError;

/// Lit le dessin `text` et le transforme en `Geometry` de taille d'Element `size` sur tous les
/// bords. Le contour extérieur est le contour fermé de plus grande aire, les autres sont des trous.
pub fn read_geometry(text: &str, size: f64) -> Result<Geometry, IoError> {
    geometry::build(reader::read_pieces(text)?, size)
}

fn error(message: String) -> IoError {
    IoError::Dxf { message }
}
