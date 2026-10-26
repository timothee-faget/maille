//! Surface STL (ASCII ou binaire) vers `Surface` : triangles soudés, fermés et orientés vers
//! l'extérieur, répartis en patches nommés. C'est l'entrée du futur mailleur 3D.
//!
//! Les unités du fichier sont reprises telles quelles : un STL n'en porte pas et le cœur ne
//! convertit rien. Le nom d'un patch deviendra le nom du Set où poser les BC.

mod reader;
mod surface;

use crate::IoError;

/// Surface fermée d'une pièce : sommets soudés et triangles dans le sens antihoraire vu de
/// l'extérieur (normale sortante).
#[derive(Debug, Clone, PartialEq)]
pub struct Surface {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
    /// Morceaux nommés de la surface. Chaque triangle est dans exactement un patch.
    pub patches: Vec<Patch>,
}

/// Morceau de surface nommé : un `solid` du fichier, ou une zone lisse limitée par des arêtes vives.
#[derive(Debug, Clone, PartialEq)]
pub struct Patch {
    pub name: String,
    /// Indices dans `Surface::triangles`, croissants.
    pub triangles: Vec<usize>,
}

impl Surface {
    /// Aire totale de la surface.
    pub fn area(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.vertices[i]);
                0.5 * norm(cross(sub(b, a), sub(c, a)))
            })
            .sum()
    }

    /// Volume enfermé par la surface (positif : les normales sortent).
    pub fn volume(&self) -> f64 {
        let o = self.vertices.first().copied().unwrap_or([0.0; 3]);
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| sub(self.vertices[i], o));
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }
}

/// Lit un fichier STL (ASCII ou binaire, reconnu à sa taille). Un fichier ASCII à plusieurs
/// `solid` donne un patch par nom ; sinon (un seul `solid`, ou binaire) la surface est découpée en
/// patches `surface_1`, `surface_2`… séparés par les arêtes dont l'angle entre les normales dépasse
/// `feature_angle` (en degrés, strictement entre 0 et 180).
///
/// Refus avec message : surface ouverte, arête partagée par plus de deux triangles, orientation
/// incohérente, triangle dégénéré, plusieurs pièces. Des normales rentrantes sur toute la surface
/// sont corrigées (le fichier n'est pas modifié, la `Surface` est retournée).
pub fn read_surface(bytes: &[u8], feature_angle: f64) -> Result<Surface, IoError> {
    if !(feature_angle > 0.0 && feature_angle < 180.0) {
        return Err(error(format!(
            "angle vif de {feature_angle}° hors de ]0, 180[ : 30 convient à la plupart des pièces"
        )));
    }
    surface::build(reader::read(bytes)?, feature_angle)
}

fn error(message: String) -> IoError {
    IoError::Stl { message }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
