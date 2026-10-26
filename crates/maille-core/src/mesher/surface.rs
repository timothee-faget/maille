//! Surface fermée d'une pièce 3D : l'entrée du mailleur 3D. `maille_io::stl::read_surface` la
//! construit et la contrôle ; ici on ne garde que sa forme.

/// Surface fermée d'une pièce : sommets soudés et triangles dans le sens antihoraire vu de
/// l'extérieur (normale sortante).
#[derive(Debug, Clone, PartialEq)]
pub struct Surface {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[usize; 3]>,
    /// Morceaux nommés de la surface. Chaque triangle est dans exactement un patch.
    pub patches: Vec<Patch>,
}

/// Morceau de surface nommé : un `solid` du fichier, ou une zone lisse limitée par des arêtes
/// vives. Son nom devient le nom du bord (le Set où poser les BC).
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
