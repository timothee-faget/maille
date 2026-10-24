//! Ce qu'on dessine, calculé sans GPU : triangles et arêtes à tracer, valeur de chaque champ aux
//! Nodes, positions déformées.

use std::collections::HashSet;

use maille_core::mesh::{Element, Mesh};
use maille_core::results::{nodal_stresses, Field, PlaneStress};
use maille_core::MailleError;

/// Part de la taille du Mesh occupée par le plus grand déplacement, dans la déformée amplifiée.
const DEFORMATION_RATIO: f64 = 0.1;

/// Champ affiché, une valeur par Node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantity {
    VonMises,
    SigmaXx,
    SigmaYy,
    TauXy,
    /// Norme du déplacement |u|.
    Displacement,
}

impl Quantity {
    pub const ALL: [Self; 5] = [
        Self::VonMises,
        Self::SigmaXx,
        Self::SigmaYy,
        Self::TauXy,
        Self::Displacement,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::VonMises => "von Mises",
            Self::SigmaXx => "σxx",
            Self::SigmaYy => "σyy",
            Self::TauXy => "τxy",
            Self::Displacement => "|u|",
        }
    }
}

pub struct Scene {
    nodes: Vec<[f64; 2]>,
    /// (ux, uy) par Node.
    displacements: Vec<[f64; 2]>,
    stresses: Vec<PlaneStress>,
    /// Triangles à remplir (indices de Nodes, antihoraires) : 1 par Tri3, 4 par Tri6.
    pub triangles: Vec<[u32; 3]>,
    /// Segments du contour des Elements, sans doublon : les arêtes internes des Tri6 découpés
    /// n'y sont pas, on voit le vrai Mesh.
    pub edges: Vec<[u32; 2]>,
}

impl Scene {
    /// Erreur si le Mesh n'est pas fait de triangles (Bar, Beam) ou si le Field vient d'un autre Mesh.
    pub fn new(mesh: &Mesh, field: &Field) -> Result<Self, MailleError> {
        let stresses = nodal_stresses(mesh, field)?;
        let mut triangles = Vec::new();
        let mut edges = HashSet::new();
        for e in &mesh.elements {
            let n: Vec<u32> = e.nodes().iter().map(|&i| i as u32).collect();
            let (tris, outline): (&[[usize; 3]], &[[usize; 2]]) = match e {
                Element::Tri3(_) => (&[[0, 1, 2]], &[[0, 1], [1, 2], [2, 0]]),
                // sommets 0, 1, 2 puis milieux 3 (0-1), 4 (1-2), 5 (2-0)
                Element::Tri6(_) => (
                    &[[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]],
                    &[[0, 3], [3, 1], [1, 4], [4, 2], [2, 5], [5, 0]],
                ),
                // déjà refusé par `nodal_stresses`
                Element::Bar(_) | Element::Beam(_) => return Err(MailleError::NoPlaneStresses),
            };
            triangles.extend(tris.iter().map(|t| t.map(|k| n[k])));
            edges.extend(
                outline
                    .iter()
                    .map(|&[a, b]| [n[a].min(n[b]), n[a].max(n[b])]),
            );
        }
        let mut edges: Vec<[u32; 2]> = edges.into_iter().collect();
        edges.sort_unstable();
        Ok(Self {
            nodes: mesh.nodes.iter().map(|p| [p.x, p.y]).collect(),
            displacements: field.displacements.as_chunks::<2>().0.to_vec(),
            stresses,
            triangles,
            edges,
        })
    }

    /// Valeur du champ à chaque Node.
    pub fn values(&self, quantity: Quantity) -> Vec<f64> {
        self.stresses
            .iter()
            .zip(&self.displacements)
            .map(|(s, u)| match quantity {
                Quantity::VonMises => s.von_mises(),
                Quantity::SigmaXx => s.sigma_xx,
                Quantity::SigmaYy => s.sigma_yy,
                Quantity::TauXy => s.tau_xy,
                Quantity::Displacement => u[0].hypot(u[1]),
            })
            .collect()
    }

    /// Coins (min, max) de la boîte englobante des positions x + scale·u (scale = 0 : non déformé).
    pub fn bounds(&self, scale: f64) -> ([f64; 2], [f64; 2]) {
        self.positions(scale).iter().fold(
            ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]),
            |(lo, hi), p| {
                (
                    [lo[0].min(p[0]), lo[1].min(p[1])],
                    [hi[0].max(p[0]), hi[1].max(p[1])],
                )
            },
        )
    }

    /// Facteur d'amplification tel que le plus grand déplacement vaille 10 % de la taille du Mesh.
    /// Nul si rien ne bouge.
    pub fn deformation_scale(&self) -> f64 {
        let (lo, hi) = self.bounds(0.0);
        let size = (hi[0] - lo[0]).max(hi[1] - lo[1]);
        let max_u = self
            .values(Quantity::Displacement)
            .into_iter()
            .fold(0.0, f64::max);
        if max_u > 0.0 {
            DEFORMATION_RATIO * size / max_u
        } else {
            0.0
        }
    }

    /// Positions x + scale·u des Nodes (scale = 0 : Mesh non déformé).
    pub fn positions(&self, scale: f64) -> Vec<[f64; 2]> {
        self.nodes
            .iter()
            .zip(&self.displacements)
            .map(|(p, u)| [p[0] + scale * u[0], p[1] + scale * u[1]])
            .collect()
    }
}

/// (min, max) des valeurs ; un champ constant donne un intervalle nul, à gérer par l'appelant.
pub fn range(values: &[f64]) -> (f64, f64) {
    values
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_core::material::{Material, Section};
    use maille_core::mesh::{Bar, Node, PlaneState, Tri3, Tri6};

    fn material() -> Material {
        Material::with_poisson(1.0, 0.3).unwrap()
    }

    /// Field fabriqué à la main : contrainte uniforme (σxx, σyy, τxy) et déplacements donnés.
    fn field(mesh: &Mesh, sigma: [f64; 3], displacements: Vec<f64>) -> Field {
        let stress = PlaneStress {
            x: 0.0,
            y: 0.0,
            sigma_xx: sigma[0],
            sigma_yy: sigma[1],
            tau_xy: sigma[2],
            sigma_zz: 0.0,
        };
        let stresses = mesh
            .elements
            .iter()
            .map(|e| match e {
                Element::Tri6(_) => vec![stress; 3],
                _ => vec![stress],
            })
            .collect();
        Field {
            reactions: vec![0.0; displacements.len()],
            displacements,
            stresses,
        }
    }

    /// Carré unité coupé en deux Tri3 qui partagent l'arête 0-2.
    fn two_tri3() -> Mesh {
        let nodes = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            .map(|(x, y)| Node::new_2d(x, y))
            .to_vec();
        let tri = |n| Element::Tri3(Tri3::new(n, material(), 1.0, PlaneState::Stress).unwrap());
        Mesh::new(nodes, vec![tri([0, 1, 2]), tri([0, 2, 3])]).unwrap()
    }

    fn one_tri6() -> Mesh {
        let nodes = [
            (0.0, 0.0),
            (2.0, 0.0),
            (0.0, 2.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        .map(|(x, y)| Node::new_2d(x, y))
        .to_vec();
        let tri6 = Tri6::new([0, 1, 2, 3, 4, 5], material(), 1.0, PlaneState::Stress).unwrap();
        Mesh::new(nodes, vec![Element::Tri6(tri6)]).unwrap()
    }

    fn signed_area(p: &[[f64; 2]], t: [u32; 3]) -> f64 {
        let [a, b, c] = t.map(|i| p[i as usize]);
        ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0
    }

    #[test]
    fn shared_edges_are_drawn_once() {
        let mesh = two_tri3();
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], vec![0.0; 8])).unwrap();
        assert_eq!(scene.triangles, vec![[0, 1, 2], [0, 2, 3]]);
        assert_eq!(scene.edges, vec![[0, 1], [0, 2], [0, 3], [1, 2], [2, 3]]);
    }

    #[test]
    fn tri6_is_split_into_four_counterclockwise_triangles_covering_it() {
        let mesh = one_tri6();
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], vec![0.0; 12])).unwrap();
        let p = scene.positions(0.0);
        assert_eq!(scene.triangles.len(), 4);
        for &t in &scene.triangles {
            assert!(
                (signed_area(&p, t) - 0.5).abs() < 1e-12,
                "{t:?} : aire 2/4, antihoraire"
            );
        }
        // contour : 3 arêtes coupées en 2 au Node milieu, pas d'arête interne (3-4, 4-5, 5-3)
        assert_eq!(scene.edges.len(), 6);
        assert!(!scene.edges.contains(&[3, 4]) && !scene.edges.contains(&[4, 5]));
    }

    #[test]
    fn values_are_the_nodal_stresses_and_displacement_norms() {
        let mesh = two_tri3();
        let u = vec![3.0, 4.0, 0.0, 0.0, -6.0, 8.0, 0.0, 1.0];
        let scene = Scene::new(&mesh, &field(&mesh, [100.0, 0.0, 50.0], u)).unwrap();
        let von_mises = (100.0_f64.powi(2) + 3.0 * 50.0_f64.powi(2)).sqrt();
        for v in scene.values(Quantity::VonMises) {
            assert!((v - von_mises).abs() < 1e-9);
        }
        assert!(scene
            .values(Quantity::SigmaXx)
            .iter()
            .all(|v| (v - 100.0).abs() < 1e-9));
        assert!(scene
            .values(Quantity::SigmaYy)
            .iter()
            .all(|v| v.abs() < 1e-9));
        assert!(scene
            .values(Quantity::TauXy)
            .iter()
            .all(|v| (v - 50.0).abs() < 1e-9));
        let norms = scene.values(Quantity::Displacement);
        for (got, expected) in norms.iter().zip([5.0, 0.0, 10.0, 1.0]) {
            assert!((got - expected).abs() < 1e-12, "{norms:?}");
        }
        let (lo, hi) = range(&norms);
        assert!(lo.abs() < 1e-12 && (hi - 10.0).abs() < 1e-12);
    }

    #[test]
    fn deformed_shape_moves_the_largest_displacement_by_a_tenth_of_the_mesh() {
        let mesh = two_tri3();
        let u = vec![0.0, 0.0, 1e-3, 0.0, 2e-3, -1e-3, 0.0, 0.0];
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], u)).unwrap();
        let scale = scene.deformation_scale();
        let (before, after) = (scene.positions(0.0), scene.positions(scale));
        let moved = before
            .iter()
            .zip(&after)
            .map(|(a, b)| (b[0] - a[0]).hypot(b[1] - a[1]))
            .fold(0.0, f64::max);
        assert!((moved - 0.1).abs() < 1e-12, "{moved}");
        let (lo, hi) = scene.bounds(0.0);
        assert!(lo.iter().all(|c| c.abs() < 1e-12));
        assert!(hi.iter().all(|c| (c - 1.0).abs() < 1e-12));
        // la boîte de la déformée englobe le Node 2 déplacé de 0,1 (u ∝ (2, −1))
        let (lo, hi) = scene.bounds(scale);
        let moved = [1.0 + 0.2 / 5_f64.sqrt(), 1.0 - 0.1 / 5_f64.sqrt()];
        assert!((hi[0] - moved[0]).abs() < 1e-12 && lo[1].abs() < 1e-12);
    }

    #[test]
    fn a_mesh_that_does_not_move_has_no_deformation() {
        let mesh = two_tri3();
        let scene = Scene::new(&mesh, &field(&mesh, [0.0; 3], vec![0.0; 8])).unwrap();
        assert!(scene.deformation_scale().abs() < 1e-12);
    }

    #[test]
    fn bars_cannot_be_shown() {
        let nodes = vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)];
        let bar = Bar {
            nodes: [0, 1],
            material: Material::new(1.0).unwrap(),
            section: Section::new(1.0).unwrap(),
        };
        let mesh = Mesh::new(nodes, vec![Element::Bar(bar)]).unwrap();
        let field = Field {
            displacements: vec![0.0; 2],
            reactions: vec![0.0; 2],
            stresses: vec![Vec::new()],
        };
        assert_eq!(
            Scene::new(&mesh, &field).err(),
            Some(MailleError::NoPlaneStresses)
        );
    }
}
