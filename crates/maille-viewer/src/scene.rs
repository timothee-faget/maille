//! Ce qu'on dessine, calculé sans GPU : triangles et arêtes à tracer, valeur de chaque champ aux
//! Nodes, positions déformées.

use std::collections::HashSet;
use std::fmt;

use maille_core::mesh::{Element, Mesh};
use maille_core::results::{nodal_plastic_strains, nodal_stresses, Field, Stress};
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
    Temperature,
    /// Déformation plastique cumulée p.
    PlasticStrain,
}

impl Quantity {
    /// Champs d'un Field statique, dans l'ordre des touches 1 à 5.
    pub const MECHANICAL: [Self; 5] = [
        Self::VonMises,
        Self::SigmaXx,
        Self::SigmaYy,
        Self::TauXy,
        Self::Displacement,
    ];

    /// Champs d'un Field élasto-plastique : les champs mécaniques, puis p (touche 6).
    pub const PLASTIC: [Self; 6] = [
        Self::VonMises,
        Self::SigmaXx,
        Self::SigmaYy,
        Self::TauXy,
        Self::Displacement,
        Self::PlasticStrain,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::VonMises => "von Mises",
            Self::SigmaXx => "σxx",
            Self::SigmaYy => "σyy",
            Self::TauXy => "τxy",
            Self::Displacement => "|u|",
            Self::Temperature => "T",
            Self::PlasticStrain => "εp (déformation plastique cumulée)",
        }
    }
}

pub struct Scene {
    nodes: Vec<[f64; 2]>,
    /// (ux, uy) par Node.
    displacements: Vec<[f64; 2]>,
    /// Vide pour un Field thermique.
    stresses: Vec<Stress>,
    /// Vide pour un Field statique.
    temperatures: Vec<f64>,
    /// Déformation plastique cumulée p aux Nodes, vide hors calcul élasto-plastique.
    plastic_strains: Vec<f64>,
    /// Triangles à remplir (indices de Nodes, antihoraires) : 1 par Tri3, 4 par Tri6.
    pub triangles: Vec<[u32; 3]>,
    /// Segments du contour des Elements, sans doublon : les arêtes internes des Tri6 découpés
    /// n'y sont pas, on voit le vrai Mesh.
    pub edges: Vec<[u32; 2]>,
}

/// Raison pour laquelle un Mesh ne peut pas être affiché.
#[derive(Debug, Clone, PartialEq)]
pub enum SceneError {
    /// Refus du cœur : Mesh de Bar ou de Beam, ou Field d'un autre Mesh.
    Core(MailleError),
    /// Mesh de tétraèdres : le viewer n'affiche que des triangles pour l'instant.
    Solid,
}

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Core(e) => e.fmt(f),
            Self::Solid => write!(
                f,
                "le viewer n'affiche que des Mesh 2D (Tri3, Tri6) pour l'instant : lire les résultats des tétraèdres dans le Field"
            ),
        }
    }
}

impl std::error::Error for SceneError {}

impl From<MailleError> for SceneError {
    fn from(e: MailleError) -> Self {
        Self::Core(e)
    }
}

impl Scene {
    /// Erreur si le Mesh n'est pas fait de triangles (Bar, Beam, tétraèdres) ou si le Field vient d'un autre Mesh.
    pub fn new(mesh: &Mesh, field: &Field) -> Result<Self, SceneError> {
        let thermal = !field.temperatures.is_empty();
        if thermal && field.temperatures.len() != mesh.nodes.len() {
            return Err(MailleError::FieldMeshMismatch.into());
        }
        let stresses = if thermal {
            Vec::new()
        } else {
            nodal_stresses(mesh, field)?
        };
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
                Element::Bar(_) | Element::Beam(_) => return Err(MailleError::NoStresses.into()),
                Element::Tet4(_) | Element::Tet10(_) => return Err(SceneError::Solid),
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
            displacements: if thermal {
                vec![[0.0; 2]; mesh.nodes.len()]
            } else {
                field.displacements.as_chunks::<2>().0.to_vec()
            },
            stresses,
            temperatures: field.temperatures.clone(),
            plastic_strains: if field.equivalent_plastic_strains.is_empty() {
                Vec::new()
            } else {
                nodal_plastic_strains(mesh, field)?
            },
            triangles,
            edges,
        })
    }

    /// Champs affichables : la température pour un Field thermique, sinon les champs mécaniques,
    /// plus p pour un Field élasto-plastique.
    pub fn quantities(&self) -> &'static [Quantity] {
        if !self.temperatures.is_empty() {
            &[Quantity::Temperature]
        } else if !self.plastic_strains.is_empty() {
            &Quantity::PLASTIC
        } else {
            &Quantity::MECHANICAL
        }
    }

    /// Vrai si le Field a des déplacements à amplifier (statique).
    pub fn has_displacements(&self) -> bool {
        self.temperatures.is_empty()
    }

    /// Valeur du champ à chaque Node (vide si le champ n'est pas dans `quantities`).
    pub fn values(&self, quantity: Quantity) -> Vec<f64> {
        let stress = |f: fn(&Stress) -> f64| self.stresses.iter().map(f).collect();
        match quantity {
            Quantity::VonMises => stress(Stress::von_mises),
            Quantity::SigmaXx => stress(|s| s.sigma_xx),
            Quantity::SigmaYy => stress(|s| s.sigma_yy),
            Quantity::TauXy => stress(|s| s.tau_xy),
            Quantity::Displacement => self
                .displacements
                .iter()
                .map(|u| u[0].hypot(u[1]))
                .collect(),
            Quantity::Temperature => self.temperatures.clone(),
            Quantity::PlasticStrain => self.plastic_strains.clone(),
        }
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
    use maille_core::mesh::{Bar, Node, PlaneState, Tet4, Tri3, Tri6};

    fn material() -> Material {
        Material::with_poisson(1.0, 0.3).unwrap()
    }

    /// Field fabriqué à la main : contrainte uniforme (σxx, σyy, τxy) et déplacements donnés.
    fn field(mesh: &Mesh, sigma: [f64; 3], displacements: Vec<f64>) -> Field {
        let stress = Stress {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            sigma_xx: sigma[0],
            sigma_yy: sigma[1],
            sigma_zz: 0.0,
            tau_xy: sigma[2],
            tau_yz: 0.0,
            tau_xz: 0.0,
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
            equivalent_plastic_strains: Vec::new(),
            temperatures: Vec::new(),
            heat_fluxes: Vec::new(),
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
    fn a_plastic_field_adds_the_plastic_strain() {
        let mesh = two_tri3();
        let elastic = Scene::new(&mesh, &field(&mesh, [100.0, 0.0, 0.0], vec![0.0; 8])).unwrap();
        assert_eq!(elastic.quantities(), &Quantity::MECHANICAL);
        // p = 0,01 dans le premier Tri3, 0 dans le second : moyenne 0,005 sur l'arête partagée 0-2
        let plastic = Field {
            equivalent_plastic_strains: vec![vec![0.01], vec![0.0]],
            ..field(&mesh, [100.0, 0.0, 0.0], vec![0.0; 8])
        };
        let scene = Scene::new(&mesh, &plastic).unwrap();
        assert_eq!(scene.quantities(), &Quantity::PLASTIC);
        let p = scene.values(Quantity::PlasticStrain);
        for (got, expected) in p.iter().zip([0.005, 0.01, 0.005, 0.0]) {
            assert!((got - expected).abs() < 1e-15, "{p:?}");
        }
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
            equivalent_plastic_strains: Vec::new(),
            temperatures: Vec::new(),
            heat_fluxes: Vec::new(),
        };
        assert_eq!(
            Scene::new(&mesh, &field).err(),
            Some(SceneError::Core(MailleError::NoStresses))
        );
    }

    #[test]
    fn a_thermal_field_shows_only_the_temperature_without_deformation() {
        let mesh = two_tri3();
        let temperatures = vec![10.0, 20.0, 30.0, 40.0];
        let field = Field {
            displacements: Vec::new(),
            reactions: vec![0.0; 4],
            stresses: Vec::new(),
            equivalent_plastic_strains: Vec::new(),
            temperatures: temperatures.clone(),
            heat_fluxes: Vec::new(),
        };
        let scene = Scene::new(&mesh, &field).unwrap();
        assert_eq!(scene.quantities(), &[Quantity::Temperature]);
        assert!(!scene.has_displacements());
        for (v, t) in scene
            .values(Quantity::Temperature)
            .iter()
            .zip(&temperatures)
        {
            assert!((v - t).abs() < 1e-12);
        }
        assert!(scene.deformation_scale().abs() < 1e-12);

        // Field thermique d'un autre Mesh (3 températures pour 4 Nodes)
        let foreign = Field {
            temperatures: vec![1.0; 3],
            ..field
        };
        assert_eq!(
            Scene::new(&mesh, &foreign).err(),
            Some(SceneError::Core(MailleError::FieldMeshMismatch))
        );
    }

    #[test]
    fn tetrahedra_cannot_be_shown() {
        let nodes = vec![
            Node::new(0.0, 0.0, 0.0),
            Node::new(1.0, 0.0, 0.0),
            Node::new(0.0, 1.0, 0.0),
            Node::new(0.0, 0.0, 1.0),
        ];
        let tet = Tet4::new([0, 1, 2, 3], material());
        let mesh = Mesh::new(nodes, vec![Element::Tet4(tet)]).unwrap();
        let field = field(&mesh, [0.0; 3], vec![0.0; 12]);
        assert_eq!(Scene::new(&mesh, &field).err(), Some(SceneError::Solid));
    }
}
