use crate::error::MailleError;
use crate::material::{Material, Section};
use crate::shape::{TET10_EDGES, TRI6_EDGES};

/// Position d'un Node. Les Bar et Beam (1D) n'utilisent que `x`, les triangles (2D) `x` et `y`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Node {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Node du plan z = 0, pour les Mesh 1D et 2D.
    pub fn new_2d(x: f64, y: f64) -> Self {
        Self::new(x, y, 0.0)
    }
}

/// Degré de liberté d'un Node : U déplacement axial (Bar), V flèche et Theta rotation (Beam),
/// Ux et Uy déplacements dans le plan (triangles), Ux, Uy et Uz (tétraèdres).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dof {
    U,
    V,
    Theta,
    Ux,
    Uy,
    Uz,
}

/// Élément barre entre deux Nodes (indices dans `Mesh::nodes`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bar {
    pub nodes: [usize; 2],
    pub material: Material,
    pub section: Section,
}

/// Poutre d'Euler-Bernoulli (flexion seule) entre deux Nodes, `nodes[0]` à gauche de `nodes[1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam {
    pub nodes: [usize; 2],
    pub material: Material,
    /// Moment quadratique I.
    pub inertia: f64,
}

impl Beam {
    pub fn new(nodes: [usize; 2], material: Material, inertia: f64) -> Result<Self, MailleError> {
        Ok(Self {
            nodes,
            material,
            inertia: crate::material::positive("I (moment quadratique)", inertia)?,
        })
    }
}

/// Hypothèse 2D : contraintes planes (plaque mince) ou déformations planes (corps long).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneState {
    Stress,
    Strain,
}

/// Triangle à 3 Nodes (déformation constante), Nodes dans le sens antihoraire.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tri3 {
    pub nodes: [usize; 3],
    pub material: Material,
    pub thickness: f64,
    pub state: PlaneState,
}

impl Tri3 {
    pub fn new(
        nodes: [usize; 3],
        material: Material,
        thickness: f64,
        state: PlaneState,
    ) -> Result<Self, MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        Ok(Self {
            nodes,
            material,
            thickness: crate::material::positive("t (épaisseur)", thickness)?,
            state,
        })
    }
}

/// Triangle quadratique à 6 Nodes : sommets 0, 1, 2 (sens antihoraire), puis milieux d'arêtes
/// 3 (0-1), 4 (1-2), 5 (2-0). Côtés droits seulement : chaque Node milieu est au milieu de son arête.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tri6 {
    pub nodes: [usize; 6],
    pub material: Material,
    pub thickness: f64,
    pub state: PlaneState,
}

impl Tri6 {
    pub fn new(
        nodes: [usize; 6],
        material: Material,
        thickness: f64,
        state: PlaneState,
    ) -> Result<Self, MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        Ok(Self {
            nodes,
            material,
            thickness: crate::material::positive("t (épaisseur)", thickness)?,
            state,
        })
    }
}

/// Tétraèdre à 4 Nodes (déformation constante). Vu depuis le Node 3, la face 0-1-2 tourne dans le
/// sens antihoraire.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tet4 {
    pub nodes: [usize; 4],
    pub material: Material,
}

impl Tet4 {
    pub fn new(nodes: [usize; 4], material: Material) -> Result<Self, MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        Ok(Self { nodes, material })
    }
}

/// Tétraèdre quadratique à 10 Nodes : sommets 0 à 3 (orientés comme le Tet4), puis milieux d'arêtes
/// 4 (0-1), 5 (1-2), 6 (2-0), 7 (0-3), 8 (1-3), 9 (2-3). Arêtes droites seulement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tet10 {
    pub nodes: [usize; 10],
    pub material: Material,
}

impl Tet10 {
    pub fn new(nodes: [usize; 10], material: Material) -> Result<Self, MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        Ok(Self { nodes, material })
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Element {
    Bar(Bar),
    Beam(Beam),
    Tri3(Tri3),
    Tri6(Tri6),
    Tet4(Tet4),
    Tet10(Tet10),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Family {
    Bar,
    Beam,
    Plane,
    Solid,
}

impl Element {
    pub fn nodes(&self) -> &[usize] {
        match self {
            Self::Bar(e) => &e.nodes,
            Self::Beam(e) => &e.nodes,
            Self::Tri3(e) => &e.nodes,
            Self::Tri6(e) => &e.nodes,
            Self::Tet4(e) => &e.nodes,
            Self::Tet10(e) => &e.nodes,
        }
    }

    fn family(&self) -> Family {
        match self {
            Self::Bar(_) => Family::Bar,
            Self::Beam(_) => Family::Beam,
            Self::Tri3(_) | Self::Tri6(_) => Family::Plane,
            Self::Tet4(_) | Self::Tet10(_) => Family::Solid,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    pub nodes: Vec<Node>,
    pub elements: Vec<Element>,
}

impl Mesh {
    pub fn new(nodes: Vec<Node>, elements: Vec<Element>) -> Result<Self, MailleError> {
        if nodes.is_empty() || elements.is_empty() {
            return Err(MailleError::EmptyMesh);
        }
        let mesh = Self { nodes, elements };
        let family = mesh.elements[0].family();
        for (i, e) in mesh.elements.iter().enumerate() {
            if e.family() != family {
                return Err(MailleError::MixedElementTypes);
            }
            for &n in e.nodes() {
                if n >= mesh.nodes.len() {
                    return Err(MailleError::UnknownNode { node: n });
                }
            }
            match e {
                Element::Bar(_) if mesh.length(i) == 0.0 => {
                    return Err(MailleError::ZeroLengthElement { element: i });
                }
                Element::Beam(b) => {
                    if mesh.length(i) == 0.0 {
                        return Err(MailleError::ZeroLengthElement { element: i });
                    }
                    if mesh.nodes[b.nodes[1]].x < mesh.nodes[b.nodes[0]].x {
                        return Err(MailleError::ReversedBeam { element: i });
                    }
                }
                Element::Tri3(t) => mesh.check_triangle(i, &t.nodes, &t.material)?,
                Element::Tri6(t) => {
                    mesh.check_triangle(i, &t.nodes, &t.material)?;
                    if !straight_edges(&mesh.coords(&t.nodes), &TRI6_EDGES, 3) {
                        return Err(MailleError::CurvedTriangle { element: i });
                    }
                }
                Element::Tet4(t) => mesh.check_tetrahedron(i, &t.nodes, &t.material)?,
                Element::Tet10(t) => {
                    mesh.check_tetrahedron(i, &t.nodes, &t.material)?;
                    if !straight_edges(&mesh.coords_3d(&t.nodes), &TET10_EDGES, 4) {
                        return Err(MailleError::CurvedTetrahedron { element: i });
                    }
                }
                _ => {}
            }
        }
        Ok(mesh)
    }

    /// Vérifie ν, l'aire et l'orientation d'un triangle ; seuls ses 3 premiers Nodes (sommets) comptent.
    fn check_triangle(
        &self,
        element: usize,
        nodes: &[usize],
        material: &Material,
    ) -> Result<(), MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        let c = self.coords(&nodes[..3]);
        let area2 =
            (c[1][0] - c[0][0]) * (c[2][1] - c[0][1]) - (c[2][0] - c[0][0]) * (c[1][1] - c[0][1]);
        let longest = TRI6_EDGES
            .iter()
            .map(|&(a, b)| distance(&c[a], &c[b]))
            .fold(0.0_f64, f64::max);
        if area2.abs() <= 1e-12 * longest * longest {
            return Err(MailleError::DegenerateTriangle { element });
        }
        if area2 < 0.0 {
            return Err(MailleError::ClockwiseTriangle { element });
        }
        Ok(())
    }

    /// Vérifie ν, le volume et l'orientation d'un tétraèdre ; seuls ses 4 premiers Nodes (sommets) comptent.
    fn check_tetrahedron(
        &self,
        element: usize,
        nodes: &[usize],
        material: &Material,
    ) -> Result<(), MailleError> {
        if material.poisson.is_none() {
            return Err(MailleError::MissingPoisson);
        }
        let c = self.coords_3d(&nodes[..4]);
        let e: [[f64; 3]; 3] =
            std::array::from_fn(|k| std::array::from_fn(|j| c[k + 1][j] - c[0][j]));
        // 6 V = e1 · (e2 × e3), positif pour l'orientation de `Tet4`
        let six_v = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1])
            - e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0])
            + e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
        let longest = TET10_EDGES
            .iter()
            .map(|&(a, b)| distance(&c[a], &c[b]))
            .fold(0.0_f64, f64::max);
        if six_v.abs() <= 1e-12 * longest.powi(3) {
            return Err(MailleError::DegenerateTetrahedron { element });
        }
        if six_v < 0.0 {
            return Err(MailleError::InvertedTetrahedron { element });
        }
        Ok(())
    }

    /// Positions (x, y, z) des Nodes donnés (indices valides garantis par `new`).
    pub fn coords_3d(&self, nodes: &[usize]) -> Vec<[f64; 3]> {
        nodes
            .iter()
            .map(|&n| [self.nodes[n].x, self.nodes[n].y, self.nodes[n].z])
            .collect()
    }

    /// Positions (x, y) des Nodes donnés (indices valides garantis par `new`).
    pub fn coords(&self, nodes: &[usize]) -> Vec<[f64; 2]> {
        nodes
            .iter()
            .map(|&n| [self.nodes[n].x, self.nodes[n].y])
            .collect()
    }

    /// Longueur de l'Element 1D `i` (indice valide garanti par `new`).
    pub fn length(&self, i: usize) -> f64 {
        let nodes = self.elements[i].nodes();
        let (a, b) = (nodes[0], nodes[1]);
        (self.nodes[b].x - self.nodes[a].x).abs()
    }

    /// 1 DOF par Node pour des Bar, 2 pour des Beam et des triangles, 3 pour des tétraèdres
    /// (`new` interdit le mélange).
    pub fn dofs_per_node(&self) -> usize {
        match self.elements[0].family() {
            Family::Bar => 1,
            Family::Beam | Family::Plane => 2,
            Family::Solid => 3,
        }
    }

    /// Indice global du DOF `dof` du Node `node`.
    pub fn dof_index(&self, node: usize, dof: Dof) -> Result<usize, MailleError> {
        if node >= self.nodes.len() {
            return Err(MailleError::UnknownNode { node });
        }
        match (self.elements[0].family(), dof) {
            (Family::Bar, Dof::U) => Ok(node),
            (Family::Beam, Dof::V) | (Family::Plane, Dof::Ux) => Ok(2 * node),
            (Family::Beam, Dof::Theta) | (Family::Plane, Dof::Uy) => Ok(2 * node + 1),
            (Family::Solid, Dof::Ux) => Ok(3 * node),
            (Family::Solid, Dof::Uy) => Ok(3 * node + 1),
            (Family::Solid, Dof::Uz) => Ok(3 * node + 2),
            _ => Err(MailleError::InvalidDof { node, dof }),
        }
    }
}

/// Distance entre deux points de même dimension.
fn distance(p: &[f64], q: &[f64]) -> f64 {
    p.iter()
        .zip(q)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Vrai si chaque Node milieu `first_mid + k` est au milieu de l'arête `edges[k]` (côtés droits).
fn straight_edges<const D: usize>(
    coords: &[[f64; D]],
    edges: &[(usize, usize)],
    first_mid: usize,
) -> bool {
    edges.iter().enumerate().all(|(k, &(a, b))| {
        let middle: [f64; D] = std::array::from_fn(|j| (coords[a][j] + coords[b][j]) / 2.0);
        distance(&coords[first_mid + k], &middle) <= 1e-9 * distance(&coords[a], &coords[b])
    })
}
