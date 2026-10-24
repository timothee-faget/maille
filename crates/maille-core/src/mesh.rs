use crate::error::MailleError;
use crate::material::{Material, Section};

/// Position d'un Node. Les Bar et Beam (1D) n'utilisent que `x` et ignorent `y`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub x: f64,
    pub y: f64,
}

/// Degré de liberté d'un Node : U déplacement axial (Bar), V flèche et Theta rotation (Beam),
/// Ux et Uy déplacements dans le plan (triangles).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dof {
    U,
    V,
    Theta,
    Ux,
    Uy,
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Element {
    Bar(Bar),
    Beam(Beam),
    Tri3(Tri3),
    Tri6(Tri6),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Family {
    Bar,
    Beam,
    Plane,
}

impl Element {
    pub fn nodes(&self) -> &[usize] {
        match self {
            Self::Bar(e) => &e.nodes,
            Self::Beam(e) => &e.nodes,
            Self::Tri3(e) => &e.nodes,
            Self::Tri6(e) => &e.nodes,
        }
    }

    fn family(&self) -> Family {
        match self {
            Self::Bar(_) => Family::Bar,
            Self::Beam(_) => Family::Beam,
            Self::Tri3(_) | Self::Tri6(_) => Family::Plane,
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
                    let c = mesh.coords(&t.nodes);
                    for (mid, (a, b)) in [(3, (0, 1)), (4, (1, 2)), (5, (2, 0))] {
                        let middle = [(c[a][0] + c[b][0]) / 2.0, (c[a][1] + c[b][1]) / 2.0];
                        let edge = (c[b][0] - c[a][0]).hypot(c[b][1] - c[a][1]);
                        let off = (c[mid][0] - middle[0]).hypot(c[mid][1] - middle[1]);
                        if off > 1e-9 * edge {
                            return Err(MailleError::CurvedTriangle { element: i });
                        }
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
        let longest = [(0, 1), (1, 2), (2, 0)]
            .iter()
            .map(|&(a, b)| (c[b][0] - c[a][0]).hypot(c[b][1] - c[a][1]))
            .fold(0.0_f64, f64::max);
        if area2.abs() <= 1e-12 * longest * longest {
            return Err(MailleError::DegenerateTriangle { element });
        }
        if area2 < 0.0 {
            return Err(MailleError::ClockwiseTriangle { element });
        }
        Ok(())
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

    /// 1 DOF par Node pour des Bar, 2 pour des Beam et des triangles (`new` interdit le mélange).
    pub fn dofs_per_node(&self) -> usize {
        match self.elements[0].family() {
            Family::Bar => 1,
            Family::Beam | Family::Plane => 2,
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
            _ => Err(MailleError::InvalidDof { node, dof }),
        }
    }
}
