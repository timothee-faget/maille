use crate::error::MailleError;
use crate::material::{Material, Section};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub x: f64,
}

/// Degré de liberté d'un Node : U déplacement axial (Bar), V flèche et Theta rotation (Beam).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dof {
    U,
    V,
    Theta,
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Element {
    Bar(Bar),
    Beam(Beam),
}

impl Element {
    pub fn nodes(&self) -> [usize; 2] {
        match self {
            Self::Bar(e) => e.nodes,
            Self::Beam(e) => e.nodes,
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
        let first_is_bar = matches!(mesh.elements[0], Element::Bar(_));
        for (i, e) in mesh.elements.iter().enumerate() {
            if matches!(e, Element::Bar(_)) != first_is_bar {
                return Err(MailleError::MixedElementTypes);
            }
            for n in e.nodes() {
                if n >= mesh.nodes.len() {
                    return Err(MailleError::UnknownNode { node: n });
                }
            }
            if mesh.length(i) == 0.0 {
                return Err(MailleError::ZeroLengthElement { element: i });
            }
            if let Element::Beam(b) = e {
                if mesh.nodes[b.nodes[1]].x < mesh.nodes[b.nodes[0]].x {
                    return Err(MailleError::ReversedBeam { element: i });
                }
            }
        }
        Ok(mesh)
    }

    /// Longueur de l'Element `i` (indice valide garanti par `new`).
    pub fn length(&self, i: usize) -> f64 {
        let [a, b] = self.elements[i].nodes();
        (self.nodes[b].x - self.nodes[a].x).abs()
    }

    /// 1 DOF par Node pour des Bar, 2 pour des Beam (`new` interdit le mélange).
    pub fn dofs_per_node(&self) -> usize {
        match self.elements[0] {
            Element::Bar(_) => 1,
            Element::Beam(_) => 2,
        }
    }

    /// Indice global du DOF `dof` du Node `node`.
    pub fn dof_index(&self, node: usize, dof: Dof) -> Result<usize, MailleError> {
        if node >= self.nodes.len() {
            return Err(MailleError::UnknownNode { node });
        }
        match (self.dofs_per_node(), dof) {
            (1, Dof::U) => Ok(node),
            (2, Dof::V) => Ok(2 * node),
            (2, Dof::Theta) => Ok(2 * node + 1),
            _ => Err(MailleError::InvalidDof { node, dof }),
        }
    }
}
