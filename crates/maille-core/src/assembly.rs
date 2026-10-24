use crate::element::{bar_stiffness, beam_stiffness};
use crate::mesh::{Element, Mesh};

/// Matrice dense carrée, stockage ligne par ligne.
#[derive(Debug, Clone, PartialEq)]
pub struct DenseMatrix {
    pub n: usize,
    pub data: Vec<f64>,
}

impl DenseMatrix {
    pub fn zeros(n: usize) -> Self {
        Self {
            n,
            data: vec![0.0; n * n],
        }
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] += v;
    }
}

/// Indices globaux des DOF de l'Element `e`, dans l'ordre de sa matrice élémentaire.
pub fn element_dofs(e: &Element) -> Vec<usize> {
    let nodes = e.nodes();
    let (a, b) = (nodes[0], nodes[1]);
    match e {
        Element::Bar(_) => vec![a, b],
        Element::Beam(_) => vec![2 * a, 2 * a + 1, 2 * b, 2 * b + 1],
    }
}

/// Assemble K global : DOF globaux = `node * dofs_per_node + local`.
pub fn assemble(mesh: &Mesh) -> DenseMatrix {
    let mut k = DenseMatrix::zeros(mesh.nodes.len() * mesh.dofs_per_node());
    for (i, e) in mesh.elements.iter().enumerate() {
        let l = mesh.length(i);
        let ke: Vec<Vec<f64>> = match e {
            Element::Bar(b) => bar_stiffness(b.material.young, b.section.area, l)
                .iter()
                .map(|r| r.to_vec())
                .collect(),
            Element::Beam(b) => beam_stiffness(b.material.young, b.inertia, l)
                .iter()
                .map(|r| r.to_vec())
                .collect(),
        };
        let dofs = element_dofs(e);
        for (a, &ga) in dofs.iter().enumerate() {
            for (b, &gb) in dofs.iter().enumerate() {
                k.add(ga, gb, ke[a][b]);
            }
        }
    }
    k
}
