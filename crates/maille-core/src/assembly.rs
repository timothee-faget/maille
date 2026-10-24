use crate::element::{bar_stiffness, beam_stiffness, plane_d, plane_stiffness};
use crate::mesh::{Element, Mesh};
use crate::quadrature::TRIANGLE_1;
use crate::shape::tri3;

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
    match e {
        Element::Bar(_) => nodes.to_vec(),
        Element::Beam(_) | Element::Tri3(_) => {
            nodes.iter().flat_map(|&n| [2 * n, 2 * n + 1]).collect()
        }
    }
}

/// Assemble K global : DOF globaux = `node * dofs_per_node + local`.
pub fn assemble(mesh: &Mesh) -> DenseMatrix {
    let mut k = DenseMatrix::zeros(mesh.nodes.len() * mesh.dofs_per_node());
    for (i, e) in mesh.elements.iter().enumerate() {
        let ke: Vec<Vec<f64>> = match e {
            Element::Bar(b) => bar_stiffness(b.material.young, b.section.area, mesh.length(i))
                .iter()
                .map(|r| r.to_vec())
                .collect(),
            Element::Beam(b) => beam_stiffness(b.material.young, b.inertia, mesh.length(i))
                .iter()
                .map(|r| r.to_vec())
                .collect(),
            Element::Tri3(t) => {
                // `Mesh::new` garantit que ν est présent
                let nu = t.material.poisson.unwrap_or(0.0);
                let d = plane_d(t.material.young, nu, t.state);
                plane_stiffness(
                    &mesh.coords(&t.nodes),
                    &TRIANGLE_1,
                    |xi, eta| tri3(xi, eta).1.to_vec(),
                    &d,
                    t.thickness,
                )
            }
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
