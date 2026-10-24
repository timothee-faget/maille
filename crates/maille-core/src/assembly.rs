use crate::element::{bar_stiffness, beam_stiffness, plane_d, plane_stiffness, plane_stresses};
use crate::material::Material;
use crate::mesh::{Element, Mesh, PlaneState, Tri3, Tri6};
use crate::quadrature::{GaussPoint, TRIANGLE_1, TRIANGLE_3};
use crate::results::PlaneStress;
use crate::shape::{tri3_vec, tri6_vec, ShapeFn};

/// Indices globaux des DOF de l'Element `e`, dans l'ordre de sa matrice élémentaire.
pub fn element_dofs(e: &Element) -> Vec<usize> {
    let nodes = e.nodes();
    match e {
        Element::Bar(_) => nodes.to_vec(),
        Element::Beam(_) | Element::Tri3(_) | Element::Tri6(_) => {
            nodes.iter().flat_map(|&n| [2 * n, 2 * n + 1]).collect()
        }
    }
}

/// Matrice de rigidité de l'Element `i`, DOF dans l'ordre de `element_dofs`.
pub fn element_stiffness(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    match &mesh.elements[i] {
        Element::Bar(b) => bar_stiffness(b.material.young, b.section.area, mesh.length(i))
            .iter()
            .map(|r| r.to_vec())
            .collect(),
        Element::Beam(b) => beam_stiffness(b.material.young, b.inertia, mesh.length(i))
            .iter()
            .map(|r| r.to_vec())
            .collect(),
        Element::Tri3(t) => PlaneKernel::tri3(t).stiffness(mesh),
        Element::Tri6(t) => PlaneKernel::tri6(t).stiffness(mesh),
    }
}

/// Ce qui distingue un triangle d'un autre pour le calcul : règle de Gauss et fonctions de forme.
struct PlaneKernel<'a> {
    nodes: &'a [usize],
    material: &'a Material,
    thickness: f64,
    state: PlaneState,
    rule: &'static [GaussPoint],
    shape: ShapeFn,
}

impl<'a> PlaneKernel<'a> {
    fn tri3(t: &'a Tri3) -> Self {
        Self {
            nodes: &t.nodes,
            material: &t.material,
            thickness: t.thickness,
            state: t.state,
            rule: &TRIANGLE_1,
            shape: tri3_vec,
        }
    }

    /// 3 points de Gauss : exact pour la rigidité d'un Tri6 à côtés droits.
    fn tri6(t: &'a Tri6) -> Self {
        Self {
            nodes: &t.nodes,
            material: &t.material,
            thickness: t.thickness,
            state: t.state,
            rule: &TRIANGLE_3,
            shape: tri6_vec,
        }
    }

    /// `Mesh::new` garantit que ν est présent.
    fn nu(&self) -> f64 {
        self.material.poisson.unwrap_or(0.0)
    }

    fn stiffness(&self, mesh: &Mesh) -> Vec<Vec<f64>> {
        let d = plane_d(self.material.young, self.nu(), self.state);
        plane_stiffness(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            &d,
            self.thickness,
        )
    }

    fn stresses(&self, mesh: &Mesh, u: &[f64]) -> Vec<PlaneStress> {
        let nu = self.nu();
        let d = plane_d(self.material.young, nu, self.state);
        let u_e: Vec<f64> = self
            .nodes
            .iter()
            .flat_map(|&n| [u[2 * n], u[2 * n + 1]])
            .collect();
        plane_stresses(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            &d,
            &u_e,
            nu,
            self.state,
        )
    }
}

/// Contraintes aux points de Gauss de l'Element `i`, à partir des déplacements globaux `u`.
/// Vide pour les Bar et Beam.
pub fn element_stresses(mesh: &Mesh, i: usize, u: &[f64]) -> Vec<PlaneStress> {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).stresses(mesh, u),
        Element::Tri6(t) => PlaneKernel::tri6(t).stresses(mesh, u),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}
