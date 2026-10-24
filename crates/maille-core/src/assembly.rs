use crate::element::{
    bar_stiffness, beam_stiffness, conduction, heat_fluxes, plane_d, plane_gauss, plane_stiffness,
    plane_stresses, solid_d, solid_gauss, solid_stiffness, solid_stresses, GaussData,
};
use crate::material::Material;
use crate::mesh::{Element, Mesh, PlaneState, Tet10, Tet4, Tri3, Tri6};
use crate::quadrature::{GaussPoint, TetGaussPoint, TETRA_1, TETRA_4, TRIANGLE_1, TRIANGLE_3};
use crate::results::{HeatFlux, Stress};
use crate::shape::{tet10, tet4, tri3_vec, tri6_vec, ShapeFn, SolidShapeFn};

/// Indices globaux des DOF de l'Element `e`, dans l'ordre de sa matrice élémentaire.
pub fn element_dofs(e: &Element) -> Vec<usize> {
    let nodes = e.nodes();
    match e {
        Element::Bar(_) => nodes.to_vec(),
        Element::Beam(_) | Element::Tri3(_) | Element::Tri6(_) => {
            nodes.iter().flat_map(|&n| [2 * n, 2 * n + 1]).collect()
        }
        Element::Tet4(_) | Element::Tet10(_) => nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect(),
    }
}

/// E du Material, que `Model::new` garantit présent en statique.
fn young(material: &Material) -> f64 {
    material.young.unwrap_or(0.0)
}

/// Matrice de rigidité de l'Element `i`, DOF dans l'ordre de `element_dofs`.
pub fn element_stiffness(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    match &mesh.elements[i] {
        Element::Bar(b) => bar_stiffness(young(&b.material), b.section.area, mesh.length(i))
            .iter()
            .map(|r| r.to_vec())
            .collect(),
        Element::Beam(b) => beam_stiffness(young(&b.material), b.inertia, mesh.length(i))
            .iter()
            .map(|r| r.to_vec())
            .collect(),
        Element::Tri3(t) => PlaneKernel::tri3(t).stiffness(mesh),
        Element::Tri6(t) => PlaneKernel::tri6(t).stiffness(mesh),
        Element::Tet4(t) => SolidKernel::tet4(t).stiffness(mesh),
        Element::Tet10(t) => SolidKernel::tet10(t).stiffness(mesh),
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

    /// `Model::new` garantit que ν est présent.
    fn nu(&self) -> f64 {
        self.material.poisson.unwrap_or(0.0)
    }

    fn stiffness(&self, mesh: &Mesh) -> Vec<Vec<f64>> {
        let d = plane_d(young(self.material), self.nu(), self.state);
        plane_stiffness(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            &d,
            self.thickness,
        )
    }

    fn gauss(&self, mesh: &Mesh) -> Vec<GaussData<2>> {
        plane_gauss(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            self.thickness,
        )
    }

    fn stresses(&self, mesh: &Mesh, u: &[f64]) -> Vec<Stress> {
        let nu = self.nu();
        let d = plane_d(young(self.material), nu, self.state);
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

/// Ce qui distingue un tétraèdre d'un autre : règle de Gauss et fonctions de forme.
struct SolidKernel<'a> {
    nodes: &'a [usize],
    material: &'a Material,
    rule: &'static [TetGaussPoint],
    shape: SolidShapeFn,
}

impl<'a> SolidKernel<'a> {
    fn tet4(t: &'a Tet4) -> Self {
        Self {
            nodes: &t.nodes,
            material: &t.material,
            rule: &TETRA_1,
            shape: tet4,
        }
    }

    /// 4 points de Gauss : exact pour la rigidité d'un Tet10 à arêtes droites.
    fn tet10(t: &'a Tet10) -> Self {
        Self {
            nodes: &t.nodes,
            material: &t.material,
            rule: &TETRA_4,
            shape: tet10,
        }
    }

    /// `Model::new` garantit que ν est présent.
    fn d(&self) -> [[f64; 6]; 6] {
        solid_d(young(self.material), self.material.poisson.unwrap_or(0.0))
    }

    fn stiffness(&self, mesh: &Mesh) -> Vec<Vec<f64>> {
        solid_stiffness(
            &mesh.coords_3d(self.nodes),
            self.rule,
            self.shape,
            &self.d(),
        )
    }

    fn gauss(&self, mesh: &Mesh) -> Vec<GaussData<3>> {
        solid_gauss(&mesh.coords_3d(self.nodes), self.rule, self.shape)
    }

    fn stresses(&self, mesh: &Mesh, u: &[f64]) -> Vec<Stress> {
        let u_e: Vec<f64> = self
            .nodes
            .iter()
            .flat_map(|&n| [u[3 * n], u[3 * n + 1], u[3 * n + 2]])
            .collect();
        solid_stresses(
            &mesh.coords_3d(self.nodes),
            self.rule,
            self.shape,
            &self.d(),
            &u_e,
        )
    }
}

/// Contraintes aux points de Gauss de l'Element `i`, à partir des déplacements globaux `u`.
/// Vide pour les Bar et Beam.
pub fn element_stresses(mesh: &Mesh, i: usize, u: &[f64]) -> Vec<Stress> {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).stresses(mesh, u),
        Element::Tri6(t) => PlaneKernel::tri6(t).stresses(mesh, u),
        Element::Tet4(t) => SolidKernel::tet4(t).stresses(mesh, u),
        Element::Tet10(t) => SolidKernel::tet10(t).stresses(mesh, u),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// k du Material, que `Model` garantit présent en thermique.
fn conductivity(material: &Material) -> f64 {
    material.conductivity.unwrap_or(0.0)
}

/// Matrice de conduction de l'Element `i`, un DOF (T) par Node dans l'ordre de ses Nodes.
/// Vide pour les Bar et Beam (refusés en thermique par `Model`).
pub fn element_conduction(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    match &mesh.elements[i] {
        Element::Tri3(t) => {
            conduction(&PlaneKernel::tri3(t).gauss(mesh), conductivity(&t.material))
        }
        Element::Tri6(t) => {
            conduction(&PlaneKernel::tri6(t).gauss(mesh), conductivity(&t.material))
        }
        Element::Tet4(t) => {
            conduction(&SolidKernel::tet4(t).gauss(mesh), conductivity(&t.material))
        }
        Element::Tet10(t) => conduction(
            &SolidKernel::tet10(t).gauss(mesh),
            conductivity(&t.material),
        ),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// Flux de chaleur aux points de Gauss de l'Element `i`, à partir des températures `t` des Nodes.
/// Vide pour les Bar et Beam.
pub fn element_heat_fluxes(mesh: &Mesh, i: usize, t: &[f64]) -> Vec<HeatFlux> {
    let e = &mesh.elements[i];
    let t_e: Vec<f64> = e.nodes().iter().map(|&n| t[n]).collect();
    let k = conductivity(e.material());
    match e {
        Element::Tri3(tri) => heat_fluxes(&PlaneKernel::tri3(tri).gauss(mesh), k, &t_e),
        Element::Tri6(tri) => heat_fluxes(&PlaneKernel::tri6(tri).gauss(mesh), k, &t_e),
        Element::Tet4(tet) => heat_fluxes(&SolidKernel::tet4(tet).gauss(mesh), k, &t_e),
        Element::Tet10(tet) => heat_fluxes(&SolidKernel::tet10(tet).gauss(mesh), k, &t_e),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}
