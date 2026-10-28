use crate::element::{
    bar_stiffness, bar_thermal_load, beam_stiffness, capacity, conduction, heat_fluxes,
    in_plane_expansion, plane_d, plane_gauss, plane_stiffness, plane_strain_response,
    plane_stresses, plane_thermal_load, solid_d, solid_gauss, solid_response, solid_stiffness,
    solid_stresses, solid_thermal_load, ElementResponse, GaussData,
};
use crate::finite_strain::{large_strain_response, SaintVenantKirchhoff};
use crate::material::Material;
use crate::mesh::{Element, Mesh, PlaneState, Tet10, Tet4, Tri3, Tri6};
use crate::plasticity::{PlasticState, VonMises};
use crate::quadrature::{
    GaussPoint, LinePoint, TetGaussPoint, LINE_2, LINE_3, TETRA_1, TETRA_14, TETRA_4, TRIANGLE_1,
    TRIANGLE_3, TRIANGLE_6,
};
use crate::results::{HeatFlux, Stress};
use crate::shape::{line2, line3, tet10, tet4, tri3_vec, tri6_vec, ShapeFn, SolidShapeFn};

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

/// Loi de von Mises du Material (E et ν garantis par `Model::new`). Sans `Plasticity`, σy est
/// infinie : la loi reste élastique.
fn von_mises(material: &Material) -> VonMises {
    let (yield_stress, hardening) = material
        .plasticity
        .map_or((f64::INFINITY, 0.0), |p| (p.yield_stress, p.hardening));
    VonMises {
        young: young(material),
        poisson: material.poisson.unwrap_or(0.0),
        yield_stress,
        hardening,
    }
}

/// Loi de Saint-Venant–Kirchhoff du Material (E et ν garantis par `Model::new`).
fn saint_venant_kirchhoff(material: &Material) -> SaintVenantKirchhoff {
    SaintVenantKirchhoff::new(young(material), material.poisson.unwrap_or(0.0))
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

    /// Le même triangle intégré avec une autre règle de Gauss.
    fn with_rule(self, rule: &'static [GaussPoint]) -> Self {
        Self { rule, ..self }
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

    /// Déplacements locaux (ux, uy) par Node, tirés des déplacements globaux `u`.
    fn displacements(&self, u: &[f64]) -> Vec<f64> {
        self.nodes
            .iter()
            .flat_map(|&n| [u[2 * n], u[2 * n + 1]])
            .collect()
    }

    /// Réponse élasto-plastique, en déformations planes (`Model` refuse les contraintes planes).
    fn response(&self, mesh: &Mesh, u: &[f64], previous: &[PlasticState]) -> ElementResponse {
        plane_strain_response(
            &self.gauss(mesh),
            &von_mises(self.material),
            &self.displacements(u),
            previous,
        )
    }

    /// Réponse en grandes déformations, en déformations planes (`Model` refuse les contraintes planes).
    fn large_response(&self, mesh: &Mesh, u: &[f64]) -> ElementResponse {
        large_strain_response(
            &self.gauss(mesh),
            &saint_venant_kirchhoff(self.material),
            &self.displacements(u),
        )
    }

    fn stresses(&self, mesh: &Mesh, u: &[f64], strain_e: &[f64]) -> Vec<Stress> {
        let nu = self.nu();
        let d = plane_d(young(self.material), nu, self.state);
        let u_e = self.displacements(u);
        plane_stresses(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            &d,
            &u_e,
            young(self.material),
            nu,
            self.state,
            strain_e,
        )
    }

    fn thermal_load(&self, mesh: &Mesh, strain_e: &[f64]) -> Vec<f64> {
        let nu = self.nu();
        plane_thermal_load(
            &mesh.coords(self.nodes),
            self.rule,
            self.shape,
            &plane_d(young(self.material), nu, self.state),
            self.thickness,
            in_plane_expansion(nu, self.state),
            strain_e,
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

    /// Le même tétraèdre intégré avec une autre règle de Gauss.
    fn with_rule(self, rule: &'static [TetGaussPoint]) -> Self {
        Self { rule, ..self }
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

    /// Déplacements locaux (ux, uy, uz) par Node, tirés des déplacements globaux `u`.
    fn displacements(&self, u: &[f64]) -> Vec<f64> {
        self.nodes
            .iter()
            .flat_map(|&n| [u[3 * n], u[3 * n + 1], u[3 * n + 2]])
            .collect()
    }

    fn response(&self, mesh: &Mesh, u: &[f64], previous: &[PlasticState]) -> ElementResponse {
        solid_response(
            &self.gauss(mesh),
            &von_mises(self.material),
            &self.displacements(u),
            previous,
        )
    }

    fn large_response(&self, mesh: &Mesh, u: &[f64]) -> ElementResponse {
        large_strain_response(
            &self.gauss(mesh),
            &saint_venant_kirchhoff(self.material),
            &self.displacements(u),
        )
    }

    fn stresses(&self, mesh: &Mesh, u: &[f64], strain_e: &[f64]) -> Vec<Stress> {
        solid_stresses(
            &mesh.coords_3d(self.nodes),
            self.rule,
            self.shape,
            &self.d(),
            &self.displacements(u),
            strain_e,
        )
    }

    fn thermal_load(&self, mesh: &Mesh, strain_e: &[f64]) -> Vec<f64> {
        solid_thermal_load(
            &mesh.coords_3d(self.nodes),
            self.rule,
            self.shape,
            &self.d(),
            strain_e,
        )
    }
}

/// Contraintes aux points de Gauss de l'Element `i`, à partir des déplacements globaux `u` ;
/// `strain_e` : α (T − T_ref) aux Nodes de l'Element (vide sans dilatation, voir
/// `element_thermal_strains`). Vide pour les Bar et Beam.
pub fn element_stresses(mesh: &Mesh, i: usize, u: &[f64], strain_e: &[f64]) -> Vec<Stress> {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).stresses(mesh, u, strain_e),
        Element::Tri6(t) => PlaneKernel::tri6(t).stresses(mesh, u, strain_e),
        Element::Tet4(t) => SolidKernel::tet4(t).stresses(mesh, u, strain_e),
        Element::Tet10(t) => SolidKernel::tet10(t).stresses(mesh, u, strain_e),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// Nombre de points de Gauss de l'Element `i` pour la rigidité et les contraintes (0 pour les Bar et
/// Beam).
pub fn element_gauss_count(mesh: &Mesh, i: usize) -> usize {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).rule.len(),
        Element::Tri6(t) => PlaneKernel::tri6(t).rule.len(),
        Element::Tet4(t) => SolidKernel::tet4(t).rule.len(),
        Element::Tet10(t) => SolidKernel::tet10(t).rule.len(),
        Element::Bar(_) | Element::Beam(_) => 0,
    }
}

/// Réponse élasto-plastique de l'Element `i` aux déplacements globaux `u`, depuis les états
/// `previous` de ses points de Gauss (DOF dans l'ordre de `element_dofs`). Triangles en
/// déformations planes et tétraèdres ; vide pour les Bar et Beam (refusées par `Model`).
pub fn element_response(
    mesh: &Mesh,
    i: usize,
    u: &[f64],
    previous: &[PlasticState],
) -> ElementResponse {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).response(mesh, u, previous),
        Element::Tri6(t) => PlaneKernel::tri6(t).response(mesh, u, previous),
        Element::Tet4(t) => SolidKernel::tet4(t).response(mesh, u, previous),
        Element::Tet10(t) => SolidKernel::tet10(t).response(mesh, u, previous),
        Element::Bar(_) | Element::Beam(_) => ElementResponse::default(),
    }
}

/// Réponse en grandes déformations (Lagrangien total, Saint-Venant–Kirchhoff) de l'Element `i` aux
/// déplacements globaux `u`, DOF dans l'ordre de `element_dofs`. Contraintes de Cauchy. Triangles en
/// déformations planes et tétraèdres ; vide pour les Bar et Beam (refusées par `Model`).
pub fn element_large_response(mesh: &Mesh, i: usize, u: &[f64]) -> ElementResponse {
    match &mesh.elements[i] {
        Element::Tri3(t) => PlaneKernel::tri3(t).large_response(mesh, u),
        Element::Tri6(t) => PlaneKernel::tri6(t).large_response(mesh, u),
        Element::Tet4(t) => SolidKernel::tet4(t).large_response(mesh, u),
        Element::Tet10(t) => SolidKernel::tet10(t).large_response(mesh, u),
        Element::Bar(_) | Element::Beam(_) => ElementResponse::default(),
    }
}

/// α (T − T_ref) aux Nodes de l'Element `i`, `t` : températures de tous les Nodes. α du Material,
/// que `Model::with_temperatures` garantit présent.
pub fn element_thermal_strains(mesh: &Mesh, i: usize, t: &[f64], reference: f64) -> Vec<f64> {
    let e = &mesh.elements[i];
    let alpha = e.material().expansion.unwrap_or(0.0);
    e.nodes()
        .iter()
        .map(|&n| alpha * (t[n] - reference))
        .collect()
}

/// Forces nodales de la dilatation de l'Element `i`, DOF dans l'ordre de `element_dofs` ;
/// `strain_e` : voir `element_thermal_strains`. Vide pour les Beam (refusées par `Model`).
pub fn element_thermal_load(mesh: &Mesh, i: usize, strain_e: &[f64]) -> Vec<f64> {
    match &mesh.elements[i] {
        Element::Bar(b) => {
            let [a, c] = b.nodes;
            let direction = (mesh.nodes[c].x - mesh.nodes[a].x).signum();
            let mean = strain_e.iter().sum::<f64>() / 2.0;
            bar_thermal_load(young(&b.material), b.section.area, mean, direction).to_vec()
        }
        Element::Tri3(t) => PlaneKernel::tri3(t).thermal_load(mesh, strain_e),
        Element::Tri6(t) => PlaneKernel::tri6(t).thermal_load(mesh, strain_e),
        Element::Tet4(t) => SolidKernel::tet4(t).thermal_load(mesh, strain_e),
        Element::Tet10(t) => SolidKernel::tet10(t).thermal_load(mesh, strain_e),
        Element::Beam(_) => Vec::new(),
    }
}

/// k du Material, que `Model` garantit présent en thermique.
fn conductivity(material: &Material) -> f64 {
    material.conductivity.unwrap_or(0.0)
}

/// Matrice de conduction de l'Element `i`, un DOF (T) par Node dans l'ordre de ses Nodes.
/// Vide pour les Bar et Beam (refusés en thermique par `Model`).
pub fn element_conduction(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    let e = &mesh.elements[i];
    let k = conductivity(e.material());
    match e {
        Element::Tri3(t) => conduction(&PlaneKernel::tri3(t).gauss(mesh), k),
        Element::Tri6(t) => conduction(&PlaneKernel::tri6(t).gauss(mesh), k),
        Element::Tet4(t) => conduction(&SolidKernel::tet4(t).gauss(mesh), k),
        Element::Tet10(t) => conduction(&SolidKernel::tet10(t).gauss(mesh), k),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// ρc du Material, que `Model` garantit présent en thermique transitoire.
fn heat_capacity(material: &Material) -> f64 {
    material.density.unwrap_or(0.0) * material.specific_heat.unwrap_or(0.0)
}

/// Matrice de capacité cohérente de l'Element `i`, un DOF (T) par Node dans l'ordre de ses Nodes.
/// Règles de Gauss qui intègrent N Nᵀ exactement : degré 2 pour Tri3 et Tet4, 4 pour Tri6 et Tet10.
/// Vide pour les Bar et Beam.
pub fn element_capacity(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    let e = &mesh.elements[i];
    let rho_c = heat_capacity(e.material());
    match e {
        Element::Tri3(t) => capacity(
            &PlaneKernel::tri3(t).with_rule(&TRIANGLE_3).gauss(mesh),
            rho_c,
        ),
        Element::Tri6(t) => capacity(
            &PlaneKernel::tri6(t).with_rule(&TRIANGLE_6).gauss(mesh),
            rho_c,
        ),
        Element::Tet4(t) => capacity(&SolidKernel::tet4(t).with_rule(&TETRA_4).gauss(mesh), rho_c),
        Element::Tet10(t) => capacity(
            &SolidKernel::tet10(t).with_rule(&TETRA_14).gauss(mesh),
            rho_c,
        ),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// Flux de chaleur aux points de Gauss d'un Element en fonction des températures globales des Nodes.
pub type FluxOperator = Box<dyn Fn(&[f64]) -> Vec<HeatFlux>>;

/// Flux de chaleur aux points de Gauss de l'Element `i` (vide pour les Bar et Beam). Les points de
/// Gauss sont calculés une seule fois : l'opérateur renvoyé s'applique à autant de champs de
/// températures qu'on veut (un par incrément en transitoire).
pub fn element_heat_flux_operator(mesh: &Mesh, i: usize) -> FluxOperator {
    let e = &mesh.elements[i];
    let (nodes, k) = (e.nodes().to_vec(), conductivity(e.material()));
    match e {
        Element::Tri3(tri) => flux_operator(PlaneKernel::tri3(tri).gauss(mesh), nodes, k),
        Element::Tri6(tri) => flux_operator(PlaneKernel::tri6(tri).gauss(mesh), nodes, k),
        Element::Tet4(tet) => flux_operator(SolidKernel::tet4(tet).gauss(mesh), nodes, k),
        Element::Tet10(tet) => flux_operator(SolidKernel::tet10(tet).gauss(mesh), nodes, k),
        Element::Bar(_) | Element::Beam(_) => Box::new(|_| Vec::new()),
    }
}

fn flux_operator<const D: usize>(
    points: Vec<GaussData<D>>,
    nodes: Vec<usize>,
    k: f64,
) -> FluxOperator {
    Box::new(move |t| {
        let t_e: Vec<f64> = nodes.iter().map(|&n| t[n]).collect();
        heat_fluxes(&points, k, &t_e)
    })
}

/// Puissances nodales équivalentes ∫ N Q dV d'une source uniforme `q` dans l'Element `i`.
/// Vide pour les Bar et Beam.
pub fn element_heat_source(mesh: &Mesh, i: usize, q: f64) -> Vec<f64> {
    match &mesh.elements[i] {
        Element::Tri3(t) => nodal_source(&PlaneKernel::tri3(t).gauss(mesh), q),
        Element::Tri6(t) => nodal_source(&PlaneKernel::tri6(t).gauss(mesh), q),
        Element::Tet4(t) => nodal_source(&SolidKernel::tet4(t).gauss(mesh), q),
        Element::Tet10(t) => nodal_source(&SolidKernel::tet10(t).gauss(mesh), q),
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// Σ_g w_g N_i q.
fn nodal_source<const D: usize>(points: &[GaussData<D>], q: f64) -> Vec<f64> {
    let mut f = vec![0.0; points.first().map_or(0, |p| p.n.len())];
    for p in points {
        for (fi, ni) in f.iter_mut().zip(&p.n) {
            *fi += p.weight * ni * q;
        }
    }
    f
}

/// Termes de convection h (T − T∞) sur une face : `matrix` = h ∫ N Nᵀ dS et `load` = h T∞ ∫ N dS,
/// indexés comme `positions` (positions des Nodes de la face dans ceux de l'Element).
#[derive(Debug, Clone, PartialEq)]
pub struct FaceConvection {
    pub positions: Vec<usize>,
    pub matrix: Vec<Vec<f64>>,
    pub load: Vec<f64>,
}

impl FaceConvection {
    /// Ajoute h ∫ N Nᵀ à la matrice `ke` de l'Element.
    pub fn add_matrix(&self, ke: &mut [Vec<f64>]) {
        for (row, &a) in self.matrix.iter().zip(&self.positions) {
            for (v, &b) in row.iter().zip(&self.positions) {
                ke[a][b] += v;
            }
        }
    }

    /// Ajoute h T∞ ∫ N au second membre global `f` ; `dofs` : DOF globaux des Nodes de l'Element.
    pub fn add_load(&self, dofs: &[usize], f: &mut [f64]) {
        for (&p, v) in self.positions.iter().zip(&self.load) {
            f[dofs[p]] += v;
        }
    }
}

/// Intégration sur une face d'Element : Nodes de la face, points de Gauss et normale sortante.
struct FaceRule {
    /// Positions des Nodes de la face dans ceux de l'Element (`Element::face_positions`).
    positions: Vec<usize>,
    /// (poids w |dS|, valeurs N) aux points de Gauss de la face. En 2D, la face est une arête et
    /// l'épaisseur t multiplie le poids.
    points: Vec<(f64, Vec<f64>)>,
    /// Normale unitaire sortante (z = 0 en 2D).
    normal: [f64; 3],
}

/// Règle d'intégration de la face `face` de l'Element `i`, ou `None` si la face n'existe pas ou si
/// l'Element est une Bar ou une Beam. La normale est orientée vers l'extérieur en la comparant à la
/// direction du centre de l'Element vers celui de la face, sans dépendre du sens des faces.
fn face_rule(mesh: &Mesh, i: usize, face: usize) -> Option<FaceRule> {
    let e = &mesh.elements[i];
    let positions = e.face_positions(face)?;
    let global: Vec<usize> = positions.iter().map(|&k| e.nodes()[k]).collect();
    let (points, normal) = match e {
        Element::Tri3(Tri3 { thickness, .. }) | Element::Tri6(Tri6 { thickness, .. }) => {
            let c = mesh.coords(&global);
            let (dx, dy) = (c[1][0] - c[0][0], c[1][1] - c[0][1]);
            let length = dx.hypot(dy);
            let linear = global.len() == 2;
            let rule: &[LinePoint] = if linear { &LINE_2 } else { &LINE_3 };
            let shape = if linear { line2 } else { line3 };
            let points = rule
                .iter()
                .map(|p| (p.weight * length * thickness, shape(p.s)))
                .collect();
            let corners = mesh.coords(&e.nodes()[..3]);
            let outward = [
                (c[0][0] + c[1][0]) / 2.0 - (corners[0][0] + corners[1][0] + corners[2][0]) / 3.0,
                (c[0][1] + c[1][1]) / 2.0 - (corners[0][1] + corners[1][1] + corners[2][1]) / 3.0,
            ];
            let n = [dy / length, -dx / length];
            let sign = if n[0] * outward[0] + n[1] * outward[1] < 0.0 {
                -1.0
            } else {
                1.0
            };
            (points, [sign * n[0], sign * n[1], 0.0])
        }
        Element::Tet4(_) | Element::Tet10(_) => {
            let c = mesh.coords_3d(&global);
            let (u, v): ([f64; 3], [f64; 3]) = (
                std::array::from_fn(|k| c[1][k] - c[0][k]),
                std::array::from_fn(|k| c[2][k] - c[0][k]),
            );
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let twice_area = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
            let (rule, shape): (&[GaussPoint], ShapeFn) = if global.len() == 3 {
                (&TRIANGLE_3, tri3_vec)
            } else {
                (&TRIANGLE_6, tri6_vec)
            };
            let points = rule
                .iter()
                .map(|g| (g.weight * twice_area, shape(g.xi, g.eta).0))
                .collect();
            let corners = mesh.coords_3d(&e.nodes()[..4]);
            let outward: [f64; 3] = std::array::from_fn(|k| {
                (c[0][k] + c[1][k] + c[2][k]) / 3.0
                    - (corners[0][k] + corners[1][k] + corners[2][k] + corners[3][k]) / 4.0
            });
            let dot: f64 = cross.iter().zip(&outward).map(|(a, b)| a * b).sum();
            let sign = if dot < 0.0 { -1.0 } else { 1.0 };
            (points, cross.map(|x| sign * x / twice_area))
        }
        Element::Bar(_) | Element::Beam(_) => return None,
    };
    Some(FaceRule {
        positions,
        points,
        normal,
    })
}

/// Convection de coefficient `film` vers la température `sink` sur la face `face` de l'Element `i`
/// (numérotation de `Element::face_nodes`). `None` si la face n'existe pas. En 2D, la face est une
/// arête et l'épaisseur t multiplie l'intégrale.
pub fn face_convection(
    mesh: &Mesh,
    i: usize,
    face: usize,
    film: f64,
    sink: f64,
) -> Option<FaceConvection> {
    let FaceRule {
        positions, points, ..
    } = face_rule(mesh, i, face)?;
    let m = positions.len();
    let mut matrix = vec![vec![0.0; m]; m];
    let mut load = vec![0.0; m];
    for (w, n) in &points {
        for a in 0..m {
            load[a] += film * sink * w * n[a];
            for b in 0..m {
                matrix[a][b] += film * w * n[a] * n[b];
            }
        }
    }
    Some(FaceConvection {
        positions,
        matrix,
        load,
    })
}

/// Forces nodales équivalentes à une pression sur une face : `forces` donne, pour chaque Node de
/// la face (dans l'ordre de `positions`), ses `dim` composantes (x, y en 2D ; x, y, z en 3D).
#[derive(Debug, Clone, PartialEq)]
pub struct FacePressure {
    pub positions: Vec<usize>,
    pub dim: usize,
    pub forces: Vec<f64>,
}

impl FacePressure {
    /// Ajoute les forces au second membre global `f` ; `dofs` : DOF globaux des Nodes de l'Element
    /// (`dim` par Node, dans l'ordre de `element_dofs`).
    pub fn add_load(&self, dofs: &[usize], f: &mut [f64]) {
        for (k, &p) in self.positions.iter().enumerate() {
            for c in 0..self.dim {
                f[dofs[self.dim * p + c]] += self.forces[self.dim * k + c];
            }
        }
    }
}

/// Pression `value` sur la face `face` de l'Element `i` : f_a = −p ∫ N_a n dS, n étant la normale
/// sortante. Une pression positive comprime ; négative, elle tire. `None` si la face n'existe pas
/// ou si l'Element est une Bar ou une Beam. En 2D, l'épaisseur t multiplie l'intégrale.
pub fn face_pressure(mesh: &Mesh, i: usize, face: usize, value: f64) -> Option<FacePressure> {
    let FaceRule {
        positions,
        points,
        normal,
    } = face_rule(mesh, i, face)?;
    let dim = if matches!(mesh.elements[i], Element::Tet4(_) | Element::Tet10(_)) {
        3
    } else {
        2
    };
    let mut forces = vec![0.0; dim * positions.len()];
    for (w, n) in &points {
        for (a, na) in n.iter().enumerate() {
            for c in 0..dim {
                forces[dim * a + c] -= value * w * na * normal[c];
            }
        }
    }
    Some(FacePressure {
        positions,
        dim,
        forces,
    })
}

#[cfg(test)]
mod capacity_tests {
    use super::*;
    use crate::mesh::Node;
    use crate::shape::{TET10_EDGES, TRI6_EDGES};

    const RHO_C: f64 = 3.0 * 0.5; // ρ = 3, c = 0,5

    fn material() -> Material {
        Material::conductive(1.0)
            .and_then(|m| m.with_density(3.0))
            .and_then(|m| m.with_specific_heat(0.5))
            .unwrap()
    }

    /// Somme de tous les monômes de degré `p` en `s` (polynôme symétrique complet h_p).
    fn complete(s: &[f64], p: u32) -> f64 {
        match s {
            [] => f64::from(u8::from(p == 0)),
            [first, rest @ ..] => (0..=p)
                .map(|k| first.powi(k as i32) * complete(rest, p - k))
                .sum(),
        }
    }

    /// ∫ (g·x)^p sur un simplexe de dimension d, de mesure `measure`, dont les sommets valent sᵢ = g·xᵢ :
    /// en coordonnées barycentriques, ∫ λ^α = measure d! α! / (|α| + d)!, d'où measure p! d! h_p(s) / (p + d)!.
    fn simplex_power(s: &[f64], p: u32, measure: f64) -> f64 {
        let fact = |n: u32| (1..=n).map(f64::from).product::<f64>();
        let d = s.len() as u32 - 1;
        measure * fact(p) * fact(d) / fact(p + d) * complete(s, p)
    }

    /// Les Nodes `v` (sommets) suivis des milieux des arêtes `edges`.
    pub(super) fn with_mids(v: &[Node], edges: &[(usize, usize)]) -> Vec<Node> {
        let mut nodes = v.to_vec();
        nodes.extend(edges.iter().map(|&(a, b)| {
            Node::new(
                (v[a].x + v[b].x) / 2.0,
                (v[a].y + v[b].y) / 2.0,
                (v[a].z + v[b].z) / 2.0,
            )
        }));
        nodes
    }

    /// C symétrique, Σ C = ρc V, et Tᵀ C T = ρc ∫ φ² dV pour φ = g·x (linéaire) ou (g·x)²
    /// (quadratique), que l'Element représente exactement.
    fn check(mesh: &Mesh, vertices: usize, volume: f64, quadratic: bool) {
        let c = element_capacity(mesh, 0);
        let g = [3.0, -2.0, 0.5];
        let s = |n: &Node| g[0] * n.x + g[1] * n.y + g[2] * n.z;
        let power = if quadratic { 2 } else { 1 };
        let t: Vec<f64> = mesh.nodes.iter().map(|n| s(n).powi(power)).collect();
        let mut total = 0.0;
        let mut energy = 0.0;
        for (i, row) in c.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!((v - c[j][i]).abs() < 1e-14 * c[0][0]);
                total += v;
                energy += t[i] * v * t[j];
            }
        }
        assert!(
            (total - RHO_C * volume).abs() < 1e-13 * volume,
            "Σ C = {total}"
        );
        let corners: Vec<f64> = mesh.nodes[..vertices].iter().map(s).collect();
        let expected = RHO_C * simplex_power(&corners, 2 * power as u32, volume);
        assert!(
            (energy - expected).abs() < 1e-12 * expected,
            "{energy} au lieu de {expected}"
        );
    }

    #[test]
    fn capacity_matrices_are_exact() {
        let m = material();
        let tri = [
            Node::new_2d(0.0, 0.0),
            Node::new_2d(2.0, 0.0),
            Node::new_2d(0.0, 3.0),
        ]; // aire 3, épaisseur 0,5
        let tri3 = Tri3::new([0, 1, 2], m, 0.5, PlaneState::Stress).unwrap();
        check(
            &Mesh::new(tri.to_vec(), vec![Element::Tri3(tri3)]).unwrap(),
            3,
            1.5,
            false,
        );
        let tri6 = Tri6::new([0, 1, 2, 3, 4, 5], m, 0.5, PlaneState::Stress).unwrap();
        check(
            &Mesh::new(with_mids(&tri, &TRI6_EDGES), vec![Element::Tri6(tri6)]).unwrap(),
            3,
            1.5,
            true,
        );

        let tet = [
            Node::new(0.0, 0.0, 0.0),
            Node::new(2.0, 0.0, 0.0),
            Node::new(0.0, 3.0, 0.0),
            Node::new(0.0, 0.0, 1.0),
        ]; // volume 1
        let tet4 = Element::Tet4(Tet4::new([0, 1, 2, 3], m));
        check(&Mesh::new(tet.to_vec(), vec![tet4]).unwrap(), 4, 1.0, false);
        let tet10 = Element::Tet10(Tet10::new(std::array::from_fn(|k| k), m));
        check(
            &Mesh::new(with_mids(&tet, &TET10_EDGES), vec![tet10]).unwrap(),
            4,
            1.0,
            true,
        );
    }
}

#[cfg(test)]
mod response_tests {
    use super::capacity_tests::with_mids;
    use super::*;
    use crate::mesh::Node;
    use crate::shape::{TET10_EDGES, TRI6_EDGES};

    /// Un Element de chaque type (triangles en déformations planes), de forme quelconque.
    fn meshes(material: Material) -> Vec<(&'static str, Mesh)> {
        let tri = [
            Node::new_2d(0.0, 0.0),
            Node::new_2d(2.0, 0.3),
            Node::new_2d(0.4, 1.8),
        ];
        let tet = [
            Node::new(1.0, 0.5, 0.0),
            Node::new(3.0, 1.0, 0.5),
            Node::new(1.5, 2.5, 0.0),
            Node::new(1.2, 1.0, 2.0),
        ];
        let strain = PlaneState::Strain;
        let tri3 = Tri3::new([0, 1, 2], material, 0.5, strain).unwrap();
        let tri6 = Tri6::new(std::array::from_fn(|k| k), material, 0.5, strain).unwrap();
        vec![
            ("Tri3", Mesh::new(tri.to_vec(), vec![Element::Tri3(tri3)])),
            (
                "Tri6",
                Mesh::new(with_mids(&tri, &TRI6_EDGES), vec![Element::Tri6(tri6)]),
            ),
            (
                "Tet4",
                Mesh::new(
                    tet.to_vec(),
                    vec![Element::Tet4(Tet4::new([0, 1, 2, 3], material))],
                ),
            ),
            (
                "Tet10",
                Mesh::new(
                    with_mids(&tet, &TET10_EDGES),
                    vec![Element::Tet10(Tet10::new(
                        std::array::from_fn(|k| k),
                        material,
                    ))],
                ),
            ),
        ]
        .into_iter()
        .map(|(name, mesh)| (name, mesh.unwrap()))
        .collect()
    }

    /// Déplacements de l'ordre de 3·10⁻³ fois la taille de l'Element, sans motif : déformations de
    /// quelques 10⁻³, au-delà de σy/E = 1,2·10⁻³.
    fn displacements(mesh: &Mesh) -> Vec<f64> {
        (0..mesh.nodes.len() * mesh.dofs_per_node())
            .map(|k| 3e-3 * (((k * 7) % 5) as f64 - 2.0 + 0.1 * k as f64))
            .collect()
    }

    fn states(mesh: &Mesh) -> Vec<PlasticState> {
        vec![PlasticState::default(); element_gauss_count(mesh, 0)]
    }

    #[test]
    fn elastic_response_is_the_linear_stiffness() {
        let elastic = Material::with_poisson(210e3, 0.3).unwrap();
        for (name, mesh) in meshes(elastic) {
            let u = displacements(&mesh);
            let r = element_response(&mesh, 0, &u, &states(&mesh));
            let k = element_stiffness(&mesh, 0);
            let scale = k[0][0];
            assert!(!r.yielding, "{name}");
            for ((row, tangent_row), f) in k.iter().zip(&r.tangent).zip(&r.internal_force) {
                let ku: f64 = row.iter().zip(&u).map(|(k, u)| k * u).sum();
                assert!(
                    (f - ku).abs() < 1e-9 * scale * 1e-2,
                    "{name} : {f} au lieu de {ku}"
                );
                for (a, b) in row.iter().zip(tangent_row) {
                    assert!((a - b).abs() < 1e-9 * scale, "{name}");
                }
            }
            for (got, expected) in r.stresses.iter().zip(element_stresses(&mesh, 0, &u, &[])) {
                for (a, b) in [
                    (got.sigma_xx, expected.sigma_xx),
                    (got.sigma_yy, expected.sigma_yy),
                    (got.sigma_zz, expected.sigma_zz),
                    (got.tau_xy, expected.tau_xy),
                    (got.tau_yz, expected.tau_yz),
                    (got.tau_xz, expected.tau_xz),
                    (got.x, expected.x),
                    (got.y, expected.y),
                    (got.z, expected.z),
                ] {
                    assert!((a - b).abs() < 1e-9 * 210e3 * 1e-2, "{name} : {got:?}");
                }
            }
        }
    }

    #[test]
    fn plastic_tangent_is_the_derivative_of_the_internal_forces() {
        let plastic = Material::with_poisson(210e3, 0.3)
            .and_then(|m| m.with_plasticity(250.0, 2000.0))
            .unwrap();
        for (name, mesh) in meshes(plastic) {
            let u = displacements(&mesh);
            let previous = states(&mesh);
            let r = element_response(&mesh, 0, &u, &previous);
            assert!(r.yielding, "{name}");
            let scale = r.tangent[0][0].abs();
            let h = 1e-9;
            for j in 0..u.len() {
                let (mut plus, mut minus) = (u.clone(), u.clone());
                plus[j] += h;
                minus[j] -= h;
                let fp = element_response(&mesh, 0, &plus, &previous).internal_force;
                let fm = element_response(&mesh, 0, &minus, &previous).internal_force;
                for (i, row) in r.tangent.iter().enumerate() {
                    let numeric = (fp[i] - fm[i]) / (2.0 * h);
                    assert!(
                        (row[j] - numeric).abs() < 1e-6 * scale,
                        "{name} : K[{i}][{j}] = {} au lieu de {numeric}",
                        row[j]
                    );
                }
            }
        }
    }

    #[test]
    fn large_tangent_is_the_derivative_of_the_internal_forces() {
        let material = Material::with_poisson(210e3, 0.3).unwrap();
        for (name, mesh) in meshes(material) {
            // déformations de l'ordre de 10 à 20 % : loin du domaine linéaire
            let u: Vec<f64> = displacements(&mesh).iter().map(|u| 40.0 * u).collect();
            let r = element_large_response(&mesh, 0, &u);
            let scale = r.tangent[0][0].abs();
            let h = 1e-6;
            for j in 0..u.len() {
                let (mut plus, mut minus) = (u.clone(), u.clone());
                plus[j] += h;
                minus[j] -= h;
                let fp = element_large_response(&mesh, 0, &plus).internal_force;
                let fm = element_large_response(&mesh, 0, &minus).internal_force;
                for (i, row) in r.tangent.iter().enumerate() {
                    let numeric = (fp[i] - fm[i]) / (2.0 * h);
                    assert!(
                        (row[j] - numeric).abs() < 1e-6 * scale,
                        "{name} : K[{i}][{j}] = {} au lieu de {numeric}",
                        row[j]
                    );
                    assert!((row[j] - r.tangent[j][i]).abs() < 1e-9 * scale, "{name}");
                }
            }
        }
    }

    #[test]
    fn large_response_tends_to_the_linear_one_for_small_displacements() {
        let material = Material::with_poisson(210e3, 0.3).unwrap();
        for (name, mesh) in meshes(material) {
            let u: Vec<f64> = displacements(&mesh).iter().map(|u| 1e-6 * u).collect();
            let r = element_large_response(&mesh, 0, &u);
            let k = element_stiffness(&mesh, 0);
            let scale = k[0][0];
            for ((row, tangent_row), f) in k.iter().zip(&r.tangent).zip(&r.internal_force) {
                let ku: f64 = row.iter().zip(&u).map(|(k, u)| k * u).sum();
                assert!((f - ku).abs() < 1e-5 * scale * 1e-5, "{name}");
                for (a, b) in row.iter().zip(tangent_row) {
                    assert!((a - b).abs() < 1e-5 * scale, "{name}");
                }
            }
        }
    }
}

#[cfg(test)]
mod pressure_tests {
    use super::*;
    use crate::mesh::Node;
    use crate::shape::TET10_EDGES;

    fn material() -> Material {
        Material::with_poisson(1.0, 0.3).unwrap()
    }

    fn close(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len(), "{a:?} contre {b:?}");
        for (x, y) in a.iter().zip(b) {
            assert!((x - y).abs() < 1e-12, "{a:?} contre {b:?}");
        }
    }

    /// Triangle rectangle de côté 2, d'épaisseur 3, sur le Mesh à un Element.
    fn triangle(quadratic: bool) -> Mesh {
        let mut nodes = vec![
            Node::new_2d(0.0, 0.0),
            Node::new_2d(2.0, 0.0),
            Node::new_2d(0.0, 2.0),
        ];
        let element = if quadratic {
            nodes.extend([
                Node::new_2d(1.0, 0.0),
                Node::new_2d(1.0, 1.0),
                Node::new_2d(0.0, 1.0),
            ]);
            Element::Tri6(
                Tri6::new([0, 1, 2, 3, 4, 5], material(), 3.0, PlaneState::Stress).unwrap(),
            )
        } else {
            Element::Tri3(Tri3::new([0, 1, 2], material(), 3.0, PlaneState::Stress).unwrap())
        };
        Mesh::new(nodes, vec![element]).unwrap()
    }

    fn tetrahedron(quadratic: bool) -> Mesh {
        let vertices = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
        let mut nodes: Vec<Node> = vertices
            .iter()
            .map(|p| Node::new(p[0], p[1], p[2]))
            .collect();
        let element = if quadratic {
            for (a, b) in TET10_EDGES {
                nodes.push(Node::new(
                    (vertices[a][0] + vertices[b][0]) / 2.0,
                    (vertices[a][1] + vertices[b][1]) / 2.0,
                    (vertices[a][2] + vertices[b][2]) / 2.0,
                ));
            }
            Element::Tet10(Tet10::new([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], material()))
        } else {
            Element::Tet4(Tet4::new([0, 1, 2, 3], material()))
        };
        Mesh::new(nodes, vec![element]).unwrap()
    }

    #[test]
    fn tri3_edge_load_is_half_the_resultant_at_each_node() {
        // arête du bas (face 0), L = 2, t = 3, p = 5 : la pression pousse vers le haut, 30 en tout
        let f = face_pressure(&triangle(false), 0, 0, 5.0).unwrap();
        assert_eq!((f.positions.clone(), f.dim), (vec![0, 1], 2));
        close(&f.forces, &[0.0, 15.0, 0.0, 15.0]);
    }

    #[test]
    fn tri6_edge_load_is_one_sixth_two_thirds_one_sixth() {
        let f = face_pressure(&triangle(true), 0, 0, 5.0).unwrap();
        close(&f.forces, &[0.0, 5.0, 0.0, 5.0, 0.0, 20.0]);
    }

    #[test]
    fn the_normal_points_out_on_every_edge() {
        // arête 1 (hypoténuse, normale (1, 1)/√2) : la force est −p L t n, L = 2√2, t = 3
        let f = face_pressure(&triangle(false), 0, 1, 2.0).unwrap();
        let half = -2.0 * 2.0 * 2.0_f64.sqrt() * 3.0 / 2.0 / 2.0_f64.sqrt();
        close(&f.forces, &[half, half, half, half]);
        // arête 2 (côté gauche, normale (−1, 0)) : poussée vers +x
        let f = face_pressure(&triangle(false), 0, 2, 2.0).unwrap();
        close(&f.forces, &[6.0, 0.0, 6.0, 0.0]);
    }

    #[test]
    fn tet4_face_load_is_a_third_of_the_resultant_at_each_node() {
        // face 0 (z = 0), A = 1/2, p = 6 : 3 vers +z, 1 par Node
        let f = face_pressure(&tetrahedron(false), 0, 0, 6.0).unwrap();
        assert_eq!(f.dim, 3);
        close(&f.forces, &[0., 0., 1., 0., 0., 1., 0., 0., 1.]);
    }

    #[test]
    fn tet10_face_load_goes_to_the_mid_nodes_only() {
        // la fonction de forme d'un sommet de Tet10 s'intègre à 0 sur la face, un milieu à A/3
        let f = face_pressure(&tetrahedron(true), 0, 0, 6.0).unwrap();
        close(
            &f.forces,
            &[
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 1.,
            ],
        );
    }

    #[test]
    fn a_slanted_face_gets_the_outward_resultant() {
        // face 2 (1-3-2, plan x + y + z = 1) : A = √3/2, n = (1, 1, 1)/√3, résultante −p A n = −3 (1, 1, 1)
        let f = face_pressure(&tetrahedron(false), 0, 2, 6.0).unwrap();
        let total: Vec<f64> = (0..3)
            .map(|c| (0..3).map(|k| f.forces[3 * k + c]).sum())
            .collect();
        close(&total, &[-3.0, -3.0, -3.0]);
    }

    #[test]
    fn missing_faces_and_one_dimensional_elements_give_none() {
        assert!(face_pressure(&triangle(false), 0, 3, 1.0).is_none());
        let bar = Mesh::new(
            vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
            vec![Element::Bar(crate::mesh::Bar {
                nodes: [0, 1],
                material: Material::new(1.0).unwrap(),
                section: crate::material::Section::new(1.0).unwrap(),
            })],
        )
        .unwrap();
        assert!(face_pressure(&bar, 0, 0, 1.0).is_none());
    }
}
