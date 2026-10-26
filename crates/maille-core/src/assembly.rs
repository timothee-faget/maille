use crate::element::{
    bar_stiffness, beam_stiffness, capacity, conduction, heat_fluxes, plane_d, plane_gauss,
    plane_stiffness, plane_stresses, solid_d, solid_gauss, solid_stiffness, solid_stresses,
    GaussData,
};
use crate::material::Material;
use crate::mesh::{Element, Mesh, PlaneState, Tet10, Tet4, Tri3, Tri6};
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
/// Règle de Gauss d'un degré au-dessus de celle de la conduction : N Nᵀ y est intégré exactement.
/// Vide pour les Bar et Beam.
pub fn element_capacity(mesh: &Mesh, i: usize) -> Vec<Vec<f64>> {
    let e = &mesh.elements[i];
    let rho_c = heat_capacity(e.material());
    match e {
        Element::Tri3(t) => {
            let kernel = PlaneKernel {
                rule: &TRIANGLE_3,
                ..PlaneKernel::tri3(t)
            };
            capacity(&kernel.gauss(mesh), rho_c)
        }
        Element::Tri6(t) => {
            let kernel = PlaneKernel {
                rule: &TRIANGLE_6,
                ..PlaneKernel::tri6(t)
            };
            capacity(&kernel.gauss(mesh), rho_c)
        }
        Element::Tet4(t) => {
            let kernel = SolidKernel {
                rule: &TETRA_4,
                ..SolidKernel::tet4(t)
            };
            capacity(&kernel.gauss(mesh), rho_c)
        }
        Element::Tet10(t) => {
            let kernel = SolidKernel {
                rule: &TETRA_14,
                ..SolidKernel::tet10(t)
            };
            capacity(&kernel.gauss(mesh), rho_c)
        }
        Element::Bar(_) | Element::Beam(_) => Vec::new(),
    }
}

/// Flux de chaleur aux points de Gauss de l'Element `i`, à partir des températures `t` des Nodes.
/// Vide pour les Bar et Beam.
pub fn element_heat_fluxes(mesh: &Mesh, i: usize, t: &[f64]) -> Vec<HeatFlux> {
    element_heat_flux_operator(mesh, i)(t)
}

/// Flux de chaleur aux points de Gauss d'un Element en fonction des températures globales des Nodes.
pub type FluxOperator = Box<dyn Fn(&[f64]) -> Vec<HeatFlux>>;

/// Comme `element_heat_fluxes`, mais les points de Gauss sont calculés une seule fois : l'opérateur
/// renvoyé s'applique à autant de champs de températures `t` qu'on veut (un par incrément en
/// transitoire).
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
    let e = &mesh.elements[i];
    let positions = e.face_positions(face)?;
    let global: Vec<usize> = positions.iter().map(|&k| e.nodes()[k]).collect();
    // (poids w |dS|, valeurs N) aux points de Gauss de la face
    let points: Vec<(f64, Vec<f64>)> = match e {
        Element::Tri3(Tri3 { thickness, .. }) | Element::Tri6(Tri6 { thickness, .. }) => {
            let c = mesh.coords(&global);
            let length = (c[1][0] - c[0][0]).hypot(c[1][1] - c[0][1]);
            let linear = global.len() == 2;
            let rule: &[LinePoint] = if linear { &LINE_2 } else { &LINE_3 };
            let shape = if linear { line2 } else { line3 };
            rule.iter()
                .map(|p| (p.weight * length * thickness, shape(p.s)))
                .collect()
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
            rule.iter()
                .map(|g| (g.weight * twice_area, shape(g.xi, g.eta).0))
                .collect()
        }
        Element::Bar(_) | Element::Beam(_) => return None,
    };
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
    fn with_mids(v: &[Node], edges: &[(usize, usize)]) -> Vec<Node> {
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
