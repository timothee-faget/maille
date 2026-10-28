//! Grandes déformations (Lagrangien total, Saint-Venant–Kirchhoff et néo-Hookéen) : objectivité,
//! cas homogènes de forme fermée, limite linéaire, convergence quadratique de Newton-Raphson,
//! pression suiveuse, découpe de l'incrément, elastica et tube épais gonflé.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{Amplitude, Bc, Kinematics, Load, Model, Pressure, Step};
use maille_core::results::{Increment, Stress};
use maille_core::MailleError;

const E: f64 = 210e3;
const NU: f64 = 0.3;
const L: f64 = 2.0;
const W: f64 = 1.0;
const D: f64 = 1.0;

fn material() -> Material {
    Material::with_poisson(E, NU).unwrap()
}

fn plane_meshes() -> Vec<(&'static str, Mesh)> {
    let (nodes, tris) = rect_grid(4, 2, L, W);
    let q = to_quadratic(&nodes, &tris);
    let strain = PlaneState::Strain;
    let tri3 = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, material(), 1.0, strain).unwrap()))
        .collect();
    let tri6 = q
        .elements
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material(), 1.0, strain).unwrap()))
        .collect();
    vec![
        ("Tri3", Mesh::new(nodes, tri3).unwrap()),
        ("Tri6", Mesh::new(q.nodes, tri6).unwrap()),
    ]
}

fn solid_meshes() -> Vec<(&'static str, Mesh)> {
    let (nodes, tets) = box_grid([3, 2, 2], [L, W, D]);
    let q = to_quadratic_tet(&nodes, &tets);
    let tet4 = tets
        .iter()
        .map(|&n| Element::Tet4(Tet4::new(n, material())))
        .collect();
    let tet10 = q
        .elements
        .iter()
        .map(|&n| Element::Tet10(Tet10::new(n, material())))
        .collect();
    vec![
        ("Tet4", Mesh::new(nodes, tet4).unwrap()),
        ("Tet10", Mesh::new(q.nodes, tet10).unwrap()),
    ]
}

fn all_meshes() -> Vec<(&'static str, Mesh)> {
    let mut meshes = plane_meshes();
    meshes.extend(solid_meshes());
    meshes
}

fn bc(node: usize, dof: Dof, value: f64) -> Bc {
    Bc { node, dof, value }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9 * L
}

fn dofs(mesh: &Mesh) -> &'static [Dof] {
    if mesh.dofs_per_node() == 2 {
        &[Dof::Ux, Dof::Uy]
    } else {
        &[Dof::Ux, Dof::Uy, Dof::Uz]
    }
}

fn large(mesh: Mesh, bcs: Vec<Bc>, loads: Vec<Load>) -> Model {
    Model::new(mesh, bcs, loads)
        .unwrap()
        .with_kinematics(Kinematics::Large)
        .unwrap()
}

fn solve(model: &Model, increments: usize) -> Vec<Increment> {
    let step = Step::quasi_static(1.0, increments).unwrap();
    model
        .solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap())
        .unwrap()
}

fn max_stress(stresses: &[Vec<Stress>]) -> f64 {
    stresses
        .iter()
        .flatten()
        .flat_map(|s| {
            [
                s.sigma_xx, s.sigma_yy, s.sigma_zz, s.tau_xy, s.tau_yz, s.tau_xz,
            ]
        })
        .fold(0.0, |m: f64, v| m.max(v.abs()))
}

/// Impose à chaque Node le déplacement `map(X) − X` : toute la pièce suit la transformation.
fn imposed_map(mesh: &Mesh, map: impl Fn(&Node) -> [f64; 3]) -> Vec<Bc> {
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        let x = map(n);
        let reference = [n.x, n.y, n.z];
        for (k, &dof) in dofs(mesh).iter().enumerate() {
            bcs.push(bc(i, dof, x[k] - reference[k]));
        }
    }
    bcs
}

#[test]
fn a_rigid_rotation_of_90_degrees_leaves_no_stress() {
    for (name, mesh) in all_meshes() {
        // (x, y, z) -> (−y, x, z) : un quart de tour autour de z, puis pour les solides un quart de
        // tour autour de x
        let solid = mesh.dofs_per_node() == 3;
        let bcs = imposed_map(&mesh, |n| {
            let (x, y, z) = (-n.y, n.x, n.z);
            if solid {
                [x, -z, y]
            } else {
                [x, y, 0.0]
            }
        });
        let model = large(mesh, bcs, Vec::new());
        let history = solve(&model, 1);
        let field = &history[1].field;
        let worst = max_stress(&field.stresses);
        assert!(worst < 1e-10 * E, "{name} : σ = {worst}");
        // les 3 Nodes imposés ne reçoivent aucune réaction
        let reaction = field.reactions.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
        assert!(reaction < 1e-10 * E, "{name} : réaction {reaction}");
    }
}

/// Traction homogène exacte de Saint-Venant–Kirchhoff, côtés latéraux libres : renvoie
/// (λ latéral, σxx de Cauchy) pour l'allongement λx, en déformations planes (`solid` faux) ou 3D.
fn uniaxial(stretch: f64, solid: bool) -> (f64, f64) {
    let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let mu = E / (2.0 * (1.0 + NU));
    let exx = (stretch * stretch - 1.0) / 2.0;
    // S latérale nulle : plan, λ (Exx + Eyy) + 2μ Eyy = 0 ; solide, λ (Exx + 2 Eyy) + 2μ Eyy = 0
    let eyy = if solid {
        -lambda * exx / (2.0 * (lambda + mu))
    } else {
        -lambda * exx / (lambda + 2.0 * mu)
    };
    let lateral = (1.0 + 2.0 * eyy).sqrt();
    let trace = exx + if solid { 2.0 * eyy } else { eyy };
    let sxx = lambda * trace + 2.0 * mu * exx;
    let j = stretch * lateral * if solid { lateral } else { 1.0 };
    (lateral, stretch * stretch * sxx / j)
}

#[test]
fn a_prescribed_uniaxial_stretch_gives_the_closed_form_stress() {
    let stretch = 1.4;
    for (name, mesh) in all_meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let (lateral, sxx) = uniaxial(stretch, solid);
        let bcs = imposed_map(&mesh, |n| [stretch * n.x, lateral * n.y, lateral * n.z]);
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        for s in history[1].field.stresses.iter().flatten() {
            assert!((s.sigma_xx - sxx).abs() < 1e-9 * sxx, "{name} : {s:?}");
            for other in [s.sigma_yy, s.tau_xy, s.tau_yz, s.tau_xz] {
                assert!(other.abs() < 1e-9 * sxx, "{name} : {s:?}");
            }
            if solid {
                assert!(s.sigma_zz.abs() < 1e-9 * sxx, "{name} : {s:?}");
            }
        }
    }
}

/// Traction par déplacement imposé en x = L, appuis glissants (ux = 0 en x = 0, uy = 0 en y = 0,
/// uz = 0 en z = 0) : la contraction latérale se fait librement.
fn free_lateral_tension(mesh: &Mesh, delta: f64) -> Vec<Bc> {
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if near(n.x, 0.0) {
            bcs.push(bc(i, Dof::Ux, 0.0));
        } else if near(n.x, L) {
            bcs.push(bc(i, Dof::Ux, delta));
        }
        if near(n.y, 0.0) {
            bcs.push(bc(i, Dof::Uy, 0.0));
        }
        if mesh.dofs_per_node() == 3 && near(n.z, 0.0) {
            bcs.push(bc(i, Dof::Uz, 0.0));
        }
    }
    bcs
}

#[test]
fn newton_finds_the_uniaxial_closed_form_with_free_lateral_sides() {
    let stretch = 1.5;
    for (name, mesh) in all_meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let (_, sxx) = uniaxial(stretch, solid);
        let bcs = free_lateral_tension(&mesh, (stretch - 1.0) * L);
        let history = solve(&large(mesh, bcs, Vec::new()), 5);
        for s in history[5].field.stresses.iter().flatten() {
            assert!((s.sigma_xx - sxx).abs() < 1e-6 * sxx, "{name} : {s:?}");
            assert!(s.sigma_yy.abs() < 1e-6 * sxx, "{name} : {s:?}");
        }
    }
}

#[test]
fn newton_converges_quadratically() {
    for (name, mesh) in all_meshes() {
        let bcs = free_lateral_tension(&mesh, 0.5 * L);
        let history = solve(&large(mesh, bcs, Vec::new()), 3);
        let mut longest = 0;
        for inc in &history[1..] {
            let r = &inc.residuals;
            longest = longest.max(r.len());
            assert!(r.len() <= 8, "{name}, t = {} : {r:?}", inc.time);
            // r[0] : déséquilibre avant d'appliquer l'incrément imposé (nul ici)
            for pair in r[1..].windows(2) {
                if pair[0] < 1e-2 {
                    assert!(
                        pair[1] <= (10.0 * pair[0] * pair[0]).max(1e-12),
                        "{name}, t = {} : {r:?}",
                        inc.time
                    );
                }
            }
        }
        assert!(
            longest >= 4,
            "{name} : le test doit exercer plusieurs itérations"
        );
    }
}

#[test]
fn a_small_load_gives_the_linear_solution() {
    let eps = 1e-6;
    for (name, mesh) in all_meshes() {
        let bcs = free_lateral_tension(&mesh, eps * L);
        let linear = Model::new(mesh.clone(), bcs.clone(), Vec::new())
            .unwrap()
            .solve()
            .unwrap();
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        let field = &history[1].field;
        let scale = E * eps;
        for (a, b) in field.displacements.iter().zip(&linear.displacements) {
            assert!(
                (a - b).abs() < 1e-5 * eps * L,
                "{name} : {a} au lieu de {b}"
            );
        }
        for (a, b) in field
            .stresses
            .iter()
            .flatten()
            .zip(linear.stresses.iter().flatten())
        {
            assert!((a.sigma_xx - b.sigma_xx).abs() < 1e-5 * scale, "{name}");
            assert!((a.sigma_yy - b.sigma_yy).abs() < 1e-5 * scale, "{name}");
        }
    }
}

#[test]
fn large_strain_inputs_are_typed_errors() {
    let tri = || {
        vec![
            Node::new_2d(0.0, 0.0),
            Node::new_2d(1.0, 0.0),
            Node::new_2d(0.0, 1.0),
        ]
    };
    let tri_mesh = |material: Material, state: PlaneState| {
        let t = Tri3::new([0, 1, 2], material, 1.0, state).unwrap();
        Mesh::new(tri(), vec![Element::Tri3(t)]).unwrap()
    };
    let bcs = || {
        vec![
            bc(0, Dof::Ux, 0.0),
            bc(0, Dof::Uy, 0.0),
            bc(1, Dof::Uy, 0.0),
        ]
    };
    let strain = PlaneState::Strain;
    let model = |m: Mesh| Model::new(m, bcs(), Vec::new()).unwrap();

    // contraintes planes
    assert_eq!(
        model(tri_mesh(material(), PlaneState::Stress)).with_kinematics(Kinematics::Large),
        Err(MailleError::UnsupportedLargeStrainElement { element: 0 })
    );
    // plasticité
    let plastic = material().with_plasticity(250.0, 2000.0).unwrap();
    assert_eq!(
        model(tri_mesh(plastic, strain)).with_kinematics(Kinematics::Large),
        Err(MailleError::LargeStrainPlasticity { element: 0 })
    );
    // Bar
    let bar = Bar {
        nodes: [0, 1],
        material: Material::new(E).unwrap(),
        section: Section::new(1.0).unwrap(),
    };
    let bar_mesh = Mesh::new(
        vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
        vec![Element::Bar(bar)],
    )
    .unwrap();
    assert_eq!(
        Model::new(bar_mesh, vec![bc(0, Dof::U, 0.0)], Vec::new())
            .unwrap()
            .with_kinematics(Kinematics::Large),
        Err(MailleError::UnsupportedLargeStrainElement { element: 0 })
    );
    // températures, dans les deux ordres
    let expandable = material().with_expansion(1e-5).unwrap();
    let hot = || {
        model(tri_mesh(expandable, strain))
            .with_temperatures(vec![10.0; 3], 0.0)
            .unwrap()
    };
    assert_eq!(
        hot().with_kinematics(Kinematics::Large),
        Err(MailleError::LargeStrainWithTemperatures)
    );
    let large_model = model(tri_mesh(expandable, strain))
        .with_kinematics(Kinematics::Large)
        .unwrap();
    assert_eq!(
        large_model.clone().with_temperatures(vec![10.0; 3], 0.0),
        Err(MailleError::LargeStrainWithTemperatures)
    );
    // solve linéaire
    assert_eq!(
        large_model.solve(),
        Err(MailleError::LargeStrainNeedsQuasiStatic)
    );
}

/// Néo-Hookéen compressible : σ_i = (μ (λ_i² − 1) + Λ ln J) / J pour un étirement principal λ_i
/// (Λ, μ : coefficients de Lamé).
mod neo_hooke {
    use super::*;

    pub fn lame() -> (f64, f64) {
        (
            E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU)),
            E / (2.0 * (1.0 + NU)),
        )
    }

    /// Contrainte de Cauchy principale pour l'étirement `stretch` et J.
    pub fn principal(stretch: f64, j: f64) -> f64 {
        let (lambda, mu) = lame();
        (mu * (stretch * stretch - 1.0) + lambda * j.ln()) / j
    }

    /// Bissection de `f` sur [low, high], où `f` change de signe.
    pub fn bisect(f: impl Fn(f64) -> f64, mut low: f64, mut high: f64) -> f64 {
        let sign = f(low).signum();
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if f(middle).signum() == sign {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    }

    /// Traction uniaxiale à côtés libres : (étirement latéral, σxx) pour l'allongement `stretch`.
    pub fn uniaxial(stretch: f64, solid: bool) -> (f64, f64) {
        let power = if solid { 2 } else { 1 };
        let j = |lateral: f64| stretch * lateral.powi(power);
        let lateral = bisect(|l| neo_hooke::principal(l, j(l)), 0.01, 10.0);
        (lateral, principal(stretch, j(lateral)))
    }

    pub fn material() -> Material {
        super::material().with_hyperelasticity(maille_core::material::Hyperelasticity::NeoHookean)
    }

    pub fn meshes() -> Vec<(&'static str, Mesh)> {
        let swap = |m: Mesh| {
            let elements = m
                .elements
                .iter()
                .map(|e| match e {
                    Element::Tri3(t) => Element::Tri3(Tri3 {
                        material: material(),
                        ..*t
                    }),
                    Element::Tri6(t) => Element::Tri6(Tri6 {
                        material: material(),
                        ..*t
                    }),
                    Element::Tet4(t) => Element::Tet4(Tet4 {
                        material: material(),
                        ..*t
                    }),
                    Element::Tet10(t) => Element::Tet10(Tet10 {
                        material: material(),
                        ..*t
                    }),
                    other => *other,
                })
                .collect();
            Mesh::new(m.nodes, elements).unwrap()
        };
        super::all_meshes()
            .into_iter()
            .map(|(name, mesh)| (name, swap(mesh)))
            .collect()
    }
}

#[test]
fn neo_hookean_simple_shear_matches_the_closed_form() {
    // F = [[1, γ], [0, 1]] : J = 1, σ = μ (b − I), b = F Fᵀ
    let gamma = 0.8;
    let (_, mu) = neo_hooke::lame();
    for (name, mesh) in neo_hooke::meshes() {
        let bcs = imposed_map(&mesh, |n| [n.x + gamma * n.y, n.y, n.z]);
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        for s in history[1].field.stresses.iter().flatten() {
            let expected = [
                (s.sigma_xx, mu * gamma * gamma),
                (s.tau_xy, mu * gamma),
                (s.sigma_yy, 0.0),
                (s.sigma_zz, 0.0),
                (s.tau_yz, 0.0),
                (s.tau_xz, 0.0),
            ];
            for (got, want) in expected {
                assert!((got - want).abs() < 1e-9 * mu, "{name} : {s:?}");
            }
        }
    }
}

#[test]
fn neo_hookean_uniaxial_tension_matches_the_closed_form() {
    let stretch = 2.0;
    for (name, mesh) in neo_hooke::meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let (lateral, sxx) = neo_hooke::uniaxial(stretch, solid);
        // étirement imposé : σyy est nulle à l'arrondi près
        let bcs = imposed_map(&mesh, |n| [stretch * n.x, lateral * n.y, lateral * n.z]);
        let history = solve(&large(mesh.clone(), bcs, Vec::new()), 4);
        for s in history[4].field.stresses.iter().flatten() {
            assert!((s.sigma_xx - sxx).abs() < 1e-8 * sxx, "{name} : {s:?}");
            assert!(s.sigma_yy.abs() < 1e-8 * sxx, "{name} : {s:?}");
        }
        // côtés libres, par Newton : même solution
        let bcs = free_lateral_tension(&mesh, (stretch - 1.0) * L);
        let history = solve(&large(mesh, bcs, Vec::new()), 8);
        for s in history[8].field.stresses.iter().flatten() {
            assert!((s.sigma_xx - sxx).abs() < 1e-6 * sxx, "{name} : {s:?}");
            assert!(s.sigma_yy.abs() < 1e-6 * sxx, "{name} : {s:?}");
        }
    }
}

#[test]
fn a_neo_hookean_rigid_rotation_leaves_no_stress() {
    for (name, mesh) in neo_hooke::meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let bcs = imposed_map(&mesh, |n| {
            let (x, y, z) = (-n.y, n.x, n.z);
            if solid {
                [x, -z, y]
            } else {
                [x, y, 0.0]
            }
        });
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        let worst = max_stress(&history[1].field.stresses);
        assert!(worst < 1e-10 * E, "{name} : σ = {worst}");
    }
}

#[test]
fn a_small_load_on_a_neo_hookean_gives_the_linear_solution() {
    let eps = 1e-6;
    for (name, mesh) in neo_hooke::meshes() {
        let bcs = free_lateral_tension(&mesh, eps * L);
        // l'hyperélasticité n'existe qu'en grandes déformations : le linéaire prend le même Material
        // en SVK
        let (_, plain) = &all_meshes().into_iter().find(|(n, _)| *n == name).unwrap();
        let linear = Model::new(plain.clone(), bcs.clone(), Vec::new())
            .unwrap()
            .solve()
            .unwrap();
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        for (a, b) in history[1]
            .field
            .displacements
            .iter()
            .zip(&linear.displacements)
        {
            assert!(
                (a - b).abs() < 1e-5 * eps * L,
                "{name} : {a} au lieu de {b}"
            );
        }
    }
}

#[test]
fn a_neo_hookean_material_needs_large_strain() {
    for (_, mesh) in neo_hooke::meshes() {
        let model = Model::new(mesh, vec![bc(0, Dof::Ux, 0.0)], Vec::new()).unwrap();
        assert_eq!(
            model.solve(),
            Err(MailleError::HyperelasticityNeedsLargeStrain { element: 0 })
        );
        let step = Step::quasi_static(1.0, 1).unwrap();
        assert_eq!(
            model.solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap()),
            Err(MailleError::HyperelasticityNeedsLargeStrain { element: 0 })
        );
    }
}

/// Appuis de symétrie (ux = 0 en x = 0, uy = 0 en y = 0, uz = 0 en z = 0) et Pressure `p` sur les faces
/// x = L, y = W (et z = D en 3D) : la pièce se comprime d'un étirement homogène.
fn hydrostatic(mesh: &Mesh, p: f64) -> (Vec<Bc>, Vec<Pressure>) {
    let solid = mesh.dofs_per_node() == 3;
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if near(n.x, 0.0) {
            bcs.push(bc(i, Dof::Ux, 0.0));
        }
        if near(n.y, 0.0) {
            bcs.push(bc(i, Dof::Uy, 0.0));
        }
        if solid && near(n.z, 0.0) {
            bcs.push(bc(i, Dof::Uz, 0.0));
        }
    }
    let mut pressures = Vec::new();
    for (element, e) in mesh.elements.iter().enumerate() {
        for face in 0..if solid { 4 } else { 3 } {
            let corners: Vec<&Node> = e
                .face_nodes(face)
                .unwrap()
                .iter()
                .map(|&n| &mesh.nodes[n])
                .collect();
            let on = |c: fn(&Node) -> f64, at: f64| corners.iter().all(|n| near(c(n), at));
            if on(|n| n.x, L) || on(|n| n.y, W) || (solid && on(|n| n.z, D)) {
                pressures.push(Pressure {
                    element,
                    face,
                    value: p,
                });
            }
        }
    }
    (bcs, pressures)
}

#[test]
fn a_following_pressure_compresses_a_neo_hookean_homogeneously() {
    // J λ³ (ou λ² en déformations planes) : σ = −p exige (μ (λ² − 1) + Λ ln J) / J = −p
    let stretch: f64 = 0.8;
    for (name, mesh) in neo_hooke::meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let j = stretch.powi(if solid { 3 } else { 2 });
        let p = -neo_hooke::principal(stretch, j);
        let (bcs, pressures) = hydrostatic(&mesh, p);
        let model = large(mesh.clone(), bcs, Vec::new())
            .with_pressures(pressures)
            .unwrap();
        let history = solve(&model, 4);
        let field = &history[4].field;
        for (n, u) in mesh
            .nodes
            .iter()
            .zip(field.displacements.chunks(mesh.dofs_per_node()))
        {
            let at = [n.x, n.y, n.z];
            for (u, x) in u.iter().zip(at) {
                assert!(
                    (u - (stretch - 1.0) * x).abs() < 1e-6 * L,
                    "{name} : u = {u}"
                );
            }
        }
        for s in field.stresses.iter().flatten() {
            assert!((s.sigma_xx + p).abs() < 1e-6 * p, "{name} : {s:?}");
            assert!((s.sigma_yy + p).abs() < 1e-6 * p, "{name} : {s:?}");
            assert!(s.tau_xy.abs() < 1e-6 * p, "{name} : {s:?}");
        }
    }
}

#[test]
fn newton_cuts_the_increment_when_it_does_not_converge() {
    // étirement de 4 en un seul incrément : trop pour Newton depuis le repos, que la coupe sauve
    let stretch = 4.0;
    for (name, mesh) in neo_hooke::meshes() {
        let solid = mesh.dofs_per_node() == 3;
        let (_, sxx) = neo_hooke::uniaxial(stretch, solid);
        let bcs = free_lateral_tension(&mesh, (stretch - 1.0) * L);
        let history = solve(&large(mesh, bcs, Vec::new()), 1);
        assert_eq!(history.len(), 2, "{name} : un seul Increment gardé");
        for s in history[1].field.stresses.iter().flatten() {
            assert!((s.sigma_xx - sxx).abs() < 1e-6 * sxx, "{name} : {s:?}");
        }
    }
}

#[test]
fn a_collapse_beyond_every_cut_is_a_typed_error() {
    // retournement de la pièce : aucun incrément, même coupé, ne mène à l'équilibre
    for (name, mesh) in neo_hooke::meshes() {
        let bcs = free_lateral_tension(&mesh, -1.2 * L);
        let step = Step::quasi_static(1.0, 1).unwrap();
        let result =
            large(mesh, bcs, Vec::new()).solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap());
        assert!(
            matches!(
                result,
                Err(MailleError::NotConverged { .. } | MailleError::Singular)
            ),
            "{name} : {result:?}"
        );
    }
}

/// Console élancée en Tri6 (L = 20, h = 1, ν = 0) sous une charge verticale P en bout : flèche de
/// l'elastica (Bisshopp et Drucker 1945). Référence calculée ici : avec φ0 la pente en bout, la
/// première intégrale de EI θ'' = −P cos θ donne, pour α = P L² / EI,
/// α = ½ (∫₀^φ0 dθ / √(sin φ0 − sin θ))², w / L = ∫ sin θ dθ / √(…) / ∫ dθ / √(…),
/// et u / L = 1 − ∫ cos θ dθ / √(…) / ∫ dθ / √(…).
mod elastica {
    /// ∫₀^φ0 g(θ) dθ / √(sin φ0 − sin θ), par θ = φ0 (1 − v²) qui lève la singularité en φ0.
    fn integral(phi0: f64, g: impl Fn(f64) -> f64) -> f64 {
        let n = 20_000;
        let h = 1.0 / n as f64;
        (0..n)
            .map(|k| {
                let v = (k as f64 + 0.5) * h;
                let theta = phi0 * (1.0 - v * v);
                2.0 * phi0 * v * g(theta) / (phi0.sin() - theta.sin()).sqrt() * h
            })
            .sum()
    }

    /// (w / L, u / L) du bout pour α = P L² / EI.
    pub fn tip(alpha: f64) -> (f64, f64) {
        let of_phi = |phi0: f64| 0.5 * integral(phi0, |_| 1.0).powi(2);
        let (mut low, mut high) = (1e-3, 1.5);
        for _ in 0..100 {
            let middle = 0.5 * (low + high);
            if of_phi(middle) < alpha {
                low = middle;
            } else {
                high = middle;
            }
        }
        let phi0 = 0.5 * (low + high);
        let total = integral(phi0, |_| 1.0);
        (
            integral(phi0, f64::sin) / total,
            1.0 - integral(phi0, f64::cos) / total,
        )
    }
}

#[test]
fn the_elastica_reference_matches_the_published_values() {
    // PL²/EI = 1, 2, 5 : w/L = 0,3017 ; 0,4935 ; 0,7138 (Mattiasson 1981)
    for (alpha, w) in [(1.0, 0.3017), (2.0, 0.4935), (5.0, 0.7138)] {
        let (got, _) = elastica::tip(alpha);
        assert!((got - w).abs() < 1e-4, "α = {alpha} : w/L = {got}");
    }
}

#[test]
fn a_slender_cantilever_follows_the_elastica() {
    let (length, depth) = (20.0, 1.0);
    let (nodes, tris) = rect_grid(40, 2, length, depth);
    let q = to_quadratic(&nodes, &tris);
    let strain = PlaneState::Strain;
    let elements = q
        .elements
        .iter()
        .map(|&n| {
            Element::Tri6(
                Tri6::new(n, Material::with_poisson(E, 0.0).unwrap(), 1.0, strain).unwrap(),
            )
        })
        .collect();
    let mesh = Mesh::new(q.nodes, elements).unwrap();
    // EI par unité d'épaisseur = E h³ / 12 (ν = 0)
    let inertia_e = E * depth.powi(3) / 12.0;
    let mut bcs = Vec::new();
    let mut tip_nodes = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if n.x.abs() < 1e-9 {
            bcs.push(bc(i, Dof::Ux, 0.0));
            bcs.push(bc(i, Dof::Uy, 0.0));
        }
        if (n.x - length).abs() < 1e-9 {
            tip_nodes.push(i);
        }
    }
    let centre = *tip_nodes
        .iter()
        .find(|&&i| (mesh.nodes[i].y - depth / 2.0).abs() < 1e-9)
        .unwrap();
    for alpha in [1.0, 2.0, 5.0] {
        // P est répartie en bout : les 5 Nodes de l'arête reçoivent la force P (aire dans l'épaisseur)
        let p = alpha * inertia_e / (length * length);
        let loads: Vec<Load> = tip_nodes
            .iter()
            .map(|&node| Load {
                node,
                dof: Dof::Uy,
                value: p / tip_nodes.len() as f64,
            })
            .collect();
        let model = large(mesh.clone(), bcs.clone(), loads);
        let history = solve(&model, 10);
        let u = &history[10].field.displacements;
        let (w, shortening) = elastica::tip(alpha);
        let (got_w, got_u) = (u[2 * centre + 1] / length, -u[2 * centre] / length);
        assert!(
            (got_w - w).abs() < 1e-2 * w,
            "α = {alpha} : w/L = {got_w} au lieu de {w}"
        );
        assert!(
            (got_u - shortening).abs() < 2e-2 * shortening,
            "α = {alpha} : u/L = {got_u} au lieu de {shortening}"
        );
    }
}

/// Tube épais néo-Hookéen (a = 10, b = 20, déformations planes), alésage élargi de δ, bord extérieur
/// libre. Il n'y a pas de solution fermée en compressible : l'équilibre radial est intégré ici. Avec
/// F = diag(λr, λθ, 1), λθ = r / R, λr = dr/dR, la première contrainte de Piola donne
/// P_rR = λθ σr, P_θΘ = λr σθ, et d(R P_rR)/dR = P_θΘ. Tir sur R P_rR en R = a, de sorte que
/// P_rR(b) = 0 ; σr vient de l'inversion de σr(λr, λθ) par bissection.
mod inflated_tube {
    use super::*;
    use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};

    pub const A: f64 = 10.0;
    pub const B: f64 = 20.0;
    pub const YOUNG: f64 = 10.0;
    pub const POISSON: f64 = 0.3;

    fn lame() -> (f64, f64) {
        (
            YOUNG * POISSON / ((1.0 + POISSON) * (1.0 - 2.0 * POISSON)),
            YOUNG / (2.0 * (1.0 + POISSON)),
        )
    }

    /// Cauchy principale pour l'étirement `l` et J = λr λθ.
    fn principal(l: f64, j: f64) -> f64 {
        let (lambda, mu) = lame();
        (mu * (l * l - 1.0) + lambda * j.ln()) / j
    }

    /// λr tel que σr(λr, λθ) = `sigma_r`, σr croissante avec λr.
    fn radial_stretch(sigma_r: f64, theta: f64) -> f64 {
        let (mut low, mut high) = (1e-3, 50.0);
        for _ in 0..100 {
            let middle = 0.5 * (low + high);
            if principal(middle, middle * theta) < sigma_r {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    }

    /// Intègre de R = a à b depuis r(a) = a + δ et R P_rR(a) = `start` : (r(b), σr(a)).
    fn shoot(bore: f64, start: f64) -> (f64, f64, f64) {
        let rhs = |radius: f64, r: f64, rho: f64| -> (f64, f64) {
            let theta = r / radius;
            let sigma_r = rho / radius / theta;
            let lr = radial_stretch(sigma_r, theta);
            let sigma_t = principal(theta, lr * theta);
            (lr, lr * sigma_t)
        };
        let n = 1000;
        let h = (B - A) / n as f64;
        let (mut radius, mut r, mut rho) = (A, A + bore, start);
        let sigma_a = start / A / ((A + bore) / A);
        for _ in 0..n {
            let (k1r, k1p) = rhs(radius, r, rho);
            let (k2r, k2p) = rhs(radius + h / 2.0, r + h / 2.0 * k1r, rho + h / 2.0 * k1p);
            let (k3r, k3p) = rhs(radius + h / 2.0, r + h / 2.0 * k2r, rho + h / 2.0 * k2p);
            let (k4r, k4p) = rhs(radius + h, r + h * k3r, rho + h * k3p);
            r += h / 6.0 * (k1r + 2.0 * k2r + 2.0 * k3r + k4r);
            rho += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
            radius += h;
        }
        (r, rho, sigma_a)
    }

    /// (r(b), pression p = −σr(a)) pour l'élargissement `bore` de l'alésage.
    pub fn reference(bore: f64) -> (f64, f64) {
        let (mut low, mut high) = (-10.0 * YOUNG * A, 0.0);
        for _ in 0..80 {
            let middle = 0.5 * (low + high);
            // R P_rR(b) croît avec la pression : on cherche le zéro
            if shoot(bore, middle).1 < 0.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        let start = 0.5 * (low + high);
        let (outer, _, sigma_a) = shoot(bore, start);
        (outer, -sigma_a)
    }

    /// (déplacement radial en b, pression) du modèle, Tri6 de taille `h`.
    pub fn solve(h: f64, bore: f64) -> (f64, f64) {
        let geometry = Geometry {
            outer: vec![
                Curve::line([A, 0.0], h, "sym_x"),
                Curve::arc([B, 0.0], [0.0, 0.0], false, h, "exterieur"),
                Curve::line([0.0, B], h, "sym_y"),
                Curve::arc([0.0, A], [0.0, 0.0], true, h, "interieur"),
            ],
            holes: vec![],
        };
        let t = triangulate(&geometry, &MeshOptions::default())
            .unwrap()
            .to_quadratic();
        let [sym_x, sym_y, inner, outer] = ["sym_x", "sym_y", "interieur", "exterieur"]
            .map(|name| t.boundary(name).unwrap().nodes.clone());
        let material = Material::with_poisson(YOUNG, POISSON)
            .unwrap()
            .with_hyperelasticity(maille_core::material::Hyperelasticity::NeoHookean);
        let elements = t
            .triangles
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material, 1.0, PlaneState::Strain).unwrap()))
            .collect();
        let mesh = Mesh::new(t.nodes, elements).unwrap();
        let mut bcs: Vec<Bc> = sym_x
            .iter()
            .map(|&node| bc(node, Dof::Uy, 0.0))
            .chain(sym_y.iter().map(|&node| bc(node, Dof::Ux, 0.0)))
            .collect();
        for &node in &inner {
            let n = mesh.nodes[node];
            let r = n.x.hypot(n.y);
            if !sym_y.contains(&node) {
                bcs.push(bc(node, Dof::Ux, bore * n.x / r));
            }
            if !sym_x.contains(&node) {
                bcs.push(bc(node, Dof::Uy, bore * n.y / r));
            }
        }
        let history = solve_model(large(mesh.clone(), bcs, Vec::new()), 6);
        let field = &history[6].field;
        // effort de cerclage du quart de tube : −Σ R_y sur y = 0 vaut p a', a' = a + δ
        let hoop: f64 = sym_x.iter().map(|&n| field.reactions[2 * n + 1]).sum();
        let mean_outer = outer
            .iter()
            .map(|&n| field.displacements[2 * n].hypot(field.displacements[2 * n + 1]))
            .sum::<f64>()
            / outer.len() as f64;
        (mean_outer, -hoop / (A + bore))
    }

    fn solve_model(model: Model, increments: usize) -> Vec<Increment> {
        super::solve(&model, increments)
    }
}

#[test]
fn an_inflated_thick_neo_hookean_tube_follows_the_radial_equilibrium() {
    let bore = 0.3 * inflated_tube::A;
    let (outer, p) = inflated_tube::reference(bore);
    let (u_outer, p_fe) = inflated_tube::solve(1.5, bore);
    assert!(
        (u_outer / (outer - inflated_tube::B) - 1.0).abs() < 1e-2,
        "u(b) = {u_outer}, référence {}",
        outer - inflated_tube::B
    );
    assert!((p_fe / p - 1.0).abs() < 1e-2, "p = {p_fe}, référence {p}");
}
