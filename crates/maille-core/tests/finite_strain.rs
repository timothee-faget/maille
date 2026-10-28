//! Grandes déformations (Lagrangien total, Saint-Venant–Kirchhoff) : objectivité, traction homogène
//! de forme fermée, limite linéaire, convergence quadratique de Newton-Raphson.

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
    // températures et pression, dans les deux ordres
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
    let pressure = || {
        vec![Pressure {
            element: 0,
            face: 0,
            value: 1.0,
        }]
    };
    assert_eq!(
        large_model.clone().with_pressures(pressure()),
        Err(MailleError::LargeStrainWithPressure)
    );
    let pressed = model(tri_mesh(expandable, strain))
        .with_pressures(pressure())
        .unwrap();
    assert_eq!(
        pressed.with_kinematics(Kinematics::Large),
        Err(MailleError::LargeStrainWithPressure)
    );
    // solve linéaire
    assert_eq!(
        large_model.solve(),
        Err(MailleError::LargeStrainNeedsQuasiStatic)
    );
}
