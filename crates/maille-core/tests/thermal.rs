//! Thermique stationnaire sur triangles et tétraèdres, comparée aux solutions analytiques.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{Analysis, Bc, Load, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

/// Conductivité (unités cohérentes au choix de l'utilisateur, ex. W/(mm·K)).
const K: f64 = 50.0;
const T: f64 = 2.0;
const L: f64 = 4.0;
const H: f64 = 1.0;

fn material() -> Material {
    Material::conductive(K).unwrap()
}

/// Mesh de triangles L × H (épaisseur T), Tri3 ou Tri6.
fn plane_mesh(nx: usize, ny: usize, quadratic: bool) -> Mesh {
    let (nodes, tris) = rect_grid(nx, ny, L, H);
    if quadratic {
        let q = to_quadratic(&nodes, &tris);
        let elements = q
            .elements
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tris
            .iter()
            .map(|&n| Element::Tri3(Tri3::new(n, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    }
}

/// Pavé de tétraèdres L × H × T, Tet4 ou Tet10.
fn solid_mesh(n: [usize; 3], quadratic: bool) -> Mesh {
    let (nodes, tets) = box_grid(n, [L, H, T]);
    if quadratic {
        let q = to_quadratic_tet(&nodes, &tets);
        let elements = q
            .elements
            .iter()
            .map(|&t| Element::Tet10(Tet10::new(t, material())))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tets
            .iter()
            .map(|&t| Element::Tet4(Tet4::new(t, material())))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    }
}

/// Les 4 Mesh de la validation : Tri3, Tri6, Tet4, Tet10. Section A = H × T dans tous les cas.
fn all_meshes() -> Vec<(&'static str, Mesh)> {
    vec![
        ("Tri3", plane_mesh(4, 2, false)),
        ("Tri6", plane_mesh(4, 2, true)),
        ("Tet4", solid_mesh([4, 1, 2], false)),
        ("Tet10", solid_mesh([4, 1, 2], true)),
    ]
}

/// BC T = `value` sur tous les Nodes du plan x = `x`.
fn fix_plane(mesh: &Mesh, x: f64, value: f64) -> Vec<Bc> {
    nodes_at(mesh, x)
        .map(|node| Bc {
            node,
            dof: Dof::T,
            value,
        })
        .collect()
}

fn nodes_at(mesh: &Mesh, x: f64) -> impl Iterator<Item = usize> + '_ {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(move |(_, p)| (p.x - x).abs() < 1e-12)
        .map(|(n, _)| n)
}

fn solve(mesh: &Mesh, bcs: Vec<Bc>, loads: Vec<Load>) -> Field {
    Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, loads)
        .unwrap()
        .solve()
        .unwrap()
}

#[test]
fn wall_with_imposed_temperatures_is_exact_with_the_right_heat_flow() {
    // T(x) = T1 + (T2 − T1) x / L, q = −k (T2 − T1) / L, puissance traversante q A
    let (t1, t2) = (100.0, 20.0);
    let q = -K * (t2 - t1) / L;
    let power = q * H * T;
    for (name, mesh) in all_meshes() {
        let mut bcs = fix_plane(&mesh, 0.0, t1);
        bcs.extend(fix_plane(&mesh, L, t2));
        let field = solve(&mesh, bcs, Vec::new());

        assert_eq!(field.temperatures.len(), mesh.nodes.len());
        assert!(field.displacements.is_empty() && field.stresses.is_empty());
        for (n, p) in mesh.nodes.iter().enumerate() {
            let exact = t1 + (t2 - t1) * p.x / L;
            assert!(
                (field.temperatures[n] - exact).abs() < 1e-10,
                "{name} Node {n}"
            );
        }
        for f in field.heat_fluxes.iter().flatten() {
            assert!((f.qx - q).abs() < 1e-9 * q, "{name} qx = {}", f.qx);
            assert!(f.qy.abs() < 1e-9 * q && f.qz.abs() < 1e-9 * q, "{name}");
        }
        // la chaleur entre par x = 0 et sort par x = L
        let sum = |x: f64| nodes_at(&mesh, x).map(|n| field.reactions[n]).sum::<f64>();
        assert!(
            (sum(0.0) - power).abs() < 1e-9 * power,
            "{name} entrée {}",
            sum(0.0)
        );
        assert!(
            (sum(L) + power).abs() < 1e-9 * power,
            "{name} sortie {}",
            sum(L)
        );
    }
}

#[test]
fn heat_injected_at_the_end_gives_t1_plus_p_x_over_k_a() {
    // T = T1 en x = 0, puissance P répartie uniformément sur l'arête x = L (parts consistantes)
    let (t1, p) = (20.0, 400.0);
    let (nx, ny) = (4, 2);
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let dy = H / ny as f64;
    for quadratic in [false, true] {
        let mesh = plane_mesh(nx, ny, quadratic);
        let mut shares: Vec<(usize, f64)> = Vec::new();
        if quadratic {
            let (nodes, tris) = rect_grid(nx, ny, L, H);
            let q = to_quadratic(&nodes, &tris);
            for j in 0..ny {
                let (a, b) = (idx(nx, j), idx(nx, j + 1));
                let segment = p * dy / H;
                shares.extend([
                    (a, segment / 6.0),
                    (q.mid(a, b), segment * 4.0 / 6.0),
                    (b, segment / 6.0),
                ]);
            }
        } else {
            for j in 0..=ny {
                let tributary = if j == 0 || j == ny { dy / 2.0 } else { dy };
                shares.push((idx(nx, j), p * tributary / H));
            }
        }
        let loads = shares
            .into_iter()
            .map(|(node, value)| Load {
                node,
                dof: Dof::T,
                value,
            })
            .collect();
        let field = solve(&mesh, fix_plane(&mesh, 0.0, t1), loads);
        for (n, node) in mesh.nodes.iter().enumerate() {
            let exact = t1 + p * node.x / (K * H * T);
            assert!(
                (field.temperatures[n] - exact).abs() < 1e-10,
                "quadratique {quadratic}, Node {n}"
            );
        }
        let entering: f64 = nodes_at(&mesh, 0.0).map(|n| field.reactions[n]).sum();
        assert!(
            (entering + p).abs() < 1e-9 * p,
            "la chaleur injectée sort par x = 0"
        );
    }
}

#[test]
fn thermal_inputs_are_typed_errors() {
    let mesh = plane_mesh(1, 1, false);
    let thermal =
        |mesh: Mesh, bcs: Vec<Bc>| Model::with_analysis(mesh, Analysis::Thermal, bcs, Vec::new());

    // DOF mécanique en thermique, DOF T en statique
    let ux = vec![Bc {
        node: 0,
        dof: Dof::Ux,
        value: 0.0,
    }];
    assert_eq!(
        thermal(mesh.clone(), ux),
        Err(MailleError::InvalidDof {
            node: 0,
            dof: Dof::Ux
        })
    );
    let t = vec![Bc {
        node: 9,
        dof: Dof::T,
        value: 0.0,
    }];
    assert_eq!(
        thermal(mesh.clone(), t),
        Err(MailleError::UnknownNode { node: 9 })
    );

    // Material sans k ; Material sans E en statique
    let steel = Material::with_poisson(210e3, 0.3).unwrap();
    let no_k = Tri3::new([0, 1, 2], steel, 1.0, PlaneState::Stress).unwrap();
    let nodes = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(1.0, 0.0),
        Node::new_2d(0.0, 1.0),
    ];
    let mesh_no_k = Mesh::new(nodes, vec![Element::Tri3(no_k)]).unwrap();
    assert_eq!(
        thermal(mesh_no_k, Vec::new()),
        Err(MailleError::MissingConductivity { element: 0 })
    );
    assert_eq!(
        Model::new(mesh, Vec::new(), Vec::new()),
        Err(MailleError::MissingYoung { element: 0 })
    );
    assert!(Material::conductive(0.0).is_err() && steel.with_conductivity(-1.0).is_err());

    // Bar : pas de thermique
    let bar = Bar {
        nodes: [0, 1],
        material: steel.with_conductivity(K).unwrap(),
        section: Section::new(1.0).unwrap(),
    };
    let bars = Mesh::new(
        vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
        vec![Element::Bar(bar)],
    )
    .unwrap();
    assert_eq!(
        thermal(bars, Vec::new()),
        Err(MailleError::UnsupportedThermalElement { element: 0 })
    );
}

#[test]
fn a_thermal_field_has_no_stresses() {
    let mesh = plane_mesh(1, 1, false);
    let field = solve(&mesh, fix_plane(&mesh, 0.0, 1.0), Vec::new());
    assert_eq!(nodal_stresses(&mesh, &field), Err(MailleError::NoStresses));
}

#[test]
fn one_material_serves_both_analyses() {
    // E, ν et k dans le même Material : le Mesh se résout en statique comme en thermique
    let m = Material::with_poisson(210e3, 0.3)
        .unwrap()
        .with_conductivity(K)
        .unwrap();
    let (nodes, tris) = rect_grid(2, 1, L, H);
    let elements = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, m, T, PlaneState::Stress).unwrap()))
        .collect();
    let mesh = Mesh::new(nodes, elements).unwrap();
    assert!(Model::new(mesh.clone(), Vec::new(), Vec::new()).is_ok());
    let field = solve(&mesh, fix_plane(&mesh, 0.0, 5.0), Vec::new());
    for t in &field.temperatures {
        assert!((t - 5.0).abs() < 1e-12, "sans flux, T uniforme");
    }
}
