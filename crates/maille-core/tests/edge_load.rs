//! Force linéique sur une arête : forces nodales équivalentes, traction (égale à la pression, sans
//! l'épaisseur), cisaillement pur, charge morte en grandes déformations, Amplitude et refus.

mod common;

use std::collections::BTreeSet;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::assembly::edge_load;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{Amplitude, Analysis, Bc, EdgeLoad, Kinematics, Model, Pressure, Step};
use maille_core::MailleError;

const E: f64 = 210e3;
const NU: f64 = 0.3;

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9 * (1.0 + a.abs().max(b.abs()))
}

/// Arêtes (élément, arête) dont tous les Nodes vérifient `on`, une seule fois par paire de sommets.
fn edges_where(mesh: &Mesh, on: impl Fn(&Node) -> bool) -> Vec<(usize, usize)> {
    let mut seen = BTreeSet::new();
    let mut found = Vec::new();
    for (i, e) in mesh.elements.iter().enumerate() {
        for k in 0..6 {
            let Some(nodes) = e.edge_nodes(k) else {
                continue;
            };
            if nodes.iter().all(|&n| on(&mesh.nodes[n]))
                && seen.insert((nodes[0].min(nodes[1]), nodes[0].max(nodes[1])))
            {
                found.push((i, k));
            }
        }
    }
    found
}

fn edge_loads(edges: &[(usize, usize)], value: [f64; 3]) -> Vec<EdgeLoad> {
    edges
        .iter()
        .map(|&(element, edge)| EdgeLoad {
            element,
            edge,
            value,
        })
        .collect()
}

fn fixed(mesh: &Mesh, dof: Dof, on: impl Fn(&Node) -> bool) -> Vec<Bc> {
    (0..mesh.nodes.len())
        .filter(|&n| on(&mesh.nodes[n]))
        .map(|node| Bc {
            node,
            dof,
            value: 0.0,
        })
        .collect()
}

fn tri3_mesh(nodes: Vec<Node>, tris: &[[usize; 3]], thickness: f64, state: PlaneState) -> Mesh {
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = tris
        .iter()
        .map(|&t| Tri3::new(t, material, thickness, state).map(Element::Tri3))
        .collect::<Result<_, _>>()
        .unwrap();
    Mesh::new(nodes, elements).unwrap()
}

// --- forces nodales équivalentes ---

#[test]
fn nodal_forces_of_a_linear_and_a_quadratic_triangle() {
    // arête oblique 3-4-5 : L = 5, q = (2, −1), l'épaisseur (7) n'y change rien
    let material = Material::with_poisson(E, NU).unwrap();
    let nodes = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(3.0, 4.0),
        Node::new_2d(0.0, 4.0),
    ];
    let linear = Mesh::new(
        nodes.clone(),
        vec![Element::Tri3(
            Tri3::new([0, 1, 2], material, 7.0, PlaneState::Stress).unwrap(),
        )],
    )
    .unwrap();
    let t = edge_load(&linear, 0, 0, [2.0, -1.0, 0.0]).unwrap();
    assert_eq!(t.positions, vec![0, 1]);
    for (f, want) in t.forces.iter().zip([5.0, -2.5, 5.0, -2.5]) {
        assert!(near(*f, want));
    }

    let quadratic = to_quadratic(&nodes, &[[0, 1, 2]]);
    let tri6 = Mesh::new(
        quadratic.nodes.clone(),
        vec![Element::Tri6(
            Tri6::new(quadratic.elements[0], material, 7.0, PlaneState::Stress).unwrap(),
        )],
    )
    .unwrap();
    let t = edge_load(&tri6, 0, 0, [2.0, -1.0, 0.0]).unwrap();
    assert_eq!(t.positions, vec![0, 1, 3]);
    // (L/6, L/6, 2L/3) · q
    let weights = [5.0 / 6.0, 5.0 / 6.0, 10.0 / 3.0];
    for (a, w) in weights.iter().enumerate() {
        assert!(near(t.forces[2 * a], 2.0 * w));
        assert!(near(t.forces[2 * a + 1], -w));
    }
}

#[test]
fn the_resultant_on_a_tetrahedron_edge_is_q_times_length() {
    let (nodes, tets) = box_grid([1, 1, 1], [1.0, 2.0, 2.0]);
    let mesh = Mesh::new(
        nodes,
        tets.iter()
            .map(|&t| Element::Tet4(Tet4::new(t, Material::with_poisson(E, NU).unwrap())))
            .collect(),
    )
    .unwrap();
    let q = [1.0, -2.0, 3.0];
    for k in 0..6 {
        let nodes = mesh.elements[0].edge_nodes(k).unwrap();
        let (a, b) = (mesh.nodes[nodes[0]], mesh.nodes[nodes[1]]);
        let length = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2) + (b.z - a.z).powi(2)).sqrt();
        let t = edge_load(&mesh, 0, k, q).unwrap();
        for (c, qc) in q.iter().enumerate() {
            let total: f64 = (0..2).map(|a| t.forces[3 * a + c]).sum();
            assert!(near(total, qc * length), "arête {k}");
            assert!(near(t.forces[c], qc * length / 2.0), "arête {k}");
        }
    }
    assert!(edge_load(&mesh, 0, 6, q).is_none());
}

// --- traction : même résultat que la pression ---

fn plate(thickness: f64) -> (Mesh, Vec<Bc>, Vec<(usize, usize)>) {
    let (nodes, tris) = rect_grid(4, 2, 4.0, 2.0);
    let mesh = tri3_mesh(nodes, &tris, thickness, PlaneState::Strain);
    let mut bcs = fixed(&mesh, Dof::Ux, |n| n.x < 1e-9);
    bcs.extend(fixed(&mesh, Dof::Uy, |n| n.y < 1e-9));
    let end = edges_where(&mesh, |n| (n.x - 4.0).abs() < 1e-9);
    (mesh, bcs, end)
}

#[test]
fn an_edge_load_pulls_like_the_pressure_without_the_thickness() {
    // σx = 100 sur une plaque d'épaisseur 2 : q = σ·t = 200 N/mm, εx = (1 − ν²) σ / E
    let (mesh, bcs, end) = plate(2.0);
    let field = Model::new(mesh.clone(), bcs.clone(), vec![])
        .unwrap()
        .with_edge_loads(edge_loads(&end, [200.0, 0.0, 0.0]))
        .unwrap()
        .solve()
        .unwrap();
    let strain = (1.0 - NU * NU) * 100.0 / E;
    for (i, n) in mesh.nodes.iter().enumerate() {
        assert!((field.displacements[2 * i] - strain * n.x).abs() < 1e-12);
    }
    let pressures = end
        .iter()
        .map(|&(element, face)| Pressure {
            element,
            face,
            value: -100.0,
        })
        .collect();
    let by_pressure = Model::new(mesh, bcs, vec![])
        .unwrap()
        .with_pressures(pressures)
        .unwrap()
        .solve()
        .unwrap();
    for (u, v) in field.displacements.iter().zip(&by_pressure.displacements) {
        assert!((u - v).abs() < 1e-12);
    }
}

// --- cisaillement pur, ce que la pression ne sait pas faire ---

#[test]
fn tangential_edge_loads_give_a_pure_shear() {
    // carré 2 × 2 en contraintes planes, τ = 10 sur les quatre côtés (t = 1) : γ = τ / G, σx = σy = 0
    let (nodes, tris) = rect_grid(2, 2, 2.0, 2.0);
    let mesh = tri3_mesh(nodes, &tris, 1.0, PlaneState::Stress);
    let tau = 10.0;
    let side = |on: fn(&Node) -> bool, q: [f64; 3]| edge_loads(&edges_where(&mesh, on), q);
    let mut loads = side(|n| n.y < 1e-9, [-tau, 0.0, 0.0]);
    loads.extend(side(|n| (n.y - 2.0).abs() < 1e-9, [tau, 0.0, 0.0]));
    loads.extend(side(|n| n.x < 1e-9, [0.0, -tau, 0.0]));
    loads.extend(side(|n| (n.x - 2.0).abs() < 1e-9, [0.0, tau, 0.0]));
    // appuis qui bloquent seulement les mouvements de corps rigide
    let bcs = vec![
        Bc {
            node: 0,
            dof: Dof::Ux,
            value: 0.0,
        },
        Bc {
            node: 0,
            dof: Dof::Uy,
            value: 0.0,
        },
        Bc {
            node: 2,
            dof: Dof::Uy,
            value: 0.0,
        },
    ];
    let field = Model::new(mesh, bcs, vec![])
        .unwrap()
        .with_edge_loads(loads)
        .unwrap()
        .solve()
        .unwrap();
    for s in field.stresses.iter().flatten() {
        assert!((s.tau_xy - tau).abs() < 1e-9, "τxy = {}", s.tau_xy);
        assert!(s.sigma_xx.abs() < 1e-9 && s.sigma_yy.abs() < 1e-9);
    }
}

// --- Amplitude ---

#[test]
fn the_edge_load_follows_the_amplitude() {
    let (mesh, bcs, end) = plate(1.0);
    let model = Model::new(mesh, bcs, vec![])
        .unwrap()
        .with_edge_loads(edge_loads(&end, [40.0, 0.0, 0.0]))
        .unwrap();
    let full = model.solve().unwrap();
    let history = model
        .solve_quasi_static(
            &Step::quasi_static(1.0, 4).unwrap(),
            &Amplitude::ramp(1.0).unwrap(),
        )
        .unwrap();
    for (k, inc) in history.iter().enumerate() {
        let lambda = k as f64 / 4.0;
        for (u, v) in inc.field.displacements.iter().zip(&full.displacements) {
            assert!((u - lambda * v).abs() < 1e-12, "incrément {k}");
        }
    }
}

// --- 3D, charge morte en grandes déformations ---

fn cube(quadratic: bool) -> (Mesh, Vec<Bc>, Vec<(usize, usize)>) {
    let (nodes, tets) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let material = Material::with_poisson(E, NU).unwrap();
    let mesh = if quadratic {
        let q = to_quadratic_tet(&nodes, &tets);
        let elements = q
            .elements
            .iter()
            .map(|&t| Element::Tet10(Tet10::new(t, material)))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tets
            .iter()
            .map(|&t| Element::Tet4(Tet4::new(t, material)))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    };
    let mut bcs = fixed(&mesh, Dof::Ux, |n| n.x < 1e-9);
    bcs.extend(fixed(&mesh, Dof::Uy, |n| n.y < 1e-9));
    bcs.extend(fixed(&mesh, Dof::Uz, |n| n.z < 1e-9));
    // les quatre côtés de la face x = 1
    let on_end = |n: &Node| (n.x - 1.0).abs() < 1e-9;
    let mut rim = edges_where(&mesh, |n| on_end(n) && n.y.abs() < 1e-9);
    rim.extend(edges_where(&mesh, |n| {
        on_end(n) && (n.y - 1.0).abs() < 1e-9
    }));
    rim.extend(edges_where(&mesh, |n| on_end(n) && n.z.abs() < 1e-9));
    rim.extend(edges_where(&mesh, |n| {
        on_end(n) && (n.z - 1.0).abs() < 1e-9
    }));
    (mesh, bcs, rim)
}

fn reaction_x(mesh: &Mesh, reactions: &[f64]) -> f64 {
    (0..mesh.nodes.len()).map(|n| reactions[3 * n]).sum()
}

#[test]
fn the_edge_load_resultant_is_q_times_the_rim_length_in_3d() {
    // périmètre 4 : la réaction vaut −q · 4, en petites déformations comme en grandes (charge morte),
    // en Tet4 comme en Tet10
    for (quadratic, kinematics) in [
        (false, Kinematics::Small),
        (false, Kinematics::Large),
        (true, Kinematics::Small),
        (true, Kinematics::Large),
    ] {
        let (mesh, bcs, rim) = cube(quadratic);
        let q = 25.0;
        let model = Model::new(mesh.clone(), bcs.clone(), vec![])
            .unwrap()
            .with_kinematics(kinematics)
            .unwrap()
            .with_edge_loads(edge_loads(&rim, [q, 0.0, 0.0]))
            .unwrap();
        let field = match kinematics {
            Kinematics::Small => model.solve().unwrap(),
            Kinematics::Large => {
                model
                    .solve_quasi_static(
                        &Step::quasi_static(1.0, 2).unwrap(),
                        &Amplitude::ramp(1.0).unwrap(),
                    )
                    .unwrap()
                    .pop()
                    .unwrap()
                    .field
            }
        };
        assert!(
            near(reaction_x(&mesh, &field.reactions), -4.0 * q),
            "quadratique : {quadratic}, {kinematics:?}"
        );
    }
}

#[test]
fn a_tet10_edge_gets_one_sixth_two_thirds_one_sixth() {
    let (nodes, tets) = box_grid([1, 1, 1], [1.0, 2.0, 2.0]);
    let q = to_quadratic_tet(&nodes, &tets);
    let material = Material::with_poisson(E, NU).unwrap();
    let mesh = Mesh::new(
        q.nodes,
        q.elements
            .iter()
            .map(|&t| Element::Tet10(Tet10::new(t, material)))
            .collect(),
    )
    .unwrap();
    let load = [1.0, -2.0, 3.0];
    for k in 0..6 {
        let nodes = mesh.elements[0].edge_nodes(k).unwrap();
        assert_eq!(nodes.len(), 3);
        let (a, b) = (mesh.nodes[nodes[0]], mesh.nodes[nodes[1]]);
        let mid = mesh.nodes[nodes[2]];
        // le Node milieu est au milieu de l'arête 0-1 de la numérotation
        assert!(near(mid.x, (a.x + b.x) / 2.0) && near(mid.z, (a.z + b.z) / 2.0));
        let length = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2) + (b.z - a.z).powi(2)).sqrt();
        let t = edge_load(&mesh, 0, k, load).unwrap();
        for (c, qc) in load.iter().enumerate() {
            for (a, w) in [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0].iter().enumerate() {
                assert!(near(t.forces[3 * a + c], qc * length * w), "arête {k}");
            }
        }
    }
}

// --- refus ---

#[test]
fn invalid_edge_loads_are_typed_errors() {
    let (mesh, bcs, end) = plate(1.0);
    let model = Model::new(mesh.clone(), bcs, vec![]).unwrap();
    let (element, edge) = end[0];
    let one = |element, edge, value: [f64; 3]| {
        vec![EdgeLoad {
            element,
            edge,
            value,
        }]
    };
    assert_eq!(
        model
            .clone()
            .with_edge_loads(one(element, 3, [1.0, 0.0, 0.0]))
            .err(),
        Some(MailleError::UnknownEdge { element, edge: 3 })
    );
    assert_eq!(
        model
            .clone()
            .with_edge_loads(one(9999, 0, [1.0, 0.0, 0.0]))
            .err(),
        Some(MailleError::UnknownElement { element: 9999 })
    );
    assert!(matches!(
        model
            .clone()
            .with_edge_loads(one(element, edge, [f64::NAN, 0.0, 0.0]))
            .err(),
        Some(MailleError::NonFinite { .. })
    ));
    assert_eq!(
        model
            .clone()
            .with_edge_loads(one(element, edge, [0.0, 0.0, 1.0]))
            .err(),
        Some(MailleError::EdgeLoadOutOfPlane { element })
    );

    let conductive = Material::conductive(1.0).unwrap();
    let thermal_elements = mesh
        .elements
        .iter()
        .map(|e| match e {
            Element::Tri3(t) => Element::Tri3(Tri3 {
                material: conductive,
                ..*t
            }),
            other => *other,
        })
        .collect();
    let thermal = Model::with_analysis(
        Mesh::new(mesh.nodes.clone(), thermal_elements).unwrap(),
        Analysis::Thermal,
        vec![],
        vec![],
    )
    .unwrap();
    assert_eq!(
        thermal
            .with_edge_loads(one(element, edge, [1.0, 0.0, 0.0]))
            .err(),
        Some(MailleError::EdgeLoadStaticOnly)
    );
}
