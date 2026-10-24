//! Validation de l'élément Tri6 : traction uniaxiale, patch tests linéaire et quadratique (flexion pure).

mod common;

use common::{rect_grid, to_quadratic, Quadratic};
use maille_core::element::plane_d;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::model::{Bc, DistributedLoad, Load, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

const E: f64 = 210e3;
const NU: f64 = 0.3;
const THICKNESS: f64 = 2.0;
const RTOL: f64 = 1e-10;

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= RTOL * scale
}

fn tri6(nodes: [usize; 6], state: PlaneState) -> Element {
    Element::Tri6(
        Tri6::new(
            nodes,
            Material::with_poisson(E, NU).unwrap(),
            THICKNESS,
            state,
        )
        .unwrap(),
    )
}

fn mesh_of(q: &Quadratic, state: PlaneState) -> Mesh {
    let elements = q.elements.iter().map(|&n| tri6(n, state)).collect();
    Mesh::new(q.nodes.clone(), elements).unwrap()
}

fn ux(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node]
}
fn uy(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node + 1]
}

#[test]
fn uniaxial_tension_is_exact_at_every_node_and_gauss_point() {
    let (nx, ny, w, h, sigma) = (3, 2, 6.0, 4.0, 100.0);
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let (nodes, tris) = rect_grid(nx, ny, w, h);
    let q = to_quadratic(&nodes, &tris);

    for state in [PlaneState::Stress, PlaneState::Strain] {
        let mesh = mesh_of(&q, state);
        // symétries : Ux = 0 sur le bord gauche, Uy = 0 sur le bord bas, milieux d'arêtes compris
        let on = |n: &Node, f: fn(&Node) -> f64| f(n).abs() < 1e-12;
        let mut bcs = Vec::new();
        for (n, node) in q.nodes.iter().enumerate() {
            if on(node, |p| p.x) {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Ux,
                    value: 0.0,
                });
            }
            if on(node, |p| p.y) {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Uy,
                    value: 0.0,
                });
            }
        }
        // forces consistantes d'une traction uniforme sur une arête quadratique : L/6, 4L/6, L/6
        let dy = h / ny as f64;
        let mut loads = Vec::new();
        for j in 0..ny {
            let (a, b) = (idx(nx, j), idx(nx, j + 1));
            let total = sigma * THICKNESS * dy;
            for (node, share) in [(a, 1.0 / 6.0), (q.mid(a, b), 4.0 / 6.0), (b, 1.0 / 6.0)] {
                loads.push(Load {
                    node,
                    dof: Dof::Ux,
                    value: total * share,
                });
            }
        }
        let f = Model::new(mesh, bcs, loads).unwrap().solve().unwrap();

        let (e_eq, nu_eq) = match state {
            PlaneState::Stress => (E, NU),
            PlaneState::Strain => (E / (1.0 - NU * NU), NU / (1.0 - NU)),
        };
        for (n, node) in q.nodes.iter().enumerate() {
            assert!(
                close(ux(&f, n), sigma * node.x / e_eq, sigma * w / e_eq),
                "{state:?} ux {n}"
            );
            assert!(
                close(uy(&f, n), -nu_eq * sigma * node.y / e_eq, sigma * h / e_eq),
                "{state:?} uy {n}"
            );
        }
        for (e, points) in f.stresses.iter().enumerate() {
            assert_eq!(points.len(), 3, "3 points de Gauss par Tri6");
            for p in points {
                assert!(close(p.sigma_xx, sigma, sigma), "{state:?} σxx élément {e}");
                assert!(close(p.sigma_yy, 0.0, sigma), "{state:?} σyy élément {e}");
                assert!(close(p.tau_xy, 0.0, sigma), "{state:?} τxy élément {e}");
            }
        }
    }
}

/// Maillage irrégulier de [0,2]² (2 Nodes intérieurs), mis en quadratique.
fn patch_quadratic() -> Quadratic {
    let nodes: Vec<Node> = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.7, 0.6),
        (1.3, 1.2),
    ]
    .map(|(x, y)| Node { x, y })
    .to_vec();
    let tris = [
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 5],
        [3, 0, 4],
        [4, 1, 5],
        [4, 5, 3],
    ];
    to_quadratic(&nodes, &tris)
}

/// Patch test : on impose `field` sur tous les Nodes du bord (sommets et milieux) ; les Nodes
/// intérieurs et les contraintes de Gauss doivent retrouver `field` et `stress` exactement.
fn run_patch(field: impl Fn(f64, f64) -> (f64, f64), stress: impl Fn(f64, f64) -> [f64; 3]) {
    let q = patch_quadratic();
    let mesh = mesh_of(&q, PlaneState::Stress);
    let boundary = |n: &Node| {
        n.x.abs() < 1e-12
            || (n.x - 2.0).abs() < 1e-12
            || n.y.abs() < 1e-12
            || (n.y - 2.0).abs() < 1e-12
    };
    let mut bcs = Vec::new();
    for (n, node) in q.nodes.iter().enumerate().filter(|(_, p)| boundary(p)) {
        let (u, v) = field(node.x, node.y);
        bcs.push(Bc {
            node: n,
            dof: Dof::Ux,
            value: u,
        });
        bcs.push(Bc {
            node: n,
            dof: Dof::Uy,
            value: v,
        });
    }
    let f = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .solve()
        .unwrap();

    let mut interior = 0;
    for (n, node) in q.nodes.iter().enumerate().filter(|(_, p)| !boundary(p)) {
        let (u, v) = field(node.x, node.y);
        assert!(close(ux(&f, n), u, 1e-3), "ux nœud intérieur {n}");
        assert!(close(uy(&f, n), v, 1e-3), "uy nœud intérieur {n}");
        interior += 1;
    }
    assert!(
        interior >= 6,
        "le test doit avoir de vrais Nodes intérieurs ({interior})"
    );

    for (e, points) in f.stresses.iter().enumerate() {
        for p in points {
            let s = stress(p.x, p.y);
            let scale = 1.0 + s.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            assert!(close(p.sigma_xx, s[0], scale), "σxx élément {e}");
            assert!(close(p.sigma_yy, s[1], scale), "σyy élément {e}");
            assert!(close(p.tau_xy, s[2], scale), "τxy élément {e}");
        }
    }

    // lissage aux Nodes : la contrainte est au plus linéaire, l'extrapolation des 3 points de Gauss
    // la retrouve exactement à chaque Node, bord compris
    let nodal = nodal_stresses(&mesh, &f).unwrap();
    assert_eq!(nodal.len(), q.nodes.len());
    for (n, p) in nodal.iter().enumerate() {
        let s = stress(q.nodes[n].x, q.nodes[n].y);
        let scale = 1.0 + s.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(close(p.x, q.nodes[n].x, 1.0) && close(p.y, q.nodes[n].y, 1.0));
        assert!(close(p.sigma_xx, s[0], scale), "σxx Node {n}");
        assert!(close(p.sigma_yy, s[1], scale), "σyy Node {n}");
        assert!(close(p.tau_xy, s[2], scale), "τxy Node {n}");
    }
}

#[test]
fn patch_test_linear_field() {
    let (a, b) = ([1e-3, 2e-3, -1.5e-3], [-5e-4, 3e-3, 2.5e-3]);
    let d = plane_d(E, NU, PlaneState::Stress);
    let eps = [a[1], b[2], a[2] + b[1]];
    let sigma: [f64; 3] = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * eps[c]).sum());
    run_patch(
        |x, y| (a[0] + a[1] * x + a[2] * y, b[0] + b[1] * x + b[2] * y),
        |_, _| sigma,
    );
}

#[test]
fn patch_test_quadratic_field_pure_bending() {
    // Flexion pure en contraintes planes : σxx = E κ y, les autres contraintes sont nulles
    // (équilibre sans force de volume). Déplacements : u = κ x y, v = −κ/2 (x² + ν y²),
    // quadratiques : le Tri6 les représente exactement.
    let kappa = 1e-3;
    run_patch(
        |x, y| (kappa * x * y, -kappa / 2.0 * (x * x + NU * y * y)),
        |_, y| [E * kappa * y, 0.0, 0.0],
    );
}

#[test]
fn gauss_points_of_a_tri6_are_centered_on_the_centroid() {
    let (nodes, tris) = rect_grid(1, 1, 2.0, 2.0);
    let q = to_quadratic(&nodes, &tris);
    let mesh = mesh_of(&q, PlaneState::Stress);
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
            node: 1,
            dof: Dof::Uy,
            value: 0.0,
        },
    ];
    let f = Model::new(
        mesh,
        bcs,
        vec![Load {
            node: 2,
            dof: Dof::Ux,
            value: 1.0,
        }],
    )
    .unwrap()
    .solve()
    .unwrap();
    // premier triangle : (0,0), (2,0), (2,2), centre de gravité (4/3, 2/3)
    let pts = &f.stresses[0];
    let (mx, my) = (
        pts.iter().map(|p| p.x).sum::<f64>() / 3.0,
        pts.iter().map(|p| p.y).sum::<f64>() / 3.0,
    );
    assert!(close(mx, 4.0 / 3.0, 1.0) && close(my, 2.0 / 3.0, 1.0));
}

#[test]
fn invalid_tri6_inputs_are_typed_errors() {
    let m = Material::with_poisson(E, NU).unwrap();
    assert_eq!(
        Tri6::new(
            [0, 1, 2, 3, 4, 5],
            Material::new(E).unwrap(),
            1.0,
            PlaneState::Stress
        ),
        Err(MailleError::MissingPoisson)
    );
    // Node milieu 3 (arête 0-1) déplacé : côté courbe refusé
    let mut nodes = vec![
        Node { x: 0.0, y: 0.0 },
        Node { x: 2.0, y: 0.0 },
        Node { x: 0.0, y: 2.0 },
        Node { x: 1.0, y: 0.0 },
        Node { x: 1.0, y: 1.0 },
        Node { x: 0.0, y: 1.0 },
    ];
    let element = Element::Tri6(Tri6::new([0, 1, 2, 3, 4, 5], m, 1.0, PlaneState::Stress).unwrap());
    assert!(Mesh::new(nodes.clone(), vec![element]).is_ok());
    nodes[3].y = 0.1;
    assert_eq!(
        Mesh::new(nodes.clone(), vec![element]),
        Err(MailleError::CurvedTriangle { element: 0 })
    );
    // sommets horaires
    nodes[3].y = 0.0;
    let reversed =
        Element::Tri6(Tri6::new([0, 2, 1, 5, 4, 3], m, 1.0, PlaneState::Stress).unwrap());
    assert_eq!(
        Mesh::new(nodes.clone(), vec![reversed]),
        Err(MailleError::ClockwiseTriangle { element: 0 })
    );
    // charge répartie refusée
    let mesh = Mesh::new(nodes, vec![element]).unwrap();
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
            node: 1,
            dof: Dof::Uy,
            value: 0.0,
        },
    ];
    assert_eq!(
        Model::new(mesh, bcs, vec![])
            .unwrap()
            .with_distributed_loads(vec![DistributedLoad {
                element: 0,
                value: 1.0
            }])
            .unwrap()
            .solve(),
        Err(MailleError::UnsupportedDistributedLoad { element: 0 })
    );
}
