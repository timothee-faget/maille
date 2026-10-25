//! Patch tests sur un maillage généré : plaque irrégulière percée d'un trou circulaire excentré,
//! triangles de Ruppert de tailles variées. Le champ imposé sur tous les Nodes du bord (pris dans
//! les bords nommés du mailleur) doit être retrouvé exactement à l'intérieur.

use maille_core::element::plane_d;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions, Triangulation};
use maille_core::model::{Bc, Model};
use maille_core::results::{nodal_stresses, Field};

const E: f64 = 210e3;
const NU: f64 = 0.3;
const RTOL: f64 = 1e-9;

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= RTOL * scale
}

/// Quadrilatère irrégulier percé d'un trou de rayon 0,3 excentré, tailles différentes par bord.
fn plate() -> Triangulation<3> {
    let outer = [[0.0, 0.0], [2.2, 0.1], [2.0, 1.9], [-0.1, 2.1]];
    let sizes = [0.4, 0.25, 0.3, 0.35];
    let g = Geometry {
        outer: (0..4)
            .map(|i| Curve::line(outer[i], sizes[i], "bord"))
            .collect(),
        holes: vec![vec![Curve::arc(
            [1.4, 0.9],
            [1.1, 0.9],
            false,
            0.08,
            "trou",
        )]],
    };
    triangulate(&g, &MeshOptions::default()).unwrap()
}

/// Nodes posés sur un bord nommé, triés.
fn boundary_nodes<const N: usize>(t: &Triangulation<N>) -> Vec<usize> {
    let mut nodes: Vec<usize> = ["bord", "trou"]
        .iter()
        .flat_map(|name| t.boundary(name).unwrap().nodes.clone())
        .collect();
    nodes.sort_unstable();
    nodes
}

fn ux(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node]
}
fn uy(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node + 1]
}

/// Impose `field` sur les Nodes du bord, résout, et vérifie déplacements, contraintes aux points de
/// Gauss et contraintes lissées aux Nodes contre `field` et `stress`.
fn run_patch(
    mesh: &Mesh,
    on_boundary: &[usize],
    field: impl Fn(f64, f64) -> (f64, f64),
    stress: impl Fn(f64, f64) -> [f64; 3],
) {
    let mut bcs = Vec::new();
    for &n in on_boundary {
        let (u, v) = field(mesh.nodes[n].x, mesh.nodes[n].y);
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

    let interior: Vec<usize> = (0..mesh.nodes.len())
        .filter(|n| on_boundary.binary_search(n).is_err())
        .collect();
    assert!(
        interior.len() >= 30,
        "le test doit avoir de vrais Nodes intérieurs ({})",
        interior.len()
    );
    let scale = 1e-3 * 3.0; // ordre de grandeur des déplacements imposés
    for &n in &interior {
        let (u, v) = field(mesh.nodes[n].x, mesh.nodes[n].y);
        assert!(close(ux(&f, n), u, scale), "ux Node intérieur {n}");
        assert!(close(uy(&f, n), v, scale), "uy Node intérieur {n}");
    }

    let check = |what: &str, p: (f64, f64), got: [f64; 3]| {
        let s = stress(p.0, p.1);
        let scale = 1.0 + s.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        for k in 0..3 {
            assert!(close(got[k], s[k], scale), "{what} : {got:?} contre {s:?}");
        }
    };
    for (e, points) in f.stresses.iter().enumerate() {
        for p in points {
            check(
                &format!("Gauss, Element {e}"),
                (p.x, p.y),
                [p.sigma_xx, p.sigma_yy, p.tau_xy],
            );
        }
    }
    let nodal = nodal_stresses(mesh, &f).unwrap();
    for (n, p) in nodal.iter().enumerate() {
        let node: &Node = &mesh.nodes[n];
        check(
            &format!("Node {n}"),
            (node.x, node.y),
            [p.sigma_xx, p.sigma_yy, p.tau_xy],
        );
    }
}

/// Champ linéaire quelconque et sa contrainte constante, en contraintes planes.
fn linear_field() -> (impl Fn(f64, f64) -> (f64, f64), [f64; 3]) {
    let (a, b) = ([1e-3, 2e-3, -1.5e-3], [-5e-4, 3e-3, 2.5e-3]);
    let d = plane_d(E, NU, PlaneState::Stress);
    let eps = [a[1], b[2], a[2] + b[1]];
    let sigma = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * eps[c]).sum());
    (
        move |x: f64, y: f64| (a[0] + a[1] * x + a[2] * y, b[0] + b[1] * x + b[2] * y),
        sigma,
    )
}

fn tri3_mesh(t: &Triangulation<3>) -> Mesh {
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = t
        .triangles
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, material, 1.0, PlaneState::Stress).unwrap()))
        .collect();
    Mesh::new(t.nodes.clone(), elements).unwrap()
}

fn tri6_mesh(t: &Triangulation<6>) -> Mesh {
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = t
        .triangles
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material, 1.0, PlaneState::Stress).unwrap()))
        .collect();
    Mesh::new(t.nodes.clone(), elements).unwrap()
}

#[test]
fn tri3_patch_test_on_generated_mesh() {
    let t = plate();
    let (field, sigma) = linear_field();
    run_patch(&tri3_mesh(&t), &boundary_nodes(&t), field, |_, _| sigma);
}

#[test]
fn tri6_patch_test_on_generated_mesh() {
    let t = plate().to_quadratic();
    let (field, sigma) = linear_field();
    run_patch(&tri6_mesh(&t), &boundary_nodes(&t), field, |_, _| sigma);
}

#[test]
fn tri6_pure_bending_on_generated_mesh() {
    // Flexion pure en contraintes planes : σxx = E κ y, u = κ x y, v = −κ/2 (x² + ν y²),
    // quadratiques : le Tri6 les représente exactement.
    let t = plate().to_quadratic();
    let kappa = 1e-3;
    run_patch(
        &tri6_mesh(&t),
        &boundary_nodes(&t),
        |x, y| (kappa * x * y, -kappa / 2.0 * (x * x + NU * y * y)),
        |_, y| [E * kappa * y, 0.0, 0.0],
    );
}

#[test]
fn quadratic_triangulation_keeps_boundaries_consistent() {
    let linear = plate();
    let (n3, faces3) = (
        linear.nodes.len(),
        linear.boundary("trou").unwrap().faces.clone(),
    );
    let t = linear.to_quadratic();
    let mesh = tri6_mesh(&t);
    // chaque arête porte un seul Node milieu : V + A Nodes, A = (3 T + arêtes de bord) / 2
    let boundary_edges: usize = ["bord", "trou"]
        .iter()
        .map(|n| t.boundary(n).unwrap().faces.len())
        .sum();
    assert_eq!(
        t.nodes.len(),
        n3 + (3 * t.triangles.len() + boundary_edges) / 2
    );
    // mêmes faces, et les 3 Nodes de chaque face (sommets puis milieu) sont sur le bord
    let trou = t.boundary("trou").unwrap();
    assert_eq!(trou.faces, faces3);
    for &(e, f) in &trou.faces {
        let nodes = mesh.elements[e].face_nodes(f).unwrap();
        assert_eq!(nodes.len(), 3);
        assert!(nodes.iter().all(|n| trou.nodes.contains(n)));
    }
    // sommets sur le cercle, milieux sur la corde (un peu à l'intérieur du cercle)
    for &n in &trou.nodes {
        let r = (t.nodes[n].x - 1.1).hypot(t.nodes[n].y - 0.9);
        assert!(r <= 0.3 + 1e-12 && r > 0.3 * (1.0 - 0.08f64.powi(2) / 0.09));
    }
}
