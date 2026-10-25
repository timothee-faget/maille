//! Mailleur 2D : bords respectés, trous vidés, bords nommés.

use std::f64::consts::PI;

use maille_core::material::Material;
use maille_core::mesh::{Element, Mesh, PlaneState, Tri3};
use maille_core::mesher::{triangulate, Curve, Geometry, Triangulation};
use maille_core::MailleError;

fn polygon(points: &[[f64; 2]], sizes: &[f64], names: &[&str]) -> Vec<Curve> {
    (0..points.len())
        .map(|i| Curve::line(points[i], sizes[i], names[i]))
        .collect()
}

fn square(size: f64) -> Vec<Curve> {
    polygon(
        &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        &[size; 4],
        &["bas", "droite", "haut", "gauche"],
    )
}

fn corners(t: &Triangulation<3>, e: usize) -> [[f64; 2]; 3] {
    t.triangles[e].map(|n| [t.nodes[n].x, t.nodes[n].y])
}

fn signed_area(ring: &[[f64; 2]]) -> f64 {
    (0..ring.len())
        .map(|i| {
            let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum::<f64>()
        / 2.0
}

/// Vérifie que chaque triangle est antihoraire et d'aire non nulle, et renvoie l'aire maillée.
fn meshed_area(t: &Triangulation<3>) -> f64 {
    (0..t.triangles.len())
        .map(|e| {
            let a = signed_area(&corners(t, e));
            assert!(a > 1e-12, "triangle {e} retourné ou plat : aire {a}");
            a
        })
        .sum()
}

fn centroid(t: &Triangulation<3>, e: usize) -> [f64; 2] {
    let c = corners(t, e);
    [0, 1].map(|j| (c[0][j] + c[1][j] + c[2][j]) / 3.0)
}

/// Les deux Nodes de chaque face d'un bord sont dans ses Nodes, et chaque Node du bord vérifie `on`.
fn check_boundary(t: &Triangulation<3>, name: &str, on: impl Fn([f64; 2]) -> bool) {
    let set = t.boundary(name).unwrap();
    for &(e, f) in &set.faces {
        let (a, b) = (t.triangles[e][f], t.triangles[e][(f + 1) % 3]);
        assert!(set.nodes.contains(&a) && set.nodes.contains(&b));
    }
    for &n in &set.nodes {
        assert!(
            on([t.nodes[n].x, t.nodes[n].y]),
            "Node {n} hors du bord {name}"
        );
    }
}

#[test]
fn square_with_named_sides() {
    let t = triangulate(&Geometry {
        outer: square(0.25),
        holes: vec![],
    })
    .unwrap();
    assert!((meshed_area(&t) - 1.0).abs() < 1e-12);
    // pas encore de raffinement : seulement les 16 points du bord, n − 2 triangles
    assert_eq!(t.nodes.len(), 16);
    assert_eq!(t.triangles.len(), 14);
    let bas = t.boundary("bas").unwrap();
    assert_eq!((bas.nodes.len(), bas.faces.len()), (5, 4));
    check_boundary(&t, "bas", |p| p[1].abs() < 1e-12);
    check_boundary(&t, "droite", |p| (p[0] - 1.0).abs() < 1e-12);
    check_boundary(&t, "haut", |p| (p[1] - 1.0).abs() < 1e-12);
    check_boundary(&t, "gauche", |p| p[0].abs() < 1e-12);
}

/// Trou mince [0,1 ; 0,9] × [0,05 ; 0,1] au-dessus d'un bas dense, côtés d'un seul morceau. Tout
/// cercle qui passe par les extrémités du long côté inférieur contient soit les coins hauts du trou,
/// soit des points du bas : ce côté n'est pas une arête de Delaunay et doit être recoupé.
#[test]
fn missing_boundary_piece_is_recovered() {
    let mut outer = square(1.0);
    outer[0].size = 0.05;
    let hole = polygon(
        &[[0.1, 0.05], [0.1, 0.1], [0.9, 0.1], [0.9, 0.05]],
        &[1.0; 4],
        &["trou"; 4],
    );
    let t = triangulate(&Geometry {
        outer,
        holes: vec![hole],
    })
    .unwrap();
    assert!((meshed_area(&t) - (1.0 - 0.8 * 0.05)).abs() < 1e-12);
    let trou = t.boundary("trou").unwrap();
    assert!(trou.faces.len() > 4, "le trou n'a pas été recoupé");
    check_boundary(&t, "trou", |p| {
        let near = |v: f64, w: f64| (v - w).abs() < 1e-12;
        ((near(p[1], 0.05) || near(p[1], 0.1)) && (0.1..=0.9).contains(&p[0]))
            || ((near(p[0], 0.1) || near(p[0], 0.9)) && (0.05..=0.1).contains(&p[1]))
    });
    for e in 0..t.triangles.len() {
        let [x, y] = centroid(&t, e);
        assert!(
            !((0.1..0.9).contains(&x) && (0.05..0.1).contains(&y)),
            "triangle {e} dans le trou"
        );
    }
}

#[test]
fn square_hole_is_carved() {
    let hole = polygon(
        &[[0.4, 0.4], [0.6, 0.4], [0.6, 0.6], [0.4, 0.6]],
        &[1.0; 4],
        &["trou"; 4],
    );
    let t = triangulate(&Geometry {
        outer: square(0.1),
        holes: vec![hole],
    })
    .unwrap();
    assert!((meshed_area(&t) - 0.96).abs() < 1e-12);
    for e in 0..t.triangles.len() {
        let [x, y] = centroid(&t, e);
        assert!(!((0.4..0.6).contains(&x) && (0.4..0.6).contains(&y)));
    }
    let on_hole = |p: [f64; 2]| {
        let edge = |v: f64| (v - 0.4).abs() < 1e-12 || (v - 0.6).abs() < 1e-12;
        (edge(p[0]) && (0.4..=0.6).contains(&p[1])) || (edge(p[1]) && (0.4..=0.6).contains(&p[0]))
    };
    check_boundary(&t, "trou", on_hole);
}

/// Quart de plaque de Kirsch : a = 1, R = 4, trou et bord extérieur en arcs.
fn kirsch(size_hole: f64, size_out: f64) -> Geometry {
    Geometry {
        outer: vec![
            Curve::line([1.0, 0.0], size_out, "sym_y"),
            Curve::arc([4.0, 0.0], [0.0, 0.0], false, size_out, "exterieur"),
            Curve::line([0.0, 4.0], size_out, "sym_x"),
            Curve::arc([0.0, 1.0], [0.0, 0.0], true, size_hole, "trou"),
        ],
        holes: vec![],
    }
}

#[test]
fn quarter_annulus_with_arcs() {
    let t = triangulate(&kirsch(0.1, 0.5)).unwrap();
    let r = |p: [f64; 2]| p[0].hypot(p[1]);
    check_boundary(&t, "trou", |p| (r(p) - 1.0).abs() < 1e-12);
    check_boundary(&t, "exterieur", |p| (r(p) - 4.0).abs() < 1e-12);
    check_boundary(&t, "sym_y", |p| p[1].abs() < 1e-12);
    check_boundary(&t, "sym_x", |p| p[0].abs() < 1e-12);
    // les cordes coupent un peu de l'aire exacte au bord extérieur et en ajoutent au trou :
    // écart borné par la flèche des cordes (h²/8R par unité de longueur d'arc)
    let exact = PI / 4.0 * (16.0 - 1.0);
    let chord_error = 0.5f64.powi(2) / 32.0 * (2.0 * PI) + 0.1f64.powi(2) / 8.0 * (PI / 2.0);
    assert!((meshed_area(&t) - exact).abs() < chord_error);
}

#[test]
fn full_circle_hole() {
    let t = triangulate(&Geometry {
        outer: square(0.1),
        holes: vec![vec![Curve::arc([0.7, 0.5], [0.5, 0.5], true, 0.05, "trou")]],
    })
    .unwrap();
    check_boundary(&t, "trou", |p| {
        ((p[0] - 0.5).hypot(p[1] - 0.5) - 0.2).abs() < 1e-12
    });
    // le polygone inscrit est plus petit que le disque : l'aire maillée dépasse 1 − π r², d'au plus
    // la flèche h²/8r fois le périmètre 2πr
    let area = meshed_area(&t);
    let disc = PI * 0.04;
    assert!(area > 1.0 - disc && area < 1.0 - disc + PI * 0.05f64.powi(2) / 4.0);
    for e in 0..t.triangles.len() {
        let c = centroid(&t, e);
        assert!((c[0] - 0.5).hypot(c[1] - 0.5) > 0.19);
    }
}

#[test]
fn generated_mesh_builds_tri3_elements() {
    let t = triangulate(&kirsch(0.2, 0.8)).unwrap();
    let material = Material::with_poisson(210e3, 0.3).unwrap();
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri3::new(n, material, 1.0, PlaneState::Stress).map(Element::Tri3))
        .collect::<Result<_, _>>()
        .unwrap();
    let mesh = Mesh::new(t.nodes.clone(), elements).unwrap();
    // les faces de bord du mailleur sont celles de Element::face_nodes
    let trou = t.boundary("trou").unwrap();
    for &(e, f) in &trou.faces {
        let nodes = mesh.elements[e].face_nodes(f).unwrap();
        assert!(nodes.iter().all(|n| trou.nodes.contains(n)));
    }
}

#[test]
fn unknown_boundary_lists_names() {
    let t = triangulate(&Geometry {
        outer: square(1.0),
        holes: vec![],
    })
    .unwrap();
    let err = t.boundary("trou").unwrap_err();
    assert!(matches!(err, MailleError::UnknownBoundary { .. }));
    assert_eq!(
        err.to_string(),
        "aucune Curve ne porte le nom « trou » : les bords nommés sont « bas », « droite », « gauche », « haut »"
    );
}
