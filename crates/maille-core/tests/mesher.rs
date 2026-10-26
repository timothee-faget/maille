//! Mailleur 2D : bords respectés, trous vidés, bords nommés.

use std::f64::consts::PI;

use maille_core::material::Material;
use maille_core::mesh::{Element, Mesh, PlaneState, Tri3};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions, Triangulation};
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

/// Les deux Nodes de chaque face d'un bord sont dans ses Nodes, chaque Node du bord vérifie `on`,
/// et les faces mises bout à bout ont la longueur `length` (le bord est couvert sans trou).
fn check_boundary(t: &Triangulation<3>, name: &str, length: f64, on: impl Fn([f64; 2]) -> bool) {
    let set = t.boundary(name).unwrap();
    let mut total = 0.0;
    for &(e, f) in &set.faces {
        let (a, b) = (t.triangles[e][f], t.triangles[e][(f + 1) % 3]);
        assert!(set.nodes.contains(&a) && set.nodes.contains(&b));
        total += (t.nodes[a].x - t.nodes[b].x).hypot(t.nodes[a].y - t.nodes[b].y);
    }
    assert!(
        (total - length).abs() < 1e-9 * length,
        "bord {name} : longueur {total}"
    );
    for &n in &set.nodes {
        assert!(
            on([t.nodes[n].x, t.nodes[n].y]),
            "Node {n} hors du bord {name}"
        );
    }
}

/// Angles des triangles en degrés, du plus petit au plus grand.
fn angles(t: &Triangulation<3>, e: usize) -> [f64; 3] {
    let c = corners(t, e);
    let mut a = [0, 1, 2].map(|k| {
        let (p, q, r) = (c[k], c[(k + 1) % 3], c[(k + 2) % 3]);
        let (u, v) = ([q[0] - p[0], q[1] - p[1]], [r[0] - p[0], r[1] - p[1]]);
        (u[0] * v[1] - u[1] * v[0])
            .atan2(u[0] * v[0] + u[1] * v[1])
            .to_degrees()
    });
    a.sort_by(f64::total_cmp);
    a
}

fn min_angle(t: &Triangulation<3>) -> f64 {
    (0..t.triangles.len())
        .map(|e| angles(t, e)[0])
        .fold(180.0, f64::min)
}

/// Plus long côté d'un triangle.
fn longest(t: &Triangulation<3>, e: usize) -> f64 {
    let c = corners(t, e);
    (0..3)
        .map(|k| (c[k][0] - c[(k + 1) % 3][0]).hypot(c[k][1] - c[(k + 1) % 3][1]))
        .fold(0.0, f64::max)
}

#[test]
fn square_with_named_sides() {
    let t = triangulate(
        &Geometry {
            outer: square(0.25),
            holes: vec![],
        },
        &MeshOptions::default(),
    )
    .unwrap();
    assert!((meshed_area(&t) - 1.0).abs() < 1e-12);
    check_boundary(&t, "bas", 1.0, |p| p[1].abs() < 1e-12);
    check_boundary(&t, "droite", 1.0, |p| (p[0] - 1.0).abs() < 1e-12);
    check_boundary(&t, "haut", 1.0, |p| (p[1] - 1.0).abs() < 1e-12);
    check_boundary(&t, "gauche", 1.0, |p| p[0].abs() < 1e-12);
    // chaque Node sert à au moins un triangle
    let mut used = vec![false; t.nodes.len()];
    for tri in &t.triangles {
        for &n in tri {
            used[n] = true;
        }
    }
    assert!(used.iter().all(|&u| u));
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
    let t = triangulate(
        &Geometry {
            outer,
            holes: vec![hole],
        },
        &MeshOptions::default(),
    )
    .unwrap();
    assert!((meshed_area(&t) - (1.0 - 0.8 * 0.05)).abs() < 1e-12);
    let trou = t.boundary("trou").unwrap();
    assert!(trou.faces.len() > 4, "le trou n'a pas été recoupé");
    check_boundary(&t, "trou", 1.7, |p| {
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
    let t = triangulate(
        &Geometry {
            outer: square(0.1),
            holes: vec![hole],
        },
        &MeshOptions::default(),
    )
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
    check_boundary(&t, "trou", 0.8, on_hole);
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
    let t = triangulate(&kirsch(0.1, 0.5), &MeshOptions::default()).unwrap();
    check_boundary(&t, "sym_y", 3.0, |p| p[1].abs() < 1e-12);
    check_boundary(&t, "sym_x", 3.0, |p| p[0].abs() < 1e-12);
    // arcs : cordes un peu plus courtes que l'arc, d'au plus h²/24R par unité de longueur
    for (name, radius, h) in [("trou", 1.0, 0.1), ("exterieur", 4.0, 0.5)] {
        let set = t.boundary(name).unwrap();
        let arc = PI / 2.0 * radius;
        let chords: f64 = set
            .faces
            .iter()
            .map(|&(e, f)| {
                let (a, b) = (t.triangles[e][f], t.triangles[e][(f + 1) % 3]);
                (t.nodes[a].x - t.nodes[b].x).hypot(t.nodes[a].y - t.nodes[b].y)
            })
            .sum();
        assert!(chords < arc && chords > arc * (1.0 - h * h / (24.0 * radius * radius)));
        for &n in &set.nodes {
            assert!((t.nodes[n].x.hypot(t.nodes[n].y) - radius).abs() < 1e-12);
        }
    }
    // les cordes coupent un peu de l'aire exacte au bord extérieur et en ajoutent au trou :
    // écart borné par la flèche des cordes (h²/8R par unité de longueur d'arc)
    let exact = PI / 4.0 * (16.0 - 1.0);
    let chord_error = 0.5f64.powi(2) / 32.0 * (2.0 * PI) + 0.1f64.powi(2) / 8.0 * (PI / 2.0);
    assert!((meshed_area(&t) - exact).abs() < chord_error);
}

#[test]
fn full_circle_hole() {
    let t = triangulate(
        &Geometry {
            outer: square(0.1),
            holes: vec![vec![Curve::arc([0.7, 0.5], [0.5, 0.5], true, 0.05, "trou")]],
        },
        &MeshOptions::default(),
    )
    .unwrap();
    let trou = t.boundary("trou").unwrap();
    for &n in &trou.nodes {
        let p = t.nodes[n];
        assert!(((p.x - 0.5).hypot(p.y - 0.5) - 0.2).abs() < 1e-12);
    }
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
    let t = triangulate(&kirsch(0.2, 0.8), &MeshOptions::default()).unwrap();
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
    let t = triangulate(
        &Geometry {
            outer: square(1.0),
            holes: vec![],
        },
        &MeshOptions::default(),
    )
    .unwrap();
    let err = t.boundary("trou").unwrap_err();
    assert!(matches!(err, MailleError::UnknownBoundary { .. }));
    assert_eq!(
        err.to_string(),
        "aucune Curve ne porte le nom « trou » : les bords nommés sont « bas », « droite », « gauche », « haut »"
    );
}

#[test]
fn minimum_angle_is_guaranteed() {
    let l = [
        [0.0, 0.0],
        [2.0, 0.0],
        [2.0, 1.0],
        [1.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
    ];
    let cases = [
        (
            "carré",
            Geometry {
                outer: square(0.1),
                holes: vec![],
            },
        ),
        (
            "L",
            Geometry {
                outer: polygon(&l, &[0.2; 6], &["bord"; 6]),
                holes: vec![],
            },
        ),
        ("Kirsch", kirsch(0.05, 0.5)),
        (
            "carré troué",
            Geometry {
                outer: square(0.2),
                holes: vec![vec![Curve::arc(
                    [0.6, 0.5],
                    [0.5, 0.5],
                    false,
                    0.02,
                    "trou",
                )]],
            },
        ),
    ];
    for (name, g) in cases {
        for min in [20.0, 25.0, 30.0] {
            let options = MeshOptions {
                min_angle: min,
                ..MeshOptions::default()
            };
            let t = triangulate(&g, &options).unwrap();
            let got = min_angle(&t);
            assert!(got >= min - 1e-9, "{name}, {min}° demandés : {got}°");
            meshed_area(&t);
        }
    }
}

#[test]
fn uniform_size_is_respected() {
    let h = 0.1;
    let t = triangulate(
        &Geometry {
            outer: square(h),
            holes: vec![],
        },
        &MeshOptions::default(),
    )
    .unwrap();
    // rayon circonscrit d'au plus h/√3, donc côtés d'au plus 2h/√3
    for e in 0..t.triangles.len() {
        assert!(longest(&t, e) <= 2.0 * h / 3f64.sqrt() * (1.0 + 1e-9));
    }
    // un triangle de rayon circonscrit h/√3 a au plus l'aire de l'équilatéral de côté h, √3/4 h² :
    // il en faut au moins 1 / (√3/4 h²) ≈ 231. Au-delà de 4 fois plus, le raffinement serait excessif.
    let n = t.triangles.len() as f64;
    let least = 1.0 / (3f64.sqrt() / 4.0 * h * h);
    assert!(n >= least && n <= 4.0 * least, "{n} triangles");
}

#[test]
fn size_grows_away_from_a_fine_boundary() {
    let t = triangulate(&kirsch(0.05, 0.5), &MeshOptions::default()).unwrap();
    let mut near_hole = 0.0_f64;
    let mut near_outer = f64::INFINITY;
    for e in 0..t.triangles.len() {
        let [x, y] = centroid(&t, e);
        let r = x.hypot(y);
        if r < 1.05 {
            near_hole = near_hole.max(longest(&t, e));
        }
        if r > 3.8 {
            near_outer = near_outer.min(longest(&t, e));
        }
    }
    // au trou, h = 0,05 ; au bord extérieur h = 0,5 (gradation 0,3 : la taille y est plafonnée)
    assert!(near_hole <= 2.0 * 0.05 / 3f64.sqrt() * 1.2, "{near_hole}");
    assert!(near_outer > 0.15, "{near_outer}");
}

#[test]
fn meshing_is_deterministic() {
    let o = MeshOptions::default();
    assert_eq!(
        triangulate(&kirsch(0.1, 0.5), &o),
        triangulate(&kirsch(0.1, 0.5), &o)
    );
}

#[test]
fn invalid_options() {
    let g = Geometry {
        outer: square(0.5),
        holes: vec![],
    };
    for min_angle in [0.0, 34.0, f64::NAN] {
        let o = MeshOptions {
            min_angle,
            ..MeshOptions::default()
        };
        assert!(matches!(
            triangulate(&g, &o),
            Err(MailleError::MinAngleOutOfRange { .. })
        ));
    }
    let o = MeshOptions {
        grading: 0.0,
        ..MeshOptions::default()
    };
    assert!(matches!(
        triangulate(&g, &o),
        Err(MailleError::InvalidProperty { .. })
    ));
}

#[test]
fn too_many_nodes_is_reported() {
    // 4 côtés de longueur 1 en morceaux de 10⁻⁷ : 4·10⁷ points de bord, au-delà du plafond
    let g = Geometry {
        outer: square(1e-7),
        holes: vec![],
    };
    let err = triangulate(&g, &MeshOptions::default()).unwrap_err();
    assert_eq!(
        err,
        MailleError::TooManyNodes {
            max: maille_core::mesher::MAX_NODES
        }
    );
    assert!(err.to_string().contains("augmenter les tailles h"));
    // plafond sur le total : 4 côtés de 4·10⁵ morceaux chacun, sous le plafond un par un
    let g = Geometry {
        outer: square(2.5e-6),
        holes: vec![],
    };
    assert!(matches!(
        triangulate(&g, &MeshOptions::default()),
        Err(MailleError::TooManyNodes { .. })
    ));
}

#[test]
fn curve_errors_name_the_curve() {
    let mut outer = square(1.0);
    outer[2].size = -1.0;
    let err = triangulate(
        &Geometry {
            outer,
            holes: vec![],
        },
        &MeshOptions::default(),
    )
    .unwrap_err();
    assert_eq!(
        err.to_string(),
        "la Curve 2 « haut » du contour extérieur a une taille h = -1 : elle doit être finie et strictement positive"
    );
}
