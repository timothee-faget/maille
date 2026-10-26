//! Import DXF : aires et périmètres exacts du dessin retrouvés sur le maillage, noms de bords par
//! calque, refus avec message.

use std::f64::consts::PI;

use maille_core::mesher::{triangulate, MeshOptions, Triangulation};
use maille_io::dxf::read_geometry;

/// Section ENTITIES à partir de blocs d'entités `(type, [(code, valeur)])`.
fn dxf(entities: &[(&str, &[(i32, &str)])]) -> String {
    let mut s = String::from("0\nSECTION\n2\nENTITIES\n");
    for (kind, groups) in entities {
        s += &format!("0\n{kind}\n");
        for (code, value) in *groups {
            s += &format!("{code}\n{value}\n");
        }
    }
    s + "0\nENDSEC\n0\nEOF\n"
}

fn line<'a>(
    layer: &'a str,
    a: [&'a str; 2],
    b: [&'a str; 2],
) -> (&'static str, Vec<(i32, &'a str)>) {
    (
        "LINE",
        vec![(8, layer), (10, a[0]), (20, a[1]), (11, b[0]), (21, b[1])],
    )
}

fn mesh_area(t: &Triangulation<3>) -> f64 {
    t.triangles
        .iter()
        .map(|&[i, j, k]| {
            let (a, b, c) = (&t.nodes[i], &t.nodes[j], &t.nodes[k]);
            0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)).abs()
        })
        .sum()
}

fn names(t: &Triangulation<3>) -> Vec<&str> {
    let mut names: Vec<&str> = t.boundary_names().collect();
    names.sort();
    names
}

fn mesh(text: &str, size: f64) -> Triangulation<3> {
    let geometry = read_geometry(text, size).unwrap();
    triangulate(&geometry, &MeshOptions::default()).unwrap()
}

fn owned(entities: Vec<(&'static str, Vec<(i32, &str)>)>) -> String {
    let refs: Vec<(&str, &[(i32, &str)])> =
        entities.iter().map(|(k, g)| (*k, g.as_slice())).collect();
    dxf(&refs)
}

#[test]
fn plate_with_circular_hole_has_the_exact_area() {
    let text = owned(vec![
        line("bas", ["0", "0"], ["100", "0"]),
        line("droite", ["100", "0"], ["100", "50"]),
        line("haut", ["100", "50"], ["0", "50"]),
        line("gauche", ["0", "50"], ["0", "0"]),
        (
            "CIRCLE",
            vec![(8, "trou"), (10, "50"), (20, "25"), (40, "10")],
        ),
    ]);
    let t = mesh(&text, 2.0);
    assert_eq!(names(&t).len(), 5);
    assert!(t.boundary("trou").is_ok());
    let exact = 5000.0 - PI * 100.0;
    assert!(
        (mesh_area(&t) - exact).abs() / exact < 1e-3,
        "aire {}",
        mesh_area(&t)
    );
}

#[test]
fn closed_polyline_with_bulges_is_a_slot() {
    // rectangle 20 x 10 fermé par deux demi-cercles de rayon 5 (arrondi 1 = demi-cercle)
    let text = owned(vec![(
        "LWPOLYLINE",
        vec![
            (8, "fente"),
            (70, "1"),
            (10, "0"),
            (20, "0"),
            (10, "20"),
            (20, "0"),
            (42, "1"),
            (10, "20"),
            (20, "10"),
            (10, "0"),
            (20, "10"),
            (42, "1"),
        ],
    )]);
    let t = mesh(&text, 0.5);
    assert_eq!(names(&t), vec!["fente"]);
    let exact = 200.0 + PI * 25.0;
    assert!(
        (mesh_area(&t) - exact).abs() / exact < 1e-3,
        "aire {}",
        mesh_area(&t)
    );
}

#[test]
fn quarter_ring_in_any_order_and_direction() {
    // entités mélangées, certaines dessinées à l'envers
    let text = owned(vec![
        (
            "ARC",
            vec![
                (8, "int"),
                (10, "0"),
                (20, "0"),
                (40, "10"),
                (50, "0"),
                (51, "90"),
            ],
        ),
        line("sym_y", ["0", "10"], ["0", "20"]),
        (
            "ARC",
            vec![
                (8, "ext"),
                (10, "0"),
                (20, "0"),
                (40, "20"),
                (50, "0"),
                (51, "90"),
            ],
        ),
        line("sym_x", ["20", "0"], ["10", "0"]),
    ]);
    let t = mesh(&text, 0.5);
    assert_eq!(names(&t), vec!["ext", "int", "sym_x", "sym_y"]);
    let exact = PI * (400.0 - 100.0) / 4.0;
    assert!(
        (mesh_area(&t) - exact).abs() / exact < 1e-3,
        "aire {}",
        mesh_area(&t)
    );
}

#[test]
fn layer_zero_gets_the_default_name() {
    let text = owned(vec![
        line("0", ["0", "0"], ["4", "0"]),
        line("0", ["4", "0"], ["4", "4"]),
        line("0", ["4", "4"], ["0", "4"]),
        line("0", ["0", "4"], ["0", "0"]),
    ]);
    assert_eq!(names(&mesh(&text, 0.5)), vec!["bord"]);
}

fn error_of(text: &str) -> String {
    read_geometry(text, 1.0).unwrap_err().to_string()
}

#[test]
fn open_contour_is_refused() {
    let text = owned(vec![
        line("a", ["0", "0"], ["4", "0"]),
        line("a", ["4", "0"], ["4", "3"]),
    ]);
    assert!(error_of(&text).contains("contour ouvert"));
}

#[test]
fn branching_contour_is_refused() {
    let text = owned(vec![
        line("a", ["0", "0"], ["4", "0"]),
        line("a", ["4", "0"], ["4", "3"]),
        line("a", ["4", "0"], ["8", "0"]),
    ]);
    assert!(error_of(&text).contains("se rejoignent"));
}

#[test]
fn unsupported_entity_is_named() {
    let text = owned(vec![("SPLINE", vec![(8, "courbe")])]);
    let message = error_of(&text);
    assert!(
        message.contains("SPLINE") && message.contains("courbe"),
        "{message}"
    );
}

#[test]
fn binary_and_garbage_are_refused() {
    assert!(error_of("AutoCAD Binary DXF\r\n\x1a\0").contains("binaire"));
    assert!(error_of("pas un dxf\nvraiment\n").contains("code de groupe"));
    assert!(error_of("0\nSECTION\n2\nHEADER\n0\nENDSEC\n0\nEOF\n").contains("ENTITIES"));
    assert!(error_of(&dxf(&[])).contains("aucune entité"));
}

#[test]
fn incomplete_entity_and_bad_number_are_refused() {
    assert!(error_of(&dxf(&[("LINE", &[(10, "0"), (20, "0")])])).contains("code 11"));
    assert!(error_of(&dxf(&[("CIRCLE", &[(10, "x"), (20, "0"), (40, "1")])])).contains("nombre"));
}

#[test]
fn separate_pieces_are_refused_by_the_mesher() {
    let square = |x: &'static str, y: &'static str, x2: &'static str, y2: &'static str| {
        vec![
            line("a", [x, y], [x2, y]),
            line("a", [x2, y], [x2, y2]),
            line("a", [x2, y2], [x, y2]),
            line("a", [x, y2], [x, y]),
        ]
    };
    let mut entities = square("0", "0", "4", "4");
    entities.extend(square("10", "0", "14", "4"));
    let geometry = read_geometry(&owned(entities), 1.0).unwrap();
    assert!(triangulate(&geometry, &MeshOptions::default()).is_err());
}
