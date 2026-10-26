//! Fichiers `.msh` : aller-retour sur les maillages du mailleur, numérotation Gmsh du Tet10,
//! orientation, refus avec message.

use maille_core::material::Material;
use maille_core::mesh::{Element, Mesh, Tet10};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::shape::{tet10, TET10_EDGES};
use maille_io::msh::{read, write, CellKind, MshFile, Set, DOMAIN_SET};
use maille_io::IoError;

fn square_geometry() -> Geometry {
    let pts = [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]];
    let names = ["bas", "droite", "haut", "gauche"];
    Geometry {
        outer: (0..4).map(|i| Curve::line(pts[i], 0.3, names[i])).collect(),
        holes: vec![],
    }
}

#[test]
fn tri3_round_trip_is_exact() {
    let t = triangulate(&square_geometry(), &MeshOptions::default()).unwrap();
    let file = MshFile::from_triangulation(&t).unwrap();
    assert_eq!(file.kind, CellKind::Tri3);
    let back = read(&write(&file).unwrap()).unwrap();
    assert_eq!(back, file);
    assert_eq!(back.sets[DOMAIN_SET].elements.len(), t.triangles.len());
    assert_eq!(back.sets["bas"].faces, t.boundary("bas").unwrap().faces);
}

#[test]
fn tri6_round_trip_is_exact() {
    let t = triangulate(&square_geometry(), &MeshOptions::default())
        .unwrap()
        .to_quadratic();
    let file = MshFile::from_triangulation(&t).unwrap();
    assert_eq!(file.kind, CellKind::Tri6);
    let back = read(&write(&file).unwrap()).unwrap();
    assert_eq!(back, file);
    // les Nodes milieux du bord sont dans le Set
    assert_eq!(back.sets["haut"].nodes, t.boundary("haut").unwrap().nodes);
}

#[test]
fn boundary_named_like_the_domain_is_refused() {
    let mut g = square_geometry();
    g.outer[0].name = DOMAIN_SET.to_string();
    let t = triangulate(&g, &MeshOptions::default()).unwrap();
    assert!(matches!(
        MshFile::from_triangulation(&t),
        Err(IoError::DuplicateSetName { .. })
    ));
}

/// Les 10 Nodes d'un tétraèdre dans l'ordre de Gmsh : sommets, puis milieux (0,1) (1,2) (2,0)
/// (3,0) (3,2) (3,1).
fn gmsh_tet10(v: [[f64; 3]; 4]) -> Vec<[f64; 3]> {
    let mid = |a: usize, b: usize| std::array::from_fn(|k| (v[a][k] + v[b][k]) / 2.0);
    let mut p = v.to_vec();
    p.extend([
        mid(0, 1),
        mid(1, 2),
        mid(2, 0),
        mid(3, 0),
        mid(3, 2),
        mid(3, 1),
    ]);
    p
}

/// Deux Tet10 écrits à la Gmsh : le second est retourné (sommets 0 et 1 permutés).
fn tet10_file() -> String {
    let a = gmsh_tet10([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
    let b = gmsh_tet10([[3., 0., 0.], [2., 0., 0.], [2., 1., 0.], [2., 0., 1.]]);
    let coords: String = a
        .iter()
        .chain(&b)
        .map(|p| format!("{} {} {}\n", p[0], p[1], p[2]))
        .collect();
    let tags: String = (1..=20).map(|i| format!("{i}\n")).collect();
    let nodes = |from: usize| {
        (from..from + 10)
            .map(|i| format!(" {i}"))
            .collect::<String>()
    };
    format!(
        "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n\
$PhysicalNames\n3\n0 3 \"coin\"\n2 2 \"base\"\n3 1 \"corps\"\n$EndPhysicalNames\n\
$Entities\n1 0 1 1\n1 0 0 0 1 3 0 0 1 0 1\n\
1 0 0 0 1 1 0 1 2\n1 0 0 0 1 1 1 1 1 0\n$EndEntities\n\
$Nodes\n1 20 1 20\n3 1 0 20\n{tags}{coords}$EndNodes\n\
$Elements\n3 4 1 4\n0 1 15 1\n1 1\n2 1 2 1\n2 1 3 2\n3 1 11 2\n3{}\n4{}\n$EndElements\n",
        nodes(1),
        nodes(11)
    )
}

#[test]
fn tet10_gmsh_numbering_and_orientation() {
    let file = read(&tet10_file()).unwrap();
    assert_eq!(file.kind, CellKind::Tet10);
    assert_eq!(file.cells.len(), 2);

    let quadratic = |p: [f64; 3]| 1.0 + p[0] + 2.0 * p[1] - p[2] + p[0] * p[1] + 3.0 * p[2] * p[2];
    for cell in &file.cells {
        // chaque milieu est au milieu de son arête, dans la numérotation de Maille
        let pt = |i: usize| {
            let n = file.nodes[cell[i]];
            [n.x, n.y, n.z]
        };
        for (k, &(a, b)) in TET10_EDGES.iter().enumerate() {
            for c in 0..3 {
                assert!((pt(4 + k)[c] - (pt(a)[c] + pt(b)[c]) / 2.0).abs() < 1e-14);
            }
        }
        // un champ quadratique est exact : interpolation aux Nodes au point (ξ, η, ζ)
        let (n, _) = tet10(0.2, 0.3, 0.1);
        let x: [f64; 3] = std::array::from_fn(|c| (0..10).map(|i| n[i] * pt(i)[c]).sum());
        let u: f64 = (0..10).map(|i| n[i] * quadratic(pt(i))).sum();
        assert!((u - quadratic(x)).abs() < 1e-12);
    }

    // orientation positive, arêtes droites : Mesh::new les vérifie
    let m = Material::with_poisson(1.0, 0.3).unwrap();
    let elements = file
        .cells
        .iter()
        .map(|c| Element::Tet10(Tet10::new(c.clone().try_into().unwrap(), m)))
        .collect();
    Mesh::new(file.nodes.clone(), elements).unwrap();
}

#[test]
fn tet10_sets() {
    let file = read(&tet10_file()).unwrap();
    assert_eq!(file.sets["corps"].elements, [0, 1]);
    assert_eq!(file.sets["coin"].nodes, [0]);
    let base = &file.sets["base"];
    assert_eq!(base.faces, [(0, 0)]);
    // face 0 : sommets 0-1-2 et milieux (0,1) (1,2) (2,0)
    let expected: Vec<usize> = [0, 1, 2, 4, 5, 6]
        .iter()
        .map(|&p| file.cells[0][p])
        .collect();
    let mut expected = expected;
    expected.sort_unstable();
    assert_eq!(base.nodes, expected);
}

#[test]
fn tet10_round_trip_is_exact() {
    let file = read(&tet10_file()).unwrap();
    assert_eq!(read(&write(&file).unwrap()).unwrap(), file);
}

#[test]
fn overlapping_sets_round_trip() {
    let t = triangulate(&square_geometry(), &MeshOptions::default()).unwrap();
    let mut file = MshFile::from_triangulation(&t).unwrap();
    file.sets.insert(
        "moitie".to_string(),
        Set {
            elements: (0..t.triangles.len() / 2).collect(),
            ..Set::default()
        },
    );
    file.sets.insert(
        "coin".to_string(),
        Set {
            nodes: vec![0, 3],
            ..Set::default()
        },
    );
    assert_eq!(read(&write(&file).unwrap()).unwrap(), file);
}

const TRIANGLE: &str = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n\
$PhysicalNames\n1\n2 1 \"plaque\"\n$EndPhysicalNames\n\
$Entities\n0 0 1 0\n1 0 0 0 1 1 0 1 1 0\n$EndEntities\n\
$Nodes\n1 3 1 3\n2 1 0 3\n1\n2\n3\n0 0 0\n1 0 0\n0 1 0\n$EndNodes\n\
$Elements\n1 1 1 1\n2 1 2 1\n1 1 2 3\n$EndElements\n";

#[test]
fn clockwise_triangle_is_flipped() {
    let file = read(&TRIANGLE.replace("1 1 2 3", "1 1 3 2")).unwrap();
    let c = &file.cells[0];
    let p = |i: usize| file.nodes[c[i]];
    let area2 = (p(1).x - p(0).x) * (p(2).y - p(0).y) - (p(2).x - p(0).x) * (p(1).y - p(0).y);
    assert!(area2 > 0.0);
    assert_eq!(file.sets["plaque"].elements, [0]);
}

#[test]
fn unnamed_group_is_named_by_its_number() {
    let text = TRIANGLE
        .replace("2 1 \"plaque\"", "")
        .replace("$PhysicalNames\n1\n", "");
    let text = text.replace("$EndPhysicalNames\n", "");
    assert_eq!(read(&text).unwrap().sets["1"].elements, [0]);
}

fn message(text: &str) -> String {
    read(text).unwrap_err().to_string()
}

#[test]
fn refuses_old_version() {
    let text = TRIANGLE.replace("4.1 0 8", "2.2 0 8");
    assert!(matches!(
        read(&text),
        Err(IoError::UnsupportedVersion { line: 2, .. })
    ));
    assert!(message(&text).contains("Mesh.MshFileVersion = 4.1"));
}

#[test]
fn refuses_binary() {
    let text = TRIANGLE.replace("4.1 0 8", "4.1 1 8");
    assert!(matches!(
        read(&text),
        Err(IoError::BinaryFormat { line: 2 })
    ));
    assert!(message(&text).contains("Mesh.Binary = 0"));
}

#[test]
fn refuses_hexahedra() {
    let text = TRIANGLE.replace("2 1 2 1\n1 1 2 3", "3 1 5 1\n1 1 2 3 1 2 3 1 2");
    assert!(matches!(
        read(&text),
        Err(IoError::UnsupportedElement { gmsh_type: 5, .. })
    ));
    assert!(message(&text).contains("hexaèdre"));
}

#[test]
fn refuses_mixed_orders() {
    let text = TRIANGLE.replace(
        "1 1 1 1\n2 1 2 1\n1 1 2 3",
        "2 2 1 2\n2 1 2 1\n1 1 2 3\n2 1 9 1\n2 1 2 3 1 2 3",
    );
    assert!(matches!(read(&text), Err(IoError::MixedCellKinds { .. })));
}

#[test]
fn refuses_file_without_cells() {
    let text = TRIANGLE.replace("1 1 1 1\n2 1 2 1\n1 1 2 3", "0 0 0 0");
    assert_eq!(read(&text), Err(IoError::NoCells));
}

#[test]
fn refuses_missing_section() {
    let text = TRIANGLE
        .replace("$Elements", "$Autre")
        .replace("$EndElements", "$EndAutre");
    assert_eq!(
        read(&text),
        Err(IoError::MissingSection { name: "Elements" })
    );
}

#[test]
fn refuses_unclosed_section() {
    let text = TRIANGLE.replace("$EndElements\n", "");
    assert!(matches!(read(&text), Err(IoError::Syntax { .. })));
}

#[test]
fn refuses_unknown_node() {
    let text = TRIANGLE.replace("1 1 2 3", "1 1 2 9");
    assert_eq!(read(&text), Err(IoError::UnknownNode { line: 25, tag: 9 }));
}

#[test]
fn refuses_non_planar_2d() {
    let text = TRIANGLE.replace("0 1 0\n$EndNodes", "0 1 0.5\n$EndNodes");
    let r = read(&text);
    assert!(matches!(r, Err(IoError::NonPlanar { .. })), "{r:?}");
}

#[test]
fn refuses_bad_number_with_line() {
    let text = TRIANGLE.replace("1 0 0\n", "1 zéro 0\n");
    assert!(matches!(read(&text), Err(IoError::Syntax { line: 19, .. })));
}

#[test]
fn refuses_duplicate_group_name() {
    // une ligne (bord) et la surface portent le même nom
    let text = TRIANGLE
        .replace("1\n2 1 \"plaque\"", "2\n2 1 \"x\"\n1 2 \"x\"")
        .replace("0 0 1 0\n", "0 1 1 0\n1 0 0 0 1 0 0 1 2 0\n")
        .replace("1 1 1 1\n2 1 2 1", "2 2 1 2\n1 1 1 1\n2 1 2\n2 1 2 1");
    assert!(matches!(read(&text), Err(IoError::DuplicateSetName { .. })));
}

#[test]
fn refuses_boundary_that_is_not_a_face() {
    // le Node 4 n'appartient à aucun triangle : la ligne 1-4 n'est la face de personne
    let text = TRIANGLE
        .replace("2 1 \"plaque\"", "1 1 \"bord\"")
        .replace(
            "0 0 1 0\n1 0 0 0 1 1 0 1 1 0",
            "0 1 0 0\n1 0 0 0 1 0 0 1 1 0",
        )
        .replace(
            "1 3 1 3\n2 1 0 3\n1\n2\n3\n",
            "1 4 1 4\n2 1 0 4\n1\n2\n3\n4\n",
        )
        .replace("0 1 0\n$EndNodes", "0 1 0\n2 2 0\n$EndNodes")
        .replace(
            "1 1 1 1\n2 1 2 1\n1 1 2 3",
            "2 2 1 2\n2 1 2 1\n1 1 2 3\n1 1 1 1\n2 1 4",
        );
    assert!(matches!(
        read(&text),
        Err(IoError::BoundaryFaceNotFound { set, .. }) if set == "bord"
    ));
}

#[test]
fn write_refuses_inconsistent_file() {
    let t = triangulate(&square_geometry(), &MeshOptions::default()).unwrap();
    let good = MshFile::from_triangulation(&t).unwrap();

    let mut file = good.clone();
    file.sets.get_mut("bas").unwrap().elements.push(0);
    assert!(matches!(write(&file), Err(IoError::Invalid { .. })));

    let mut file = good.clone();
    file.cells[0][0] = file.nodes.len();
    assert!(matches!(write(&file), Err(IoError::Invalid { .. })));

    let mut file = good.clone();
    file.sets.insert(
        "a\"b".to_string(),
        Set {
            nodes: vec![0],
            ..Set::default()
        },
    );
    assert!(matches!(write(&file), Err(IoError::Invalid { .. })));

    let mut file = good;
    file.sets.insert("vide".to_string(), Set::default());
    assert!(matches!(write(&file), Err(IoError::Invalid { .. })));
}
