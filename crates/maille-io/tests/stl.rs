//! Import STL : volume et aire exacts du solide retrouvés sur la surface, patches par solid ou par
//! angle vif, refus avec message.

use std::f64::consts::PI;

use maille_io::stl::{read_surface, Surface};
use maille_io::IoError;

type Facet = [[f64; 3]; 3];

/// Boîte a × b × c dont un coin est l'origine : 12 triangles sortants, dans l'ordre des faces
/// -x, +x, -y, +y, -z, +z (deux triangles chacune).
fn cuboid(a: f64, b: f64, c: f64, origin: [f64; 3]) -> Vec<Vec<Facet>> {
    let p = |x: f64, y: f64, z: f64| [origin[0] + x * a, origin[1] + y * b, origin[2] + z * c];
    let quad = |q: [[f64; 3]; 4]| vec![[q[0], q[1], q[2]], [q[0], q[2], q[3]]];
    vec![
        quad([p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)]),
        quad([p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)]),
        quad([p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)]),
        quad([p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)]),
        quad([p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)]),
        quad([p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)]),
    ]
}

fn flat(faces: &[Vec<Facet>]) -> Vec<Facet> {
    faces.iter().flatten().copied().collect()
}

fn ascii(solids: &[(&str, &[Facet])]) -> Vec<u8> {
    let mut s = String::new();
    for (name, facets) in solids {
        s += &format!("solid {name}\n");
        for f in *facets {
            s += "  facet normal 0 0 0\n    outer loop\n";
            for v in f {
                s += &format!("      vertex {:e} {:e} {:e}\n", v[0], v[1], v[2]);
            }
            s += "    endloop\n  endfacet\n";
        }
        s += &format!("endsolid {name}\n");
    }
    s.into_bytes()
}

fn binary(facets: &[Facet]) -> Vec<u8> {
    let mut b = vec![0u8; 80];
    b[..5].copy_from_slice(b"solid");
    b.extend((facets.len() as u32).to_le_bytes());
    for f in facets {
        b.extend([0u8; 12]);
        for v in f {
            for x in v {
                b.extend((*x as f32).to_le_bytes());
            }
        }
        b.extend([0u8; 2]);
    }
    b
}

fn message(r: Result<Surface, IoError>) -> String {
    r.expect_err("le fichier devait être refusé").to_string()
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * b.abs()
}

#[test]
fn boite_ascii_volume_aire_et_six_patches() {
    let s = read_surface(
        &ascii(&[("boite", &flat(&cuboid(2., 3., 4., [0.; 3])))]),
        30.0,
    )
    .unwrap();
    assert_eq!(s.vertices.len(), 8);
    assert_eq!(s.triangles.len(), 12);
    assert!(near(s.volume(), 24.0), "{}", s.volume());
    assert!(near(s.area(), 52.0), "{}", s.area());
    let names: Vec<_> = s.patches.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "surface_1",
            "surface_2",
            "surface_3",
            "surface_4",
            "surface_5",
            "surface_6"
        ]
    );
    assert!(s.patches.iter().all(|p| p.triangles.len() == 2));
}

#[test]
fn binaire_donne_la_meme_surface() {
    let facets = flat(&cuboid(2., 3., 4., [0.; 3]));
    let a = read_surface(&ascii(&[("b", &facets)]), 30.0).unwrap();
    let b = read_surface(&binary(&facets), 30.0).unwrap();
    assert_eq!(a, b);
}

#[test]
fn solids_nommes_donnent_les_patches() {
    let faces = cuboid(1., 1., 1., [0.; 3]);
    let names = ["gauche", "droite", "devant", "derriere", "bas", "haut"];
    let solids: Vec<(&str, &[Facet])> = names
        .iter()
        .zip(&faces)
        .map(|(n, f)| (*n, f.as_slice()))
        .collect();
    let s = read_surface(&ascii(&solids), 30.0).unwrap();
    let got: Vec<_> = s.patches.iter().map(|p| p.name.as_str()).collect();
    assert_eq!(got, names);
    assert!(near(s.volume(), 1.0));
}

#[test]
fn normales_rentrantes_corrigees() {
    let inward: Vec<Facet> = flat(&cuboid(1., 2., 3., [0.; 3]))
        .into_iter()
        .map(|[a, b, c]| [a, c, b])
        .collect();
    let s = read_surface(&ascii(&[("b", &inward)]), 30.0).unwrap();
    assert!(near(s.volume(), 6.0), "{}", s.volume());
}

/// Cylindre facetté de n côtés (rayon 1, hauteur 2) : volume ½ n sin(2π/n) · h, trois patches.
#[test]
fn cylindre_trois_patches_par_angle_vif() {
    let n = 64;
    let at = |i: usize, z: f64| {
        let t = 2.0 * PI * i as f64 / n as f64;
        [t.cos(), t.sin(), z]
    };
    let mut facets = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        facets.push([at(i, 0.), at(j, 0.), at(j, 2.)]);
        facets.push([at(i, 0.), at(j, 2.), at(i, 2.)]);
        facets.push([[0., 0., 0.], at(j, 0.), at(i, 0.)]);
        facets.push([[0., 0., 2.], at(i, 2.), at(j, 2.)]);
    }
    let s = read_surface(&ascii(&[("cyl", &facets)]), 30.0).unwrap();
    let exact = 0.5 * n as f64 * (2.0 * PI / n as f64).sin() * 2.0;
    assert!(near(s.volume(), exact), "{} contre {exact}", s.volume());
    assert_eq!(s.patches.len(), 3);
    let mut sizes: Vec<_> = s.patches.iter().map(|p| p.triangles.len()).collect();
    sizes.sort_unstable();
    assert_eq!(sizes, [n, n, 2 * n]);
}

#[test]
fn sommets_voisins_soudes() {
    let mut facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    facets[3][0][0] += 1e-9;
    let s = read_surface(&ascii(&[("b", &facets)]), 30.0).unwrap();
    assert_eq!(s.vertices.len(), 8);
}

#[test]
fn surface_ouverte_refusee() {
    let mut facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    facets.pop();
    let m = message(read_surface(&ascii(&[("b", &facets)]), 30.0));
    assert!(m.contains("surface ouverte"), "{m}");
}

#[test]
fn orientation_incoherente_refusee() {
    let mut facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    facets[0].swap(1, 2);
    let m = message(read_surface(&ascii(&[("b", &facets)]), 30.0));
    assert!(m.contains("orientation incohérente"), "{m}");
}

#[test]
fn deux_pieces_refusees() {
    let mut facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    facets.extend(flat(&cuboid(1., 1., 1., [5., 0., 0.])));
    let m = message(read_surface(&ascii(&[("b", &facets)]), 30.0));
    assert!(m.contains("2 pièces"), "{m}");
}

#[test]
fn facette_degeneree_refusee() {
    let mut facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    facets.push([[0., 0., 0.], [0.5, 0., 0.], [1., 0., 0.]]);
    let m = message(read_surface(&ascii(&[("b", &facets)]), 30.0));
    assert!(m.contains("dégénérée"), "{m}");
}

#[test]
fn fichiers_illisibles_refuses() {
    assert!(message(read_surface(b"bonjour", 30.0)).contains("ni STL binaire"));
    let mut tronque = ascii(&[("b", &flat(&cuboid(1., 1., 1., [0.; 3])))]);
    tronque.truncate(tronque.len() - 40);
    let m = message(read_surface(&tronque, 30.0));
    assert!(m.contains("attendu"), "{m}");
    let facets = flat(&cuboid(1., 1., 1., [0.; 3]));
    assert!(message(read_surface(&ascii(&[("b", &facets)]), 0.0)).contains("angle vif"));
    assert!(message(read_surface(b"solid vide\nendsolid vide\n", 30.0)).contains("aucune facette"));
}
