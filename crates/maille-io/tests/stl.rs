//! Import STL : volume et aire exacts du solide retrouvés sur la surface, patches par solid ou par
//! angle vif, refus avec message.

use std::f64::consts::PI;

use std::collections::HashMap;

use maille_core::mesher::{tetrahedralize, TetOptions, Tetrahedralization};
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

/// Cylindre facetté de n côtés (rayon 1, hauteur 2).
fn cylinder(n: usize) -> Vec<Facet> {
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
    facets
}

/// Volume ½ n sin(2π/n) · h, trois patches.
#[test]
fn cylindre_trois_patches_par_angle_vif() {
    let n = 64;
    let facets = cylinder(n);
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

/// Prisme en L (équerre de 2 × 2 dont on retire un carré, épaisseur 1) : volume 3, arête rentrante.
fn l_prism() -> Vec<Facet> {
    let mut facets = Vec::new();
    for (x, y) in [(0., 0.), (1., 0.), (0., 1.)] {
        for z in [0., 1.] {
            let q = |dx: f64, dy: f64| [x + dx, y + dy, z];
            let (t1, t2) = (
                [q(0., 0.), q(1., 0.), q(1., 1.)],
                [q(0., 0.), q(1., 1.), q(0., 1.)],
            );
            if z == 1. {
                facets.extend([t1, t2]);
            } else {
                facets.extend([[t1[0], t1[2], t1[1]], [t2[0], t2[2], t2[1]]]);
            }
        }
    }
    let ring = [
        (0., 0.),
        (1., 0.),
        (2., 0.),
        (2., 1.),
        (1., 1.),
        (1., 2.),
        (0., 2.),
        (0., 1.),
    ];
    for k in 0..ring.len() {
        let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
        let at = |c: (f64, f64), z: f64| [c.0, c.1, z];
        facets.push([at(p, 0.), at(q, 0.), at(q, 1.)]);
        facets.push([at(p, 0.), at(q, 1.), at(p, 1.)]);
    }
    facets
}

fn tet_volume(t: &Tetrahedralization<4>, tet: [usize; 4]) -> f64 {
    let p = tet.map(|i| [t.nodes[i].x, t.nodes[i].y, t.nodes[i].z]);
    let (u, v, w) = (
        [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]],
        [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]],
        [p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2]],
    );
    (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0]))
        / 6.0
}

/// Volume des tétraèdres égal à celui du solide, tétraèdres positifs, faces de bord = faces à un
/// seul tétraèdre, aire du bord égale à celle de la surface.
fn check_tetrahedralization(s: &Surface, t: &Tetrahedralization<4>) {
    let mut volume = 0.0;
    let mut count: HashMap<[usize; 3], usize> = HashMap::new();
    for &tet in &t.tetrahedra {
        let v = tet_volume(t, tet);
        assert!(v > 0.0, "tétraèdre {tet:?} non positif : {v}");
        volume += v;
        for face in [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]] {
            let mut f = face.map(|k| tet[k]);
            f.sort_unstable();
            *count.entry(f).or_default() += 1;
        }
    }
    assert!(near(volume, s.volume()), "{volume} contre {}", s.volume());
    assert!(count.values().all(|&c| c <= 2));
    let on_surface = count.values().filter(|&&c| c == 1).count();
    let (mut faces, mut area) = (0, 0.0);
    for name in t.boundary_names() {
        let set = t.boundary(name).unwrap();
        for &(e, f) in &set.faces {
            let face = [[0, 1, 2], [0, 3, 1], [1, 3, 2], [2, 3, 0]][f].map(|k| t.tetrahedra[e][k]);
            let [a, b, c] = face.map(|i| [t.nodes[i].x, t.nodes[i].y, t.nodes[i].z]);
            let (u, v) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
            );
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            area += 0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            faces += 1;
        }
    }
    assert_eq!(
        faces, on_surface,
        "faces de bord contre faces à un tétraèdre"
    );
    assert!(
        near(area, s.area()),
        "aire du bord {area} contre {}",
        s.area()
    );
}

#[test]
fn boite_tetraedres_conformes() {
    let s = read_surface(
        &ascii(&[("boite", &flat(&cuboid(2., 3., 4., [0.; 3])))]),
        30.0,
    )
    .unwrap();
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    check_tetrahedralization(&s, &t);
    let names: Vec<_> = t.boundary_names().collect();
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
}

#[test]
fn solids_nommes_donnent_les_bords() {
    let faces = cuboid(1., 1., 1., [0.; 3]);
    let names = ["gauche", "droite", "devant", "derriere", "bas", "haut"];
    let solids: Vec<(&str, &[Facet])> = names
        .iter()
        .zip(&faces)
        .map(|(n, f)| (*n, f.as_slice()))
        .collect();
    let s = read_surface(&ascii(&solids), 30.0).unwrap();
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    check_tetrahedralization(&s, &t);
    assert!(t
        .boundary("gauche")
        .unwrap()
        .nodes
        .iter()
        .all(|&i| t.nodes[i].x == 0.0));
    assert!(t
        .boundary("haut")
        .unwrap()
        .nodes
        .iter()
        .all(|&i| t.nodes[i].z == 1.0));
    assert!(t.boundary("dessous").is_err());
}

#[test]
fn cylindre_tetraedres_conformes() {
    let s = read_surface(&ascii(&[("cyl", &cylinder(64))]), 30.0).unwrap();
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    check_tetrahedralization(&s, &t);
    assert_eq!(t.boundary_names().count(), 3);
}

#[test]
fn prisme_en_l_non_convexe() {
    let s = read_surface(&ascii(&[("l", &l_prism())]), 30.0).unwrap();
    assert!(near(s.volume(), 3.0));
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    check_tetrahedralization(&s, &t);
}

/// 256 sommets cosphériques : le Delaunay a des égalités partout, la surface reste récupérée.
#[test]
fn cylindre_fin_de_128_cotes() {
    let s = read_surface(&ascii(&[("cyl", &cylinder(128))]), 30.0).unwrap();
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    check_tetrahedralization(&s, &t);
}

/// Éventail de 200 triangles de 1,8° au centre : hors de portée, refusé avec un message.
#[test]
fn eventail_trop_effile_refuse() {
    let s = read_surface(&ascii(&[("cyl", &cylinder(200))]), 30.0).unwrap();
    let m = tetrahedralize(&s, &TetOptions::default())
        .unwrap_err()
        .to_string();
    assert!(m.contains("ne se laisse pas récupérer"), "{m}");
}

#[test]
fn surface_inutilisable_refusee() {
    let mut s = read_surface(&ascii(&[("b", &flat(&cuboid(1., 1., 1., [0.; 3])))]), 30.0).unwrap();
    s.triangles[0][0] = 99;
    let m = tetrahedralize(&s, &TetOptions::default())
        .unwrap_err()
        .to_string();
    assert!(m.contains("sommet qui n'existe pas"), "{m}");
}

#[test]
fn tetraedres_quadratiques_gardent_volume_et_bords() {
    let faces = cuboid(1., 2., 3., [0.; 3]);
    let names = ["gauche", "droite", "devant", "derriere", "bas", "haut"];
    let solids: Vec<(&str, &[Facet])> = names
        .iter()
        .zip(&faces)
        .map(|(n, f)| (*n, f.as_slice()))
        .collect();
    let s = read_surface(&ascii(&solids), 30.0).unwrap();
    let linear = tetrahedralize(&s, &TetOptions::default()).unwrap();
    let t = linear.clone().to_quadratic();
    let edges: std::collections::HashSet<_> = linear
        .tetrahedra
        .iter()
        .flat_map(|t| {
            (0..4).flat_map(move |i| (i + 1..4).map(move |j| (t[i].min(t[j]), t[i].max(t[j]))))
        })
        .collect();
    assert_eq!(t.nodes.len(), linear.nodes.len() + edges.len());
    // sommets inchangés, milieux au milieu des arêtes (ordre de Tet10 : 0-1, 1-2, 2-0, 0-3, 1-3, 2-3)
    for (lin, quad) in linear.tetrahedra.iter().zip(&t.tetrahedra) {
        assert_eq!(quad[..4], lin[..]);
        for (k, (a, b)) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
            .into_iter()
            .enumerate()
        {
            let (p, q, m) = (t.nodes[lin[a]], t.nodes[lin[b]], t.nodes[quad[4 + k]]);
            assert!(((p.x + q.x) / 2.0 - m.x).abs() < 1e-12);
            assert!(((p.y + q.y) / 2.0 - m.y).abs() < 1e-12);
            assert!(((p.z + q.z) / 2.0 - m.z).abs() < 1e-12);
        }
    }
    // les Nodes du bord « gauche » (x = 0) comprennent les milieux des arêtes de ses faces
    let gauche = t.boundary("gauche").unwrap();
    assert!(gauche.nodes.iter().all(|&i| t.nodes[i].x == 0.0));
    let on_x0 = (0..t.nodes.len()).filter(|&i| t.nodes[i].x == 0.0).count();
    assert_eq!(gauche.nodes.len(), on_x0);
}

#[test]
fn fichier_msh_des_tetraedres_se_relit() {
    use maille_io::msh::{self, CellKind, MshFile, DOMAIN_SET};
    let s = read_surface(&ascii(&[("b", &flat(&cuboid(1., 1., 1., [0.; 3])))]), 30.0).unwrap();
    let t = tetrahedralize(&s, &TetOptions::default()).unwrap();
    let file = MshFile::from_tetrahedralization(&t).unwrap();
    let back = msh::read(&msh::write(&file).unwrap()).unwrap();
    assert_eq!(back.kind, CellKind::Tet4);
    assert_eq!(back.cells.len(), t.tetrahedra.len());
    assert_eq!(back.sets.len(), 7);
    assert_eq!(back.sets[DOMAIN_SET].elements.len(), t.tetrahedra.len());
    for name in t.boundary_names() {
        assert_eq!(
            back.sets[name].faces.len(),
            t.boundary(name).unwrap().faces.len(),
            "{name}"
        );
    }
}
