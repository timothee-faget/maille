//! Mailleur 3D à taille h : surface remaillée puis intérieur raffiné. Volume exact, qualité
//! mesurée, champ affine exact (test de patch) et tube de Lamé épais en Tet10, avec convergence
//! en h.

use std::collections::HashMap;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Tet10, Tet4};
use maille_core::mesher::{
    tetrahedralize, tetrahedralize_brep, Brep, BrepEdge, BrepFace, EdgeCurve, FaceSurface, Patch,
    Surface, TetOptions, Tetrahedralization,
};
use maille_core::model::{Bc, Model, Pressure};
use maille_core::MailleError;

type Triangle = [[f64; 3]; 3];

/// Surface d'une soupe de triangles nommés, sommets soudés à l'identique.
fn surface(patches: &[(&str, Vec<Triangle>)]) -> Surface {
    let mut index: HashMap<[u64; 3], usize> = HashMap::new();
    let mut s = Surface {
        vertices: Vec::new(),
        triangles: Vec::new(),
        patches: Vec::new(),
    };
    for (name, tris) in patches {
        let mut patch = Patch {
            name: name.to_string(),
            triangles: Vec::new(),
        };
        for tri in tris {
            let t = tri.map(|p| {
                *index
                    .entry(p.map(|x| (x + 0.0).to_bits()))
                    .or_insert_with(|| {
                        s.vertices.push(p);
                        s.vertices.len() - 1
                    })
            });
            patch.triangles.push(s.triangles.len());
            s.triangles.push(t);
        }
        s.patches.push(patch);
    }
    s
}

/// Quadrilatère p0 p1 p2 p3, antihoraire vu de dehors, en deux triangles.
fn quad(p: [[f64; 3]; 4]) -> Vec<Triangle> {
    vec![[p[0], p[1], p[2]], [p[0], p[2], p[3]]]
}

/// Cube [0, 1]³, une face par Patch (`x0`, `x1`, `y0`, `y1`, `z0`, `z1`), deux triangles par face.
fn cube() -> Surface {
    let p = |x: f64, y: f64, z: f64| [x, y, z];
    surface(&[
        (
            "x0",
            quad([p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)]),
        ),
        (
            "x1",
            quad([p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)]),
        ),
        (
            "y0",
            quad([p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)]),
        ),
        (
            "y1",
            quad([p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)]),
        ),
        (
            "z0",
            quad([p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)]),
        ),
        (
            "z1",
            quad([p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)]),
        ),
    ])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn corners(t: &Tetrahedralization<4>, e: [usize; 4]) -> [[f64; 3]; 4] {
    e.map(|i| {
        let n = t.nodes[i];
        [n.x, n.y, n.z]
    })
}

fn volume(t: &Tetrahedralization<4>) -> f64 {
    t.tetrahedra
        .iter()
        .map(|&e| {
            let [a, b, c, d] = corners(t, e);
            dot(sub(b, a), cross(sub(c, a), sub(d, a))) / 6.0
        })
        .sum()
}

/// Rayon circonscrit et plus petite arête d'un tétraèdre.
fn radius_and_shortest([a, b, c, d]: [[f64; 3]; 4]) -> (f64, f64) {
    let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
    let (x, y, z) = (cross(v, w), cross(w, u), cross(u, v));
    let det = 2.0 * dot(u, x);
    let o: [f64; 3] =
        std::array::from_fn(|k| (dot(u, u) * x[k] + dot(v, v) * y[k] + dot(w, w) * z[k]) / det);
    let p = [a, b, c, d];
    let mut shortest = f64::INFINITY;
    for i in 0..4 {
        for j in i + 1..4 {
            shortest = shortest.min(dot(sub(p[i], p[j]), sub(p[i], p[j])).sqrt());
        }
    }
    (dot(o, o).sqrt(), shortest)
}

/// Qualité : part des tétraèdres de rapport rayon / arête au plus 2, rapport maximal, et rayon
/// circonscrit maximal rapporté à celui du tétraèdre régulier de côté h.
fn quality(t: &Tetrahedralization<4>, h: f64) -> (f64, f64, f64) {
    let (mut good, mut worst, mut biggest) = (0, 0.0_f64, 0.0_f64);
    for &e in &t.tetrahedra {
        let (r, l) = radius_and_shortest(corners(t, e));
        if r <= 2.0 * l * (1.0 + 1e-9) {
            good += 1;
        }
        worst = worst.max(r / l);
        biggest = biggest.max(r / (h * 6f64.sqrt() / 4.0));
    }
    (good as f64 / t.tetrahedra.len() as f64, worst, biggest)
}

fn sized(h: f64) -> TetOptions {
    TetOptions {
        size: Some(h),
        ..TetOptions::default()
    }
}

#[test]
fn cube_a_la_taille_h() {
    let s = cube();
    for h in [0.5, 0.25, 0.125] {
        let t = tetrahedralize(&s, &sized(h)).unwrap();
        assert!((volume(&t) - 1.0).abs() < 1e-12, "volume à h = {h}");
        // mesuré : 100 % de rapport ≤ 2 (au plus 1,2) ; rayon au plus 1,09 fois la cible, les
        // centres qui empiéteraient sur la surface étant refusés
        let (good, worst, biggest) = quality(&t, h);
        assert!(good >= 0.99 && worst < 4.0, "h = {h} : {good}, {worst}");
        assert!(biggest < 1.25, "h = {h} : rayon {biggest}");
        // les Patches gardent leur aire : 2 × (1/h)² triangles au moins par face
        for name in ["x0", "x1", "y0", "y1", "z0", "z1"] {
            assert!(t.boundary(name).unwrap().faces.len() as f64 >= 2.0 / (h * h));
        }
    }
}

/// Test de patch : u = (εx + γ y, ...) affine imposé sur tout le bord, retrouvé à l'intérieur.
#[test]
fn champ_affine_exact_en_tet4() {
    let t = tetrahedralize(&cube(), &sized(0.25)).unwrap();
    let material = Material::with_poisson(210e3, 0.3).unwrap();
    let elements = t
        .tetrahedra
        .iter()
        .map(|&n| Element::Tet4(Tet4::new(n, material)))
        .collect();
    let mesh = Mesh::new(t.nodes.clone(), elements).unwrap();
    let exact = |n: usize| {
        let p = mesh.nodes[n];
        [
            1e-3 * p.x + 2e-4 * p.y,
            -3e-4 * p.y + 1e-4 * p.z,
            5e-4 * p.z - 2e-4 * p.x,
        ]
    };
    let mut on_boundary: Vec<usize> = t
        .boundary_names()
        .flat_map(|name| t.boundary(name).unwrap().nodes.clone())
        .collect();
    on_boundary.sort_unstable();
    on_boundary.dedup();
    let mut bcs = Vec::new();
    for &node in &on_boundary {
        for (k, dof) in [Dof::Ux, Dof::Uy, Dof::Uz].into_iter().enumerate() {
            bcs.push(Bc {
                node,
                dof,
                value: exact(node)[k],
            });
        }
    }
    let field = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .solve()
        .unwrap();
    let inner = mesh.nodes.len() - on_boundary.len();
    assert!(inner > 0, "des Nodes intérieurs");
    for n in 0..mesh.nodes.len() {
        for (k, dof) in [Dof::Ux, Dof::Uy, Dof::Uz].into_iter().enumerate() {
            let got = field.displacements[mesh.dof_index(n, dof).unwrap()];
            assert!((got - exact(n)[k]).abs() < 1e-10 * 1e-3, "Node {n}");
        }
    }
}

// --- quart de tube de Lamé épais, déformations planes (uz = 0 aux deux bouts) ---

const E: f64 = 210e3;
const NU: f64 = 0.3;
const INNER: f64 = 10.0;
const OUTER: f64 = 20.0;
const HEIGHT: f64 = 5.0;
const PRESSURE: f64 = 100.0;

/// Quart de tube facetté en `m` côtés : Patches `interieur`, `exterieur`, `sym_x` (y = 0),
/// `sym_y` (x = 0), `bas` (z = 0) et `haut` (z = HEIGHT).
fn tube(m: usize) -> Surface {
    let at = |r: f64, i: usize, z: f64| {
        if i == m {
            return [0.0, r, z];
        }
        let a = std::f64::consts::FRAC_PI_2 * i as f64 / m as f64;
        [r * a.cos(), r * a.sin(), z]
    };
    let (mut inner, mut outer, mut bottom, mut top) = (vec![], vec![], vec![], vec![]);
    for i in 0..m {
        let j = i + 1;
        inner.extend(quad([
            at(INNER, i, 0.0),
            at(INNER, i, HEIGHT),
            at(INNER, j, HEIGHT),
            at(INNER, j, 0.0),
        ]));
        outer.extend(quad([
            at(OUTER, i, 0.0),
            at(OUTER, j, 0.0),
            at(OUTER, j, HEIGHT),
            at(OUTER, i, HEIGHT),
        ]));
        bottom.extend(quad([
            at(INNER, i, 0.0),
            at(INNER, j, 0.0),
            at(OUTER, j, 0.0),
            at(OUTER, i, 0.0),
        ]));
        top.extend(quad([
            at(INNER, i, HEIGHT),
            at(OUTER, i, HEIGHT),
            at(OUTER, j, HEIGHT),
            at(INNER, j, HEIGHT),
        ]));
    }
    surface(&[
        ("interieur", inner),
        ("exterieur", outer),
        (
            "sym_x",
            quad([
                at(INNER, 0, 0.0),
                at(OUTER, 0, 0.0),
                at(OUTER, 0, HEIGHT),
                at(INNER, 0, HEIGHT),
            ]),
        ),
        (
            "sym_y",
            quad([
                at(INNER, m, 0.0),
                at(INNER, m, HEIGHT),
                at(OUTER, m, HEIGHT),
                at(OUTER, m, 0.0),
            ]),
        ),
        ("bas", bottom),
        ("haut", top),
    ])
}

/// Déplacement radial de Lamé en déformations planes sous la pression intérieure.
fn lame_ur(r: f64) -> f64 {
    let span = OUTER * OUTER - INNER * INNER;
    let a = PRESSURE * INNER * INNER / span;
    let b = PRESSURE * INNER * INNER * OUTER * OUTER / span;
    (1.0 + NU) / E * ((1.0 - 2.0 * NU) * a * r + b / r)
}

/// Erreur maximale du déplacement radial aux Nodes (rapportée à u_r intérieur) et plus grand |uz|,
/// quart de tube en Tet10 de taille h.
fn tube_errors(m: usize, h: f64) -> (f64, f64, usize) {
    lame_errors(tetrahedralize(&tube(m), &sized(h)).unwrap())
}

/// Erreurs de Lamé (comme `tube_errors`) d'un quart de tube déjà maillé, passé en Tet10.
fn lame_errors(t: Tetrahedralization<4>) -> (f64, f64, usize) {
    let t = t.to_quadratic();
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = t
        .tetrahedra
        .iter()
        .map(|&n| Element::Tet10(Tet10::new(n, material)))
        .collect();
    let mesh = Mesh::new(t.nodes.clone(), elements).unwrap();
    let mut bcs = Vec::new();
    for (name, dof) in [
        ("sym_x", Dof::Uy),
        ("sym_y", Dof::Ux),
        ("bas", Dof::Uz),
        ("haut", Dof::Uz),
    ] {
        bcs.extend(t.boundary(name).unwrap().nodes.iter().map(|&node| Bc {
            node,
            dof,
            value: 0.0,
        }));
    }
    let loads: Vec<Pressure> = t
        .boundary("interieur")
        .unwrap()
        .faces
        .iter()
        .map(|&(element, face)| Pressure {
            element,
            face,
            value: PRESSURE,
        })
        .collect();
    let field = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .with_pressures(loads)
        .unwrap()
        .solve()
        .unwrap();
    let u = |n: usize, dof: Dof| field.displacements[mesh.dof_index(n, dof).unwrap()];
    let (mut worst, mut uz) = (0.0_f64, 0.0_f64);
    for (n, p) in mesh.nodes.iter().enumerate() {
        let r = p.x.hypot(p.y);
        let ur = (u(n, Dof::Ux) * p.x + u(n, Dof::Uy) * p.y) / r;
        worst = worst.max((ur - lame_ur(r)).abs() / lame_ur(INNER));
        uz = uz.max(u(n, Dof::Uz).abs() / lame_ur(INNER));
    }
    (worst, uz, mesh.nodes.len())
}

/// Le tube est facetté : l'erreur vient de la flèche des cordes (1 − cos(π / 4m) du rayon), pas
/// des Tet10. Mesuré à h = 1 : 1,8e-3 à 16 côtés, 4,3e-4 à 32 (1,4 fois la flèche, divisée par 4
/// quand les côtés doublent).
#[test]
fn tube_de_lame_en_tet10() {
    let sag = |m: usize| 1.0 - (std::f64::consts::PI / (4.0 * m as f64)).cos();
    let mut errors = Vec::new();
    for m in [16, 32] {
        let (worst, uz, _) = tube_errors(m, 1.0);
        assert!(worst < 2.0 * sag(m), "{m} côtés : erreur {worst:e}");
        assert!(uz < 1e-3, "{m} côtés : uz {uz:e}");
        errors.push(worst);
    }
    assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
}

/// Le même quart de tube en B-rep exact : arcs de cercle, faces cylindriques et planes.
fn tube_brep() -> Brep {
    // sommets a = (Ri, 0), b = (Ro, 0), c = (0, Ro), d = (0, Ri) en bas, puis en haut
    let mut vertices = Vec::new();
    for z in [0.0, HEIGHT] {
        for (x, y) in [(INNER, 0.0), (OUTER, 0.0), (0.0, OUTER), (0.0, INNER)] {
            vertices.push([x, y, z]);
        }
    }
    let mut edges = Vec::new();
    let mut edge = |start: usize, end: usize, curve: EdgeCurve| {
        edges.push(BrepEdge { start, end, curve });
        edges.len() - 1
    };
    let line = |a: usize, b: usize| EdgeCurve::Line {
        from: vertices[a],
        to: vertices[b],
    };
    let arc = |r: f64, z: f64| EdgeCurve::Arc {
        center: [0.0, 0.0, z],
        u: [1.0, 0.0, 0.0],
        v: [0.0, 1.0, 0.0],
        radii: [r, r],
        start_angle: 0.0,
        sweep: std::f64::consts::FRAC_PI_2,
    };
    let (mut inner, mut outer, mut sym_x, mut sym_y) = ([0; 2], [0; 2], [0; 2], [0; 2]);
    for k in 0..2 {
        let (a, b, c, d) = (4 * k, 4 * k + 1, 4 * k + 2, 4 * k + 3);
        let z = k as f64 * HEIGHT;
        inner[k] = edge(a, d, arc(INNER, z));
        outer[k] = edge(b, c, arc(OUTER, z));
        sym_x[k] = edge(a, b, line(a, b));
        sym_y[k] = edge(d, c, line(d, c));
    }
    let up: Vec<usize> = (0..4).map(|i| edge(i, i + 4, line(i, i + 4))).collect();
    let cylinder = |radius: f64| FaceSurface::Cylinder {
        origin: [0.0; 3],
        axis: [0.0, 0.0, 1.0],
        radius,
    };
    let plane = |origin: [f64; 3], normal: [f64; 3]| FaceSurface::Plane { origin, normal };
    // vu de dehors, chaque face est à gauche de son contour (`true` : arête à rebours)
    let faces = [
        (
            "bas",
            plane([0.0; 3], [0.0, 0.0, -1.0]),
            [
                (inner[0], false),
                (sym_y[0], false),
                (outer[0], true),
                (sym_x[0], true),
            ],
        ),
        (
            "haut",
            plane([0.0, 0.0, HEIGHT], [0.0, 0.0, 1.0]),
            [
                (sym_x[1], false),
                (outer[1], false),
                (sym_y[1], true),
                (inner[1], true),
            ],
        ),
        (
            "sym_x",
            plane([0.0; 3], [0.0, -1.0, 0.0]),
            [
                (sym_x[0], false),
                (up[1], false),
                (sym_x[1], true),
                (up[0], true),
            ],
        ),
        (
            "sym_y",
            plane([0.0; 3], [-1.0, 0.0, 0.0]),
            [
                (up[3], false),
                (sym_y[1], false),
                (up[2], true),
                (sym_y[0], true),
            ],
        ),
        (
            "exterieur",
            cylinder(OUTER),
            [
                (outer[0], false),
                (up[2], false),
                (outer[1], true),
                (up[1], true),
            ],
        ),
        (
            "interieur",
            cylinder(INNER),
            [
                (up[0], false),
                (inner[1], false),
                (up[3], true),
                (inner[0], true),
            ],
        ),
    ];
    Brep {
        vertices,
        edges,
        faces: faces
            .into_iter()
            .map(|(name, surface, lp)| BrepFace {
                name: name.to_string(),
                surface,
                loops: vec![lp.to_vec()],
            })
            .collect(),
    }
}

/// Les sommets des faces courbes sont sur les cylindres : l'erreur suit la flèche des cordes de
/// longueur h sur le rayon intérieur, (h / Ri)² / 8. Mesuré : 4,5e-3 à h = 2 (flèche 5e-3),
/// 1,27e-3 à h = 1 (flèche 1,25e-3).
#[test]
fn tube_de_lame_lu_en_b_rep() {
    let mut errors = Vec::new();
    for h in [2.0, 1.0] {
        let (worst, uz, _) = lame_errors(tetrahedralize_brep(&tube_brep(), &sized(h)).unwrap());
        let sag = (h / INNER).powi(2) / 8.0;
        assert!(worst < 1.5 * sag, "h = {h} : erreur {worst:e}");
        assert!(uz < 1e-3, "h = {h} : uz {uz:e}");
        errors.push(worst);
    }
    assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
}

#[test]
fn b_rep_sans_taille_refusee() {
    let m = tetrahedralize_brep(&tube_brep(), &TetOptions::default())
        .unwrap_err()
        .to_string();
    assert!(m.contains("taille h"), "{m}");
}

#[test]
fn options_invalides_refusees() {
    let s = cube();
    for bad in [
        TetOptions {
            size: Some(0.0),
            ..TetOptions::default()
        },
        TetOptions {
            size: Some(f64::NAN),
            ..TetOptions::default()
        },
        TetOptions {
            size: Some(0.5),
            grading: -1.0,
        },
    ] {
        assert!(matches!(
            tetrahedralize(&s, &bad),
            Err(MailleError::InvalidProperty { .. })
        ));
    }
}
