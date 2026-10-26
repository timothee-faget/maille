//! Pression sur une face : tube de Lamé en pression directe (intérieure et extérieure), cube sous
//! pression hydrostatique en Tet4 et Tet10, signe, Amplitude et refus.

mod common;

use common::{box_grid, rect_grid, to_quadratic_tet};
use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Analysis, Bc, Model, Pressure, Step};
use maille_core::results::Stress;
use maille_core::MailleError;

/// Faces (élément, face) dont tous les Nodes vérifient `on`.
fn faces_where(mesh: &Mesh, on: impl Fn(&Node) -> bool) -> Vec<(usize, usize)> {
    let mut found = Vec::new();
    for (i, e) in mesh.elements.iter().enumerate() {
        for f in 0..4 {
            if let Some(nodes) = e.face_nodes(f) {
                if nodes.iter().all(|&n| on(&mesh.nodes[n])) {
                    found.push((i, f));
                }
            }
        }
    }
    found
}

fn pressures(faces: &[(usize, usize)], value: f64) -> Vec<Pressure> {
    faces
        .iter()
        .map(|&(element, face)| Pressure {
            element,
            face,
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

// --- tube de Lamé en déformations planes ---

const E: f64 = 210e3;
const NU: f64 = 0.3;
const INNER: f64 = 10.0;
const OUTER: f64 = 20.0;

/// (σr, σθ, σz) de Lamé en r, pressions pi et po (positives en compression).
fn lame(r: f64, pi: f64, po: f64) -> [f64; 3] {
    let span = OUTER * OUTER - INNER * INNER;
    let a = (pi * INNER * INNER - po * OUTER * OUTER) / span;
    let b = (pi - po) * INNER * INNER * OUTER * OUTER / span;
    [a - b / (r * r), a + b / (r * r), 2.0 * NU * a]
}

fn polar(s: &Stress) -> [f64; 3] {
    let r = s.x.hypot(s.y);
    let (c, sn) = (s.x / r, s.y / r);
    [
        s.sigma_xx * c * c + s.sigma_yy * sn * sn + 2.0 * s.tau_xy * c * sn,
        s.sigma_xx * sn * sn + s.sigma_yy * c * c - 2.0 * s.tau_xy * c * sn,
        s.sigma_zz,
    ]
}

/// Quart de tube en Tri6 de taille `size` sous pression : appuis glissants en x = 0 et y = 0.
fn tube(pi: f64, po: f64, size: f64) -> (Mesh, Vec<Vec<Stress>>) {
    let geometry = Geometry {
        outer: vec![
            Curve::line([INNER, 0.0], size, "sym_x"),
            Curve::arc([OUTER, 0.0], [0.0, 0.0], false, size, "exterieur"),
            Curve::line([0.0, OUTER], size, "sym_y"),
            Curve::arc([0.0, INNER], [0.0, 0.0], true, size, "interieur"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())
        .unwrap()
        .to_quadratic();
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, 1.0, PlaneState::Strain).map(Element::Tri6))
        .collect::<Result<_, _>>()
        .unwrap();
    let mesh = Mesh::new(t.nodes.clone(), elements).unwrap();
    let mut bcs = Vec::new();
    for (name, dof) in [("sym_x", Dof::Uy), ("sym_y", Dof::Ux)] {
        bcs.extend(t.boundary(name).unwrap().nodes.iter().map(|&node| Bc {
            node,
            dof,
            value: 0.0,
        }));
    }
    let mut loads = pressures(&t.boundary("interieur").unwrap().faces, pi);
    loads.extend(pressures(&t.boundary("exterieur").unwrap().faces, po));
    let field = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .with_pressures(loads)
        .unwrap()
        .solve()
        .unwrap();
    (mesh, field.stresses)
}

/// Erreur maximale et erreur RMS des contraintes aux points de Gauss, relatives à `scale`.
fn lame_errors(stresses: &[Vec<Stress>], pi: f64, po: f64, scale: f64) -> (f64, f64) {
    let (mut worst, mut sum, mut count) = (0.0_f64, 0.0, 0);
    for s in stresses.iter().flatten() {
        let (got, want) = (polar(s), lame(s.x.hypot(s.y), pi, po));
        for k in 0..3 {
            let error = (got[k] - want[k]).abs() / scale;
            worst = worst.max(error);
            sum += error * error;
            count += 1;
        }
    }
    (worst, (sum / count as f64).sqrt())
}

#[test]
fn tube_under_internal_pressure_follows_lame() {
    let (_, stresses) = tube(100.0, 0.0, 1.0);
    let (worst, rms) = lame_errors(&stresses, 100.0, 0.0, lame(INNER, 100.0, 0.0)[1].abs());
    assert!(
        worst < 5e-2 && rms < 3e-3,
        "pire {worst:.2e}, RMS {rms:.2e}"
    );
}

#[test]
fn tube_under_external_pressure_follows_lame() {
    let (_, stresses) = tube(0.0, 80.0, 1.0);
    let (worst, rms) = lame_errors(&stresses, 0.0, 80.0, lame(INNER, 0.0, 80.0)[1].abs());
    assert!(
        worst < 5e-2 && rms < 3e-3,
        "pire {worst:.2e}, RMS {rms:.2e}"
    );
}

#[test]
fn tube_under_both_pressures_follows_lame() {
    let (_, stresses) = tube(100.0, 30.0, 1.0);
    let (worst, rms) = lame_errors(&stresses, 100.0, 30.0, lame(INNER, 100.0, 30.0)[1].abs());
    assert!(
        worst < 5e-2 && rms < 3e-3,
        "pire {worst:.2e}, RMS {rms:.2e}"
    );
}

#[test]
fn the_lame_error_shrinks_with_the_mesh_size() {
    // l'écart aux points de Gauss est de la discrétisation (bord en cordes, Delaunay irrégulier), pas
    // un biais : il décroît quand on raffine
    let scale = lame(INNER, 100.0, 0.0)[1].abs();
    let rms = |size| lame_errors(&tube(100.0, 0.0, size).1, 100.0, 0.0, scale).1;
    let (coarse, fine) = (rms(1.0), rms(0.5));
    assert!(fine < coarse / 2.0, "RMS {coarse:.2e} puis {fine:.2e}");
}

// --- cube sous pression hydrostatique ---

fn cube(quadratic: bool) -> Mesh {
    let (nodes, tets) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let material = Material::with_poisson(E, NU).unwrap();
    if quadratic {
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
    }
}

#[test]
fn cube_under_hydrostatic_pressure_shrinks_uniformly() {
    // σ = −p I : ε = −p (1 − 2ν) / E dans les trois directions, u = ε x, exact en Tet4 comme en Tet10
    for quadratic in [false, true] {
        let mesh = cube(quadratic);
        let p = 50.0;
        let tol = 1e-9 * p / E;
        let mut bcs = fixed(&mesh, Dof::Ux, |n| n.x < tol);
        bcs.extend(fixed(&mesh, Dof::Uy, |n| n.y < tol));
        bcs.extend(fixed(&mesh, Dof::Uz, |n| n.z < tol));
        let mut faces = faces_where(&mesh, |n| (n.x - 1.0).abs() < tol);
        faces.extend(faces_where(&mesh, |n| (n.y - 1.0).abs() < tol));
        faces.extend(faces_where(&mesh, |n| (n.z - 1.0).abs() < tol));
        let field = Model::new(mesh.clone(), bcs, vec![])
            .unwrap()
            .with_pressures(pressures(&faces, p))
            .unwrap()
            .solve()
            .unwrap();
        let strain = -p * (1.0 - 2.0 * NU) / E;
        for (i, n) in mesh.nodes.iter().enumerate() {
            for (c, x) in [n.x, n.y, n.z].into_iter().enumerate() {
                let u = field.displacements[3 * i + c];
                assert!((u - strain * x).abs() < 1e-12, "quadratic {quadratic}");
            }
        }
        for s in field.stresses.iter().flatten() {
            assert!((s.sigma_xx + p).abs() < 1e-6 && (s.tau_xy).abs() < 1e-6);
        }
    }
}

// --- signe, Amplitude, refus ---

/// Plaque rectangulaire 4 × 2 en déformations planes, appuis glissants, extrémité x = 4 chargée.
fn plate() -> (Mesh, Vec<Bc>, Vec<(usize, usize)>) {
    let (nodes, tris) = rect_grid(4, 2, 4.0, 2.0);
    let material = Material::with_poisson(E, NU).unwrap();
    let elements = tris
        .iter()
        .map(|&t| Tri3::new(t, material, 1.0, PlaneState::Strain).map(Element::Tri3))
        .collect::<Result<_, _>>()
        .unwrap();
    let mesh = Mesh::new(nodes, elements).unwrap();
    let mut bcs = fixed(&mesh, Dof::Ux, |n| n.x < 1e-9);
    bcs.extend(fixed(&mesh, Dof::Uy, |n| n.y < 1e-9));
    let end = faces_where(&mesh, |n| (n.x - 4.0).abs() < 1e-9);
    (mesh, bcs, end)
}

#[test]
fn a_negative_pressure_pulls() {
    // traction −p = 100 en bout : σx = 100, εx = (1 − ν²) σ / E en déformations planes (εz = 0)
    let (mesh, bcs, end) = plate();
    let field = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .with_pressures(pressures(&end, -100.0))
        .unwrap()
        .solve()
        .unwrap();
    let strain = (1.0 - NU * NU) * 100.0 / E;
    for (i, n) in mesh.nodes.iter().enumerate() {
        let u = field.displacements[2 * i];
        assert!((u - strain * n.x).abs() < 1e-12, "x = {}", n.x);
    }
}

#[test]
fn the_pressure_follows_the_amplitude() {
    let (mesh, bcs, end) = plate();
    let model = Model::new(mesh, bcs, vec![])
        .unwrap()
        .with_pressures(pressures(&end, 40.0))
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

#[test]
fn invalid_pressures_are_typed_errors() {
    let (mesh, bcs, end) = plate();
    let model = Model::new(mesh.clone(), bcs, vec![]).unwrap();
    let (element, face) = end[0];
    let one = |element, face, value| {
        vec![Pressure {
            element,
            face,
            value,
        }]
    };
    assert_eq!(
        model.clone().with_pressures(one(element, 3, 1.0)).err(),
        Some(MailleError::UnknownFace { element, face: 3 })
    );
    assert_eq!(
        model.clone().with_pressures(one(9999, 0, 1.0)).err(),
        Some(MailleError::UnknownElement { element: 9999 })
    );
    assert!(matches!(
        model
            .clone()
            .with_pressures(one(element, face, f64::NAN))
            .err(),
        Some(MailleError::NonFinite { .. })
    ));

    let thermal_material = Material::conductive(1.0).unwrap();
    let thermal_elements = mesh
        .elements
        .iter()
        .map(|e| match e {
            Element::Tri3(t) => Element::Tri3(Tri3 {
                material: thermal_material,
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
        thermal.with_pressures(one(element, face, 1.0)).err(),
        Some(MailleError::PressureStaticOnly)
    );

    let bar = Mesh::new(
        vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
        vec![Element::Bar(Bar {
            nodes: [0, 1],
            material: Material::new(E).unwrap(),
            section: Section::new(1.0).unwrap(),
        })],
    )
    .unwrap();
    assert_eq!(
        Model::new(bar, vec![], vec![])
            .unwrap()
            .with_pressures(one(0, 0, 1.0))
            .err(),
        Some(MailleError::UnsupportedPressureElement { element: 0 })
    );
}
