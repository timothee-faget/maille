//! Thermique stationnaire sur triangles et tétraèdres, comparée aux solutions analytiques.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Analysis, Bc, Convection, HeatSource, Load, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

/// Conductivité (unités cohérentes au choix de l'utilisateur, ex. W/(mm·K)).
const K: f64 = 50.0;
const T: f64 = 2.0;
const L: f64 = 4.0;
const H: f64 = 1.0;

fn material() -> Material {
    Material::conductive(K).unwrap()
}

/// Mesh de triangles L × H (épaisseur T), Tri3 ou Tri6.
fn plane_mesh(nx: usize, ny: usize, quadratic: bool) -> Mesh {
    let (nodes, tris) = rect_grid(nx, ny, L, H);
    if quadratic {
        let q = to_quadratic(&nodes, &tris);
        let elements = q
            .elements
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tris
            .iter()
            .map(|&n| Element::Tri3(Tri3::new(n, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    }
}

/// Pavé de tétraèdres L × H × T, Tet4 ou Tet10.
fn solid_mesh(n: [usize; 3], quadratic: bool) -> Mesh {
    let (nodes, tets) = box_grid(n, [L, H, T]);
    if quadratic {
        let q = to_quadratic_tet(&nodes, &tets);
        let elements = q
            .elements
            .iter()
            .map(|&t| Element::Tet10(Tet10::new(t, material())))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tets
            .iter()
            .map(|&t| Element::Tet4(Tet4::new(t, material())))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    }
}

/// Les 4 Mesh de la validation : Tri3, Tri6, Tet4, Tet10. Section A = H × T dans tous les cas.
fn all_meshes() -> Vec<(&'static str, Mesh)> {
    vec![
        ("Tri3", plane_mesh(4, 2, false)),
        ("Tri6", plane_mesh(4, 2, true)),
        ("Tet4", solid_mesh([4, 1, 2], false)),
        ("Tet10", solid_mesh([4, 1, 2], true)),
    ]
}

/// BC T = `value` sur tous les Nodes du plan x = `x`.
fn fix_plane(mesh: &Mesh, x: f64, value: f64) -> Vec<Bc> {
    nodes_at(mesh, x)
        .map(|node| Bc {
            node,
            dof: Dof::T,
            value,
        })
        .collect()
}

fn nodes_at(mesh: &Mesh, x: f64) -> impl Iterator<Item = usize> + '_ {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(move |(_, p)| (p.x - x).abs() < 1e-12)
        .map(|(n, _)| n)
}

fn solve(mesh: &Mesh, bcs: Vec<Bc>, loads: Vec<Load>) -> Field {
    Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, loads)
        .unwrap()
        .solve()
        .unwrap()
}

#[test]
fn wall_with_imposed_temperatures_is_exact_with_the_right_heat_flow() {
    // T(x) = T1 + (T2 − T1) x / L, q = −k (T2 − T1) / L, puissance traversante q A
    let (t1, t2) = (100.0, 20.0);
    let q = -K * (t2 - t1) / L;
    let power = q * H * T;
    for (name, mesh) in all_meshes() {
        let mut bcs = fix_plane(&mesh, 0.0, t1);
        bcs.extend(fix_plane(&mesh, L, t2));
        let field = solve(&mesh, bcs, Vec::new());

        assert_eq!(field.temperatures.len(), mesh.nodes.len());
        assert!(field.displacements.is_empty() && field.stresses.is_empty());
        for (n, p) in mesh.nodes.iter().enumerate() {
            let exact = t1 + (t2 - t1) * p.x / L;
            assert!(
                (field.temperatures[n] - exact).abs() < 1e-10,
                "{name} Node {n}"
            );
        }
        for f in field.heat_fluxes.iter().flatten() {
            assert!((f.qx - q).abs() < 1e-9 * q, "{name} qx = {}", f.qx);
            assert!(f.qy.abs() < 1e-9 * q && f.qz.abs() < 1e-9 * q, "{name}");
        }
        // la chaleur entre par x = 0 et sort par x = L
        let sum = |x: f64| nodes_at(&mesh, x).map(|n| field.reactions[n]).sum::<f64>();
        assert!(
            (sum(0.0) - power).abs() < 1e-9 * power,
            "{name} entrée {}",
            sum(0.0)
        );
        assert!(
            (sum(L) + power).abs() < 1e-9 * power,
            "{name} sortie {}",
            sum(L)
        );
    }
}

#[test]
fn heat_injected_at_the_end_gives_t1_plus_p_x_over_k_a() {
    // T = T1 en x = 0, puissance P répartie uniformément sur l'arête x = L (parts consistantes)
    let (t1, p) = (20.0, 400.0);
    let (nx, ny) = (4, 2);
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let dy = H / ny as f64;
    for quadratic in [false, true] {
        let mesh = plane_mesh(nx, ny, quadratic);
        let mut shares: Vec<(usize, f64)> = Vec::new();
        if quadratic {
            let (nodes, tris) = rect_grid(nx, ny, L, H);
            let q = to_quadratic(&nodes, &tris);
            for j in 0..ny {
                let (a, b) = (idx(nx, j), idx(nx, j + 1));
                let segment = p * dy / H;
                shares.extend([
                    (a, segment / 6.0),
                    (q.mid(a, b), segment * 4.0 / 6.0),
                    (b, segment / 6.0),
                ]);
            }
        } else {
            for j in 0..=ny {
                let tributary = if j == 0 || j == ny { dy / 2.0 } else { dy };
                shares.push((idx(nx, j), p * tributary / H));
            }
        }
        let loads = shares
            .into_iter()
            .map(|(node, value)| Load {
                node,
                dof: Dof::T,
                value,
            })
            .collect();
        let field = solve(&mesh, fix_plane(&mesh, 0.0, t1), loads);
        for (n, node) in mesh.nodes.iter().enumerate() {
            let exact = t1 + p * node.x / (K * H * T);
            assert!(
                (field.temperatures[n] - exact).abs() < 1e-10,
                "quadratique {quadratic}, Node {n}"
            );
        }
        let entering: f64 = nodes_at(&mesh, 0.0).map(|n| field.reactions[n]).sum();
        assert!(
            (entering + p).abs() < 1e-9 * p,
            "la chaleur injectée sort par x = 0"
        );
    }
}

#[test]
fn thermal_inputs_are_typed_errors() {
    let mesh = plane_mesh(1, 1, false);
    let thermal =
        |mesh: Mesh, bcs: Vec<Bc>| Model::with_analysis(mesh, Analysis::Thermal, bcs, Vec::new());

    // DOF mécanique en thermique, DOF T en statique
    let ux = vec![Bc {
        node: 0,
        dof: Dof::Ux,
        value: 0.0,
    }];
    assert_eq!(
        thermal(mesh.clone(), ux),
        Err(MailleError::InvalidDof {
            node: 0,
            dof: Dof::Ux
        })
    );
    let t = vec![Bc {
        node: 9,
        dof: Dof::T,
        value: 0.0,
    }];
    assert_eq!(
        thermal(mesh.clone(), t),
        Err(MailleError::UnknownNode { node: 9 })
    );

    // Material sans k ; Material sans E en statique
    let steel = Material::with_poisson(210e3, 0.3).unwrap();
    let no_k = Tri3::new([0, 1, 2], steel, 1.0, PlaneState::Stress).unwrap();
    let nodes = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(1.0, 0.0),
        Node::new_2d(0.0, 1.0),
    ];
    let mesh_no_k = Mesh::new(nodes, vec![Element::Tri3(no_k)]).unwrap();
    assert_eq!(
        thermal(mesh_no_k, Vec::new()),
        Err(MailleError::MissingConductivity { element: 0 })
    );
    assert_eq!(
        Model::new(mesh, Vec::new(), Vec::new()),
        Err(MailleError::MissingYoung { element: 0 })
    );
    assert!(Material::conductive(0.0).is_err() && steel.with_conductivity(-1.0).is_err());

    // Bar : pas de thermique
    let bar = Bar {
        nodes: [0, 1],
        material: steel.with_conductivity(K).unwrap(),
        section: Section::new(1.0).unwrap(),
    };
    let bars = Mesh::new(
        vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
        vec![Element::Bar(bar)],
    )
    .unwrap();
    assert_eq!(
        thermal(bars, Vec::new()),
        Err(MailleError::UnsupportedThermalElement { element: 0 })
    );
}

#[test]
fn a_thermal_field_has_no_stresses() {
    let mesh = plane_mesh(1, 1, false);
    let field = solve(&mesh, fix_plane(&mesh, 0.0, 1.0), Vec::new());
    assert_eq!(nodal_stresses(&mesh, &field), Err(MailleError::NoStresses));
}

#[test]
fn one_material_serves_both_analyses() {
    // E, ν et k dans le même Material : le Mesh se résout en statique comme en thermique
    let m = Material::with_poisson(210e3, 0.3)
        .unwrap()
        .with_conductivity(K)
        .unwrap();
    let (nodes, tris) = rect_grid(2, 1, L, H);
    let elements = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, m, T, PlaneState::Stress).unwrap()))
        .collect();
    let mesh = Mesh::new(nodes, elements).unwrap();
    assert!(Model::new(mesh.clone(), Vec::new(), Vec::new()).is_ok());
    let field = solve(&mesh, fix_plane(&mesh, 0.0, 5.0), Vec::new());
    for t in &field.temperatures {
        assert!((t - 5.0).abs() < 1e-12, "sans flux, T uniforme");
    }
}

/// Dalle de largeur L (T = 0 en x = 0 et x = L), source Q uniforme : T = Q x (L − x) / 2k et
/// qx = Q (x − L/2). Vérifie le bilan (toute la chaleur produite, Q V, sort par les deux faces) et
/// renvoie les écarts maximaux de T aux Nodes et de qx aux points de Gauss, rapportés à leurs
/// maximums Q L² / 8k et Q L / 2.
fn slab_errors(mesh: &Mesh) -> (f64, f64) {
    let q = 1e3;
    let sources = (0..mesh.elements.len())
        .map(|element| HeatSource { element, value: q })
        .collect();
    let mut bcs = fix_plane(mesh, 0.0, 0.0);
    bcs.extend(fix_plane(mesh, L, 0.0));
    let field = Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, Vec::new())
        .unwrap()
        .with_heat_sources(sources)
        .unwrap()
        .solve()
        .unwrap();
    let produced = q * L * H * T;
    let out: f64 = field.reactions.iter().sum();
    assert!(
        (out + produced).abs() < 1e-9 * produced,
        "bilan {out} contre {produced}"
    );
    let temperature = mesh
        .nodes
        .iter()
        .zip(&field.temperatures)
        .map(|(p, t)| (t - q * p.x * (L - p.x) / (2.0 * K)).abs())
        .fold(0.0, f64::max);
    let flux = field
        .heat_fluxes
        .iter()
        .flatten()
        .map(|f| (f.qx - q * (f.x - L / 2.0)).abs())
        .fold(0.0, f64::max);
    (temperature / (q * L * L / (8.0 * K)), flux / (q * L / 2.0))
}

#[test]
fn quadratic_elements_are_exact_for_a_uniform_heat_source() {
    // T quadratique et q linéaire : dans l'espace du Tri6 et du Tet10
    for mesh in [plane_mesh(2, 1, true), solid_mesh([2, 1, 1], true)] {
        let (t, q) = slab_errors(&mesh);
        assert!(t < 1e-10 && q < 1e-10, "écarts T {t}, q {q}");
    }
}

#[test]
fn linear_elements_are_exact_at_the_nodes_and_their_flux_converges() {
    // problème qui ne dépend que de x avec des charges consistantes : les Elements linéaires sont
    // exacts aux Nodes (comme en 1D), mais leur flux est constant par Element (erreur en h)
    let tri3 = [
        plane_mesh(4, 2, false),
        plane_mesh(8, 2, false),
        plane_mesh(16, 2, false),
    ];
    let tet4 = [
        solid_mesh([4, 1, 1], false),
        solid_mesh([8, 1, 1], false),
        solid_mesh([16, 1, 1], false),
    ];
    for meshes in [tri3, tet4] {
        let errors: Vec<(f64, f64)> = meshes.iter().map(slab_errors).collect();
        for &(t, q) in &errors {
            assert!(t < 1e-10, "{errors:?}");
            assert!(
                q > 1e-3,
                "un flux constant par Element ne peut pas être exact : {errors:?}"
            );
        }
        // ordre 1 : l'erreur de flux est divisée par ~2 à chaque raffinement
        for w in errors.windows(2) {
            let ratio = w[0].1 / w[1].1;
            assert!((1.8..2.2).contains(&ratio), "{errors:?}");
        }
    }
}

#[test]
fn heat_sources_need_a_thermal_model_and_a_known_element() {
    let mesh = plane_mesh(1, 1, false);
    let source = |element| {
        vec![HeatSource {
            element,
            value: 1.0,
        }]
    };
    let steel = Material::with_poisson(210e3, 0.3).unwrap();
    let (nodes, tris) = rect_grid(1, 1, L, H);
    let elements = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, steel, T, PlaneState::Stress).unwrap()))
        .collect();
    let static_model =
        Model::new(Mesh::new(nodes, elements).unwrap(), Vec::new(), Vec::new()).unwrap();
    assert_eq!(
        static_model.with_heat_sources(source(0)),
        Err(MailleError::ThermalOnly)
    );
    let thermal = Model::with_analysis(mesh, Analysis::Thermal, Vec::new(), Vec::new()).unwrap();
    assert_eq!(
        thermal.with_heat_sources(source(2)),
        Err(MailleError::UnknownElement { element: 2 })
    );
}

/// (Element, face) dont tous les Nodes vérifient `on`.
fn faces_where(mesh: &Mesh, on: impl Fn(&Node) -> bool) -> Vec<(usize, usize)> {
    let mut faces = Vec::new();
    for (i, e) in mesh.elements.iter().enumerate() {
        for face in 0..4 {
            if let Some(nodes) = e.face_nodes(face) {
                if nodes.iter().all(|&n| on(&mesh.nodes[n])) {
                    faces.push((i, face));
                }
            }
        }
    }
    faces
}

fn convective(faces: &[(usize, usize)], film: f64, sink: f64) -> Vec<Convection> {
    faces
        .iter()
        .map(|&(element, face)| Convection {
            element,
            face,
            film,
            sink,
        })
        .collect()
}

fn solve_with_convection(mesh: &Mesh, bcs: Vec<Bc>, convection: Vec<Convection>) -> Field {
    Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, Vec::new())
        .unwrap()
        .with_convection(convection)
        .unwrap()
        .solve()
        .unwrap()
}

#[test]
fn wall_cooled_by_convection_is_exact() {
    // T = T1 en x = 0, convection h (T − T∞) sur x = L : T = T1 + (T∞ − T1) h x / (k + h L)
    let (t1, h, sink) = (20.0, 30.0, 300.0);
    let power = H * T * K * h * (sink - t1) / (K + h * L);
    for (name, mesh) in all_meshes() {
        let faces = faces_where(&mesh, |p| (p.x - L).abs() < 1e-12);
        assert!(!faces.is_empty(), "{name}");
        let field = solve_with_convection(
            &mesh,
            fix_plane(&mesh, 0.0, t1),
            convective(&faces, h, sink),
        );
        for (n, p) in mesh.nodes.iter().enumerate() {
            let exact = t1 + (sink - t1) * h * p.x / (K + h * L);
            assert!(
                (field.temperatures[n] - exact).abs() < 1e-10,
                "{name} Node {n}"
            );
        }
        // la chaleur venue du fluide sort par la face à T imposée
        let out: f64 = nodes_at(&mesh, 0.0).map(|n| field.reactions[n]).sum();
        assert!(
            (out + power).abs() < 1e-9 * power,
            "{name} : {out} au lieu de {}",
            -power
        );
    }
}

/// Ailette de longueur `length`, encastrée à T0 en x = 0, refroidie par convection sur ses faces
/// latérales (`lateral`), bout isolé. Théorie des ailettes (Biot petit) :
/// θ / θ0 = cosh m (L − x) / cosh m L, m² = h P / (k A). Renvoie l'écart maximal aux Nodes,
/// rapporté à θ0.
fn fin_error(mesh: &Mesh, length: f64, m: f64, h: f64, lateral: impl Fn(&Node) -> bool) -> f64 {
    let (t0, sink) = (100.0, 20.0);
    let faces = faces_where(mesh, lateral);
    let field = solve_with_convection(mesh, fix_plane(mesh, 0.0, t0), convective(&faces, h, sink));
    mesh.nodes
        .iter()
        .zip(&field.temperatures)
        .map(|(p, t)| {
            let exact = sink + (t0 - sink) * (m * (length - p.x)).cosh() / (m * length).cosh();
            (t - exact).abs() / (t0 - sink)
        })
        .fold(0.0, f64::max)
}

#[test]
fn thin_fins_match_fin_theory() {
    // h H / 2k = 10⁻³ : la température est quasi uniforme dans l'épaisseur, la théorie 1D s'applique
    let (h, length) = (0.1, 16.0);

    // 2D : bande length × H d'épaisseur T, convection sur y = 0 et y = H (P = 2 T, A = H T)
    let (nodes, tris) = rect_grid(32, 2, length, H);
    let q = to_quadratic(&nodes, &tris);
    let elements = q
        .elements
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material(), T, PlaneState::Stress).unwrap()))
        .collect();
    let strip = Mesh::new(q.nodes, elements).unwrap();
    let m = (2.0 * h / (K * H)).sqrt();
    let e = fin_error(&strip, length, m, h, |p| {
        p.y.abs() < 1e-12 || (p.y - H).abs() < 1e-12
    });
    assert!(e < 0.002, "Tri6 : écart {e}");

    // 3D : barre length × H × T, convection sur les 4 faces latérales (P = 2 (H + T), A = H T)
    let (nodes, tets) = box_grid([16, 1, 1], [length, H, T]);
    let q = to_quadratic_tet(&nodes, &tets);
    let elements = q
        .elements
        .iter()
        .map(|&t| Element::Tet10(Tet10::new(t, material())))
        .collect();
    let bar = Mesh::new(q.nodes, elements).unwrap();
    let m = (h * 2.0 * (H + T) / (K * H * T)).sqrt();
    let side = |c: f64, l: f64| c.abs() < 1e-12 || (c - l).abs() < 1e-12;
    let e = fin_error(&bar, length, m, h, |p| side(p.y, H) || side(p.z, T));
    assert!(e < 0.002, "Tet10 : écart {e}");
}

#[test]
fn faces_follow_the_abaqus_numbering() {
    let mesh = solid_mesh([1, 1, 1], true);
    let e = &mesh.elements[0];
    let n = e.nodes();
    // face 1 = 0-3-1, puis milieux 0-3 (Node 7), 3-1 (8), 1-0 (4)
    assert_eq!(
        e.face_nodes(1),
        Some(vec![n[0], n[3], n[1], n[7], n[8], n[4]])
    );
    assert_eq!(e.face_nodes(4), None);
    let tri = &plane_mesh(1, 1, true).elements[0];
    let t = tri.nodes();
    assert_eq!(tri.face_nodes(2), Some(vec![t[2], t[0], t[5]]));
    assert_eq!(tri.face_nodes(3), None);
}

#[test]
fn convection_inputs_are_typed_errors() {
    let mesh = plane_mesh(1, 1, false);
    let thermal =
        || Model::with_analysis(mesh.clone(), Analysis::Thermal, Vec::new(), Vec::new()).unwrap();
    let c = |element, face, film| {
        vec![Convection {
            element,
            face,
            film,
            sink: 20.0,
        }]
    };
    assert_eq!(
        thermal().with_convection(c(0, 3, 10.0)),
        Err(MailleError::UnknownFace {
            element: 0,
            face: 3
        })
    );
    assert_eq!(
        thermal().with_convection(c(5, 0, 10.0)),
        Err(MailleError::UnknownElement { element: 5 })
    );
    assert!(matches!(
        thermal().with_convection(c(0, 0, 0.0)),
        Err(MailleError::InvalidProperty { .. })
    ));
}

/// Cylindre creux rᵢ < r < rₑ (quart, symétries adiabatiques naturelles), T(r) exacte imposée sur
/// les deux bords : T = Tᵢ + (Tₑ − Tᵢ) ln(r/rᵢ) / ln(rₑ/rᵢ). Grille polaire `n` × `n`, rayons
/// régulièrement espacés (en progression géométrique, la solution serait linéaire dans l'indice
/// radial et le Tri3 exact aux Nodes : le test ne départagerait rien). Renvoie l'écart maximal aux Nodes rapporté à Tᵢ − Tₑ, et la puissance
/// qui entre par le bord intérieur rapportée à la valeur exacte (π/2) k t (Tᵢ − Tₑ) / ln(rₑ/rᵢ).
fn hollow_cylinder(n: usize, quadratic: bool) -> (f64, f64) {
    let (ri, ro, ti, to) = (10.0, 40.0, 200.0, 20.0);
    let exact = |p: &Node| ti + (to - ti) * (p.x.hypot(p.y) / ri).ln() / (ro / ri).ln();
    let (grid, tris) = rect_grid(n, n, 1.0, 1.0);
    let corners: Vec<Node> = grid
        .iter()
        .map(|g| {
            let (r, t) = (ri + (ro - ri) * g.x, std::f64::consts::FRAC_PI_2 * g.y);
            Node::new_2d(r * t.cos(), r * t.sin())
        })
        .collect();
    let idx = |i: usize, j: usize| i + j * (n + 1);
    // bords r = rᵢ (i = 0) et r = rₑ (i = n) : sommets, et milieux de leurs arêtes en Tri6
    let mut inner: Vec<usize> = (0..=n).map(|j| idx(0, j)).collect();
    let mut outer: Vec<usize> = (0..=n).map(|j| idx(n, j)).collect();
    let (nodes, elements): (Vec<Node>, Vec<Element>) = if quadratic {
        let q = to_quadratic(&corners, &tris);
        inner.extend((0..n).map(|j| q.mid(idx(0, j), idx(0, j + 1))));
        outer.extend((0..n).map(|j| q.mid(idx(n, j), idx(n, j + 1))));
        let elements = q
            .elements
            .iter()
            .map(|&e| Element::Tri6(Tri6::new(e, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        (q.nodes, elements)
    } else {
        let elements = tris
            .iter()
            .map(|&e| Element::Tri3(Tri3::new(e, material(), T, PlaneState::Stress).unwrap()))
            .collect();
        (corners, elements)
    };
    let bcs = inner
        .iter()
        .chain(&outer)
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: exact(&nodes[node]),
        })
        .collect();
    let mesh = Mesh::new(nodes, elements).unwrap();
    let field = solve(&mesh, bcs, Vec::new());
    let error = mesh
        .nodes
        .iter()
        .zip(&field.temperatures)
        .map(|(p, t)| (t - exact(p)).abs())
        .fold(0.0, f64::max);
    let power = std::f64::consts::FRAC_PI_2 * K * T * (ti - to) / (ro / ri).ln();
    let entering: f64 = inner.iter().map(|&n| field.reactions[n]).sum();
    (error / (ti - to), entering / power)
}

/// Même cylindre creux, maillé par le mailleur avec la taille `h` partout. Renvoie les mêmes
/// grandeurs que `hollow_cylinder`.
fn hollow_cylinder_meshed(h: f64, quadratic: bool) -> (f64, f64) {
    let (ri, ro, ti, to) = (10.0, 40.0, 200.0, 20.0);
    let exact = |p: &Node| ti + (to - ti) * (p.x.hypot(p.y) / ri).ln() / (ro / ri).ln();
    let g = Geometry {
        outer: vec![
            Curve::line([ri, 0.0], h, "sym"),
            Curve::arc([ro, 0.0], [0.0, 0.0], false, h, "exterieur"),
            Curve::line([0.0, ro], h, "sym"),
            Curve::arc([0.0, ri], [0.0, 0.0], true, h, "interieur"),
        ],
        holes: vec![],
    };
    let t3 = triangulate(&g, &MeshOptions::default()).unwrap();
    let (nodes, elements, inner, outer): (Vec<Node>, Vec<Element>, Vec<usize>, Vec<usize>) =
        if quadratic {
            let t = t3.to_quadratic();
            let elements = t
                .triangles
                .iter()
                .map(|&e| Element::Tri6(Tri6::new(e, material(), T, PlaneState::Stress).unwrap()))
                .collect();
            let (i, o) = (
                t.boundary("interieur").unwrap(),
                t.boundary("exterieur").unwrap(),
            );
            (t.nodes.clone(), elements, i.nodes.clone(), o.nodes.clone())
        } else {
            let elements = t3
                .triangles
                .iter()
                .map(|&e| Element::Tri3(Tri3::new(e, material(), T, PlaneState::Stress).unwrap()))
                .collect();
            let (i, o) = (
                t3.boundary("interieur").unwrap(),
                t3.boundary("exterieur").unwrap(),
            );
            (t3.nodes.clone(), elements, i.nodes.clone(), o.nodes.clone())
        };
    let bcs = inner
        .iter()
        .chain(&outer)
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: exact(&nodes[node]),
        })
        .collect();
    let mesh = Mesh::new(nodes, elements).unwrap();
    let field = solve(&mesh, bcs, Vec::new());
    let error = mesh
        .nodes
        .iter()
        .zip(&field.temperatures)
        .map(|(p, t)| (t - exact(p)).abs())
        .fold(0.0, f64::max);
    let power = std::f64::consts::FRAC_PI_2 * K * T * (ti - to) / (ro / ri).ln();
    let entering: f64 = inner.iter().map(|&n| field.reactions[n]).sum();
    (error / (ti - to), entering / power)
}

#[test]
fn hollow_cylinder_on_generated_mesh() {
    // Maillage non structuré : l'erreur ne baisse pas régulièrement d'un maillage au suivant, on
    // mesure donc l'ordre sur un facteur 4 en h. Tri3 : rapport > 4 (ordre > 1, ~2 attendu) ;
    // Tri6 : rapport > 16 (ordre > 2).
    for (name, quadratic, min_ratio) in [("Tri3", false, 4.0), ("Tri6", true, 16.0)] {
        let (coarse, _) = hollow_cylinder_meshed(4.0, quadratic);
        let (fine, _) = hollow_cylinder_meshed(1.0, quadratic);
        assert!(
            coarse / fine > min_ratio,
            "{name} : {coarse:.2e} puis {fine:.2e}"
        );
    }
    // Tri6, h = 2 : au moins aussi près que la grille polaire 8 × 8
    let (error, power) = hollow_cylinder_meshed(2.0, true);
    assert!(
        error < 2e-4 && (power - 1.0).abs() < 1e-4,
        "Tri6 : {error:.2e}, puissance {power:.6}"
    );
}

#[test]
fn hollow_cylinder_converges_to_the_logarithmic_profile() {
    let (e8, p8) = hollow_cylinder(8, false);
    let (e16, p16) = hollow_cylinder(16, false);
    let (q8, pq8) = hollow_cylinder(8, true);
    let (q16, _) = hollow_cylinder(16, true);
    // Tri3 : ordre 2 en T et en puissance (erreur divisée par ~4 quand h est divisé par 2)
    assert!(
        (3.5..4.5).contains(&(e8 / e16)),
        "Tri3 : {e8:.2e} puis {e16:.2e}"
    );
    let (d8, d16) = ((p8 - 1.0).abs(), (p16 - 1.0).abs());
    assert!(
        (3.5..4.5).contains(&(d8 / d16)),
        "Tri3 : {p8:.5} puis {p16:.5}"
    );
    // Tri6 : ordre plus élevé, et déjà 1e-4 près en 8 × 8
    assert!(q8 / q16 > 8.0, "Tri6 : {q8:.2e} puis {q16:.2e}");
    assert!(
        q8 < 2e-4 && (pq8 - 1.0).abs() < 1e-4,
        "Tri6 : {q8:.2e}, puissance {pq8:.6}"
    );
}

#[test]
fn any_finite_fluid_temperature_is_accepted() {
    let mesh = plane_mesh(1, 1, false);
    let thermal =
        || Model::with_analysis(mesh.clone(), Analysis::Thermal, Vec::new(), Vec::new()).unwrap();
    let c = |sink| {
        vec![Convection {
            element: 0,
            face: 0,
            film: 10.0,
            sink,
        }]
    };
    assert!(thermal().with_convection(c(-40.0)).is_ok());
    assert!(thermal().with_convection(c(0.0)).is_ok());
    assert!(matches!(
        thermal().with_convection(c(f64::NAN)),
        Err(MailleError::NonFinite { .. })
    ));
}
