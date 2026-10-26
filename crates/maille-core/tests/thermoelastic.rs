//! Couplage thermo-mécanique séquentiel : dilatation α (T − T_ref) d'un champ de températures
//! imposé à une analyse statique, comparée aux solutions analytiques.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Beam, Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions, Triangulation};
use maille_core::model::{Analysis, Bc, Model};
use maille_core::results::{Field, Stress};
use maille_core::MailleError;

const E: f64 = 200e3;
const NU: f64 = 0.3;
const ALPHA: f64 = 1.2e-5;
const L: f64 = 4.0;
const H: f64 = 1.0;
const THICKNESS: f64 = 0.5;
const T_REF: f64 = 20.0;

fn material() -> Material {
    Material::with_poisson(E, NU)
        .and_then(|m| m.with_expansion(ALPHA))
        .unwrap()
}

/// Les Mesh de la validation : triangles dans les deux états plans, tétraèdres. L × H (× THICKNESS
/// en 3D), 4 × 2 (× 1).
fn meshes() -> Vec<(&'static str, Mesh)> {
    let (nodes, tris) = rect_grid(4, 2, L, H);
    let q = to_quadratic(&nodes, &tris);
    let plane = |state: PlaneState, quadratic: bool| {
        let elements = if quadratic {
            q.elements
                .iter()
                .map(|&n| Element::Tri6(Tri6::new(n, material(), THICKNESS, state).unwrap()))
                .collect()
        } else {
            tris.iter()
                .map(|&n| Element::Tri3(Tri3::new(n, material(), THICKNESS, state).unwrap()))
                .collect()
        };
        let nodes = if quadratic {
            q.nodes.clone()
        } else {
            nodes.clone()
        };
        Mesh::new(nodes, elements).unwrap()
    };
    let (box_nodes, tets) = box_grid([4, 2, 1], [L, H, THICKNESS]);
    let qt = to_quadratic_tet(&box_nodes, &tets);
    vec![
        ("Tri3 CP", plane(PlaneState::Stress, false)),
        ("Tri6 CP", plane(PlaneState::Stress, true)),
        ("Tri3 DP", plane(PlaneState::Strain, false)),
        ("Tri6 DP", plane(PlaneState::Strain, true)),
        (
            "Tet4",
            Mesh::new(
                box_nodes.clone(),
                tets.iter()
                    .map(|&t| Element::Tet4(Tet4::new(t, material())))
                    .collect(),
            )
            .unwrap(),
        ),
        (
            "Tet10",
            Mesh::new(
                qt.nodes.clone(),
                qt.elements
                    .iter()
                    .map(|&t| Element::Tet10(Tet10::new(t, material())))
                    .collect(),
            )
            .unwrap(),
        ),
    ]
}

fn solid(mesh: &Mesh) -> bool {
    matches!(mesh.elements[0], Element::Tet4(_) | Element::Tet10(_))
}

fn plane_strain(mesh: &Mesh) -> bool {
    matches!(
        mesh.elements[0],
        Element::Tri3(Tri3 {
            state: PlaneState::Strain,
            ..
        }) | Element::Tri6(Tri6 {
            state: PlaneState::Strain,
            ..
        })
    )
}

/// Dilatation libre ε = a(x) I, avec a(x) = s (α ou α(1 + ν) en déformations planes) fois
/// ΔT(x) = d0 + g·x : u = s (d0 x + (g·x) x − |x|² g / 2). Le terme antisymétrique (g·x) x − x (g·x)
/// du gradient est une rotation : la déformation est bien a(x) I.
fn free_expansion(p: &Node, scale: f64, d0: f64, g: [f64; 3]) -> [f64; 3] {
    let x = [p.x, p.y, p.z];
    let gx: f64 = (0..3).map(|k| g[k] * x[k]).sum();
    let r2: f64 = x.iter().map(|v| v * v).sum();
    std::array::from_fn(|k| scale * (d0 * x[k] + gx * x[k] - r2 * g[k] / 2.0))
}

/// BC qui ne bloquent que les mouvements de corps rigide, avec la valeur du champ exact `u` :
/// (0, 0, 0) bloqué, (L, 0, 0) en y (et z), en 3D (0, H, 0) en z.
fn rigid_body_bcs(mesh: &Mesh, u: impl Fn(&Node) -> [f64; 3]) -> Vec<Bc> {
    let find = |x: f64, y: f64| {
        mesh.nodes
            .iter()
            .position(|p| (p.x - x).abs() < 1e-12 && (p.y - y).abs() < 1e-12 && p.z.abs() < 1e-12)
            .unwrap()
    };
    let (o, a, b) = (find(0.0, 0.0), find(L, 0.0), find(0.0, H));
    let bc = |node: usize, dof: Dof, k: usize| Bc {
        node,
        dof,
        value: u(&mesh.nodes[node])[k],
    };
    let mut bcs = vec![bc(o, Dof::Ux, 0), bc(o, Dof::Uy, 1), bc(a, Dof::Uy, 1)];
    if solid(mesh) {
        bcs.extend([bc(o, Dof::Uz, 2), bc(a, Dof::Uz, 2), bc(b, Dof::Uz, 2)]);
    }
    bcs
}

fn heated(mesh: &Mesh, bcs: Vec<Bc>, t: impl Fn(&Node) -> f64) -> Field {
    let temperatures = mesh.nodes.iter().map(t).collect();
    Model::new(mesh.clone(), bcs, Vec::new())
        .and_then(|m| m.with_temperatures(temperatures, T_REF))
        .and_then(|m| m.solve())
        .unwrap()
}

fn stresses(field: &Field) -> impl Iterator<Item = &Stress> {
    field.stresses.iter().flatten()
}

/// Composantes de contrainte autres que σzz.
fn others_than_zz(s: &Stress) -> [f64; 5] {
    [s.sigma_xx, s.sigma_yy, s.tau_xy, s.tau_yz, s.tau_xz]
}

#[test]
fn free_uniform_heating_expands_without_stress() {
    let dt = 80.0;
    for (name, mesh) in meshes() {
        let scale = if plane_strain(&mesh) {
            ALPHA * (1.0 + NU)
        } else {
            ALPHA
        };
        let exact = |p: &Node| free_expansion(p, scale, dt, [0.0; 3]);
        let field = heated(&mesh, rigid_body_bcs(&mesh, exact), |_| T_REF + dt);
        let d = if solid(&mesh) { 3 } else { 2 };
        for (n, p) in mesh.nodes.iter().enumerate() {
            for k in 0..d {
                let got = field.displacements[d * n + k];
                assert!(
                    (got - exact(p)[k]).abs() < 1e-12 * L,
                    "{name} Node {n} : {got} au lieu de {}",
                    exact(p)[k]
                );
            }
        }
        // déformations planes : εzz bloquée, σzz = −E α ΔT ; ailleurs, aucune contrainte
        let zz = if plane_strain(&mesh) {
            -E * ALPHA * dt
        } else {
            0.0
        };
        for s in stresses(&field) {
            assert!(
                others_than_zz(s).iter().all(|v| v.abs() < 1e-8),
                "{name} {s:?}"
            );
            assert!(
                (s.sigma_zz - zz).abs() < 1e-8,
                "{name} σzz = {}",
                s.sigma_zz
            );
        }
        assert!(field.reactions.iter().all(|r| r.abs() < 1e-8), "{name}");
    }
}

#[test]
fn clamped_uniform_heating_gives_the_blocked_stress() {
    // tous les Nodes bloqués : contraintes planes σxx = σyy = −E α ΔT / (1 − ν) ; déformations
    // planes et 3D, état hydrostatique −E α ΔT / (1 − 2ν)
    let dt = 80.0;
    for (name, mesh) in meshes() {
        let d = if solid(&mesh) { 3 } else { 2 };
        let dofs = [Dof::Ux, Dof::Uy, Dof::Uz];
        let bcs = (0..mesh.nodes.len())
            .flat_map(|node| {
                dofs[..d].iter().map(move |&dof| Bc {
                    node,
                    dof,
                    value: 0.0,
                })
            })
            .collect();
        let field = heated(&mesh, bcs, |_| T_REF + dt);
        let plane_stress = !solid(&mesh) && !plane_strain(&mesh);
        let (normal, zz) = if plane_stress {
            (-E * ALPHA * dt / (1.0 - NU), 0.0)
        } else {
            let p = -E * ALPHA * dt / (1.0 - 2.0 * NU);
            (p, p)
        };
        for s in stresses(&field) {
            assert!(
                (s.sigma_xx - normal).abs() < 1e-9 * normal.abs()
                    && (s.sigma_yy - normal).abs() < 1e-9 * normal.abs()
                    && (s.sigma_zz - zz).abs() < 1e-9 * normal.abs(),
                "{name} {s:?}"
            );
            assert!(s.tau_xy.abs() < 1e-9 * normal.abs(), "{name}");
        }
        // les réactions équilibrent la charge thermique : leur somme est nulle dans chaque direction
        for k in 0..d {
            let sum: f64 = field.reactions.iter().skip(k).step_by(d).sum();
            assert!(sum.abs() < 1e-8 * normal.abs(), "{name} somme {sum}");
        }
    }
}

#[test]
fn free_linear_temperature_is_stress_free_with_quadratic_elements() {
    // ΔT linéaire : déformation compatible, aucune contrainte. Le déplacement est quadratique :
    // exact en Tri6 et Tet10 (pas en Tri3 et Tet4, comme en flexion)
    let (d0, g) = (30.0, [15.0, -40.0, 25.0]);
    for (name, mesh) in meshes()
        .into_iter()
        .filter(|(name, _)| name.starts_with("Tri6") || *name == "Tet10")
    {
        let g = if solid(&mesh) { g } else { [g[0], g[1], 0.0] };
        let scale = if plane_strain(&mesh) {
            ALPHA * (1.0 + NU)
        } else {
            ALPHA
        };
        let exact = |p: &Node| free_expansion(p, scale, d0, g);
        let temperature = |p: &Node| T_REF + d0 + g[0] * p.x + g[1] * p.y + g[2] * p.z;
        let field = heated(&mesh, rigid_body_bcs(&mesh, exact), temperature);
        let d = if solid(&mesh) { 3 } else { 2 };
        for (n, p) in mesh.nodes.iter().enumerate() {
            for k in 0..d {
                let got = field.displacements[d * n + k];
                assert!((got - exact(p)[k]).abs() < 1e-12 * L, "{name} Node {n}");
            }
        }
        for s in stresses(&field) {
            assert!(
                others_than_zz(s).iter().all(|v| v.abs() < 1e-7),
                "{name} {s:?}"
            );
            let zz = if plane_strain(&mesh) {
                -E * ALPHA * (temperature(&Node::new(s.x, s.y, s.z)) - T_REF)
            } else {
                0.0
            };
            assert!(
                (s.sigma_zz - zz).abs() < 1e-7,
                "{name} σzz = {}",
                s.sigma_zz
            );
        }
    }
}

fn bar_mesh(n: usize, reversed: bool) -> Mesh {
    let material = Material::new(E)
        .and_then(|m| m.with_expansion(ALPHA))
        .unwrap();
    let nodes = (0..=n)
        .map(|i| Node::new_2d(L * i as f64 / n as f64, 0.0))
        .collect();
    let elements = (0..n)
        .map(|i| {
            let nodes = if reversed { [i + 1, i] } else { [i, i + 1] };
            Element::Bar(Bar {
                nodes,
                material,
                section: Section::new(2.0).unwrap(),
            })
        })
        .collect();
    Mesh::new(nodes, elements).unwrap()
}

#[test]
fn heated_bars_expand_freely_or_push_on_their_supports() {
    let (dt, area) = (80.0, 2.0);
    let fix = |node: usize| Bc {
        node,
        dof: Dof::U,
        value: 0.0,
    };
    for reversed in [false, true] {
        // libre : u = α ΔT x
        let mesh = bar_mesh(4, reversed);
        let field = heated(&mesh, vec![fix(0)], |_| T_REF + dt);
        for (n, p) in mesh.nodes.iter().enumerate() {
            let exact = ALPHA * dt * p.x;
            assert!((field.displacements[n] - exact).abs() < 1e-14, "Node {n}");
        }
        // T linéaire de T_REF à T_REF + 2 ΔT : u(L) = α ΔT L (allongement moyen)
        let field = heated(&mesh, vec![fix(0)], |p| T_REF + 2.0 * dt * p.x / L);
        assert!((field.displacements[4] - ALPHA * dt * L).abs() < 1e-14);
        // bridée aux deux bouts : les appuis poussent la barre vers l'intérieur, ±E A α ΔT
        let field = heated(&mesh, vec![fix(0), fix(4)], |_| T_REF + dt);
        let push = E * area * ALPHA * dt;
        assert!(
            (field.reactions[0] - push).abs() < 1e-9 * push,
            "réversée {reversed}"
        );
        assert!(
            (field.reactions[4] + push).abs() < 1e-9 * push,
            "réversée {reversed}"
        );
    }
}

/// Cylindre creux rᵢ < r < rₑ, températures imposées Tᵢ et Tₑ (profil en ln r), déformations planes
/// (εz = 0), surfaces libres : Timoshenko et Goodier, avec I(r) = ∫ᵣᵢʳ ΔT r dr,
/// σr = αE / ((1 − ν) r²) [(r² − rᵢ²)/(rₑ² − rᵢ²) I(rₑ) − I(r)],
/// σθ = αE / ((1 − ν) r²) [(r² + rᵢ²)/(rₑ² − rᵢ²) I(rₑ) + I(r) − ΔT r²],
/// σz = ν (σr + σθ) − E α ΔT.
mod cylinder {
    use super::*;

    pub const INNER: f64 = 10.0;
    pub const OUTER: f64 = 20.0;
    pub const T_INNER: f64 = 120.0;
    pub const T_OUTER: f64 = T_REF;
    const K: f64 = 0.05;

    fn beta() -> f64 {
        (T_OUTER - T_INNER) / (OUTER / INNER).ln()
    }

    /// ΔT(r) = T(r) − T_REF.
    pub fn delta(r: f64) -> f64 {
        T_INNER - T_REF + beta() * (r / INNER).ln()
    }

    /// I(r) = ∫ᵣᵢʳ ΔT r dr, avec ∫ ln(r/rᵢ) r dr = r²/2 ln(r/rᵢ) − r²/4.
    fn integral(r: f64) -> f64 {
        (T_INNER - T_REF) * (r * r - INNER * INNER) / 2.0
            + beta() * (r * r / 2.0 * (r / INNER).ln() - r * r / 4.0 + INNER * INNER / 4.0)
    }

    /// (σr, σθ, σz) exactes en r.
    pub fn exact(r: f64) -> [f64; 3] {
        let c = ALPHA * E / ((1.0 - NU) * r * r);
        let span = OUTER * OUTER - INNER * INNER;
        let sr = c * ((r * r - INNER * INNER) / span * integral(OUTER) - integral(r));
        let st =
            c * ((r * r + INNER * INNER) / span * integral(OUTER) + integral(r) - delta(r) * r * r);
        [sr, st, NU * (sr + st) - E * ALPHA * delta(r)]
    }

    /// Quart de section maillé (taille h), calcul thermique puis statique sur le même Mesh.
    /// Renvoie le Mesh et le Field mécanique.
    pub fn solve(h: f64, quadratic: bool) -> (Mesh, Field) {
        let geometry = Geometry {
            outer: vec![
                Curve::line([INNER, 0.0], h, "sym_x"),
                Curve::arc([OUTER, 0.0], [0.0, 0.0], false, h, "exterieur"),
                Curve::line([0.0, OUTER], h, "sym_y"),
                Curve::arc([0.0, INNER], [0.0, 0.0], true, h, "interieur"),
            ],
            holes: vec![],
        };
        let t = triangulate(&geometry, &MeshOptions::default()).unwrap();
        let material = material().with_conductivity(K).unwrap();
        let (mesh, [sym_x, outer, sym_y, inner]) = if quadratic {
            build(t.to_quadratic(), |n| {
                Element::Tri6(Tri6::new(n, material, 1.0, PlaneState::Strain).unwrap())
            })
        } else {
            build(t, |n| {
                Element::Tri3(Tri3::new(n, material, 1.0, PlaneState::Strain).unwrap())
            })
        };

        // thermique : solution exacte imposée sur les deux rayons (Nodes pas tous sur le cercle)
        let bcs = inner
            .iter()
            .chain(&outer)
            .map(|&node| Bc {
                node,
                dof: Dof::T,
                value: T_REF + delta(mesh.nodes[node].x.hypot(mesh.nodes[node].y)),
            })
            .collect();
        let thermal = Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, Vec::new())
            .unwrap()
            .solve()
            .unwrap();

        // mécanique : symétries, surfaces libres, températures du calcul thermique
        let roller = |nodes: &[usize], dof: Dof| {
            nodes
                .iter()
                .map(|&node| Bc {
                    node,
                    dof,
                    value: 0.0,
                })
                .collect::<Vec<_>>()
        };
        let mut bcs = roller(&sym_x, Dof::Uy);
        bcs.extend(roller(&sym_y, Dof::Ux));
        let field = Model::new(mesh.clone(), bcs, Vec::new())
            .and_then(|m| m.with_temperatures(thermal.temperatures, T_REF))
            .and_then(|m| m.solve())
            .unwrap();
        (mesh, field)
    }

    /// Mesh des triangles de `t` et Nodes des bords sym_x, exterieur, sym_y, interieur.
    fn build<const N: usize>(
        t: Triangulation<N>,
        element: impl Fn([usize; N]) -> Element,
    ) -> (Mesh, [Vec<usize>; 4]) {
        let sets = ["sym_x", "exterieur", "sym_y", "interieur"]
            .map(|name| t.boundary(name).unwrap().nodes.clone());
        let elements = t.triangles.iter().map(|&n| element(n)).collect();
        (Mesh::new(t.nodes, elements).unwrap(), sets)
    }

    /// Erreur RMS sur (σr, σθ, σz) aux points de Gauss, rapportée au plus grand |σθ| exact.
    pub fn error(field: &Field) -> f64 {
        let scale = exact(INNER)[1].abs().max(exact(OUTER)[1].abs());
        let (mut sum, mut count) = (0.0, 0);
        for s in field.stresses.iter().flatten() {
            let r = s.x.hypot(s.y);
            let (c, sn) = (s.x / r, s.y / r);
            let sr = s.sigma_xx * c * c + s.sigma_yy * sn * sn + 2.0 * s.tau_xy * c * sn;
            let st = s.sigma_xx * sn * sn + s.sigma_yy * c * c - 2.0 * s.tau_xy * c * sn;
            let e = exact(r);
            sum += (sr - e[0]).powi(2) + (st - e[1]).powi(2) + (s.sigma_zz - e[2]).powi(2);
            count += 3;
        }
        (sum / count as f64).sqrt() / scale
    }
}

#[test]
fn thick_cylinder_thermal_stresses_follow_timoshenko() {
    // surfaces libres : σr nulle en rᵢ et rₑ
    assert!(cylinder::exact(cylinder::INNER)[0].abs() < 1e-9);
    assert!(cylinder::exact(cylinder::OUTER)[0].abs() < 1e-9);
    let (_, tri6) = cylinder::solve(1.0, true);
    let e6 = cylinder::error(&tri6);
    assert!(e6 < 2e-3, "Tri6 h = 1 : {e6:.2e}");
    // Tri3 : contrainte constante par Element, convergence d'ordre 1 en h
    let (e_coarse, e_fine) = (
        cylinder::error(&cylinder::solve(2.0, false).1),
        cylinder::error(&cylinder::solve(1.0, false).1),
    );
    let ratio = e_coarse / e_fine;
    assert!(
        (1.6..2.5).contains(&ratio),
        "Tri3 : {e_coarse:.2e} puis {e_fine:.2e}"
    );
}

#[test]
fn temperature_inputs_are_typed_errors() {
    let mesh = meshes().remove(0).1;
    let n = mesh.nodes.len();
    let model = || Model::new(mesh.clone(), Vec::new(), Vec::new()).unwrap();
    assert_eq!(
        model().with_temperatures(vec![0.0; 2], T_REF),
        Err(MailleError::TemperatureCount {
            expected: n,
            got: 2
        })
    );
    let mut nan = vec![0.0; n];
    nan[3] = f64::NAN;
    assert_eq!(
        model().with_temperatures(nan, T_REF),
        Err(MailleError::NonFiniteTemperature { node: 3 })
    );
    assert!(matches!(
        model().with_temperatures(vec![0.0; n], f64::INFINITY),
        Err(MailleError::NonFinite { .. })
    ));

    let thermal = Material::conductive(1.0).unwrap();
    let tri = Element::Tri3(Tri3::new([0, 1, 2], thermal, 1.0, PlaneState::Stress).unwrap());
    let three = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(1.0, 0.0),
        Node::new_2d(0.0, 1.0),
    ];
    let thermal_model = Model::with_analysis(
        Mesh::new(three.clone(), vec![tri]).unwrap(),
        Analysis::Thermal,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(
        thermal_model.with_temperatures(vec![0.0; 3], T_REF),
        Err(MailleError::StaticOnly)
    );

    // E et ν sans α
    let elastic = Material::with_poisson(E, NU).unwrap();
    let tri = Element::Tri3(Tri3::new([0, 1, 2], elastic, 1.0, PlaneState::Stress).unwrap());
    let no_alpha =
        Model::new(Mesh::new(three, vec![tri]).unwrap(), Vec::new(), Vec::new()).unwrap();
    assert_eq!(
        no_alpha.with_temperatures(vec![0.0; 3], T_REF),
        Err(MailleError::MissingExpansion { element: 0 })
    );

    // Beam : pas de dilatation
    let m = Material::new(E)
        .and_then(|m| m.with_expansion(ALPHA))
        .unwrap();
    let beam = Element::Beam(Beam::new([0, 1], m, 1.0).unwrap());
    let beams = Model::new(
        Mesh::new(
            vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
            vec![beam],
        )
        .unwrap(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(
        beams.with_temperatures(vec![0.0; 2], T_REF),
        Err(MailleError::UnsupportedThermalStrain { element: 0 })
    );
}
