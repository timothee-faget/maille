//! Élasto-plasticité de von Mises : Step quasi-statique, Amplitude, Newton-Raphson, comparés à la
//! loi au point de Gauss et aux solutions analytiques de la traction.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{Amplitude, Analysis, Bc, Load, Model, Step};
use maille_core::plasticity::{PlasticState, VonMises};
use maille_core::results::Increment;
use maille_core::MailleError;

const E: f64 = 210e3;
const NU: f64 = 0.3;
const SY: f64 = 250.0;
const H: f64 = 2000.0;
const L: f64 = 2.0;
const W: f64 = 1.0;
const D: f64 = 1.0;

fn plastic(hardening: f64) -> Material {
    Material::with_poisson(E, NU)
        .and_then(|m| m.with_plasticity(SY, hardening))
        .unwrap()
}

fn plane_meshes(material: Material) -> Vec<(&'static str, Mesh)> {
    let (nodes, tris) = rect_grid(4, 2, L, W);
    let q = to_quadratic(&nodes, &tris);
    let strain = PlaneState::Strain;
    let tri3 = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, material, 1.0, strain).unwrap()))
        .collect();
    let tri6 = q
        .elements
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material, 1.0, strain).unwrap()))
        .collect();
    vec![
        ("Tri3", Mesh::new(nodes, tri3).unwrap()),
        ("Tri6", Mesh::new(q.nodes, tri6).unwrap()),
    ]
}

fn solid_meshes(material: Material) -> Vec<(&'static str, Mesh)> {
    let (nodes, tets) = box_grid([3, 2, 2], [L, W, D]);
    let q = to_quadratic_tet(&nodes, &tets);
    let tet4 = tets
        .iter()
        .map(|&n| Element::Tet4(Tet4::new(n, material)))
        .collect();
    let tet10 = q
        .elements
        .iter()
        .map(|&n| Element::Tet10(Tet10::new(n, material)))
        .collect();
    vec![
        ("Tet4", Mesh::new(nodes, tet4).unwrap()),
        ("Tet10", Mesh::new(q.nodes, tet10).unwrap()),
    ]
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9 * L
}

fn dofs(mesh: &Mesh) -> &'static [Dof] {
    if mesh.dofs_per_node() == 2 {
        &[Dof::Ux, Dof::Uy]
    } else {
        &[Dof::Ux, Dof::Uy, Dof::Uz]
    }
}

/// Charge jusqu'à t = 1, puis décharge jusqu'à λ = `last` à t = 2, en 8 incréments.
fn load_unload(model: &Model, last: f64) -> Vec<Increment> {
    let step = Step::quasi_static(2.0, 8).unwrap();
    let amplitude = Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (2.0, last)]).unwrap();
    model.solve_quasi_static(&step, &amplitude).unwrap()
}

#[test]
fn amplitude_is_piecewise_linear_then_constant() {
    let a = Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (3.0, -1.0)]).unwrap();
    for (t, expected) in [
        (0.0, 0.0),
        (0.25, 0.25),
        (1.0, 1.0),
        (2.0, 0.0),
        (2.5, -0.5),
        (3.0, -1.0),
        (10.0, -1.0),
    ] {
        assert!((a.at(t) - expected).abs() < 1e-15, "λ({t}) = {}", a.at(t));
    }
    let ramp = Amplitude::ramp(4.0).unwrap();
    assert!((ramp.at(1.0) - 0.25).abs() < 1e-15 && (ramp.at(5.0) - 1.0).abs() < 1e-15);
}

#[test]
fn amplitude_inputs_are_typed_errors() {
    for (points, bad) in [
        (vec![], 0),
        (vec![(0.0, 0.5), (1.0, 1.0)], 0),
        (vec![(0.1, 0.0), (1.0, 1.0)], 0),
        (vec![(0.0, 0.0), (1.0, 1.0), (1.0, 2.0)], 2),
        (vec![(0.0, 0.0), (1.0, f64::NAN)], 1),
        (vec![(0.0, 0.0), (f64::INFINITY, 1.0)], 1),
    ] {
        assert_eq!(
            Amplitude::new(points),
            Err(MailleError::InvalidAmplitude { point: bad })
        );
    }
    assert!(matches!(
        Amplitude::ramp(0.0),
        Err(MailleError::InvalidProperty { .. })
    ));
}

/// Traction selon x d'un pavé L × W × D : symétries sur x = 0, y = 0, z = 0 et déplacement δ imposé
/// sur x = L. État uniaxial exact, que les tétraèdres représentent.
fn tension_bcs(mesh: &Mesh, delta: f64) -> Vec<Bc> {
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        for (c, dof) in [(n.x, Dof::Ux), (n.y, Dof::Uy), (n.z, Dof::Uz)] {
            if near(c, 0.0) {
                bcs.push(Bc {
                    node: i,
                    dof,
                    value: 0.0,
                });
            }
        }
        if near(n.x, L) {
            bcs.push(Bc {
                node: i,
                dof: Dof::Ux,
                value: delta,
            });
        }
    }
    bcs
}

/// Contrainte de traction uniaxiale σ(ε) en charge monotone : élastique, puis σy + Eₜ (ε − σy/E)
/// avec Eₜ = E H / (E + H) ; et la déformation plastique p = ε − σ/E.
fn uniaxial(strain: f64, hardening: f64) -> (f64, f64) {
    let sigma = if strain <= SY / E {
        E * strain
    } else {
        SY + E * hardening / (E + hardening) * (strain - SY / E)
    };
    (sigma, strain - sigma / E)
}

#[test]
fn uniaxial_tension_loads_and_unloads_like_the_closed_form() {
    let strain_max = 3.0 * SY / E;
    for hardening in [H, 0.0] {
        let (sigma_max, p_max) = uniaxial(strain_max, hardening);
        // décharge jusqu'à force nulle : il reste la déformation plastique
        let last = p_max / strain_max;
        for (name, mesh) in solid_meshes(plastic(hardening)) {
            let bcs = tension_bcs(&mesh, strain_max * L);
            let model = Model::new(mesh.clone(), bcs, Vec::new()).unwrap();
            let history = load_unload(&model, last);
            assert_eq!(history.len(), 9);
            let corner = mesh
                .nodes
                .iter()
                .position(|n| near(n.x, L) && near(n.y, W) && near(n.z, D))
                .unwrap();
            for inc in &history {
                let t = inc.time;
                let lambda = if t <= 1.0 {
                    t
                } else {
                    1.0 + (last - 1.0) * (t - 1.0)
                };
                let strain = lambda * strain_max;
                let (sigma, p) = if t <= 1.0 {
                    uniaxial(strain, hardening)
                } else {
                    (sigma_max - E * (strain_max - strain), p_max)
                };
                let label = format!("{name}, H = {hardening}, t = {t}");
                let force: f64 = mesh
                    .nodes
                    .iter()
                    .enumerate()
                    .filter(|(_, n)| near(n.x, L))
                    .map(|(i, _)| inc.field.reactions[3 * i])
                    .sum();
                assert!(
                    (force - sigma * W * D).abs() < 1e-6 * SY,
                    "{label} : F = {force}"
                );
                for (points, ps) in inc
                    .field
                    .stresses
                    .iter()
                    .zip(&inc.field.equivalent_plastic_strains)
                {
                    for (s, got_p) in points.iter().zip(ps) {
                        assert!((s.sigma_xx - sigma).abs() < 1e-6 * SY, "{label} : {s:?}");
                        for other in [s.sigma_yy, s.sigma_zz, s.tau_xy, s.tau_yz, s.tau_xz] {
                            assert!(other.abs() < 1e-6 * SY, "{label} : {s:?}");
                        }
                        assert!(
                            (got_p - p).abs() < 1e-9,
                            "{label} : p = {got_p} au lieu de {p}"
                        );
                    }
                }
                // contraction latérale : élastique −ν σ/E, plastique −p/2 (volume conservé)
                let uy = inc.field.displacements[3 * corner + 1];
                let expected = -(NU * sigma / E + p / 2.0) * W;
                assert!(
                    (uy - expected).abs() < 1e-6 * strain_max * L,
                    "{label} : uy = {uy}"
                );
            }
            // la fin de la décharge ramène la force à zéro
            let last_force: f64 = history[8].field.reactions.iter().step_by(3).sum();
            assert!(last_force.abs() < 1e-6 * SY, "{name} : {last_force}");
        }
    }
}

/// Déplacement u = ε x d'une déformation uniforme (Voigt de `solid_d`, glissements d'ingénieur).
fn uniform(strain: &[f64; 6], n: &Node) -> [f64; 3] {
    [
        strain[0] * n.x + strain[5] / 2.0 * n.y + strain[4] / 2.0 * n.z,
        strain[5] / 2.0 * n.x + strain[1] * n.y + strain[3] / 2.0 * n.z,
        strain[4] / 2.0 * n.x + strain[3] / 2.0 * n.y + strain[2] * n.z,
    ]
}

#[test]
fn homogeneous_strain_history_matches_the_material_point() {
    // charge, puis décharge au-delà de zéro : plastification inverse
    let plane = [4e-3, -1e-3, 0.0, 0.0, 0.0, 2e-3];
    let solid = [4e-3, -1e-3, -1.5e-3, 1e-3, -0.5e-3, 2e-3];
    let law = VonMises {
        young: E,
        poisson: NU,
        yield_stress: SY,
        hardening: H,
    };
    let cases = plane_meshes(plastic(H))
        .into_iter()
        .map(|(name, mesh)| (name, mesh, plane))
        .chain(
            solid_meshes(plastic(H))
                .into_iter()
                .map(|(name, mesh)| (name, mesh, solid)),
        );
    for (name, mesh, strain) in cases {
        // tout le bord suit le champ uniforme ; les Nodes intérieurs sont libres
        let on_boundary = |n: &Node| {
            near(n.x, 0.0)
                || near(n.x, L)
                || near(n.y, 0.0)
                || near(n.y, W)
                || (mesh.dofs_per_node() == 3 && (near(n.z, 0.0) || near(n.z, D)))
        };
        let mut bcs = Vec::new();
        for (i, n) in mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| on_boundary(n))
        {
            for (&dof, value) in dofs(&mesh).iter().zip(uniform(&strain, n)) {
                bcs.push(Bc {
                    node: i,
                    dof,
                    value,
                });
            }
        }
        let model = Model::new(mesh.clone(), bcs, Vec::new()).unwrap();
        let history = load_unload(&model, -0.5);

        let mut state = PlasticState::default();
        for inc in &history[1..] {
            let t = inc.time;
            let lambda = if t <= 1.0 { t } else { 1.0 - 1.5 * (t - 1.0) };
            let point = law.update(&strain.map(|e| lambda * e), &state);
            state = point.state;
            let s = point.stress;
            for (points, ps) in inc
                .field
                .stresses
                .iter()
                .zip(&inc.field.equivalent_plastic_strains)
            {
                for (got, got_p) in points.iter().zip(ps) {
                    let got_s = [
                        got.sigma_xx,
                        got.sigma_yy,
                        got.sigma_zz,
                        got.tau_yz,
                        got.tau_xz,
                        got.tau_xy,
                    ];
                    for (a, b) in got_s.iter().zip(s) {
                        assert!(
                            (a - b).abs() < 1e-6 * SY,
                            "{name}, t = {t} : {got_s:?} au lieu de {s:?}"
                        );
                    }
                    assert!((got_p - state.p).abs() < 1e-9, "{name}, t = {t}");
                }
            }
        }
        assert!(state.p > 0.0, "{name} : le chemin doit plastifier");
    }
}

#[test]
fn an_elastic_model_gives_the_linear_solution_in_one_iteration() {
    let elastic = Material::with_poisson(E, NU).unwrap();
    for (name, mesh) in solid_meshes(elastic) {
        // traction imposée en force sur x = L, pour que F intervienne aussi
        let bcs: Vec<Bc> = tension_bcs(&mesh, 0.0)
            .into_iter()
            .filter(|bc| !near(mesh.nodes[bc.node].x, L))
            .collect();
        let loads: Vec<Load> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| near(n.x, L) && near(n.y, W))
            .map(|(i, _)| Load {
                node: i,
                dof: Dof::Ux,
                value: 10.0,
            })
            .collect();
        let model = Model::new(mesh, bcs, loads).unwrap();
        let linear = model.solve().unwrap();
        let step = Step::quasi_static(1.0, 2).unwrap();
        let history = model
            .solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap())
            .unwrap();
        for inc in &history[1..] {
            assert_eq!(inc.residuals.len(), 2, "{name} : {:?}", inc.residuals);
            for (u, v) in inc.field.displacements.iter().zip(&linear.displacements) {
                assert!((u - inc.time * v).abs() < 1e-9 * 1e-3, "{name}");
            }
        }
        let last = &history[2].field;
        for (a, b) in last.reactions.iter().zip(&linear.reactions) {
            assert!((a - b).abs() < 1e-6, "{name}");
        }
        for (pa, pb) in last.stresses.iter().zip(&linear.stresses) {
            for (a, b) in pa.iter().zip(pb) {
                assert!((a.von_mises() - b.von_mises()).abs() < 1e-6, "{name}");
            }
        }
        assert!(last
            .equivalent_plastic_strains
            .iter()
            .flatten()
            .all(|&p| p == 0.0));
    }
}

#[test]
fn newton_converges_quadratically() {
    // bande encastrée en x = 0, fléchie par un déplacement imposé en x = L : zone plastique non
    // uniforme qui grandit à chaque incrément
    for (name, mesh) in plane_meshes(plastic(H)) {
        let mut bcs = Vec::new();
        for (i, n) in mesh.nodes.iter().enumerate() {
            if near(n.x, 0.0) {
                bcs.push(Bc {
                    node: i,
                    dof: Dof::Ux,
                    value: 0.0,
                });
                bcs.push(Bc {
                    node: i,
                    dof: Dof::Uy,
                    value: 0.0,
                });
            } else if near(n.x, L) {
                bcs.push(Bc {
                    node: i,
                    dof: Dof::Uy,
                    value: 0.02,
                });
            }
        }
        let model = Model::new(mesh, bcs, Vec::new()).unwrap();
        let step = Step::quasi_static(1.0, 5).unwrap();
        let history = model
            .solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap())
            .unwrap();
        let mut longest = 0;
        for inc in &history[1..] {
            let r = &inc.residuals;
            longest = longest.max(r.len());
            assert!(r.len() <= 8, "{name}, t = {} : {r:?}", inc.time);
            // r[0] : déséquilibre avant d'appliquer l'incrément imposé (nul ici) ; ensuite, une fois
            // proche de la solution, l'erreur est élevée au carré à chaque itération
            for pair in r[1..].windows(2) {
                if pair[0] < 1e-2 {
                    assert!(
                        pair[1] <= (10.0 * pair[0] * pair[0]).max(1e-12),
                        "{name}, t = {} : {r:?}",
                        inc.time
                    );
                }
            }
        }
        assert!(
            longest >= 4,
            "{name} : le test doit exercer plusieurs itérations"
        );
        let last = &history[5].field;
        assert!(last
            .equivalent_plastic_strains
            .iter()
            .flatten()
            .any(|&p| p > 0.0));
    }
}

#[test]
fn a_force_beyond_the_limit_load_is_a_typed_error() {
    // traction en déformations planes, parfaitement plastique : σxx limite = 2 σy / √3
    let (nodes, tris) = rect_grid(4, 2, L, W);
    let elements = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, plastic(0.0), 1.0, PlaneState::Strain).unwrap()))
        .collect();
    let mesh = Mesh::new(nodes, elements).unwrap();
    let mut bcs = Vec::new();
    let mut loads = Vec::new();
    let total = 1.2 * 2.0 * SY / 3.0_f64.sqrt() * W;
    for (i, n) in mesh.nodes.iter().enumerate() {
        if near(n.x, 0.0) {
            bcs.push(Bc {
                node: i,
                dof: Dof::Ux,
                value: 0.0,
            });
        }
        if near(n.x, 0.0) && near(n.y, 0.0) {
            bcs.push(Bc {
                node: i,
                dof: Dof::Uy,
                value: 0.0,
            });
        }
        if near(n.x, L) {
            // 3 Nodes sur le bord : poids ¼, ½, ¼ d'une traction uniforme
            let share = if near(n.y, 0.0) || near(n.y, W) {
                0.25
            } else {
                0.5
            };
            loads.push(Load {
                node: i,
                dof: Dof::Ux,
                value: share * total,
            });
        }
    }
    let model = Model::new(mesh, bcs, loads).unwrap();
    let step = Step::quasi_static(1.0, 6).unwrap();
    let result = model.solve_quasi_static(&step, &Amplitude::ramp(1.0).unwrap());
    assert!(
        matches!(
            result,
            Err(MailleError::PlasticCollapse { .. } | MailleError::NotConverged { .. })
        ),
        "{result:?}"
    );
}

#[test]
fn plasticity_inputs_are_typed_errors() {
    let tri = || {
        vec![
            Node::new_2d(0.0, 0.0),
            Node::new_2d(1.0, 0.0),
            Node::new_2d(0.0, 1.0),
        ]
    };
    let mesh = |material: Material, state: PlaneState| {
        let t = Tri3::new([0, 1, 2], material, 1.0, state).unwrap();
        Mesh::new(tri(), vec![Element::Tri3(t)]).unwrap()
    };
    let bcs = || {
        vec![
            Bc {
                node: 0,
                dof: Dof::Ux,
                value: 0.0,
            },
            Bc {
                node: 0,
                dof: Dof::Uy,
                value: 0.0,
            },
            Bc {
                node: 1,
                dof: Dof::Uy,
                value: 0.0,
            },
        ]
    };
    let step = Step::quasi_static(1.0, 1).unwrap();
    let ramp = Amplitude::ramp(1.0).unwrap();

    // contraintes planes : pas de version plastique
    assert_eq!(
        Model::new(mesh(plastic(H), PlaneState::Stress), bcs(), Vec::new()),
        Err(MailleError::UnsupportedPlasticElement { element: 0 })
    );
    let elastic = Material::with_poisson(E, NU)
        .and_then(|m| m.with_expansion(1e-5))
        .unwrap();
    let plane_stress = Model::new(mesh(elastic, PlaneState::Stress), bcs(), Vec::new()).unwrap();
    assert_eq!(
        plane_stress.solve_quasi_static(&step, &ramp),
        Err(MailleError::UnsupportedPlasticElement { element: 0 })
    );

    let model = Model::new(mesh(plastic(H), PlaneState::Strain), bcs(), Vec::new()).unwrap();
    assert_eq!(
        model.solve(),
        Err(MailleError::PlasticityNeedsQuasiStatic { element: 0 })
    );
    assert_eq!(
        model.clone().with_temperatures(vec![20.0; 3], 20.0),
        Err(MailleError::QuasiStaticWithTemperatures)
    );
    let heated = Model::new(mesh(elastic, PlaneState::Strain), bcs(), Vec::new())
        .and_then(|m| m.with_temperatures(vec![30.0; 3], 20.0))
        .unwrap();
    assert_eq!(
        heated.solve_quasi_static(&step, &ramp),
        Err(MailleError::QuasiStaticWithTemperatures)
    );

    let conductive = Material::conductive(50.0).unwrap();
    let thermal = Model::with_analysis(
        mesh(conductive, PlaneState::Strain),
        Analysis::Thermal,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(
        thermal.solve_quasi_static(&step, &ramp),
        Err(MailleError::StaticOnly)
    );
}
