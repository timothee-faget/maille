//! Thermique transitoire (θ-méthode) sur triangles et tétraèdres, comparée aux solutions analytiques.

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{Analysis, Bc, Convection, HeatSource, Model, Step};
use maille_core::results::Increment;
use maille_core::MailleError;

/// Unités cohérentes au choix de l'utilisateur. Diffusivité α = k / (ρc) = 1.
const K: f64 = 1.0;
const RHO: f64 = 2.0;
const C: f64 = 0.5;
const L: f64 = 1.0;
const H: f64 = 0.25;
const T: f64 = 0.25;

fn material(k: f64) -> Material {
    Material::conductive(k)
        .and_then(|m| m.with_density(RHO))
        .and_then(|m| m.with_specific_heat(C))
        .unwrap()
}

/// Mesh L × H × T (nx Elements en x, 1 en y et z), avec la conductivité `k`.
fn mesh(name: &str, nx: usize, k: f64) -> Mesh {
    let m = material(k);
    match name {
        "Tri3" | "Tri6" => {
            let (nodes, tris) = rect_grid(nx, 1, L, H);
            if name == "Tri6" {
                let q = to_quadratic(&nodes, &tris);
                let elements = q
                    .elements
                    .iter()
                    .map(|&n| Element::Tri6(Tri6::new(n, m, T, PlaneState::Stress).unwrap()))
                    .collect();
                Mesh::new(q.nodes, elements).unwrap()
            } else {
                let elements = tris
                    .iter()
                    .map(|&n| Element::Tri3(Tri3::new(n, m, T, PlaneState::Stress).unwrap()))
                    .collect();
                Mesh::new(nodes, elements).unwrap()
            }
        }
        _ => {
            let (nodes, tets) = box_grid([nx, 1, 1], [L, H, T]);
            if name == "Tet10" {
                let q = to_quadratic_tet(&nodes, &tets);
                let elements = q
                    .elements
                    .iter()
                    .map(|&t| Element::Tet10(Tet10::new(t, m)))
                    .collect();
                Mesh::new(q.nodes, elements).unwrap()
            } else {
                let elements = tets
                    .iter()
                    .map(|&t| Element::Tet4(Tet4::new(t, m)))
                    .collect();
                Mesh::new(nodes, elements).unwrap()
            }
        }
    }
}

const ELEMENTS: [&str; 4] = ["Tri3", "Tri6", "Tet4", "Tet10"];

/// Volume du Mesh : L H T dans tous les cas (épaisseur T en 2D).
const VOLUME: f64 = L * H * T;

fn fix_plane(mesh: &Mesh, x: f64, value: f64) -> Vec<Bc> {
    mesh.nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| (p.x - x).abs() < 1e-12)
        .map(|(node, _)| Bc {
            node,
            dof: Dof::T,
            value,
        })
        .collect()
}

/// Convection sur toutes les faces du plan x = `x`.
fn convection_at(mesh: &Mesh, x: f64, film: f64, sink: f64) -> Vec<Convection> {
    let mut convection = Vec::new();
    for (element, e) in mesh.elements.iter().enumerate() {
        for face in 0..4 {
            let on = |n: &usize| (mesh.nodes[*n].x - x).abs() < 1e-12;
            if e.face_nodes(face).is_some_and(|nodes| nodes.iter().all(on)) {
                convection.push(Convection {
                    element,
                    face,
                    film,
                    sink,
                });
            }
        }
    }
    convection
}

fn thermal(mesh: &Mesh, bcs: Vec<Bc>) -> Model {
    Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, Vec::new()).unwrap()
}

/// Moyenne volumique de T, ∫ T dV / V = Σ Tᵢ ρc ∫ Nᵢ dV / (ρc V), `weights` donnant ρc ∫ Nᵢ dV.
fn mean(increment: &Increment, weights: &[f64]) -> f64 {
    let t = &increment.field.temperatures;
    t.iter().zip(weights).map(|(t, w)| t * w).sum::<f64>() / (RHO * C * VOLUME)
}

/// Σ_j C_ij = ρc ∫ Nᵢ dV, calculé ici sans le cœur. ∫ Nᵢ dV vaut V/3 aux sommets d'un Tri3 et V/4
/// d'un Tet4 ; 0 aux sommets et V/3 aux milieux d'un Tri6 ; −V/20 aux sommets et V/5 aux milieux
/// d'un Tet10.
fn weights(mesh: &Mesh) -> Vec<f64> {
    let mut w = vec![0.0; mesh.nodes.len()];
    for e in &mesh.elements {
        let nodes = e.nodes();
        let edge = |k: usize| {
            let (a, b) = (&mesh.nodes[nodes[0]], &mesh.nodes[nodes[k]]);
            [b.x - a.x, b.y - a.y, b.z - a.z]
        };
        let (u, v, s) = (edge(1), edge(2), edge(3.min(nodes.len() - 1)));
        let area = (u[0] * v[1] - u[1] * v[0]).abs() / 2.0;
        let tet = (u[0] * (v[1] * s[2] - v[2] * s[1]) - u[1] * (v[0] * s[2] - v[2] * s[0])
            + u[2] * (v[0] * s[1] - v[1] * s[0]))
            .abs()
            / 6.0;
        // (volume, nombre de sommets, part d'un sommet, part d'un milieu)
        let (volume, vertices, at_vertex, at_middle) = match e {
            Element::Tri3(t) => (area * t.thickness, 3, 1.0 / 3.0, 0.0),
            Element::Tri6(t) => (area * t.thickness, 3, 0.0, 1.0 / 3.0),
            Element::Tet4(_) => (tet, 4, 0.25, 0.0),
            Element::Tet10(_) => (tet, 4, -0.05, 0.2),
            Element::Bar(_) | Element::Beam(_) => unreachable!("pas de Bar ni de Beam ici"),
        };
        for (k, &n) in nodes.iter().enumerate() {
            w[n] += RHO * C * volume * if k < vertices { at_vertex } else { at_middle };
        }
    }
    w
}

#[test]
fn increments_start_from_the_initial_field_at_regular_times() {
    let m = mesh("Tri3", 4, K);
    let initial: Vec<f64> = m.nodes.iter().map(|p| 10.0 * p.x).collect();
    let step = Step::transient(0.5, 4).unwrap();
    assert!((step.increment() - 0.125).abs() < 1e-15 && (step.theta() - 1.0).abs() < 1e-15);
    let history = thermal(&m, fix_plane(&m, 0.0, 30.0))
        .solve_step(&step, &initial)
        .unwrap();
    assert_eq!(history.len(), 5);
    for (s, inc) in history.iter().enumerate() {
        assert!((inc.time - 0.125 * s as f64).abs() < 1e-15);
        assert_eq!(inc.field.temperatures.len(), m.nodes.len());
        assert_eq!(inc.field.heat_fluxes.len(), m.elements.len());
    }
    // t = 0 : la condition initiale, sauf aux Nodes imposés, qui prennent la BC dès le départ
    let first = &history[0].field;
    for (n, p) in m.nodes.iter().enumerate() {
        let expected = if p.x.abs() < 1e-12 { 30.0 } else { 10.0 * p.x };
        assert!((first.temperatures[n] - expected).abs() < 1e-15, "Node {n}");
    }
    assert!(first.reactions.iter().all(|r| r.abs() < 1e-15));
}

#[test]
fn starting_from_the_steady_state_stays_there() {
    // mur T1 / T2 avec source : le stationnaire est un point fixe de la θ-méthode
    for name in ELEMENTS {
        let m = mesh(name, 4, K);
        let mut bcs = fix_plane(&m, 0.0, 100.0);
        bcs.extend(fix_plane(&m, L, 20.0));
        let sources: Vec<HeatSource> = (0..m.elements.len())
            .map(|element| HeatSource {
                element,
                value: 300.0,
            })
            .collect();
        let model = thermal(&m, bcs).with_heat_sources(sources).unwrap();
        let steady = model.solve().unwrap();
        let step = Step::transient(1.0, 5).unwrap().with_theta(0.5).unwrap();
        let history = model.solve_step(&step, &steady.temperatures).unwrap();
        for inc in &history[1..] {
            for (a, b) in inc.field.temperatures.iter().zip(&steady.temperatures) {
                assert!((a - b).abs() < 1e-10, "{name} t = {}", inc.time);
            }
            for (a, b) in inc.field.reactions.iter().zip(&steady.reactions) {
                assert!((a - b).abs() < 1e-9, "{name} réaction {a} au lieu de {b}");
            }
        }
    }
}

#[test]
fn long_time_reaches_the_steady_state() {
    for name in ELEMENTS {
        let m = mesh(name, 4, K);
        let mut bcs = fix_plane(&m, 0.0, 100.0);
        bcs.extend(fix_plane(&m, L, 20.0));
        let sources: Vec<HeatSource> = (0..m.elements.len())
            .map(|element| HeatSource {
                element,
                value: 300.0,
            })
            .collect();
        let model = thermal(&m, bcs).with_heat_sources(sources).unwrap();
        let steady = model.solve().unwrap();
        // temps de diffusion L²/α = 1 : 40 fois plus long, Euler implicite à grands pas
        let history = model
            .solve_step(
                &Step::transient(40.0, 40).unwrap(),
                &vec![0.0; m.nodes.len()],
            )
            .unwrap();
        let last = &history[40].field;
        for (a, b) in last.temperatures.iter().zip(&steady.temperatures) {
            assert!((a - b).abs() < 1e-8, "{name} : {a} au lieu de {b}");
        }
    }
}

/// Mur de Fourier : T = 0 au départ, T = Ts sur x = 0 et x = L pour t > 0 :
/// T(x, t) = Ts [1 − Σ_{n impair} 4/(nπ) sin(nπx/L) exp(−α (nπ/L)² t)].
fn fourier(x: f64, t: f64, ts: f64) -> f64 {
    let alpha = K / (RHO * C);
    let pi = std::f64::consts::PI;
    let sum: f64 = (0..500)
        .map(|i| {
            let n = (2 * i + 1) as f64;
            4.0 / (n * pi) * (n * pi * x / L).sin() * (-alpha * (n * pi / L).powi(2) * t).exp()
        })
        .sum();
    ts * (1.0 - sum)
}

/// Écart maximal aux Nodes avec la série de Fourier, pour le mur trempé (`name`, `nx` Elements en x,
/// Step de 0,1 L²/α en `increments` pas d'Euler implicite), en t = 0,02, 0,05 et 0,1.
fn quench_errors(name: &str, nx: usize, increments: usize) -> [f64; 3] {
    let ts = 100.0;
    let m = mesh(name, nx, K);
    let mut bcs = fix_plane(&m, 0.0, ts);
    bcs.extend(fix_plane(&m, L, ts));
    let step = Step::transient(0.1, increments).unwrap();
    let history = thermal(&m, bcs)
        .solve_step(&step, &vec![0.0; m.nodes.len()])
        .unwrap();
    [increments / 5, increments / 2, increments].map(|s| {
        let inc = &history[s];
        m.nodes
            .iter()
            .zip(&inc.field.temperatures)
            .map(|(p, t)| (t - fourier(p.x, inc.time, ts)).abs() / ts)
            .fold(0.0, f64::max)
    })
}

#[test]
fn quenched_slab_follows_the_fourier_series() {
    // Euler implicite : un départ discontinu fait osciller Crank-Nicolson. De t = 0,02 (front
    // encore raide) à t = 0,1 (profil presque sinusoïdal), écart sous 0,4 % de Ts.
    for (name, nx, increments) in [
        ("Tri3", 40, 400),
        ("Tri6", 10, 400),
        ("Tet4", 40, 400),
        ("Tet10", 10, 400),
    ] {
        let errors = quench_errors(name, nx, increments);
        assert!(errors.iter().all(|&e| e < 4e-3), "{name} : {errors:?}");
    }
}

#[test]
fn quenched_slab_converges_in_time_and_space() {
    // Tri6 et Tet10 : erreur en espace négligeable, Euler d'ordre 1 (Δt divisé par 4, erreur ÷ 4)
    for name in ["Tri6", "Tet10"] {
        let ratio = quench_errors(name, 20, 100)[2] / quench_errors(name, 20, 400)[2];
        assert!((3.4..4.4).contains(&ratio), "{name} : rapport {ratio:.2}");
    }
    // Tri3 et Tet4 : l'erreur en espace domine ; h divisé par 2 la divise par 3 à 3,5 (ordre
    // proche de 2 : l'erreur en temps n'est pas nulle)
    for name in ["Tri3", "Tet4"] {
        let ratio = quench_errors(name, 10, 400)[2] / quench_errors(name, 20, 400)[2];
        assert!((2.8..4.4).contains(&ratio), "{name} : rapport {ratio:.2}");
    }
}

#[test]
fn small_biot_cooling_follows_the_discrete_amplification_factor() {
    // conductivité énorme (Biot 10⁻⁷) : T reste uniforme, ρcV dT/dt = −h A (T − T∞) sur la face
    // x = L d'aire A = H T, d'où τ = ρcV/(hA) ; un pas θ multiplie T − T∞ par
    // g = (1 − (1 − θ) r) / (1 + θ r), r = Δt/τ
    let (film, sink, t0) = (0.4, 20.0, 100.0);
    let tau = RHO * C * VOLUME / (film * H * T);
    for name in ELEMENTS {
        let m = mesh(name, 2, 1e7 * film * L);
        let w = weights(&m);
        let model = thermal(&m, Vec::new())
            .with_convection(convection_at(&m, L, film, sink))
            .unwrap();
        for theta in [1.0, 0.75, 0.5] {
            let step = Step::transient(tau, 8).unwrap().with_theta(theta).unwrap();
            let r = step.increment() / tau;
            let g = (1.0 - (1.0 - theta) * r) / (1.0 + theta * r);
            let history = model.solve_step(&step, &vec![t0; m.nodes.len()]).unwrap();
            for (s, inc) in history.iter().enumerate() {
                let expected = sink + (t0 - sink) * g.powi(s as i32);
                let got = mean(inc, &w);
                assert!(
                    (got - expected).abs() < 1e-5 * (t0 - sink),
                    "{name} θ = {theta}, incrément {s} : {got} au lieu de {expected}"
                );
                let spread = inc
                    .field
                    .temperatures
                    .iter()
                    .map(|t| (t - got).abs())
                    .fold(0.0, f64::max);
                assert!(spread < 1e-5 * (t0 - sink), "{name} : T non uniforme");
            }
        }
    }
}

#[test]
fn time_order_is_one_for_euler_and_two_for_crank_nicolson() {
    // même refroidissement, comparé à exp(−t/τ) en t = τ avec Δt = τ/10 puis τ/20
    let (film, sink, t0) = (0.4, 20.0, 100.0);
    let tau = RHO * C * VOLUME / (film * H * T);
    let m = mesh("Tri3", 2, 1e7 * film * L);
    let w = weights(&m);
    let model = thermal(&m, Vec::new())
        .with_convection(convection_at(&m, L, film, sink))
        .unwrap();
    let exact = sink + (t0 - sink) * (-1.0_f64).exp();
    let error = |theta: f64, n: usize| {
        let step = Step::transient(tau, n).unwrap().with_theta(theta).unwrap();
        let history = model.solve_step(&step, &vec![t0; m.nodes.len()]).unwrap();
        (mean(&history[n], &w) - exact).abs()
    };
    let euler = error(1.0, 10) / error(1.0, 20);
    let crank = error(0.5, 10) / error(0.5, 20);
    assert!((1.9..2.1).contains(&euler), "Euler : rapport {euler:.3}");
    assert!(
        (3.9..4.1).contains(&crank),
        "Crank-Nicolson : rapport {crank:.3}"
    );
}

#[test]
fn energy_balance_is_exact() {
    // ∫ ρc (T_fin − T_0) dV = Σ_incréments Δt (puissance des BC + puissance des sources) :
    // la conduction ne crée ni ne détruit de chaleur, quel que soit θ
    let q = 50.0;
    for name in ELEMENTS {
        let m = mesh(name, 4, K);
        let w = weights(&m);
        let sources: Vec<HeatSource> = (0..m.elements.len())
            .map(|element| HeatSource { element, value: q })
            .collect();
        let model = thermal(&m, fix_plane(&m, 0.0, 80.0))
            .with_heat_sources(sources)
            .unwrap();
        let step = Step::transient(0.3, 12).unwrap().with_theta(0.5).unwrap();
        let initial: Vec<f64> = m.nodes.iter().map(|p| 20.0 + 5.0 * p.x).collect();
        let history = model.solve_step(&step, &initial).unwrap();
        let stored = RHO * C * VOLUME * (mean(&history[12], &w) - mean(&history[0], &w));
        let supplied: f64 = history[1..]
            .iter()
            .map(|inc| step.increment() * (inc.field.reactions.iter().sum::<f64>() + q * VOLUME))
            .sum();
        assert!(
            (stored - supplied).abs() < 1e-10 * stored.abs(),
            "{name} : stocké {stored}, apporté {supplied}"
        );
    }
}

#[test]
fn step_inputs_are_typed_errors() {
    assert_eq!(
        Step::transient(0.0, 10),
        Err(MailleError::InvalidProperty {
            name: "durée de la Step",
            value: 0.0
        })
    );
    assert_eq!(Step::transient(1.0, 0), Err(MailleError::NoIncrements));
    let step = Step::transient(1.0, 10).unwrap();
    assert_eq!(
        step.with_theta(0.4),
        Err(MailleError::ThetaOutOfRange { theta: 0.4 })
    );
    assert!(step.with_theta(f64::NAN).is_err());

    let m = mesh("Tri3", 1, K);
    let n = m.nodes.len();
    let model = thermal(&m, Vec::new());
    assert_eq!(
        model.solve_step(&step, &[0.0]),
        Err(MailleError::InitialTemperatureCount {
            expected: n,
            got: 1
        })
    );
    let mut nan = vec![0.0; n];
    nan[2] = f64::NAN;
    assert_eq!(
        model.solve_step(&step, &nan),
        Err(MailleError::NonFiniteInitialTemperature { node: 2 })
    );

    // Model statique, Material sans ρ ou sans c
    let elastic = Material::with_poisson(210e3, 0.3).unwrap();
    let tri =
        |material| Element::Tri3(Tri3::new([0, 1, 2], material, T, PlaneState::Stress).unwrap());
    let nodes = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(1.0, 0.0),
        Node::new_2d(0.0, 1.0),
    ];
    let one = |material| Mesh::new(nodes.clone(), vec![tri(material)]).unwrap();
    let static_model = Model::new(one(elastic), Vec::new(), Vec::new()).unwrap();
    assert_eq!(
        static_model.solve_step(&step, &[0.0; 3]),
        Err(MailleError::ThermalOnly)
    );
    let no_density = Material::conductive(K)
        .and_then(|m| m.with_specific_heat(C))
        .unwrap();
    assert_eq!(
        thermal(&one(no_density), Vec::new()).solve_step(&step, &[0.0; 3]),
        Err(MailleError::MissingDensity { element: 0 })
    );
    let no_heat = Material::conductive(K)
        .and_then(|m| m.with_density(RHO))
        .unwrap();
    assert_eq!(
        thermal(&one(no_heat), Vec::new()).solve_step(&step, &[0.0; 3]),
        Err(MailleError::MissingSpecificHeat { element: 0 })
    );
}
