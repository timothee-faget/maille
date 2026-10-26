//! Solveur itératif (gradients conjugués préconditionnés par IC(0)) comparé au solveur direct sur
//! les trois chemins de calcul : statique, transitoire et quasi-statique (Newton-Raphson).

mod common;

use common::{box_grid, to_quadratic_tet};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tet10, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Analysis, Bc, Load, Model, Step};
use maille_core::solve::SolverKind;
use maille_core::MailleError;

const E: f64 = 210e3;

fn bc(node: usize, dof: Dof, value: f64) -> Bc {
    Bc { node, dof, value }
}

/// Plus grand écart |a − b| rapporté au plus grand |b| (écart absolu si b est nul).
fn relative_gap(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = b
        .iter()
        .fold(0.0_f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    let scale = if scale < 1e-300 { 1.0 } else { scale };
    a.iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
        / scale
}

fn tet10(n: [usize; 3], size: [f64; 3], material: Material) -> Mesh {
    let (nodes, tets) = box_grid(n, size);
    let q = to_quadratic_tet(&nodes, &tets);
    let elements = q
        .elements
        .iter()
        .map(|&t| Element::Tet10(Tet10::new(t, material)))
        .collect();
    Mesh::new(q.nodes, elements).unwrap()
}

/// Console 10 × 1 × 1 encastrée en x = 0, effort selon z sur x = L.
fn cantilever() -> Model {
    let mesh = tet10(
        [10, 2, 2],
        [10.0, 1.0, 1.0],
        Material::with_poisson(E, 0.3).unwrap(),
    );
    let mut bcs = Vec::new();
    let mut loads = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        if p.x.abs() < 1e-12 {
            bcs.extend([Dof::Ux, Dof::Uy, Dof::Uz].map(|dof| bc(node, dof, 0.0)));
        } else if (p.x - 10.0).abs() < 1e-12 {
            loads.push(Load {
                node,
                dof: Dof::Uz,
                value: 1.0,
            });
        }
    }
    Model::new(mesh, bcs, loads).unwrap()
}

#[test]
fn iterative_static_solution_matches_the_direct_one() {
    let model = cantilever();
    let direct = model
        .clone()
        .with_solver(SolverKind::Direct)
        .solve()
        .unwrap();
    let iterative = model.with_solver(SolverKind::Iterative).solve().unwrap();
    assert!(relative_gap(&iterative.displacements, &direct.displacements) < 1e-8);
    assert!(relative_gap(&iterative.reactions, &direct.reactions) < 1e-8);
    let von_mises = |f: &maille_core::results::Field| -> Vec<f64> {
        f.stresses.iter().flatten().map(|s| s.sigma_xx).collect()
    };
    assert!(relative_gap(&von_mises(&iterative), &von_mises(&direct)) < 1e-7);
}

#[test]
fn auto_picks_the_direct_solver_for_small_models() {
    // 3 000 DOF environ, sous le seuil : Auto doit donner exactement le résultat du direct
    let model = cantilever();
    let auto = model.clone().solve().unwrap();
    let direct = model.with_solver(SolverKind::Direct).solve().unwrap();
    assert_eq!(auto.displacements, direct.displacements);
}

#[test]
fn iterative_transient_matches_the_direct_one() {
    let material = Material::conductive(1.0)
        .and_then(|m| m.with_density(2.0))
        .and_then(|m| m.with_specific_heat(0.5))
        .unwrap();
    let mesh = tet10([8, 2, 2], [4.0, 1.0, 1.0], material);
    let bcs = (0..mesh.nodes.len())
        .filter(|&node| mesh.nodes[node].x.abs() < 1e-12)
        .map(|node| bc(node, Dof::T, 1.0))
        .collect();
    let initial = vec![0.0; mesh.nodes.len()];
    let model = Model::with_analysis(mesh, Analysis::Thermal, bcs, Vec::new()).unwrap();
    let step = Step::transient(1.0, 10).unwrap();
    let direct = model
        .clone()
        .with_solver(SolverKind::Direct)
        .solve_step(&step, &initial)
        .unwrap();
    let iterative = model
        .with_solver(SolverKind::Iterative)
        .solve_step(&step, &initial)
        .unwrap();
    for (d, i) in direct.iter().zip(&iterative) {
        assert!(relative_gap(&i.field.temperatures, &d.field.temperatures) < 1e-8);
    }
}

/// Quart de tube a = 10, b = 20 en Tri6 plastique, déplacement radial imposé à l'alésage.
fn plastic_tube() -> Model {
    let geometry = Geometry {
        outer: vec![
            Curve::line([10.0, 0.0], 2.0, "sym_x"),
            Curve::arc([20.0, 0.0], [0.0, 0.0], false, 2.0, "exterieur"),
            Curve::line([0.0, 20.0], 2.0, "sym_y"),
            Curve::arc([0.0, 10.0], [0.0, 0.0], true, 2.0, "interieur"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())
        .unwrap()
        .to_quadratic();
    let [sym_x, sym_y, inner] =
        ["sym_x", "sym_y", "interieur"].map(|name| t.boundary(name).unwrap().nodes.clone());
    let material = Material::with_poisson(E, 0.49)
        .and_then(|m| m.with_plasticity(250.0, 0.0))
        .unwrap();
    let elements = t
        .triangles
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material, 1.0, PlaneState::Strain).unwrap()))
        .collect();
    let mesh = Mesh::new(t.nodes, elements).unwrap();
    let mut bcs: Vec<Bc> = sym_x
        .iter()
        .map(|&node| bc(node, Dof::Uy, 0.0))
        .chain(sym_y.iter().map(|&node| bc(node, Dof::Ux, 0.0)))
        .collect();
    for &node in &inner {
        let p = mesh.nodes[node];
        let r = p.x.hypot(p.y);
        if !sym_y.contains(&node) {
            bcs.push(bc(node, Dof::Ux, 0.03 * p.x / r));
        }
        if !sym_x.contains(&node) {
            bcs.push(bc(node, Dof::Uy, 0.03 * p.y / r));
        }
    }
    Model::new(mesh, bcs, Vec::new()).unwrap()
}

#[test]
fn iterative_newton_raphson_matches_the_direct_one() {
    let model = plastic_tube();
    let (step, ramp) = (
        Step::quasi_static(1.0, 6).unwrap(),
        Amplitude::ramp(1.0).unwrap(),
    );
    let direct = model
        .clone()
        .with_solver(SolverKind::Direct)
        .solve_quasi_static(&step, &ramp)
        .unwrap();
    let iterative = model
        .with_solver(SolverKind::Iterative)
        .solve_quasi_static(&step, &ramp)
        .unwrap();
    assert_eq!(direct.len(), iterative.len());
    for (d, i) in direct.iter().zip(&iterative) {
        assert!(relative_gap(&i.field.displacements, &d.field.displacements) < 1e-7);
        assert_eq!(
            i.residuals.len(),
            d.residuals.len(),
            "même nombre d'itérations de Newton"
        );
    }
    // le tube a bien plastifié : sinon le test ne dirait rien de la tangente
    let last = &iterative.last().unwrap().field;
    assert!(last
        .equivalent_plastic_strains
        .iter()
        .flatten()
        .any(|&p| p > 1e-4));
}

#[test]
fn iterative_solver_reports_an_unconstrained_structure() {
    // console sans encastrement, chargée : pas de solution, le solveur ne doit pas rendre de bruit
    let mesh = tet10(
        [4, 1, 1],
        [4.0, 1.0, 1.0],
        Material::with_poisson(E, 0.3).unwrap(),
    );
    let load = Load {
        node: 0,
        dof: Dof::Uz,
        value: 1.0,
    };
    let result = Model::new(mesh, vec![], vec![load])
        .unwrap()
        .with_solver(SolverKind::Iterative)
        .solve();
    assert!(
        matches!(
            result,
            Err(MailleError::Singular | MailleError::IterativeNotConverged { .. })
        ),
        "{:?}",
        result.map(|f| f.displacements.len())
    );
}

/// Console 10 × 1 × 1 d'environ 10⁶ DOF, hors de portée du Cholesky en ligne de ciel, comparée à
/// la flèche de la théorie des poutres. Lancer : `cargo test --release --test iterative -- --ignored
/// --nocapture`.
#[test]
#[ignore = "long : environ 10⁶ DOF"]
fn million_dof_cantilever_matches_beam_theory() {
    let mesh = tet10(
        [100, 20, 20],
        [10.0, 1.0, 1.0],
        Material::with_poisson(E, 0.3).unwrap(),
    );
    let dofs = 3 * mesh.nodes.len();
    let mut bcs = Vec::new();
    let mut tip = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        if p.x.abs() < 1e-12 {
            bcs.extend([Dof::Ux, Dof::Uy, Dof::Uz].map(|dof| bc(node, dof, 0.0)));
        } else if (p.x - 10.0).abs() < 1e-12 {
            tip.push(node);
        }
    }
    // effort total 1 réparti au prorata égal sur les Nodes de la section (la flexion ne dépend que
    // de la résultante, à l'effet de Saint-Venant près)
    let loads: Vec<Load> = tip
        .iter()
        .map(|&node| Load {
            node,
            dof: Dof::Uz,
            value: 1.0 / tip.len() as f64,
        })
        .collect();
    let model = Model::new(mesh.clone(), bcs, loads).unwrap();
    let start = std::time::Instant::now();
    let field = model.with_solver(SolverKind::Iterative).solve().unwrap();
    eprintln!("{dofs} DOF résolus en {:?}", start.elapsed());

    // flèche moyenne de la section d'extrémité : P L³ / (3 E I), I = 1/12, plus le cisaillement
    let mean: f64 = tip
        .iter()
        .map(|&n| field.displacements[3 * n + 2])
        .sum::<f64>()
        / tip.len() as f64;
    let bending = 1.0 * 1000.0 / (3.0 * E / 12.0);
    assert!(
        (mean / bending - 1.0).abs() < 0.05,
        "flèche {mean:e}, poutre {bending:e}"
    );
}
