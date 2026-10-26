//! Banc de performance (M8) : temps de `solve`, `solve_step` et `solve_quasi_static` sur des cas
//! types, et du solveur creux seul (Cuthill-McKee inverse, Cholesky) sur une grille 3D.
//! `cargo bench -p maille-core -- <filtre>` pour n'en lancer qu'une partie.

mod common;

use std::time::Duration;

use common::{box_grid, to_quadratic_tet};
use criterion::measurement::WallTime;
use criterion::{
    criterion_group, criterion_main, BatchSize, BenchmarkGroup, BenchmarkId, Criterion, Throughput,
};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tet10, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Analysis, Bc, Load, Model, Step};
use maille_core::solve::{reverse_cuthill_mckee, SkylineMatrix, SolverKind};

/// Unités cohérentes (N, mm, MPa).
const E: f64 = 210e3;
const NU: f64 = 0.3;
const SY: f64 = 250.0;

fn bc(node: usize, dof: Dof, value: f64) -> Bc {
    Bc { node, dof, value }
}

fn tet10_mesh(n: [usize; 3], size: [f64; 3], material: Material) -> Mesh {
    let (nodes, tets) = box_grid(n, size);
    let q = to_quadratic_tet(&nodes, &tets);
    let elements = q
        .elements
        .iter()
        .map(|&t| Element::Tet10(Tet10::new(t, material)))
        .collect();
    Mesh::new(q.nodes, elements).unwrap()
}

/// Console 10 × 1 × 1 en Tet10, encastrée en x = 0, effort selon z réparti sur les Nodes de x = L.
fn cantilever(n: [usize; 3]) -> (Model, usize) {
    let mesh = tet10_mesh(n, [10.0, 1.0, 1.0], Material::with_poisson(E, NU).unwrap());
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
    let nodes = mesh.nodes.len();
    (Model::new(mesh, bcs, loads).unwrap(), nodes)
}

/// Quart de tube a = 10, b = 20 en Tri6 (déformations planes) de taille h, avec un déplacement
/// radial `bore` à l'alésage ; plastique si `plastic`.
fn tube(h: f64, bore: f64, plastic: bool) -> (Model, usize) {
    let geometry = Geometry {
        outer: vec![
            Curve::line([10.0, 0.0], h, "sym_x"),
            Curve::arc([20.0, 0.0], [0.0, 0.0], false, h, "exterieur"),
            Curve::line([0.0, 20.0], h, "sym_y"),
            Curve::arc([0.0, 10.0], [0.0, 0.0], true, h, "interieur"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())
        .unwrap()
        .to_quadratic();
    let [sym_x, sym_y, inner] =
        ["sym_x", "sym_y", "interieur"].map(|name| t.boundary(name).unwrap().nodes.clone());
    let mut material = Material::with_poisson(E, 0.49).unwrap();
    if plastic {
        material = material.with_plasticity(SY, 0.0).unwrap();
    }
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
            bcs.push(bc(node, Dof::Ux, bore * p.x / r));
        }
        if !sym_x.contains(&node) {
            bcs.push(bc(node, Dof::Uy, bore * p.y / r));
        }
    }
    let nodes = mesh.nodes.len();
    (Model::new(mesh, bcs, Vec::new()).unwrap(), nodes)
}

/// Bloc thermique 4 × 1 × 1 en Tet10, T = 1 imposée en x = 0, départ à T = 0.
fn heated_block(k: usize) -> (Model, Vec<f64>, usize) {
    let material = Material::conductive(1.0)
        .and_then(|m| m.with_density(2.0))
        .and_then(|m| m.with_specific_heat(0.5))
        .unwrap();
    let mesh = tet10_mesh([4 * k, k, k], [4.0, 1.0, 1.0], material);
    let bcs = (0..mesh.nodes.len())
        .filter(|&node| mesh.nodes[node].x.abs() < 1e-12)
        .map(|node| bc(node, Dof::T, 1.0))
        .collect();
    let initial = vec![0.0; mesh.nodes.len()];
    let nodes = mesh.nodes.len();
    let model = Model::with_analysis(mesh, Analysis::Thermal, bcs, Vec::new()).unwrap();
    (model, initial, nodes)
}

fn slow<'a>(c: &'a mut Criterion, name: &str) -> BenchmarkGroup<'a, WallTime> {
    let mut group = c.benchmark_group(name);
    group
        .sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(10));
    group
}

fn bench_cantilever(c: &mut Criterion) {
    let mut group = slow(c, "cantilever_tet10_solve");
    for n in [[10, 2, 2], [20, 4, 4], [30, 6, 6], [40, 8, 8]] {
        let (model, nodes) = cantilever(n);
        let dofs = 3 * nodes;
        group.throughput(Throughput::Elements(dofs as u64));
        for (name, kind) in [
            ("direct", SolverKind::Direct),
            ("iterative", SolverKind::Iterative),
        ] {
            let model = model.clone().with_solver(kind);
            group.bench_function(BenchmarkId::new(name, dofs), |b| {
                b.iter(|| model.solve().unwrap())
            });
        }
    }
    group.finish();
}

fn bench_tube(c: &mut Criterion) {
    let mut group = slow(c, "tube_tri6");
    for h in [2.0, 1.0, 0.5] {
        let (elastic, nodes) = tube(h, 0.01, false);
        let (plastic, _) = tube(h, 0.03, true);
        let dofs = 2 * nodes;
        group.throughput(Throughput::Elements(dofs as u64));
        group.bench_function(BenchmarkId::new("elastic_solve", dofs), |b| {
            b.iter(|| elastic.solve().unwrap())
        });
        let step = Step::quasi_static(1.0, 8).unwrap();
        let ramp = Amplitude::ramp(1.0).unwrap();
        group.bench_function(BenchmarkId::new("plastic_quasi_static", dofs), |b| {
            b.iter(|| plastic.solve_quasi_static(&step, &ramp).unwrap())
        });
    }
    group.finish();
}

fn bench_transient(c: &mut Criterion) {
    let mut group = slow(c, "transient_tet10");
    let step = Step::transient(1.0, 20).unwrap();
    for k in [2, 4, 6] {
        let (model, initial, dofs) = heated_block(k);
        group.throughput(Throughput::Elements(dofs as u64));
        group.bench_function(BenchmarkId::new("dofs_20_steps", dofs), |b| {
            b.iter(|| model.solve_step(&step, &initial).unwrap())
        });
    }
    group.finish();
}

/// Graphe de la grille n × n × n à 7 points (un DOF par sommet), numérotée avec x qui varie le plus
/// vite, c'est-à-dire un bloc dont le Cuthill-McKee inverse doit retrouver la largeur de bande n².
fn grid_graph(n: usize) -> Vec<Vec<usize>> {
    let id = |i: usize, j: usize, k: usize| i + n * (j + n * k);
    let mut adjacency = vec![Vec::new(); n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                for (di, dj, dk) in [(1, 0, 0), (0, 1, 0), (0, 0, 1)] {
                    let (a, b) = (id(i, j, k), id(i + di, j + dj, k + dk));
                    if i + di < n && j + dj < n && k + dk < n {
                        adjacency[a].push(b);
                        adjacency[b].push(a);
                    }
                }
            }
        }
    }
    adjacency
}

/// Matrice du laplacien de `adjacency` (4 + degré sur la diagonale, −1 hors diagonale) dans le
/// profil de la numérotation `order`.
fn laplacian(adjacency: &[Vec<usize>], order: &[usize]) -> SkylineMatrix {
    let mut rank = vec![0; order.len()];
    for (new, &old) in order.iter().enumerate() {
        rank[old] = new;
    }
    let mut first: Vec<usize> = (0..rank.len()).collect();
    for (v, neighbours) in adjacency.iter().enumerate() {
        first[rank[v]] = neighbours
            .iter()
            .map(|&w| rank[w])
            .fold(rank[v], usize::min);
    }
    let mut a = SkylineMatrix::zeros(first);
    for (v, neighbours) in adjacency.iter().enumerate() {
        a.add(rank[v], rank[v], 4.0 + neighbours.len() as f64);
        for &w in neighbours.iter().filter(|&&w| w < v) {
            a.add(rank[v], rank[w], -1.0);
        }
    }
    a
}

fn bench_sparse_solver(c: &mut Criterion) {
    let mut group = slow(c, "sparse_solver_grid3d");
    for n in [10, 16, 22] {
        let adjacency = grid_graph(n);
        let order = reverse_cuthill_mckee(&adjacency);
        let matrix = laplacian(&adjacency, &order);
        let dofs = n * n * n;
        println!(
            "grille {n}³ : {dofs} DOF, {} termes stockés",
            matrix.stored()
        );
        group.throughput(Throughput::Elements(dofs as u64));
        group.bench_function(BenchmarkId::new("rcm", dofs), |b| {
            b.iter(|| reverse_cuthill_mckee(&adjacency))
        });
        group.bench_function(BenchmarkId::new("factorize", dofs), |b| {
            b.iter_batched(
                || matrix.clone(),
                |m| m.factorize().unwrap(),
                BatchSize::LargeInput,
            )
        });
        let factor = matrix.factorize().unwrap();
        let rhs = vec![1.0; dofs];
        group.bench_function(BenchmarkId::new("solve", dofs), |b| {
            b.iter(|| factor.solve(&rhs))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_cantilever,
    bench_tube,
    bench_transient,
    bench_sparse_solver
);
criterion_main!(benches);
