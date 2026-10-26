//! Les exemples de `examples/` : chacun est résolu par le binaire `maille` et comparé à sa solution
//! analytique, pour qu'ils restent justes. Un test `#[ignore]` régénère les `.msh` (le mailleur
//! n'a pas de commande avant M9) : `cargo test -p maille-cli --test examples -- --ignored`.

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use maille_core::mesh::Node;
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_io::msh::{write, CellKind, MshFile, Set, DOMAIN_SET};
use maille_io::vtk::{read_pvd, read_results, VtuResults};

fn examples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples")
}

// --- génération des maillages ---

fn polygon(points: &[[f64; 2]], size: f64, names: &[&str]) -> Geometry {
    Geometry {
        outer: (0..points.len())
            .map(|i| Curve::line(points[i], size, names[i]))
            .collect(),
        holes: vec![],
    }
}

fn mesh_text(geometry: &Geometry, quadratic: bool) -> String {
    let t = triangulate(geometry, &MeshOptions::default()).unwrap();
    let file = if quadratic {
        MshFile::from_triangulation(&t.to_quadratic())
    } else {
        MshFile::from_triangulation(&t)
    }
    .unwrap();
    write(&file).unwrap()
}

/// Cube unité en 6 tétraèdres Tet4 ; Sets de Nodes x0, x1, y0, z0.
fn cube() -> MshFile {
    let nodes: Vec<Node> = (0..8)
        .map(|k| Node::new((k & 1) as f64, ((k >> 1) & 1) as f64, ((k >> 2) & 1) as f64))
        .collect();
    let orders = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let cells = orders
        .iter()
        .map(|o| {
            let (mut p, mut cell) = (0usize, vec![0usize]);
            for &axis in o {
                p |= 1 << axis;
                cell.push(p);
            }
            cell
        })
        .collect();
    let on = |axis: usize, bit: usize| Set {
        nodes: (0..8).filter(|k| (k >> axis) & 1 == bit).collect(),
        ..Set::default()
    };
    let mut sets = BTreeMap::new();
    for (name, axis, bit) in [("x0", 0, 0), ("x1", 0, 1), ("y0", 1, 0), ("z0", 2, 0)] {
        sets.insert(name.to_string(), on(axis, bit));
    }
    sets.insert(
        DOMAIN_SET.to_string(),
        Set {
            elements: (0..6).collect(),
            ..Set::default()
        },
    );
    MshFile {
        nodes,
        kind: CellKind::Tet4,
        cells,
        sets,
    }
}

#[test]
#[ignore = "régénère les .msh d'examples/"]
fn regenerate_example_meshes() {
    let dir = examples_dir();
    let plate = polygon(
        &[[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]],
        0.5,
        &["bas", "droite", "haut", "gauche"],
    );
    let quench = polygon(
        &[[0.0, 0.0], [20.0, 0.0], [20.0, 10.0], [0.0, 10.0]],
        1.0,
        &["isole", "trempe", "isole", "trempe"],
    );
    let tube = Geometry {
        outer: vec![
            Curve::line([10.0, 0.0], 1.0, "sym_x"),
            Curve::arc([20.0, 0.0], [0.0, 0.0], false, 1.0, "exterieur"),
            Curve::line([0.0, 20.0], 1.0, "sym_y"),
            Curve::arc([0.0, 10.0], [0.0, 0.0], true, 1.0, "interieur"),
        ],
        holes: vec![],
    };
    fs::write(dir.join("plaque.msh"), mesh_text(&plate, true)).unwrap();
    fs::write(dir.join("trempe.msh"), mesh_text(&quench, true)).unwrap();
    fs::write(dir.join("tube.msh"), mesh_text(&tube, true)).unwrap();
    fs::write(dir.join("cube.msh"), write(&cube()).unwrap()).unwrap();
}

// --- exécution du binaire sur une copie des exemples ---

struct Workdir(PathBuf);

impl Workdir {
    /// Copie des `.toml` et `.msh` d'`examples/` dans un dossier temporaire : les résultats n'y
    /// salissent pas le dépôt.
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("maille-{test}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for entry in fs::read_dir(examples_dir()).unwrap() {
            let path = entry.unwrap().path();
            if matches!(
                path.extension().and_then(|e| e.to_str()),
                Some("toml" | "msh")
            ) {
                fs::copy(&path, dir.join(path.file_name().unwrap())).unwrap();
            }
        }
        Self(dir)
    }

    fn maille(&self, args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_maille"))
            .args(args)
            .current_dir(&self.0)
            .output()
            .unwrap()
    }

    fn solve(&self, model: &str) {
        let out = self.maille(&["solve", model]);
        assert!(
            out.status.success(),
            "{model} : {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }

    /// Résultats de l'incrément `k` du `.pvd` (le dernier si `None`), avec son temps.
    fn results(&self, pvd: &str, k: Option<usize>) -> (f64, VtuResults) {
        let entries = read_pvd(&fs::read_to_string(self.0.join(pvd)).unwrap()).unwrap();
        let entry = &entries[k.unwrap_or(entries.len() - 1)];
        let text = fs::read_to_string(self.0.join(&entry.file)).unwrap();
        (entry.time, read_results(&text).unwrap())
    }
}

impl Drop for Workdir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn array<'a>(r: &'a VtuResults, name: &str) -> &'a [f64] {
    &r.arrays.iter().find(|a| a.name == name).unwrap().values
}

#[test]
fn traction_follows_the_analytic_field() {
    let w = Workdir::new("traction");
    w.solve("traction.toml");
    let (_, r) = w.results("traction.pvd", None);
    assert_eq!(r.units.as_deref(), Some("N, mm, MPa"));
    let (u, s) = (array(&r, "displacement"), array(&r, "stress"));
    for (i, p) in r.nodes.iter().enumerate() {
        assert!((u[3 * i] - 0.001 * p.x).abs() < 1e-9);
        assert!((u[3 * i + 1] + 0.3 * 0.001 * p.y).abs() < 1e-9);
        assert!((s[6 * i] - 210.0).abs() < 1e-6, "σx = {}", s[6 * i]);
        assert!(s[6 * i + 1].abs() < 1e-6);
    }
}

#[test]
fn cube_follows_the_analytic_field() {
    let w = Workdir::new("cube");
    w.solve("cube.toml");
    let (_, r) = w.results("cube.pvd", None);
    assert_eq!(r.kind, CellKind::Tet4);
    let u = array(&r, "displacement");
    for (i, p) in r.nodes.iter().enumerate() {
        assert!((u[3 * i] - 0.01 * p.x).abs() < 1e-12);
        assert!((u[3 * i + 1] + 0.25 * 0.01 * p.y).abs() < 1e-12);
        assert!((u[3 * i + 2] + 0.25 * 0.01 * p.z).abs() < 1e-12);
    }
}

#[test]
fn quench_follows_the_fourier_series() {
    let w = Workdir::new("trempe");
    w.solve("trempe.toml");
    let (alpha, width) = (0.05 / (7.85e-6 * 460.0), 20.0);
    for k in [20, 100, 200] {
        let (time, r) = w.results("trempe.pvd", Some(k));
        let t = array(&r, "temperature");
        for (p, temp) in r.nodes.iter().zip(t) {
            let sum: f64 = (0..500)
                .map(|i| {
                    let n = (2 * i + 1) as f64;
                    let kn = n * PI / width;
                    4.0 / (n * PI) * (kn * p.x).sin() * (-alpha * kn * kn * time).exp()
                })
                .sum();
            let exact = 20.0 + 830.0 * sum;
            assert!(
                (temp - exact).abs() < 8.3,
                "t = {time} : écart {}",
                temp - exact
            );
        }
    }
}

#[test]
fn tube_thermal_stresses_follow_timoshenko() {
    let w = Workdir::new("tube");
    // le couplage lit les résultats du calcul thermique : l'ordre compte
    let missing = w.maille(&["solve", "tube_contraintes.toml"]);
    assert!(!missing.status.success());
    assert!(String::from_utf8_lossy(&missing.stderr).contains("tube_thermique.pvd"));
    w.solve("tube_thermique.toml");
    w.solve("tube_contraintes.toml");

    let (_, r) = w.results("tube_contraintes.pvd", None);
    let s = array(&r, "stress");
    let (e, nu, alpha, ri, ro): (f64, f64, f64, f64, f64) = (210e3, 0.3, 1.2e-5, 10.0, 20.0);
    let beta = (20.0 - 120.0) / (ro / ri).ln();
    let delta = |r: f64| 100.0 + beta * (r / ri).ln();
    let integral = |r: f64| {
        100.0 * (r * r - ri * ri) / 2.0
            + beta * (r * r / 2.0 * (r / ri).ln() - r * r / 4.0 + ri * ri / 4.0)
    };
    let exact = |r: f64| {
        let c = alpha * e / ((1.0 - nu) * r * r);
        let span = ro * ro - ri * ri;
        let sr = c * ((r * r - ri * ri) / span * integral(ro) - integral(r));
        let st = c * ((r * r + ri * ri) / span * integral(ro) + integral(r) - delta(r) * r * r);
        [sr, st, nu * (sr + st) - e * alpha * delta(r)]
    };
    let scale = exact(ri)[1].abs();
    let (mut sum, mut count) = (0.0, 0);
    for (i, p) in r.nodes.iter().enumerate() {
        let radius = p.x.hypot(p.y);
        let (c, sn) = (p.x / radius, p.y / radius);
        let (xx, yy, xy) = (s[6 * i], s[6 * i + 1], s[6 * i + 3]);
        let got = [
            xx * c * c + yy * sn * sn + 2.0 * xy * c * sn,
            xx * sn * sn + yy * c * c - 2.0 * xy * c * sn,
            s[6 * i + 2],
        ];
        for (g, x) in got.iter().zip(exact(radius)) {
            sum += (g - x) * (g - x);
            count += 1;
            assert!(
                (g - x).abs() < 5e-2 * scale,
                "r = {radius} : {g} contre {x}"
            );
        }
    }
    let rms = (sum / count as f64).sqrt() / scale;
    assert!(rms < 3e-3, "erreur RMS {rms:.2e}");
}

#[test]
fn commands_report_errors_and_usage() {
    let w = Workdir::new("erreurs");
    let help = w.maille(&[]);
    assert!(help.status.success());
    assert!(String::from_utf8_lossy(&help.stdout).contains("maille solve"));

    let fail = |args: &[&str], expected: &str| {
        let out = w.maille(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr.contains(expected), "{args:?} : {stderr}");
    };
    fail(&["solve", "absent.toml"], "absent.toml");
    fail(&["solve"], "usage");
    fail(&["view"], "usage");
    fail(&["view", "cube.toml"], ".pvd ou .vtu");
    fail(&["inconnu"], "commande inconnue");
    // un résultat 3D ne s'affiche pas dans le viewer 2D
    w.solve("cube.toml");
    fail(&["view", "cube.pvd"], "ParaView");
    fail(&["view", "cube.pvd", "5"], "incrément 5");
    // une clé mal écrite dans le modèle
    let text = fs::read_to_string(w.0.join("cube.toml")).unwrap();
    fs::write(w.0.join("faute.toml"), text.replace("young", "youngg")).unwrap();
    fail(&["solve", "faute.toml"], "youngg");
}

#[test]
fn output_prefix_is_honoured() {
    let w = Workdir::new("sortie");
    fs::create_dir(w.0.join("out")).unwrap();
    let out = w.maille(&["solve", "cube.toml", "-o", "out/essai"]);
    assert!(out.status.success());
    assert!(w.0.join("out/essai.vtu").exists() && w.0.join("out/essai.pvd").exists());
    assert!(!w.0.join("cube.pvd").exists());
}
