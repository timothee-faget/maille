//! Les exemples de `examples/` : chacun est résolu par le binaire `maille` et comparé à sa solution
//! analytique, pour qu'ils restent justes. Un test `#[ignore]` régénère les `.msh` 2D depuis les
//! `.dxf` avec `maille mesh` : `cargo test -p maille-cli --test examples -- --ignored`.

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use maille_core::mesh::Node;
use maille_io::msh::{write, CellKind, MshFile, Set, DOMAIN_SET};
use maille_io::vtk::{read_pvd, read_results, VtuResults};

fn examples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples")
}

// --- génération des maillages ---

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
    for (name, size) in [("plaque", "0.5"), ("trempe", "1"), ("tube", "1")] {
        let out = Command::new(env!("CARGO_BIN_EXE_maille"))
            .args([
                "mesh",
                &format!("{name}.dxf"),
                "--size",
                size,
                "--quadratic",
            ])
            .current_dir(&dir)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{name} : {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    fs::write(dir.join("cube.msh"), write(&cube()).unwrap()).unwrap();
}

// --- exécution du binaire sur une copie des exemples ---

struct Workdir(PathBuf);

impl Workdir {
    /// Copie des `.toml`, `.msh`, `.dxf`, `.stl` et `.step` d'`examples/` dans un dossier temporaire : les résultats n'y
    /// salissent pas le dépôt.
    fn new(test: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("maille-{test}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for entry in fs::read_dir(examples_dir()).unwrap() {
            let path = entry.unwrap().path();
            if matches!(
                path.extension().and_then(|e| e.to_str()),
                Some("toml" | "msh" | "dxf" | "stl" | "step")
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
fn dxf_drawing_meshed_by_the_cli_follows_the_analytic_field() {
    let w = Workdir::new("dxf");
    let out = w.maille(&[
        "mesh",
        "plaque.dxf",
        "--size",
        "0.5",
        "--quadratic",
        "-o",
        "dessin.msh",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout).contains("bas, droite, gauche, haut"));
    // le modèle de la traction, sur le maillage issu du dessin
    let model = fs::read_to_string(w.0.join("traction.toml")).unwrap();
    fs::write(
        w.0.join("dessin.toml"),
        model.replace("plaque.msh", "dessin.msh"),
    )
    .unwrap();
    w.solve("dessin.toml");
    let (_, r) = w.results("dessin.pvd", None);
    let (u, s) = (array(&r, "displacement"), array(&r, "stress"));
    for (i, p) in r.nodes.iter().enumerate() {
        assert!((u[3 * i] - 0.001 * p.x).abs() < 1e-9);
        assert!((s[6 * i] - 210.0).abs() < 1e-6, "σx = {}", s[6 * i]);
    }
}

/// `cube.stl` (6 solids nommés comme les Sets de `cube.toml`) meshé par la CLI, à la finesse du STL
/// puis à la taille h : l'élasticité linéaire retrouve exactement le champ affine, en Tet4 comme en
/// Tet10.
#[test]
fn stl_surface_meshed_by_the_cli_follows_the_analytic_field() {
    let w = Workdir::new("stl");
    let cases = [
        (false, None, CellKind::Tet4),
        (true, None, CellKind::Tet10),
        (false, Some("0.25"), CellKind::Tet4),
        (true, Some("0.5"), CellKind::Tet10),
    ];
    for (quadratic, size, kind) in cases {
        let mut args = vec!["mesh", "cube.stl", "-o", "piece.msh"];
        if quadratic {
            args.push("--quadratic");
        }
        if let Some(h) = size {
            args.extend(["--size", h]);
        }
        let out = w.maille(&args);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(String::from_utf8_lossy(&out.stdout).contains("bords x0, x1, y0, y1, z0, z1"));
        let model = fs::read_to_string(w.0.join("cube.toml")).unwrap();
        fs::write(
            w.0.join("piece.toml"),
            model.replace("cube.msh", "piece.msh"),
        )
        .unwrap();
        w.solve("piece.toml");
        let (_, r) = w.results("piece.pvd", None);
        assert_eq!(r.kind, kind);
        let u = array(&r, "displacement");
        for (i, p) in r.nodes.iter().enumerate() {
            assert!((u[3 * i] - 0.01 * p.x).abs() < 1e-10, "ux en {p:?}");
            assert!((u[3 * i + 1] + 0.25 * 0.01 * p.y).abs() < 1e-10);
            assert!((u[3 * i + 2] + 0.25 * 0.01 * p.z).abs() < 1e-10);
        }
    }
}

/// Pièce STEP maillée par la CLI : cylindre de rayon 2 et de hauteur 3 à couture, comme l'écrit
/// OCCT (`tests/cylindre.step`). Volume des tétraèdres : 12π par défaut, de moins que la flèche
/// relative des cordes (h / r)² / 8 ; sans `--size`, refus.
#[test]
fn step_piece_meshed_by_the_cli() {
    let w = Workdir::new("step");
    fs::write(w.0.join("cylindre.step"), include_str!("cylindre.step")).unwrap();
    let out = w.maille(&["mesh", "cylindre.step", "--size", "0.5"]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(stdout.contains("unité de longueur mm"), "{stdout}");
    assert!(stdout.contains("bords cote, fond, haut"), "{stdout}");
    let file =
        maille_io::msh::read(&fs::read_to_string(w.0.join("cylindre.msh")).unwrap()).unwrap();
    let p = |i: usize| [file.nodes[i].x, file.nodes[i].y, file.nodes[i].z];
    let volume: f64 = file
        .cells
        .iter()
        .map(|c| {
            let [a, b, c, d] = [c[0], c[1], c[2], c[3]].map(p);
            let (u, v, w) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                [d[0] - a[0], d[1] - a[1], d[2] - a[2]],
            );
            (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]))
                / 6.0
        })
        .sum();
    let e = (12.0 * PI - volume) / (12.0 * PI);
    assert!(e > 0.0 && e < (0.5f64 / 2.0).powi(2) / 8.0, "écart {e:e}");
    let out = w.maille(&["mesh", "cylindre.step"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--size manque"));
}

/// `maille view --png` écrit une image PNG d'un résultat 3D, avec coupe et déformée. Sans carte
/// graphique (machine d'intégration continue), le test est sauté.
/// `tube3d.stl` (quart de tube facetté en 16 côtés) maillé à la taille h en Tet10 par la CLI, puis
/// `tube3d.toml` : pressions intérieure et extérieure, bouts bloqués en z. Lamé en déformations
/// planes, à l'erreur de facettage près (flèche des cordes : 1,2e-3 du rayon). Mesuré : u_r à
/// 1,85e-3 près, contraintes à 4,4e-3 près en RMS.
#[test]
fn tube3d_meshed_from_its_stl_follows_lame() {
    let w = Workdir::new("tube3d");
    let out = w.maille(&["mesh", "tube3d.stl", "--size", "2", "--quadratic"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    w.solve("tube3d.toml");
    let (_, r) = w.results("tube3d.pvd", None);
    assert_eq!(r.kind, CellKind::Tet10);
    let (u, s) = (array(&r, "displacement"), array(&r, "stress"));
    let (pi, po, ri, ro, e, nu): (f64, f64, f64, f64, f64, f64) =
        (100.0, 20.0, 10.0, 20.0, 210e3, 0.3);
    let span = ro * ro - ri * ri;
    let a = (pi * ri * ri - po * ro * ro) / span;
    let b = (pi - po) * ri * ri * ro * ro / span;
    let ur = |r: f64| (1.0 + nu) / e * ((1.0 - 2.0 * nu) * a * r + b / r);
    let scale = a + b / (ri * ri);
    let (mut worst_u, mut sum, mut count) = (0.0_f64, 0.0, 0);
    for (i, p) in r.nodes.iter().enumerate() {
        let radius = p.x.hypot(p.y);
        let (c, sn) = (p.x / radius, p.y / radius);
        let got_u = u[3 * i] * c + u[3 * i + 1] * sn;
        worst_u = worst_u.max((got_u - ur(radius)).abs() / ur(ri));
        let (xx, yy, xy) = (s[6 * i], s[6 * i + 1], s[6 * i + 3]);
        let got = [
            xx * c * c + yy * sn * sn + 2.0 * xy * c * sn,
            xx * sn * sn + yy * c * c - 2.0 * xy * c * sn,
            s[6 * i + 2],
        ];
        let want = [
            a - b / (radius * radius),
            a + b / (radius * radius),
            2.0 * nu * a,
        ];
        for (g, x) in got.iter().zip(want) {
            sum += (g - x) * (g - x);
            count += 1;
        }
    }
    let rms = (sum / count as f64).sqrt() / scale;
    assert!(worst_u < 2.5e-3, "u_r : {worst_u:.2e}");
    assert!(rms < 6e-3, "σ RMS : {rms:.2e}");
}

/// `sphere.step` (huitième de sphère creuse, faces sphériques exactes) maillé à la taille h en
/// Tet10 par la CLI, puis `sphere.toml` : pressions intérieure et extérieure, plans de coupe en
/// appui glissant. Lamé pour la sphère, à l'erreur de géométrie près : les sommets sont sur les
/// sphères, les triangles en dessous. Mesuré : u_r à 1,32e-2 près, σ à 4,2e-3 près en RMS à h = 2 ;
/// 3,6e-3 et 1,2e-3 à h = 1 (ordre 2, environ 2,6 fois la flèche (h / Ri)² / 8 d'un cylindre : la
/// sphère est courbée dans les deux sens).
#[test]
fn sphere_meshed_from_its_step_follows_lame() {
    let w = Workdir::new("sphere");
    let out = w.maille(&["mesh", "sphere.step", "--size", "2", "--quadratic"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    w.solve("sphere.toml");
    let (_, r) = w.results("sphere.pvd", None);
    assert_eq!(r.kind, CellKind::Tet10);
    let (u, s) = (array(&r, "displacement"), array(&r, "stress"));
    let (pi, po, ri, ro, e, nu): (f64, f64, f64, f64, f64, f64) =
        (100.0, 20.0, 10.0, 20.0, 210e3, 0.3);
    let span = ro.powi(3) - ri.powi(3);
    let c = (pi * ri.powi(3) - po * ro.powi(3)) / span;
    let d = (pi - po) * ri.powi(3) * ro.powi(3) / span;
    let ur = |r: f64| (1.0 - 2.0 * nu) / e * c * r + (1.0 + nu) / (2.0 * e) * d / (r * r);
    let scale = pi;
    let (mut worst_u, mut sum, mut count) = (0.0_f64, 0.0, 0);
    for (i, p) in r.nodes.iter().enumerate() {
        let radius = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        let n = [p.x / radius, p.y / radius, p.z / radius];
        let got_u = (0..3).map(|k| u[3 * i + k] * n[k]).sum::<f64>();
        worst_u = worst_u.max((got_u - ur(radius)).abs() / ur(ri));
        // σ dans l'ordre xx, yy, zz, xy, yz, xz : σ_rr = nᵀ σ n, σ_θθ = (tr σ − σ_rr) / 2
        let t = &s[6 * i..6 * i + 6];
        let sigma = [[t[0], t[3], t[5]], [t[3], t[1], t[4]], [t[5], t[4], t[2]]];
        let rr: f64 = (0..3)
            .flat_map(|a| (0..3).map(move |b| (a, b)))
            .map(|(a, b)| n[a] * sigma[a][b] * n[b])
            .sum();
        let hoop = (t[0] + t[1] + t[2] - rr) / 2.0;
        let r3 = radius.powi(3);
        for (g, x) in [(rr, c - d / r3), (hoop, c + d / (2.0 * r3))] {
            sum += (g - x) * (g - x);
            count += 1;
        }
    }
    let rms = (sum / count as f64).sqrt() / scale;
    assert!(worst_u < 1.6e-2, "u_r : {worst_u:.2e}");
    assert!(rms < 6e-3, "σ RMS : {rms:.2e}");
}

#[test]
fn view_writes_a_png_of_a_3d_result() {
    let w = Workdir::new("png");
    let out = w.maille(&["mesh", "cube.stl", "-o", "piece.msh"]);
    assert!(out.status.success());
    let model = fs::read_to_string(w.0.join("cube.toml")).unwrap();
    fs::write(
        w.0.join("piece.toml"),
        model.replace("cube.msh", "piece.msh"),
    )
    .unwrap();
    w.solve("piece.toml");
    let out = w.maille(&[
        "view",
        "piece.pvd",
        "--png",
        "vue.png",
        "--size",
        "300x200",
        "--angles",
        "30,20",
        "--field",
        "8",
        "--deformed",
        "--cut",
        "x=0.5",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    if stderr.contains("aucune carte graphique") {
        eprintln!("test sauté : {stderr}");
        return;
    }
    assert!(out.status.success(), "{stderr}");
    let png = fs::read(w.0.join("vue.png")).unwrap();
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    // IHDR : largeur et hauteur demandées
    assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 300);
    assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 200);
}

#[test]
fn view_command_reports_errors() {
    let w = Workdir::new("view-erreurs");
    let fail = |args: &[&str], expected: &str| {
        let out = w.maille(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr.contains(expected), "{args:?} : {stderr}");
    };
    fail(&["view"], "usage");
    // les options sont lues avant d'ouvrir le fichier
    fail(
        &["view", "a.pvd", "--png", "x.png", "--size", "grand"],
        "largeurxhauteur",
    );
    fail(
        &["view", "a.pvd", "--png", "x.png", "--size", "0x10"],
        "1 à 8192",
    );
    fail(
        &["view", "a.pvd", "--png", "x.png", "--angles", "30"],
        "lacet,tangage",
    );
    fail(
        &["view", "a.pvd", "--png", "x.png", "--field", "0"],
        "numéro de la touche",
    );
    fail(
        &["view", "a.pvd", "--png", "x.png", "--cut", "w=0.5"],
        "axe « w » inconnu",
    );
    fail(
        &["view", "a.pvd", "--png", "x.png", "--cut", "x"],
        "axe=position",
    );
}

#[test]
fn stl_mesh_command_reports_errors() {
    let w = Workdir::new("stl-erreurs");
    let fail = |args: &[&str], expected: &str| {
        let out = w.maille(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr.contains(expected), "{args:?} : {stderr}");
    };
    fail(&["mesh", "cube.stl", "--size", "0"], "taille h");
    fail(&["mesh", "cube.stl", "--size", "fin"], "attend un nombre");
    fail(
        &["mesh", "plaque.dxf", "--size", "1", "--feature-angle", "20"],
        "ne sert qu'aux fichiers STL",
    );
    fail(
        &["mesh", "cube.stl", "--feature-angle", "beaucoup"],
        "attend un nombre",
    );
    fail(&["mesh", "cube.stl", "--feature-angle", "0"], "angle vif");
    fail(
        &["mesh", "cube.toml"],
        "extension « .toml » non prise en charge",
    );
    fail(&["mesh", "absent.stl"], "absent.stl");
    // une face manque : surface ouverte
    let text = fs::read_to_string(w.0.join("cube.stl")).unwrap();
    let open = text.replacen("facet normal", "facette normal", 1);
    fs::write(w.0.join("casse.stl"), open).unwrap();
    fail(&["mesh", "casse.stl"], "STL");
}

#[test]
fn mesh_command_reports_errors() {
    let w = Workdir::new("mesh-erreurs");
    let fail = |args: &[&str], expected: &str| {
        let out = w.maille(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr.contains(expected), "{args:?} : {stderr}");
    };
    fail(&["mesh"], "usage");
    fail(&["mesh", "plaque.dxf"], "--size manque");
    fail(
        &["mesh", "plaque.dxf", "--size", "grand"],
        "--size attend un nombre",
    );
    fail(
        &["mesh", "plaque.dxf", "--size", "-1"],
        "strictement positive",
    );
    fail(&["mesh", "absent.dxf", "--size", "1"], "absent.dxf");
    // un côté manque : contour ouvert
    let text = fs::read_to_string(w.0.join("plaque.dxf")).unwrap();
    let open = &text[..text.rfind("0\nLINE").unwrap()];
    fs::write(w.0.join("ouvert.dxf"), format!("{open}0\nENDSEC\n0\nEOF\n")).unwrap();
    fail(&["mesh", "ouvert.dxf", "--size", "1"], "contour ouvert");
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
fn tube_under_pressure_follows_lame() {
    let w = Workdir::new("tube-pression");
    w.solve("tube_pression.toml");
    let (_, r) = w.results("tube_pression.pvd", None);
    let s = array(&r, "stress");
    let (pi, po, ri, ro, nu): (f64, f64, f64, f64, f64) = (100.0, 20.0, 10.0, 20.0, 0.3);
    let span = ro * ro - ri * ri;
    let a = (pi * ri * ri - po * ro * ro) / span;
    let b = (pi - po) * ri * ri * ro * ro / span;
    let scale = a + b / (ri * ri);
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
        let want = [
            a - b / (radius * radius),
            a + b / (radius * radius),
            2.0 * nu * a,
        ];
        for (g, x) in got.iter().zip(want) {
            sum += (g - x) * (g - x);
            count += 1;
            assert!(
                (g - x).abs() < 5e-2 * scale,
                "r = {radius} : {g} contre {x}"
            );
        }
    }
    let rms = (sum / count as f64).sqrt() / scale;
    assert!(rms < 4e-3, "erreur RMS {rms:.2e}");
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
    // le cube n'a qu'un incrément
    w.solve("cube.toml");
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
