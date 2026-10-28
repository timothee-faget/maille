//! Fichier de modèle TOML : les calculs passent par `.msh` et `.toml` et retrouvent la solution
//! analytique et le résultat du cœur appelé directement ; un test par message d'erreur.
//!
//! Les BC et Load du TOML ont une valeur par Set : les cas à valeur variable le long d'un bord
//! (Lamé, Kirsch en déplacements imposés) ne s'écrivent pas encore, ils restent dans le CLI.

use std::f64::consts::PI;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Bc, Model, Step};
use maille_core::MailleError;
use maille_io::model::{build, parse, Loaded, Run};
use maille_io::msh::{read, write, CellKind, MshFile, Set, DOMAIN_SET};
use maille_io::IoError;

/// Rectangle w × h, bords bas, droite, haut, gauche, passé par un vrai fichier `.msh`.
fn rectangle(w: f64, h: f64, size: f64, quadratic: bool) -> MshFile {
    let pts = [[0.0, 0.0], [w, 0.0], [w, h], [0.0, h]];
    let names = ["bas", "droite", "haut", "gauche"];
    let geometry = Geometry {
        outer: (0..4)
            .map(|i| Curve::line(pts[i], size, names[i]))
            .collect(),
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default()).unwrap();
    let file = if quadratic {
        MshFile::from_triangulation(&t.to_quadratic())
    } else {
        MshFile::from_triangulation(&t)
    }
    .unwrap();
    read(&write(&file).unwrap()).unwrap()
}

fn load(toml: &str, mesh: &MshFile) -> Loaded {
    build(&parse(toml).unwrap(), mesh, None).unwrap()
}

fn error(toml: &str, mesh: &MshFile) -> IoError {
    build(&parse(toml).unwrap(), mesh, None).unwrap_err()
}

const UNIAXIAL: &str = r#"
units = "N, mm, MPa"
mesh = "plaque.msh"
analysis = "static"

[materials.acier]
young = 210e3
poisson = 0.3

[[sections]]
set = "domaine"
material = "acier"
thickness = 2.0
plane_state = "stress"

[[bcs]]
set = "gauche"
dof = "ux"
value = 0.0

[[bcs]]
set = "droite"
dof = "ux"
value = 0.004

[[bcs]]
set = "bas"
dof = "uy"
value = 0.0
"#;

#[test]
fn uniaxial_traction_matches_the_analytic_field_and_the_core() {
    let mesh = rectangle(4.0, 2.0, 0.7, true);
    let loaded = load(UNIAXIAL, &mesh);
    assert_eq!(loaded.units, "N, mm, MPa");
    assert_eq!(loaded.run, Run::Linear);
    let field = loaded.model.solve().unwrap();

    // ε = 0,004 / 4 : ux = ε x, uy = −ν ε y, σx = E ε uniforme
    let (eps, nu, young, thickness) = (0.001, 0.3, 210e3, 2.0);
    for (i, p) in mesh.nodes.iter().enumerate() {
        assert!((field.displacements[2 * i] - eps * p.x).abs() < 1e-12);
        assert!((field.displacements[2 * i + 1] + nu * eps * p.y).abs() < 1e-12);
    }
    // effort de réaction sur x = 4 : σ × aire de la section (hauteur 2 × épaisseur 2)
    let right = &mesh.sets["droite"].nodes;
    let force: f64 = right.iter().map(|&n| field.reactions[2 * n]).sum();
    assert!((force - young * eps * 2.0 * thickness).abs() < 1e-6 * force);

    // même résultat que le cœur construit à la main
    let material = Material::with_poisson(young, nu).unwrap();
    let elements = mesh
        .cells
        .iter()
        .map(|c| {
            Tri6::new(
                c.clone().try_into().unwrap(),
                material,
                thickness,
                PlaneState::Stress,
            )
            .map(Element::Tri6)
        })
        .collect::<Result<_, _>>()
        .unwrap();
    let mut bcs = Vec::new();
    for (set, dof, value) in [
        ("gauche", Dof::Ux, 0.0),
        ("droite", Dof::Ux, 0.004),
        ("bas", Dof::Uy, 0.0),
    ] {
        bcs.extend(
            mesh.sets[set]
                .nodes
                .iter()
                .map(|&node| Bc { node, dof, value }),
        );
    }
    let direct = Model::new(
        Mesh::new(mesh.nodes.clone(), elements).unwrap(),
        bcs,
        vec![],
    )
    .unwrap()
    .solve()
    .unwrap();
    for (a, b) in field.displacements.iter().zip(&direct.displacements) {
        assert!((a - b).abs() < 1e-12);
    }
}

/// Cube unité en 6 tétraèdres (un par ordre des axes), Sets de Nodes sur les faces x = 0, x = 1,
/// y = 0, z = 0.
fn cube() -> MshFile {
    let mut nodes = Vec::new();
    for k in 0..8 {
        nodes.push(Node::new(
            (k & 1) as f64,
            ((k >> 1) & 1) as f64,
            ((k >> 2) & 1) as f64,
        ));
    }
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
    let mut sets = std::collections::BTreeMap::new();
    sets.insert("x0".to_string(), on(0, 0));
    sets.insert("x1".to_string(), on(0, 1));
    sets.insert("y0".to_string(), on(1, 0));
    sets.insert("z0".to_string(), on(2, 0));
    sets.insert(
        DOMAIN_SET.to_string(),
        Set {
            elements: (0..6).collect(),
            ..Set::default()
        },
    );
    let file = MshFile {
        nodes,
        kind: CellKind::Tet4,
        cells,
        sets,
    };
    read(&write(&file).unwrap()).unwrap()
}

#[test]
fn cube_in_traction_matches_the_analytic_field() {
    let mesh = cube();
    let toml = r#"
units = "N, mm, MPa"
mesh = "cube.msh"
analysis = "static"
[materials.acier]
young = 200e3
poisson = 0.25
[[sections]]
set = "domaine"
material = "acier"
[[bcs]]
set = "x0"
dof = "ux"
value = 0.0
[[bcs]]
set = "x1"
dof = "ux"
value = 0.01
[[bcs]]
set = "y0"
dof = "uy"
value = 0.0
[[bcs]]
set = "z0"
dof = "uz"
value = 0.0
"#;
    let field = load(toml, &mesh).model.solve().unwrap();
    // ε = 0,01 : ux = ε x, uy = −ν ε y, uz = −ν ε z
    for (i, p) in mesh.nodes.iter().enumerate() {
        assert!((field.displacements[3 * i] - 0.01 * p.x).abs() < 1e-12);
        assert!((field.displacements[3 * i + 1] + 0.25 * 0.01 * p.y).abs() < 1e-12);
        assert!((field.displacements[3 * i + 2] + 0.25 * 0.01 * p.z).abs() < 1e-12);
    }
}

const THERMAL: &str = r#"
units = "W, mm, K"
mesh = "plaque.msh"
analysis = "thermal"

[materials.m]
conductivity = 2.0

[[sections]]
set = "domaine"
material = "m"

[[bcs]]
set = "gauche"
dof = "t"
value = 100.0
"#;

#[test]
fn steady_conduction_with_convection_is_linear() {
    // T(x) = T0 − (T0 − T∞) h x / (k + h L) = 100 − 10 x pour h = 0,5, T∞ = 20, L = 4, k = 2
    let mesh = rectangle(4.0, 2.0, 0.7, true);
    let toml = format!("{THERMAL}\n[[convections]]\nset = \"droite\"\nfilm = 0.5\nsink = 20.0\n");
    let field = load(&toml, &mesh).model.solve().unwrap();
    for (p, t) in mesh.nodes.iter().zip(&field.temperatures) {
        assert!((t - (100.0 - 10.0 * p.x)).abs() < 1e-9, "x = {}", p.x);
    }
}

#[test]
fn steady_conduction_with_a_heat_source_is_parabolic() {
    // T = 0 aux deux extrémités, source Q : T = Q x (L − x) / (2 k) = 3 x (4 − x) / 4
    let mesh = rectangle(4.0, 2.0, 0.7, true);
    let toml = format!(
        "{}\n[[bcs]]\nset = \"droite\"\ndof = \"t\"\nvalue = 0.0\n\n[[heat_sources]]\nset = \"domaine\"\nvalue = 3.0\n",
        THERMAL.replace("value = 100.0", "value = 0.0")
    );
    let field = load(&toml, &mesh).model.solve().unwrap();
    for (p, t) in mesh.nodes.iter().zip(&field.temperatures) {
        assert!((t - 0.75 * p.x * (4.0 - p.x)).abs() < 1e-9, "x = {}", p.x);
    }
}

#[test]
fn plate_quench_follows_the_fourier_series() {
    // même cas que `maille-cli trempe` : plaque d'acier de 20 mm, de 850 °C à 20 °C en surface
    let (width, k, rho, c) = (20.0, 0.05, 7.85e-6, 460.0);
    let (t_start, t_surface) = (850.0, 20.0);
    let pts = [[0.0, 0.0], [width, 0.0], [width, 10.0], [0.0, 10.0]];
    let names = ["isole", "trempe", "isole", "trempe"];
    let geometry = Geometry {
        outer: (0..4).map(|i| Curve::line(pts[i], 1.0, names[i])).collect(),
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())
        .unwrap()
        .to_quadratic();
    let mesh = read(&write(&MshFile::from_triangulation(&t).unwrap()).unwrap()).unwrap();
    let toml = r#"
units = "W, mm, K, s, kg"
mesh = "trempe.msh"
analysis = "thermal"
[materials.acier]
conductivity = 0.05
density = 7.85e-6
specific_heat = 460.0
[[sections]]
set = "domaine"
material = "acier"
[[bcs]]
set = "trempe"
dof = "t"
value = 20.0
[step]
kind = "transient"
duration = 10.0
increments = 200
initial_temperature = 850.0
"#;
    let Run::Transient { step, initial } = load(toml, &mesh).run else {
        panic!("calcul transitoire attendu");
    };
    assert_eq!(initial, vec![850.0; mesh.nodes.len()]);
    let history = load(toml, &mesh).model.solve_step(&step, &initial).unwrap();
    let alpha = k / (rho * c);
    let exact = |x: f64, time: f64| {
        let sum: f64 = (0..500)
            .map(|i| {
                let n = (2 * i + 1) as f64;
                let kn = n * PI / width;
                4.0 / (n * PI) * (kn * x).sin() * (-alpha * kn * kn * time).exp()
            })
            .sum();
        t_surface + (t_start - t_surface) * sum
    };
    for s in [20, 100, 200] {
        let inc = &history[s];
        for (p, temp) in mesh.nodes.iter().zip(&inc.field.temperatures) {
            let e = (temp - exact(p.x, inc.time)).abs();
            assert!(
                e < 1e-2 * (t_start - t_surface),
                "t = {} : écart {e}",
                inc.time
            );
        }
    }
}

#[test]
fn quasi_static_step_matches_the_core() {
    let mesh = rectangle(4.0, 2.0, 0.8, false);
    let toml = r#"
units = "N, mm, MPa"
mesh = "plaque.msh"
analysis = "static"
[materials.acier]
young = 210e3
poisson = 0.49
yield_stress = 250.0
hardening = 1000.0
[[sections]]
set = "domaine"
material = "acier"
plane_state = "strain"
[[bcs]]
set = "gauche"
dof = "ux"
value = 0.0
[[bcs]]
set = "droite"
dof = "ux"
value = 0.02
[[bcs]]
set = "bas"
dof = "uy"
value = 0.0
[step]
kind = "quasi_static"
duration = 2.0
increments = 8
amplitude = [[0, 0], [1, 1], [2, 0.5]]
"#;
    let Run::QuasiStatic { step, amplitude } = load(toml, &mesh).run else {
        panic!("calcul quasi-statique attendu");
    };
    assert_eq!(step, Step::quasi_static(2.0, 8).unwrap());
    assert_eq!(
        amplitude,
        Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.5)]).unwrap()
    );
    let history = load(toml, &mesh)
        .model
        .solve_quasi_static(&step, &amplitude)
        .unwrap();

    let material = Material::with_poisson(210e3, 0.49)
        .unwrap()
        .with_plasticity(250.0, 1000.0)
        .unwrap();
    let elements = mesh
        .cells
        .iter()
        .map(|c| {
            Tri3::new(
                c.clone().try_into().unwrap(),
                material,
                1.0,
                PlaneState::Strain,
            )
            .map(Element::Tri3)
        })
        .collect::<Result<_, _>>()
        .unwrap();
    let mut bcs = Vec::new();
    for (set, dof, value) in [
        ("gauche", Dof::Ux, 0.0),
        ("droite", Dof::Ux, 0.02),
        ("bas", Dof::Uy, 0.0),
    ] {
        bcs.extend(
            mesh.sets[set]
                .nodes
                .iter()
                .map(|&node| Bc { node, dof, value }),
        );
    }
    let direct = Model::new(
        Mesh::new(mesh.nodes.clone(), elements).unwrap(),
        bcs,
        vec![],
    )
    .unwrap()
    .solve_quasi_static(&step, &amplitude)
    .unwrap();
    assert_eq!(history.len(), direct.len());
    for (a, b) in history.iter().zip(&direct) {
        for (x, y) in a.field.displacements.iter().zip(&b.field.displacements) {
            assert!((x - y).abs() < 1e-12);
        }
    }
    let plastic = history[4].field.equivalent_plastic_strains.iter().flatten();
    assert!(
        plastic.into_iter().any(|&p| p > 0.0),
        "le calcul doit plastifier"
    );
}

#[test]
fn explicit_core_error_is_wrapped() {
    // Material sans ν sur des triangles : le cœur le refuse
    let toml = UNIAXIAL.replace("poisson = 0.3\n", "");
    let e = error(&toml, &rectangle(4.0, 2.0, 1.0, false));
    assert_eq!(
        e,
        IoError::Model(MailleError::MissingPoisson { element: 0 })
    );
}

fn plate() -> MshFile {
    rectangle(4.0, 2.0, 1.0, false)
}

#[test]
fn syntax_errors_give_line_and_column() {
    let e = parse("units = \"N\"\nmesh = \n").unwrap_err();
    assert!(
        matches!(&e, IoError::Toml(m) if m.contains("line 2")),
        "{e}"
    );
}

#[test]
fn unknown_key_is_refused() {
    let e = parse(&UNIAXIAL.replace("young = 210e3", "youngg = 210e3")).unwrap_err();
    assert!(
        matches!(&e, IoError::Toml(m) if m.contains("youngg")),
        "{e}"
    );
}

#[test]
fn units_are_required() {
    let e = parse(&UNIAXIAL.replace("units = \"N, mm, MPa\"\n", "")).unwrap_err();
    assert!(matches!(&e, IoError::Toml(m) if m.contains("units")), "{e}");
    let e = error(&UNIAXIAL.replace("N, mm, MPa", "  "), &plate());
    assert!(e.to_string().contains("units"), "{e}");
}

#[test]
fn unknown_set_lists_the_known_ones() {
    let e = error(&UNIAXIAL.replace("\"gauche\"", "\"gauch\""), &plate());
    assert!(
        matches!(&e, IoError::UnknownSet { key, name, .. } if key == "bcs[0]" && name == "gauch")
    );
    let text = e.to_string();
    assert!(
        text.contains("« gauche »") && text.contains("« domaine »"),
        "{text}"
    );
}

#[test]
fn wrong_set_kind_is_refused() {
    let as_section = UNIAXIAL.replace("set = \"domaine\"", "set = \"bas\"");
    assert!(matches!(
        error(&as_section, &plate()),
        IoError::WrongSetKind { .. }
    ));
    let as_bc = UNIAXIAL.replace("set = \"gauche\"", "set = \"domaine\"");
    assert!(matches!(
        error(&as_bc, &plate()),
        IoError::WrongSetKind { .. }
    ));
    let convection = format!(
        "{}\n[[convections]]\nset = \"domaine\"\nfilm = 1.0\nsink = 0.0\n",
        THERMAL
    );
    assert!(matches!(
        error(&convection, &plate()),
        IoError::WrongSetKind { .. }
    ));
}

#[test]
fn unknown_material_lists_the_defined_ones() {
    let e = error(
        &UNIAXIAL.replace("material = \"acier\"", "material = \"alu\""),
        &plate(),
    );
    assert!(matches!(&e, IoError::UnknownMaterial { name, .. } if name == "alu"));
    assert!(e.to_string().contains("« acier »"));
}

#[test]
fn sections_must_cover_every_element_once() {
    let mut mesh = plate();
    let n = mesh.cells.len();
    mesh.sets
        .get_mut(DOMAIN_SET)
        .unwrap()
        .elements
        .truncate(n - 1);
    assert_eq!(
        error(UNIAXIAL, &mesh),
        IoError::UncoveredElements {
            count: 1,
            first: n - 1
        }
    );
    let twice = format!(
        "{UNIAXIAL}\n[[sections]]\nset = \"domaine\"\nmaterial = \"acier\"\nplane_state = \"stress\"\n"
    );
    assert!(matches!(
        error(&twice, &plate()),
        IoError::OverlappingSections {
            element: 0,
            first: 0,
            second: 1
        }
    ));
}

#[test]
fn dof_must_fit_analysis_and_dimension() {
    let uz = UNIAXIAL.replace("dof = \"uy\"", "dof = \"uz\"");
    assert!(error(&uz, &plate()).to_string().contains("2D"));
    let t = UNIAXIAL.replace("dof = \"uy\"", "dof = \"t\"");
    assert!(error(&t, &plate()).to_string().contains("thermique"));
    let ux_thermal = THERMAL.replace("dof = \"t\"", "dof = \"ux\"");
    assert!(error(&ux_thermal, &plate()).to_string().contains("`t`"));
}

#[test]
fn section_keys_must_fit_the_elements() {
    let missing = UNIAXIAL.replace("plane_state = \"stress\"\n", "");
    assert!(error(&missing, &plate())
        .to_string()
        .contains("plane_state"));
    let toml = r#"
units = "N, mm, MPa"
mesh = "cube.msh"
analysis = "static"
[materials.m]
young = 1.0
poisson = 0.3
[[sections]]
set = "domaine"
material = "m"
thickness = 1.0
"#;
    assert!(error(toml, &cube()).to_string().contains("tétraèdres"));
}

#[test]
fn incomplete_materials_are_refused() {
    let none = UNIAXIAL.replace("young = 210e3\npoisson = 0.3", "density = 1.0");
    assert!(error(&none, &plate()).to_string().contains("young"));
    let alone = UNIAXIAL.replace("young = 210e3\n", "");
    assert!(error(&alone, &plate()).to_string().contains("poisson"));
    let hardening = UNIAXIAL.replace("poisson = 0.3", "poisson = 0.3\nhardening = 5.0");
    assert!(error(&hardening, &plate())
        .to_string()
        .contains("yield_stress"));
}

#[test]
fn temperatures_section_and_argument_go_together() {
    let toml = format!(
        "{UNIAXIAL}\n[temperatures]\nfrom = \"trempe.pvd\"\nincrement = \"last\"\nreference = 20.0\n"
    );
    let mesh = plate();
    let file = parse(&toml).unwrap();
    let e = build(&file, &mesh, None).unwrap_err();
    assert!(e.to_string().contains("trempe.pvd"), "{e}");
    let e = build(&parse(UNIAXIAL).unwrap(), &mesh, Some(&[0.0])).unwrap_err();
    assert!(e.to_string().contains("[temperatures]"), "{e}");
}

#[test]
fn step_keys_must_fit_the_kind() {
    let transient = |extra: &str, toml: &str| {
        format!("{toml}\n[step]\nkind = \"transient\"\nduration = 1.0\nincrements = 4\n{extra}")
    };
    let no_t0 = transient("", THERMAL);
    assert!(error(&no_t0, &plate())
        .to_string()
        .contains("initial_temperature"));
    let in_static = transient("initial_temperature = 0.0\n", UNIAXIAL);
    assert!(error(&in_static, &plate()).to_string().contains("thermal"));
    let amplitude = transient("initial_temperature = 0.0\namplitude = [[0, 0]]\n", THERMAL);
    assert!(error(&amplitude, &plate())
        .to_string()
        .contains("amplitude"));

    let quasi = |extra: &str, toml: &str| {
        format!("{toml}\n[step]\nkind = \"quasi_static\"\nduration = 1.0\nincrements = 4\n{extra}")
    };
    let in_thermal = quasi("", THERMAL);
    assert!(error(&in_thermal, &plate()).to_string().contains("static"));
    let theta = quasi("theta = 0.5\n", UNIAXIAL);
    assert!(error(&theta, &plate()).to_string().contains("transient"));
    // une Amplitude qui ne part pas de (0, 0) : refus du cœur
    let bad = quasi("amplitude = [[0, 1]]\n", UNIAXIAL);
    assert!(matches!(
        error(&bad, &plate()),
        IoError::Model(MailleError::InvalidAmplitude { point: 0 })
    ));
}

#[test]
fn pressure_on_a_boundary_pulls_or_pushes() {
    // traction −100 sur x = 4, déformations planes : σx = 100, εx = (1 − ν²) σ / E
    let mesh = rectangle(4.0, 2.0, 0.7, true);
    let toml = UNIAXIAL
        .replace("plane_state = \"stress\"", "plane_state = \"strain\"")
        .replace("thickness = 2.0\n", "")
        .replace("value = 0.004", "value = 0.0")
        .replace("[[bcs]]\nset = \"droite\"\ndof = \"ux\"\nvalue = 0.0\n", "")
        + "\n[[pressures]]\nset = \"droite\"\nvalue = -100.0\n";
    let field = load(&toml, &mesh).model.solve().unwrap();
    let strain = (1.0 - 0.3 * 0.3) * 100.0 / 210e3;
    for (i, p) in mesh.nodes.iter().enumerate() {
        assert!(
            (field.displacements[2 * i] - strain * p.x).abs() < 1e-12,
            "x = {}",
            p.x
        );
    }
}

#[test]
fn pressure_needs_a_boundary_and_a_static_analysis() {
    let mesh = plate();
    let on_domain = format!("{UNIAXIAL}\n[[pressures]]\nset = \"domaine\"\nvalue = 1.0\n");
    assert!(matches!(
        error(&on_domain, &mesh),
        IoError::WrongSetKind { .. }
    ));
    let unknown = format!("{UNIAXIAL}\n[[pressures]]\nset = \"dessus\"\nvalue = 1.0\n");
    assert!(
        matches!(&error(&unknown, &mesh), IoError::UnknownSet { key, .. } if key == "pressures[0]")
    );
    let thermal = format!("{THERMAL}\n[[pressures]]\nset = \"droite\"\nvalue = 1.0\n");
    assert_eq!(
        error(&thermal, &mesh),
        IoError::Model(MailleError::PressureStaticOnly)
    );
}

#[test]
fn large_strain_keys_build_a_quasi_static_model() {
    let mesh = rectangle(4.0, 2.0, 0.8, false);
    let toml = |kinematics: &str, material: &str, step: &str| {
        format!(
            r#"
units = "N, mm, MPa"
mesh = "plaque.msh"
analysis = "static"
{kinematics}
[materials.caoutchouc]
young = 10.0
poisson = 0.3
{material}
[[sections]]
set = "domaine"
material = "caoutchouc"
plane_state = "strain"
[[bcs]]
set = "gauche"
dof = "ux"
value = 0.0
[[bcs]]
set = "bas"
dof = "uy"
value = 0.0
[[bcs]]
set = "droite"
dof = "ux"
value = 2.0
{step}
"#
        )
    };
    let step = "[step]\nkind = \"quasi_static\"\nduration = 1.0\nincrements = 4\n";

    // Saint-Venant–Kirchhoff, puis néo-Hookéen : le même modèle, deux lois, deux contraintes
    let svk = load(&toml("kinematics = \"large\"", "", step), &mesh);
    let neo = load(
        &toml(
            "kinematics = \"large\"",
            "hyperelasticity = \"neo_hookean\"",
            step,
        ),
        &mesh,
    );
    assert!(svk.large_strain && neo.large_strain);
    let stress = |loaded: &Loaded| match &loaded.run {
        Run::QuasiStatic { step, amplitude } => {
            let history = loaded.model.solve_quasi_static(step, amplitude).unwrap();
            history[4].field.stresses[0][0].sigma_xx
        }
        other => panic!("{other:?}"),
    };
    assert!((stress(&svk) - stress(&neo)).abs() > 1e-3);

    // petites déformations par défaut, ou écrites
    assert!(!load(&toml("", "", step), &mesh).large_strain);
    assert!(!load(&toml("kinematics = \"small\"", "", step), &mesh).large_strain);

    // les cas refusés disent quoi corriger
    let no_step = error(&toml("kinematics = \"large\"", "", ""), &mesh).to_string();
    assert!(no_step.contains("quasi_static"), "{no_step}");
    let law = error(&toml("", "hyperelasticity = \"neo_hookean\"", step), &mesh).to_string();
    assert!(law.contains("caoutchouc") && law.contains("large"), "{law}");
    let bad = parse(&toml("kinematics = \"enorme\"", "", step)).unwrap_err();
    assert!(bad.to_string().contains("enorme"), "{bad}");
    let bad = parse(&toml("", "hyperelasticity = \"mooney\"", step)).unwrap_err();
    assert!(bad.to_string().contains("mooney"), "{bad}");
}

#[test]
fn large_strain_needs_a_static_analysis() {
    let toml = format!(
        "{}\n[step]\nkind = \"quasi_static\"\nduration = 1.0\nincrements = 2\n",
        THERMAL.replace(
            "analysis = \"thermal\"",
            "analysis = \"thermal\"\nkinematics = \"large\""
        )
    );
    let msg = error(&toml, &plate()).to_string();
    assert!(msg.contains("static"), "{msg}");
}
