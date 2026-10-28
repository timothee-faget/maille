//! Résultats VTK : aller-retour exact, types de cellules, refus avec message, et couplage
//! thermo-mécanique entièrement par fichiers (`.msh`, `.toml`, `.vtu`, `.pvd`) contre Timoshenko et
//! Goodier.

use std::collections::BTreeMap;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Analysis, Bc, Model};
use maille_core::results::Stress;
use maille_io::model::{build, parse, Run};
use maille_io::msh::{read, write, CellKind, MshFile};
use maille_io::vtk::{
    read_point_data, read_pvd, read_results, result_arrays, write_pvd, write_vtu, Entry, PointArray,
};
use maille_io::IoError;

fn vtk_message(e: IoError) -> String {
    assert!(matches!(e, IoError::Vtk { .. }), "{e:?}");
    e.to_string()
}

fn rectangle(quadratic: bool) -> MshFile {
    let pts = [[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [0.0, 2.0]];
    let names = ["bas", "droite", "haut", "gauche"];
    let geometry = Geometry {
        outer: (0..4).map(|i| Curve::line(pts[i], 1.0, names[i])).collect(),
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default()).unwrap();
    if quadratic {
        MshFile::from_triangulation(&t.to_quadratic()).unwrap()
    } else {
        MshFile::from_triangulation(&t).unwrap()
    }
}

/// Un seul Element, Nodes aux sommets d'un tétraèdre et (Tet10) aux milieux de ses arêtes.
fn one_tet(kind: CellKind) -> MshFile {
    let v = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let mut nodes: Vec<Node> = v.iter().map(|p| Node::new(p[0], p[1], p[2])).collect();
    if kind == CellKind::Tet10 {
        for (a, b) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
            nodes.push(Node::new(
                (v[a][0] + v[b][0]) / 2.0,
                (v[a][1] + v[b][1]) / 2.0,
                (v[a][2] + v[b][2]) / 2.0,
            ));
        }
    }
    MshFile {
        cells: vec![(0..nodes.len()).collect()],
        nodes,
        kind,
        sets: BTreeMap::new(),
    }
}

/// Valeurs numériques du `DataArray` de nom `name` de tout le fichier.
fn cell_array(text: &str, name: &str) -> Vec<usize> {
    let tag = text.find(&format!("Name=\"{name}\"")).unwrap();
    let body = &text[tag..];
    let body = &body[body.find('>').unwrap() + 1..];
    body[..body.find("</DataArray>").unwrap()]
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect()
}

#[test]
fn point_data_round_trips_to_the_last_bit() {
    let mesh = rectangle(true);
    let n = mesh.nodes.len();
    let awkward = [
        0.1 + 0.2,
        1e-300,
        -2.5e17,
        1.0 / 3.0,
        f64::MIN_POSITIVE,
        123_456_789.123_456_79,
    ];
    let scalar = PointArray {
        name: "temperature".to_string(),
        components: 1,
        values: (0..n)
            .map(|i| awkward[i % awkward.len()] * (1.0 + i as f64))
            .collect(),
    };
    let vector = PointArray {
        name: "displacement".to_string(),
        components: 3,
        values: (0..3 * n)
            .map(|i| awkward[i % awkward.len()] - i as f64)
            .collect(),
    };
    let text = write_vtu(
        &mesh,
        "N, mm, MPa",
        false,
        &[scalar.clone(), vector.clone()],
    )
    .unwrap();
    assert_eq!(read_point_data(&text, "temperature").unwrap(), scalar);
    assert_eq!(read_point_data(&text, "displacement").unwrap(), vector);
    assert!(text.contains("<!-- units: N, mm, MPa -->"));
}

#[test]
fn cell_types_and_connectivity_follow_vtk() {
    for (mesh, code, width) in [
        (rectangle(false), 5, 3),
        (rectangle(true), 22, 6),
        (one_tet(CellKind::Tet4), 10, 4),
        (one_tet(CellKind::Tet10), 24, 10),
    ] {
        let text = write_vtu(&mesh, "u", false, &[]).unwrap();
        let cells = mesh.cells.len();
        assert_eq!(cell_array(&text, "types"), vec![code; cells]);
        assert_eq!(
            cell_array(&text, "offsets"),
            (1..=cells).map(|k| k * width).collect::<Vec<_>>()
        );
        assert_eq!(
            cell_array(&text, "connectivity"),
            mesh.cells.concat(),
            "type {code}"
        );
        assert!(text.contains(&format!("NumberOfCells=\"{cells}\"")));
    }
}

#[test]
fn whole_file_round_trips_for_every_cell_kind() {
    for mesh in [
        rectangle(false),
        rectangle(true),
        one_tet(CellKind::Tet4),
        one_tet(CellKind::Tet10),
    ] {
        let t = PointArray {
            name: "temperature".to_string(),
            components: 1,
            values: (0..mesh.nodes.len()).map(|i| 0.1 * i as f64).collect(),
        };
        let text = write_vtu(&mesh, "W, mm, K", false, std::slice::from_ref(&t)).unwrap();
        let back = read_results(&text).unwrap();
        assert_eq!(back.kind, mesh.kind);
        assert_eq!(back.cells, mesh.cells);
        assert_eq!(back.nodes, mesh.nodes);
        assert_eq!(back.arrays, vec![t]);
        assert_eq!(back.units.as_deref(), Some("W, mm, K"));
    }
}

#[test]
fn whole_file_errors() {
    let mesh = rectangle(false);
    let text = write_vtu(&mesh, "u", false, &[]).unwrap();
    let odd = text.replace(
        "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n5",
        "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n9",
    );
    assert!(vtk_message(read_results(&odd).unwrap_err()).contains("type de cellule VTK 9"));
    let marker = "Name=\"connectivity\" format=\"ascii\">\n";
    let at = text.find(marker).unwrap() + marker.len();
    let first = mesh.cells[0][0].to_string().len();
    let outside = format!("{}999999{}", &text[..at], &text[at + first..]);
    assert!(vtk_message(read_results(&outside).unwrap_err()).contains("connectivité"));
    assert!(vtk_message(read_results("<VTKFile/>").unwrap_err()).contains("Points"));
}

#[test]
fn write_refuses_inconsistent_arrays() {
    let mesh = rectangle(false);
    let n = mesh.nodes.len();
    let good = PointArray {
        name: "t".to_string(),
        components: 1,
        values: vec![1.0; n],
    };
    let short = PointArray {
        values: vec![1.0; n - 1],
        ..good.clone()
    };
    assert!(vtk_message(write_vtu(&mesh, "u", false, &[short]).unwrap_err()).contains("valeurs"));
    let mut values = vec![1.0; n];
    values[3] = f64::NAN;
    let nan = PointArray {
        values,
        ..good.clone()
    };
    assert!(vtk_message(write_vtu(&mesh, "u", false, &[nan]).unwrap_err()).contains("Node 3"));
    let named = PointArray {
        name: "a\"b".to_string(),
        ..good
    };
    assert!(vtk_message(write_vtu(&mesh, "u", false, &[named]).unwrap_err()).contains("nom"));
    assert!(vtk_message(write_vtu(&mesh, "a -- b", false, &[]).unwrap_err()).contains("units"));
}

#[test]
fn read_refuses_binary_and_names_what_exists() {
    let mesh = rectangle(false);
    let t = PointArray {
        name: "temperature".to_string(),
        components: 1,
        values: vec![1.0; mesh.nodes.len()],
    };
    let text = write_vtu(&mesh, "u", false, &[t]).unwrap();
    let missing = vtk_message(read_point_data(&text, "pression").unwrap_err());
    assert!(
        missing.contains("« pression »") && missing.contains("« temperature »"),
        "{missing}"
    );
    let binary = text.replacen("format=\"ascii\"", "format=\"binary\"", 1);
    assert!(vtk_message(read_point_data(&binary, "temperature").unwrap_err()).contains("ASCII"));
    let none = write_vtu(&mesh, "u", false, &[]).unwrap();
    assert!(vtk_message(read_point_data(&none, "temperature").unwrap_err()).contains("aucune"));
    assert!(vtk_message(read_point_data("<VTKFile/>", "t").unwrap_err()).contains("PointData"));
}

#[test]
fn pvd_round_trips_and_increments_are_picked() {
    let entries = vec![
        Entry {
            time: 0.0,
            file: "r_0000.vtu".to_string(),
        },
        Entry {
            time: 0.05,
            file: "r_0001.vtu".to_string(),
        },
        Entry {
            time: 10.0,
            file: "r_0002.vtu".to_string(),
        },
    ];
    let text = write_pvd(&entries).unwrap();
    assert_eq!(read_pvd(&text).unwrap(), entries);

    let pick = |increment: &str| {
        let toml = format!(
            "units = \"N\"\nmesh = \"m.msh\"\nanalysis = \"static\"\n[temperatures]\nfrom = \"r.pvd\"\nreference = 20.0\n{increment}\n"
        );
        let file = parse(&toml).unwrap();
        file.temperatures
            .unwrap()
            .pick(&entries)
            .map(|e| e.file.clone())
    };
    assert_eq!(pick("").unwrap(), "r_0002.vtu");
    assert_eq!(pick("increment = \"last\"").unwrap(), "r_0002.vtu");
    assert_eq!(pick("increment = 1").unwrap(), "r_0001.vtu");
    assert!(pick("increment = 3")
        .unwrap_err()
        .to_string()
        .contains("numérotés depuis 0"));
    assert!(pick("increment = \"premier\"")
        .unwrap_err()
        .to_string()
        .contains("last"));
}

#[test]
fn pvd_errors() {
    assert!(vtk_message(write_pvd(&[]).unwrap_err()).contains("sans fichier"));
    assert!(
        vtk_message(read_pvd("<VTKFile><Collection/></VTKFile>").unwrap_err()).contains("DataSet")
    );
    let bad = "<DataSet timestep=\"x\" file=\"a.vtu\"/>";
    assert!(vtk_message(read_pvd(bad).unwrap_err()).contains("timestep"));
}

#[test]
fn result_arrays_list_what_each_analysis_writes() {
    let mesh = rectangle(false);
    let toml = |analysis: &str, extra: &str| {
        format!(
            "units = \"u\"\nmesh = \"m.msh\"\nanalysis = \"{analysis}\"\n[materials.m]\n{extra}\n[[sections]]\nset = \"domaine\"\nmaterial = \"m\"\nplane_state = \"strain\"\n[[bcs]]\nset = \"gauche\"\ndof = \"{}\"\nvalue = 0.0\n[[bcs]]\nset = \"droite\"\ndof = \"{}\"\nvalue = 1.0\n",
            if analysis == "thermal" { "t" } else { "ux" },
            if analysis == "thermal" { "t" } else { "ux" }
        ) + if analysis == "static" {
            "[[bcs]]\nset = \"bas\"\ndof = \"uy\"\nvalue = 0.0\n"
        } else {
            ""
        }
    };
    let names = |arrays: Vec<PointArray>| arrays.into_iter().map(|a| a.name).collect::<Vec<_>>();

    let thermal = build(
        &parse(&toml("thermal", "conductivity = 1.0")).unwrap(),
        &mesh,
        None,
    )
    .unwrap();
    let field = thermal.model.solve().unwrap();
    assert_eq!(
        names(result_arrays(&thermal.mesh, &field).unwrap()),
        ["temperature", "reaction"]
    );

    let elastic = build(
        &parse(&toml("static", "young = 100.0\npoisson = 0.3")).unwrap(),
        &mesh,
        None,
    )
    .unwrap();
    let field = elastic.model.solve().unwrap();
    let arrays = result_arrays(&elastic.mesh, &field).unwrap();
    assert_eq!(
        names(arrays.clone()),
        ["displacement", "reaction", "stress", "von_mises"]
    );
    assert_eq!(arrays[0].components, 3);
    // allongement ε = 1/4 imposé : ux = x / 4 sur chaque Node
    for (p, u) in mesh.nodes.iter().zip(arrays[0].values.chunks(3)) {
        assert!((u[0] - p.x / 4.0).abs() < 1e-12 && u[2] == 0.0);
    }

    let plastic = toml("static", "young = 100.0\npoisson = 0.3\nyield_stress = 1.0")
        + "[step]\nkind = \"quasi_static\"\nduration = 1.0\nincrements = 2\n";
    let loaded = build(&parse(&plastic).unwrap(), &mesh, None).unwrap();
    let Run::QuasiStatic { step, amplitude } = loaded.run.clone() else {
        panic!("quasi-statique attendu");
    };
    let history = loaded.model.solve_quasi_static(&step, &amplitude).unwrap();
    let arrays = result_arrays(&loaded.mesh, &history[2].field).unwrap();
    assert_eq!(arrays.last().unwrap().name, "plastic_strain");
    assert!(arrays.last().unwrap().values.iter().any(|&p| p > 0.0));
}

// --- couplage par fichiers : tube épais chauffé de l'intérieur ---

const YOUNG: f64 = 210e3;
const POISSON: f64 = 0.3;
const EXPANSION: f64 = 1.2e-5;
const CONDUCTIVITY: f64 = 0.05;
const INNER: f64 = 10.0;
const OUTER: f64 = 20.0;
const T_INNER: f64 = 120.0;
const T_OUTER: f64 = 20.0;
const T_REF: f64 = 20.0;

fn quarter_tube() -> MshFile {
    let size = 1.0;
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
    read(&write(&MshFile::from_triangulation(&t).unwrap()).unwrap()).unwrap()
}

fn beta() -> f64 {
    (T_OUTER - T_INNER) / (OUTER / INNER).ln()
}

fn delta(r: f64) -> f64 {
    T_INNER - T_REF + beta() * (r / INNER).ln()
}

fn integral(r: f64) -> f64 {
    (T_INNER - T_REF) * (r * r - INNER * INNER) / 2.0
        + beta() * (r * r / 2.0 * (r / INNER).ln() - r * r / 4.0 + INNER * INNER / 4.0)
}

/// (σr, σθ, σz) de Timoshenko et Goodier en r.
fn exact(r: f64) -> [f64; 3] {
    let c = EXPANSION * YOUNG / ((1.0 - POISSON) * r * r);
    let span = OUTER * OUTER - INNER * INNER;
    let sr = c * ((r * r - INNER * INNER) / span * integral(OUTER) - integral(r));
    let st =
        c * ((r * r + INNER * INNER) / span * integral(OUTER) + integral(r) - delta(r) * r * r);
    [sr, st, POISSON * (sr + st) - YOUNG * EXPANSION * delta(r)]
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

const THERMAL_TUBE: &str = r#"
units = "W, mm, K"
mesh = "tube.msh"
analysis = "thermal"

[materials.acier]
conductivity = 0.05

[[sections]]
set = "domaine"
material = "acier"

[[bcs]]
set = "interieur"
dof = "t"
value = 120.0

[[bcs]]
set = "exterieur"
dof = "t"
value = 20.0
"#;

const STATIC_TUBE: &str = r#"
units = "N, mm, MPa"
mesh = "tube.msh"
analysis = "static"

[materials.acier]
young = 210e3
poisson = 0.3
expansion = 1.2e-5

[[sections]]
set = "domaine"
material = "acier"
plane_state = "strain"

[[bcs]]
set = "sym_x"
dof = "uy"
value = 0.0

[[bcs]]
set = "sym_y"
dof = "ux"
value = 0.0

[temperatures]
from = "tube.pvd"
increment = "last"
reference = 20.0
"#;

#[test]
fn thermal_then_static_through_files_follows_timoshenko() {
    let mesh = quarter_tube();

    // 1. calcul thermique, résultats écrits
    let thermal = build(&parse(THERMAL_TUBE).unwrap(), &mesh, None).unwrap();
    assert_eq!(thermal.run, Run::Linear);
    let field = thermal.model.solve().unwrap();
    let arrays = result_arrays(&thermal.mesh, &field).unwrap();
    let vtu = write_vtu(&mesh, &thermal.units, false, &arrays).unwrap();
    let pvd = write_pvd(&[Entry {
        time: 0.0,
        file: "tube_0000.vtu".to_string(),
    }])
    .unwrap();

    // 2. le modèle statique relit ces fichiers
    let file = parse(STATIC_TUBE).unwrap();
    let def = file.temperatures.as_ref().unwrap();
    let entries = read_pvd(&pvd).unwrap();
    assert_eq!(def.pick(&entries).unwrap().file, "tube_0000.vtu");
    let temperatures = read_point_data(&vtu, "temperature").unwrap();
    assert_eq!(temperatures.values, field.temperatures);
    let loaded = build(&file, &mesh, Some(&temperatures.values)).unwrap();
    let mechanical = loaded.model.solve().unwrap();

    // 3. même résultat que le cœur appelé directement
    let material = Material::with_poisson(YOUNG, POISSON)
        .unwrap()
        .with_expansion(EXPANSION)
        .unwrap()
        .with_conductivity(CONDUCTIVITY)
        .unwrap();
    let elements = mesh
        .cells
        .iter()
        .map(|c| {
            Tri6::new(
                c.clone().try_into().unwrap(),
                material,
                1.0,
                PlaneState::Strain,
            )
            .map(Element::Tri6)
        })
        .collect::<Result<_, _>>()
        .unwrap();
    let core_mesh = Mesh::new(mesh.nodes.clone(), elements).unwrap();
    let bc = |set: &str, dof, value| -> Vec<Bc> {
        mesh.sets[set]
            .nodes
            .iter()
            .map(|&node| Bc { node, dof, value })
            .collect()
    };
    let direct_thermal = Model::with_analysis(
        core_mesh.clone(),
        Analysis::Thermal,
        [
            bc("interieur", Dof::T, T_INNER),
            bc("exterieur", Dof::T, T_OUTER),
        ]
        .concat(),
        vec![],
    )
    .unwrap()
    .solve()
    .unwrap();
    let direct = Model::new(
        core_mesh,
        [bc("sym_x", Dof::Uy, 0.0), bc("sym_y", Dof::Ux, 0.0)].concat(),
        vec![],
    )
    .unwrap()
    .with_temperatures(direct_thermal.temperatures, T_REF)
    .unwrap()
    .solve()
    .unwrap();
    for (a, b) in mechanical.displacements.iter().zip(&direct.displacements) {
        assert!((a - b).abs() < 1e-12);
    }

    // 4. Timoshenko et Goodier, mêmes seuils que `maille-cli cylindre-thermo`
    let scale = exact(INNER)[1].abs();
    let (mut sum, mut count) = (0.0, 0);
    for s in mechanical.stresses.iter().flatten() {
        let (got, want) = (polar(s), exact(s.x.hypot(s.y)));
        for k in 0..3 {
            let error = (got[k] - want[k]).abs();
            assert!(error < 2e-2 * scale, "({}, {}) composante {k}", s.x, s.y);
            sum += error * error;
            count += 1;
        }
    }
    let rms = (sum / count as f64).sqrt() / scale;
    assert!(rms < 2e-3, "erreur RMS {rms:.2e}");
}

#[test]
fn coupling_with_a_foreign_mesh_is_refused() {
    let mesh = quarter_tube();
    let file = parse(STATIC_TUBE).unwrap();
    let e = build(&file, &mesh, Some(&[20.0; 3])).unwrap_err();
    assert!(matches!(e, IoError::Model(_)), "{e:?}");
    assert!(e.to_string().contains("3"), "{e}");
}

#[test]
fn large_strain_is_marked_in_the_file_and_read_back() {
    let mesh = rectangle(false);
    let plain = write_vtu(&mesh, "u", false, &[]).unwrap();
    assert!(!read_results(&plain).unwrap().large_strain);
    let large = write_vtu(&mesh, "u", true, &[]).unwrap();
    assert!(read_results(&large).unwrap().large_strain);
    // le commentaire n'empêche pas de relire les unités
    assert_eq!(read_results(&large).unwrap().units.as_deref(), Some("u"));
}
