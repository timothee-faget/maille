//! Bibliothèque de Materials : la livrée est valide, la conversion d'unités garde la physique (la
//! diffusivité, la plaque en traction), les erreurs disent quoi corriger.

use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_io::materials::{Library, UnitSystem};
use maille_io::model::{build, parse};
use maille_io::msh::{read, write, MshFile};
use maille_io::IoError;

/// Écart relatif.
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}

#[test]
fn builtin_is_valid() {
    let library = Library::builtin();
    let names: Vec<&str> = library.names().collect();
    assert!(names.contains(&"acier") && names.contains(&"aluminium"));
    for name in names {
        let m = library.get(name).unwrap();
        assert!(!m.description.is_empty(), "{name} sans description");
        assert!(m.young.is_some() && m.poisson.is_some() && m.density.is_some());
    }
}

#[test]
fn steel_in_the_usual_systems() {
    let library = Library::builtin();
    let si = library.material("acier", UnitSystem::SI).unwrap();
    assert_eq!(si.young, Some(210e9));
    assert_eq!(si.density, Some(7850.0));

    // mm, N, s : MPa, tonnes/mm³, mJ/(t·K), N/(s·K)
    let mm = library.material("acier", UnitSystem::MM_N_S).unwrap();
    assert!(close(mm.young.unwrap(), 210e3));
    // la plasticité se demande à part : le Material converti reste élastique
    assert_eq!(mm.yield_stress, None);
    let steel = library.get("acier").unwrap();
    assert!(close(
        steel.yield_stress_in(UnitSystem::MM_N_S).unwrap(),
        235.0
    ));
    assert_eq!(steel.yield_stress_in(UnitSystem::SI), Some(235e6));
    assert!(close(mm.density.unwrap(), 7.85e-9));
    assert!(close(mm.specific_heat.unwrap(), 4.6e8));
    assert!(close(mm.conductivity.unwrap(), 45.0));
    assert_eq!(mm.poisson, Some(0.3));
    assert_eq!(mm.expansion, Some(1.2e-5));
    assert_eq!(mm.hardening, None);
}

#[test]
fn thermal_diffusivity_keeps_its_dimension_in_every_system() {
    // a = k / (ρ c) est une longueur² / temps : en SI puis dans un système, a_système = a_SI T / L²
    let library = Library::builtin();
    for text in ["m, N, s", "mm, N, s", "cm, kN, s", "in, lbf, s", "m, kN, s"] {
        let units = UnitSystem::parse(text).unwrap();
        let (length, time) = match text {
            "m, N, s" | "m, kN, s" => (1.0, 1.0),
            "mm, N, s" => (1e-3, 1.0),
            "cm, kN, s" => (1e-2, 1.0),
            _ => (0.0254, 1.0),
        };
        for name in library.names() {
            let si = library.material(name, UnitSystem::SI).unwrap();
            let m = library.material(name, units).unwrap();
            let a_si = si.conductivity.unwrap() / (si.density.unwrap() * si.specific_heat.unwrap());
            let a = m.conductivity.unwrap() / (m.density.unwrap() * m.specific_heat.unwrap());
            assert!(
                close(a, a_si * time / (length * length)),
                "{name} en {text}"
            );
        }
    }
}

#[test]
fn units_are_parsed_and_refused_with_a_message() {
    let units = UnitSystem::parse(" mm , N , s ").unwrap();
    assert_eq!(units, UnitSystem::MM_N_S);

    let err = UnitSystem::parse("mm, N").unwrap_err();
    assert!(err.to_string().contains("longueur, force, temps"), "{err}");
    let err = UnitSystem::parse("pouce, N, s").unwrap_err();
    assert!(err.to_string().contains("pouce") && err.to_string().contains("m, cm, mm, in"));
    assert!(UnitSystem::new(0.0, 1.0, 1.0).is_err());
    assert!(UnitSystem::new(1.0, f64::NAN, 1.0).is_err());
}

#[test]
fn user_library_overrides_and_unknown_names_list_the_known_ones() {
    let user = Library::parse(
        r#"
[materials.acier]
description = "Mon acier"
young = 200e9
poisson = 0.3

[materials.maison]
young = 1e9
poisson = 0.4
"#,
    )
    .unwrap();
    let library = Library::builtin().merged(user);
    let steel = library.material("acier", UnitSystem::SI).unwrap();
    assert_eq!(steel.young, Some(200e9));
    assert_eq!(steel.density, None);
    assert!(library.get("maison").is_some() && library.get("cuivre").is_some());

    match library.material("mithril", UnitSystem::SI).unwrap_err() {
        IoError::UnknownMaterial { name, known, .. } => {
            assert_eq!(name, "mithril");
            assert!(known.contains(&"maison".to_string()));
        }
        e => panic!("{e}"),
    }
}

#[test]
fn absurd_values_and_unknown_keys_are_refused() {
    let err = Library::parse("[materials.x]\nyoung = -1.0\n").unwrap_err();
    assert!(err.to_string().contains("« x » : young = -1"), "{err}");
    let err = Library::parse("[materials.x]\npoisson = 0.5\n").unwrap_err();
    assert!(err.to_string().contains("poisson"), "{err}");
    let err = Library::parse("[materials.x]\nyoungg = 1.0\n").unwrap_err();
    assert!(matches!(err, IoError::Toml(_)), "{err}");
}

/// Plaque `w` × `h` en traction (ε = 0,1 %) avec le Material `acier` de la bibliothèque, dans le
/// système `units` ; longueurs et épaisseur en unités de ce système (`scale` = taille d'une unité
/// du système de référence mm). Rend la contrainte σx au premier triangle.
fn plate_stress(units: UnitSystem, scale: f64) -> f64 {
    let (w, h) = (4.0 * scale, 2.0 * scale);
    let pts = [[0.0, 0.0], [w, 0.0], [w, h], [0.0, h]];
    let names = ["bas", "droite", "haut", "gauche"];
    let geometry = Geometry {
        outer: (0..4)
            .map(|i| Curve::line(pts[i], 0.7 * scale, names[i]))
            .collect(),
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default()).unwrap();
    let mesh =
        read(&write(&MshFile::from_triangulation(&t.to_quadratic()).unwrap()).unwrap()).unwrap();

    let mut file = parse(&format!(
        r#"
units = "cohérent"
mesh = "plaque.msh"
analysis = "static"

[[sections]]
set = "domaine"
material = "acier"
thickness = {thickness}
plane_state = "stress"

[[bcs]]
set = "gauche"
dof = "ux"
value = 0.0

[[bcs]]
set = "droite"
dof = "ux"
value = {shift}

[[bcs]]
set = "bas"
dof = "uy"
value = 0.0
"#,
        thickness = 2.0 * scale,
        shift = 0.004 * scale,
    ))
    .unwrap();
    file.materials.insert(
        "acier".to_string(),
        Library::builtin().material("acier", units).unwrap(),
    );
    let field = build(&file, &mesh, None).unwrap().model.solve().unwrap();

    // réaction sur x = w, divisée par l'aire de la section : σx moyenne
    let force: f64 = mesh.sets["droite"]
        .nodes
        .iter()
        .map(|&n| field.reactions[2 * n])
        .sum();
    force / (h * 2.0 * scale)
}

#[test]
fn tensile_plate_gives_the_same_stress_in_si_and_in_mm() {
    // σx = E ε = 210 MPa, quel que soit le système d'unités
    let mm = plate_stress(UnitSystem::MM_N_S, 1.0);
    assert!((mm - 210.0).abs() < 1e-6 * 210.0, "{mm} MPa");
    let si = plate_stress(UnitSystem::SI, 1e-3);
    assert!((si - 210e6).abs() < 1e-6 * 210e6, "{si} Pa");
    assert!((si / 1e6 - mm).abs() < 1e-9 * mm);
}
