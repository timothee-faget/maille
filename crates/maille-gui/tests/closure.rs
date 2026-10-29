//! Clôture de M13 : le parcours de l'interface (DXF, maillage, modèle, calcul) donne le même `.vtu` que
//! la chaîne de la ligne de commande (`maille mesh`, modèle écrit à la main, `maille solve`).

use std::fs;
use std::path::{Path, PathBuf};

use maille_core::model::Flow;
use maille_gui::mesh_space::{self, Imported};
use maille_gui::model_space::{diagnose, Draft, Severity};
use maille_io::meshing::{self, Recipe, Source};
use maille_io::model::{self, DofDef, ValueDef};
use maille_io::{msh, solving};

fn examples() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples")
}

fn sandbox(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("maille-closure-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn solve(model: &Path) -> Vec<PathBuf> {
    solving::solve_file(model, None, &mut |_| Flow::Continue, &mut |_| {}).unwrap()
}

/// La plaque en traction de `examples/traction.toml`, de bout en bout.
#[test]
fn the_interface_chain_gives_the_same_vtu_as_the_command_line_chain() {
    let dxf = examples().join("plaque.dxf");
    let size = 0.5;

    // chaîne de l'interface : importer, mailler, remplir le modèle, enregistrer, calculer
    let gui = sandbox("gui");
    let Imported::Cad(cad) = mesh_space::import(&dxf).unwrap() else {
        panic!("plaque.dxf est un dessin");
    };
    let built = mesh_space::mesh(&cad, size, false, &Default::default()).unwrap();
    let mut draft = Draft::default();
    draft.material = Some(draft.steel());
    draft.thickness = 2.0;
    let bc = |set: &str, value| ValueDef {
        set: set.into(),
        dof: DofDef::Ux,
        value,
    };
    draft.bcs = vec![
        bc("gauche", 0.0),
        bc("droite", 0.004),
        ValueDef {
            set: "bas".into(),
            dof: DofDef::Uy,
            value: 0.0,
        },
    ];
    let diagnostics = diagnose(&draft, &built.file, &[]);
    assert!(
        diagnostics.iter().all(|d| d.severity != Severity::Error),
        "{diagnostics:?}"
    );
    fs::write(gui.join("plaque.msh"), msh::write(&built.file).unwrap()).unwrap();
    fs::write(
        gui.join("plaque_modele.toml"),
        draft.to_toml("plaque.msh", &built.file, None).unwrap(),
    )
    .unwrap();
    let from_gui = solve(&gui.join("plaque_modele.toml"));

    // chaîne de la ligne de commande : `maille mesh plaque.dxf --size 0.5`, puis le modèle de l'exemple
    let cli = sandbox("cli");
    let recipe = Recipe {
        size: Some(size),
        ..Recipe::default()
    };
    let meshed =
        meshing::mesh(Source::of(&dxf).unwrap(), &fs::read(&dxf).unwrap(), &recipe).unwrap();
    fs::write(cli.join("plaque.msh"), msh::write(&meshed.file).unwrap()).unwrap();
    fs::copy(examples().join("traction.toml"), cli.join("traction.toml")).unwrap();
    let from_cli = solve(&cli.join("traction.toml"));

    // le maillage est le même, octet pour octet, donc les résultats aussi
    assert_eq!(
        fs::read(gui.join("plaque.msh")).unwrap(),
        fs::read(cli.join("plaque.msh")).unwrap()
    );
    let vtu = |files: &[PathBuf]| -> Vec<u8> {
        let path = files
            .iter()
            .find(|p| p.extension().is_some_and(|e| e == "vtu"));
        fs::read(path.expect("un .vtu")).unwrap()
    };
    assert_eq!(vtu(&from_gui), vtu(&from_cli));

    // et le modèle écrit par l'interface est celui de l'exemple (à l'`analysis` près, qu'il déduit)
    let mut example =
        model::parse(&fs::read_to_string(cli.join("traction.toml")).unwrap()).unwrap();
    example.analysis = None;
    let written =
        model::parse(&fs::read_to_string(gui.join("plaque_modele.toml")).unwrap()).unwrap();
    assert_eq!(written, example);
}
