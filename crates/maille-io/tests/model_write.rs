//! Écriture du fichier de modèle : `parse` → `write` → `parse` redonne le même `ModelFile` pour
//! tous les exemples de `examples/`, et le texte écrit se relit sans clé perdue.

use std::path::PathBuf;

use maille_io::model::{parse, write, ModelFile};

fn examples() -> Vec<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples");
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "toml"))
        .collect();
    files.sort();
    files
}

#[test]
fn round_trip_on_every_example() {
    let files = examples();
    assert!(files.len() >= 10, "exemples introuvables : {files:?}");
    for path in files {
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        let original = parse(&std::fs::read_to_string(&path).unwrap()).unwrap();
        let text = write(&original).unwrap_or_else(|e| panic!("{name} : écriture : {e}"));
        let again = parse(&text).unwrap_or_else(|e| panic!("{name} : relecture : {e}\n{text}"));
        assert_eq!(original, again, "{name} : aller-retour\n{text}");
    }
}

/// Un modèle qui remplit toutes les sections se relit à l'identique, et une clé absente reste absente.
#[test]
fn writes_only_the_keys_that_are_set() {
    let original = parse(
        r#"
units = "N, mm, MPa"
mesh = "a.msh"
analysis = "static"

[materials.acier]
young = 210e3
poisson = 0.3

[[bcs]]
set = "gauche"
dof = "ux"
value = 0.0
"#,
    )
    .unwrap();
    let text = write(&original).unwrap();
    for absent in [
        "kinematics",
        "loads",
        "contacts",
        "step",
        "temperatures",
        "density",
    ] {
        assert!(
            !text.contains(absent),
            "`{absent}` ne devait pas être écrit :\n{text}"
        );
    }
    let again: ModelFile = parse(&text).unwrap();
    assert_eq!(original, again);
}
