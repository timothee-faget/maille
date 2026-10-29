//! Calcul déduit du contenu (`analysis` facultative) et découpe du calcul couplé : un test par cas
//! et par message d'erreur.

use maille_io::model::{parse, physics, split, Physics, StepKind};
use maille_io::IoError;

const HEAD: &str = "units = \"N, mm, MPa, K\"\nmesh = \"a.msh\"\n";

fn file(body: &str) -> maille_io::model::ModelFile {
    parse(&format!("{HEAD}{body}")).unwrap()
}

const THERMAL_BC: &str = "[[bcs]]\nset = \"in\"\ndof = \"t\"\nvalue = 120.0\n";
const MECHANICAL_BC: &str = "[[bcs]]\nset = \"sym\"\ndof = \"ux\"\nvalue = 0.0\n";
const COUPLING: &str = "[temperatures]\nreference = 20.0\n";

fn message(body: &str) -> String {
    match physics(&file(body)) {
        Err(IoError::InvalidModel { problem }) => problem,
        other => panic!("une erreur de modèle était attendue : {other:?}"),
    }
}

#[test]
fn content_decides_the_physics() {
    assert_eq!(physics(&file(THERMAL_BC)).unwrap(), Physics::Thermal);
    assert_eq!(
        physics(&file(
            "[[convections]]\nset = \"o\"\nfilm = 1.0\nsink = 0.0\n"
        ))
        .unwrap(),
        Physics::Thermal
    );
    assert_eq!(physics(&file(MECHANICAL_BC)).unwrap(), Physics::Static);
    assert_eq!(
        physics(&file("[[pressures]]\nset = \"o\"\nvalue = 1.0\n")).unwrap(),
        Physics::Static
    );
    let both = format!("{THERMAL_BC}{MECHANICAL_BC}{COUPLING}");
    assert_eq!(physics(&file(&both)).unwrap(), Physics::Coupled);
}

#[test]
fn analysis_written_wins_over_the_content() {
    let m = file(&format!("analysis = \"thermal\"\n{MECHANICAL_BC}"));
    assert_eq!(physics(&m).unwrap(), Physics::Thermal);
}

#[test]
fn legacy_two_file_coupling_stays_static() {
    let m = file(&format!(
        "{MECHANICAL_BC}[temperatures]\nfrom = \"t.pvd\"\nreference = 20.0\n"
    ));
    assert_eq!(physics(&m).unwrap(), Physics::Static);
}

#[test]
fn diagnostics_say_what_to_fix() {
    assert!(message("").contains("ne dit pas quel calcul"));
    assert!(message(&format!("{THERMAL_BC}{MECHANICAL_BC}")).contains("[temperatures]"));
    assert!(message(&format!("{MECHANICAL_BC}{COUPLING}")).contains("partie thermique"));
    assert!(message(&format!(
        "{THERMAL_BC}[temperatures]\nfrom = \"t.pvd\"\nreference = 20.0\n"
    ))
    .contains("propre partie thermique"));
    let m = file("analysis = \"static\"\n[thermal_step]\nduration = 1.0\nincrements = 2\ninitial_temperature = 20.0\n");
    assert!(
        matches!(physics(&m), Err(IoError::InvalidModel { problem }) if problem.contains("thermal_step"))
    );
}

#[test]
fn split_gives_one_file_per_physics() {
    let m = file(&format!(
        "{THERMAL_BC}{MECHANICAL_BC}{COUPLING}\
         [thermal_step]\nduration = 10.0\nincrements = 5\ninitial_temperature = 20.0\n\
         [step]\nkind = \"quasi_static\"\nduration = 1.0\nincrements = 2\n"
    ));
    let (thermal, mechanical) = split(&m).unwrap();
    assert_eq!(physics(&thermal).unwrap(), Physics::Thermal);
    assert_eq!(physics(&mechanical).unwrap(), Physics::Static);
    assert_eq!(thermal.bcs.len(), 1);
    assert_eq!(mechanical.bcs.len(), 1);
    let t = thermal.step.unwrap();
    assert_eq!((t.kind, t.increments), (StepKind::Transient, 5));
    assert_eq!(t.initial_temperature, Some(20.0));
    assert_eq!(mechanical.step.unwrap().kind, StepKind::QuasiStatic);
    assert!(mechanical.temperatures.is_some() && thermal.temperatures.is_none());
}

#[test]
fn split_refuses_a_transient_step_in_a_coupled_model() {
    let m = file(&format!(
        "{THERMAL_BC}{MECHANICAL_BC}{COUPLING}\
         [step]\nkind = \"transient\"\nduration = 1.0\nincrements = 2\ninitial_temperature = 20.0\n"
    ));
    assert!(
        matches!(split(&m), Err(IoError::InvalidModel { problem }) if problem.contains("[thermal_step]"))
    );
}
