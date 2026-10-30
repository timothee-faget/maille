//! Quel calcul le fichier demande-t-il ? Avec `analysis`, ce qu'elle dit. Sans elle, le contenu décide :
//! des BC ou des chargements thermiques (DOF `t`, `[[heat_sources]]`, `[[convections]]`) donnent un
//! calcul thermique, des BC ou des efforts mécaniques (`ux`, `uy`, `uz`, `[[pressures]]`,
//! `[[edge_loads]]`, `[[contacts]]`, `kinematics`, `[temperatures]`) un calcul statique, et les deux un calcul couplé :
//! `split` coupe alors le fichier en sa partie thermique et sa partie mécanique.

use super::schema::*;
use crate::IoError;

fn invalid(problem: &str) -> IoError {
    IoError::InvalidModel {
        problem: problem.to_string(),
    }
}

/// Calcul déduit d'un fichier de modèle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Physics {
    Static,
    Thermal,
    /// Thermique puis mécanique : `split` donne les deux parties.
    Coupled,
}

/// Physique demandée par le fichier : `analysis` si elle est écrite, le contenu sinon.
pub fn physics(file: &ModelFile) -> Result<Physics, IoError> {
    if let Some(analysis) = file.analysis {
        if file.thermal_step.is_some() {
            return Err(invalid(
                "`[thermal_step]` ne sert qu'au calcul couplé : retirer `analysis`, ou écrire `[step]` pour un calcul thermique seul",
            ));
        }
        return Ok(match analysis {
            AnalysisKind::Static => Physics::Static,
            AnalysisKind::Thermal => Physics::Thermal,
        });
    }
    let thermal_dof = |d: &ValueDef| d.dof == DofDef::T;
    let thermal = file.bcs.iter().any(thermal_dof)
        || file.loads.iter().any(thermal_dof)
        || !file.heat_sources.is_empty()
        || !file.convections.is_empty()
        || file.thermal_step.is_some()
        || file
            .step
            .as_ref()
            .is_some_and(|s| s.kind == StepKind::Transient);
    let mechanical = file.bcs.iter().any(|d| !thermal_dof(d))
        || file.loads.iter().any(|d| !thermal_dof(d))
        || !file.pressures.is_empty()
        || !file.edge_loads.is_empty()
        || !file.contacts.is_empty()
        || file.kinematics.is_some()
        || file.temperatures.is_some()
        || file
            .step
            .as_ref()
            .is_some_and(|s| s.kind == StepKind::QuasiStatic);
    let from = file.temperatures.as_ref().map(|t| t.from.is_some());
    match (thermal, mechanical, from) {
        (false, false, _) => Err(invalid(
            "le modèle ne dit pas quel calcul faire : ajouter une BC (`t` pour le thermique, `ux`, `uy` ou `uz` pour le statique), ou écrire `analysis`",
        )),
        (true, _, Some(true)) => Err(invalid(
            "[temperatures] `from` lit les températures d'un autre calcul, mais ce modèle a sa propre partie thermique : retirer `from` pour les coupler, ou retirer les BC et chargements thermiques",
        )),
        (false, true, Some(false)) => Err(invalid(
            "[temperatures] sans `from` demande une partie thermique dans ce modèle : ajouter une BC `t`, ou donner `from` (le `.pvd` d'un calcul thermique)",
        )),
        (true, true, None) => Err(invalid(
            "ce modèle a une partie thermique et une partie mécanique : ajouter `[temperatures]` avec `reference` pour les coupler (le thermique, puis le mécanique dilaté)",
        )),
        (true, false, _) => Ok(Physics::Thermal),
        (false, true, _) => Ok(Physics::Static),
        (true, true, Some(_)) => Ok(Physics::Coupled),
    }
}

/// Analyse d'un calcul à une seule physique : refuse le calcul couplé, que `split` doit couper.
pub(super) fn analysis_of(file: &ModelFile) -> Result<AnalysisKind, IoError> {
    match physics(file)? {
        Physics::Static => Ok(AnalysisKind::Static),
        Physics::Thermal => Ok(AnalysisKind::Thermal),
        Physics::Coupled => Err(invalid(
            "calcul couplé : le couper avec `model::split`, puis construire chaque partie",
        )),
    }
}

/// Coupe un calcul couplé en `(thermique, mécanique)`, deux fichiers à une seule physique : le
/// thermique avec ses BC en `t`, ses sources, sa convection et `[thermal_step]` pour Step ; le
/// mécanique avec le reste et `[temperatures]` (sans `from`, que l'appelant remplace par les
/// températures du calcul thermique).
pub fn split(file: &ModelFile) -> Result<(ModelFile, ModelFile), IoError> {
    if physics(file)? != Physics::Coupled {
        return Err(invalid(
            "ce modèle n'est pas couplé : le construire directement",
        ));
    }
    let (thermal_step, step) = match (&file.thermal_step, &file.step) {
        (_, Some(s)) if s.kind == StepKind::Transient => {
            return Err(invalid(
                "calcul couplé : le transitoire thermique s'écrit `[thermal_step]` ; `[step]` est celui du calcul mécanique (kind = \"quasi_static\")",
            ))
        }
        (t, s) => (t.clone(), s.clone()),
    };
    let on_t = |d: &&ValueDef| d.dof == DofDef::T;
    let off_t = |d: &&ValueDef| d.dof != DofDef::T;
    let thermal = ModelFile {
        analysis: Some(AnalysisKind::Thermal),
        kinematics: None,
        bcs: file.bcs.iter().filter(on_t).cloned().collect(),
        loads: file.loads.iter().filter(on_t).cloned().collect(),
        pressures: vec![],
        edge_loads: vec![],
        contacts: vec![],
        temperatures: None,
        step: thermal_step.map(|t| StepDef {
            kind: StepKind::Transient,
            duration: t.duration,
            increments: t.increments,
            theta: t.theta,
            initial_temperature: Some(t.initial_temperature),
            amplitude: None,
        }),
        thermal_step: None,
        ..file.clone()
    };
    let mechanical = ModelFile {
        analysis: Some(AnalysisKind::Static),
        bcs: file.bcs.iter().filter(off_t).cloned().collect(),
        loads: file.loads.iter().filter(off_t).cloned().collect(),
        heat_sources: vec![],
        convections: vec![],
        step,
        thermal_step: None,
        ..file.clone()
    };
    Ok((thermal, mechanical))
}
