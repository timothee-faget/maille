//! Construction du `Model` du cœur à partir du fichier de modèle et du maillage.

use std::collections::BTreeMap;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::model::{
    Amplitude, Analysis, Bc, Convection, HeatSource, Load, Model, Pressure, Step,
};

use super::schema::*;
use super::{Loaded, Run};
use crate::msh::{CellKind, MshFile, Set};
use crate::IoError;

fn invalid(problem: String) -> IoError {
    IoError::InvalidModel { problem }
}

/// Construit le `Model` et le calcul demandé. Les erreurs disent quel Set ou quelle clé corriger.
///
/// `temperatures` : une valeur par Node, lue par le CLI dans le fichier de `[temperatures]` ; elle
/// est obligatoire si le modèle a cette section, et refusée sinon.
pub fn build(
    file: &ModelFile,
    mesh: &MshFile,
    temperatures: Option<&[f64]>,
) -> Result<Loaded, IoError> {
    if file.units.trim().is_empty() {
        return Err(invalid(
            "`units` est vide : écrire le système d'unités cohérent utilisé, par exemple \"N, mm, MPa\""
                .to_string(),
        ));
    }
    match (&file.temperatures, temperatures) {
        (Some(def), None) => {
            return Err(invalid(format!(
                "[temperatures] demande les températures de {} : les lire (vtk::read_point_data) et les passer à build",
                def.from
            )))
        }
        (None, Some(_)) => {
            return Err(invalid(
                "des températures sont fournies mais le modèle n'a pas de [temperatures] : ajouter la section ou ne rien fournir"
                    .to_string(),
            ))
        }
        _ => {}
    }
    let analysis = match file.analysis {
        AnalysisKind::Static => Analysis::Static,
        AnalysisKind::Thermal => Analysis::Thermal,
    };

    let materials = file
        .materials
        .iter()
        .map(|(name, def)| Ok((name.as_str(), material(name, def)?)))
        .collect::<Result<BTreeMap<_, _>, IoError>>()?;
    let core_mesh = Mesh::new(
        mesh.nodes.clone(),
        elements(file, mesh, analysis, &materials)?,
    )?;

    let dof = |key: &str, d: DofDef| node_dof(key, d, analysis, mesh.kind);
    let mut bcs = Vec::new();
    for (i, def) in file.bcs.iter().enumerate() {
        let key = format!("bcs[{i}]");
        let dof = dof(&key, def.dof)?;
        for &node in &node_set(mesh, &key, &def.set)?.nodes {
            bcs.push(Bc {
                node,
                dof,
                value: def.value,
            });
        }
    }
    let mut loads = Vec::new();
    for (i, def) in file.loads.iter().enumerate() {
        let key = format!("loads[{i}]");
        let dof = dof(&key, def.dof)?;
        for &node in &node_set(mesh, &key, &def.set)?.nodes {
            loads.push(Load {
                node,
                dof,
                value: def.value,
            });
        }
    }
    let mut model = Model::with_analysis(core_mesh.clone(), analysis, bcs, loads)?;
    if let (Some(def), Some(t)) = (&file.temperatures, temperatures) {
        model = model.with_temperatures(t.to_vec(), def.reference)?;
    }

    if !file.heat_sources.is_empty() {
        let mut sources = Vec::new();
        for (i, def) in file.heat_sources.iter().enumerate() {
            let set = find(mesh, &format!("heat_sources[{i}]"), &def.set)?;
            expect(!set.elements.is_empty(), &def.set, "d'Elements")?;
            sources.extend(set.elements.iter().map(|&element| HeatSource {
                element,
                value: def.value,
            }));
        }
        model = model.with_heat_sources(sources)?;
    }
    if !file.convections.is_empty() {
        let mut convection = Vec::new();
        for (i, def) in file.convections.iter().enumerate() {
            let set = find(mesh, &format!("convections[{i}]"), &def.set)?;
            expect(!set.faces.is_empty(), &def.set, "de faces (un bord)")?;
            convection.extend(set.faces.iter().map(|&(element, face)| Convection {
                element,
                face,
                film: def.film,
                sink: def.sink,
            }));
        }
        model = model.with_convection(convection)?;
    }

    if !file.pressures.is_empty() {
        let mut pressures = Vec::new();
        for (i, def) in file.pressures.iter().enumerate() {
            let set = find(mesh, &format!("pressures[{i}]"), &def.set)?;
            expect(!set.faces.is_empty(), &def.set, "de faces (un bord)")?;
            pressures.extend(set.faces.iter().map(|&(element, face)| Pressure {
                element,
                face,
                value: def.value,
            }));
        }
        model = model.with_pressures(pressures)?;
    }

    let run = match &file.step {
        None => Run::Linear,
        Some(def) => step(def, file.analysis, mesh.nodes.len())?,
    };
    Ok(Loaded {
        units: file.units.clone(),
        mesh: core_mesh,
        model,
        run,
    })
}

fn find<'a>(mesh: &'a MshFile, key: &str, name: &str) -> Result<&'a Set, IoError> {
    mesh.sets.get(name).ok_or_else(|| IoError::UnknownSet {
        key: key.to_string(),
        name: name.to_string(),
        known: mesh.sets.keys().cloned().collect(),
    })
}

fn expect(ok: bool, name: &str, expected: &'static str) -> Result<(), IoError> {
    if ok {
        Ok(())
    } else {
        Err(IoError::WrongSetKind {
            name: name.to_string(),
            expected,
        })
    }
}

/// Set utilisable pour des BC et des Load : il a des Nodes.
fn node_set<'a>(mesh: &'a MshFile, key: &str, name: &str) -> Result<&'a Set, IoError> {
    let set = find(mesh, key, name)?;
    expect(
        !set.nodes.is_empty(),
        name,
        "de Nodes (un bord ou un point)",
    )?;
    Ok(set)
}

fn node_dof(key: &str, dof: DofDef, analysis: Analysis, kind: CellKind) -> Result<Dof, IoError> {
    match (analysis, dof) {
        (Analysis::Static, DofDef::Ux) => Ok(Dof::Ux),
        (Analysis::Static, DofDef::Uy) => Ok(Dof::Uy),
        (Analysis::Static, DofDef::Uz) if kind.dim() == 3 => Ok(Dof::Uz),
        (Analysis::Static, DofDef::Uz) => Err(invalid(format!(
            "{key} : `uz` n'existe pas dans un maillage 2D, utiliser `ux` ou `uy`"
        ))),
        (Analysis::Thermal, DofDef::T) => Ok(Dof::T),
        (Analysis::Static, DofDef::T) => Err(invalid(format!(
            "{key} : `t` est la température, réservée à l'analyse thermique : en statique, utiliser `ux`, `uy` ou `uz`"
        ))),
        (Analysis::Thermal, _) => Err(invalid(format!(
            "{key} : en analyse thermique, le seul DOF est `t`"
        ))),
    }
}

fn material(name: &str, def: &MaterialDef) -> Result<Material, IoError> {
    let mut m = match (def.young, def.poisson, def.conductivity) {
        (Some(e), Some(nu), _) => Material::with_poisson(e, nu)?,
        (Some(e), None, _) => Material::new(e)?,
        (None, Some(_), _) => {
            return Err(invalid(format!(
                "matériau « {name} » : `poisson` sans `young`"
            )))
        }
        (None, None, Some(k)) => Material::conductive(k)?,
        (None, None, None) => {
            return Err(invalid(format!(
                "matériau « {name} » : donner `young` (statique) ou `conductivity` (thermique)"
            )))
        }
    };
    if let (Some(_), Some(k)) = (def.young, def.conductivity) {
        m = m.with_conductivity(k)?;
    }
    if let Some(v) = def.density {
        m = m.with_density(v)?;
    }
    if let Some(v) = def.specific_heat {
        m = m.with_specific_heat(v)?;
    }
    if let Some(v) = def.expansion {
        m = m.with_expansion(v)?;
    }
    match (def.yield_stress, def.hardening) {
        (Some(sy), h) => m = m.with_plasticity(sy, h.unwrap_or(0.0))?,
        (None, Some(_)) => {
            return Err(invalid(format!(
                "matériau « {name} » : `hardening` sans `yield_stress`"
            )))
        }
        (None, None) => {}
    }
    Ok(m)
}

fn array<const N: usize>(nodes: &[usize]) -> Result<[usize; N], IoError> {
    nodes
        .try_into()
        .map_err(|_| invalid(format!("un Element a {} Nodes, {N} attendus", nodes.len())))
}

/// Un Element par cellule, avec le Material de sa section.
fn elements(
    file: &ModelFile,
    mesh: &MshFile,
    analysis: Analysis,
    materials: &BTreeMap<&str, Material>,
) -> Result<Vec<Element>, IoError> {
    let triangles = mesh.kind.dim() == 2;
    let mut section_of: Vec<Option<usize>> = vec![None; mesh.cells.len()];
    for (i, def) in file.sections.iter().enumerate() {
        let key = format!("sections[{i}]");
        let set = find(mesh, &key, &def.set)?;
        expect(
            !set.elements.is_empty(),
            &def.set,
            "d'Elements (le domaine)",
        )?;
        if !materials.contains_key(def.material.as_str()) {
            return Err(IoError::UnknownMaterial {
                key,
                name: def.material.clone(),
                known: materials.keys().map(|k| k.to_string()).collect(),
            });
        }
        if !triangles && (def.thickness.is_some() || def.plane_state.is_some()) {
            return Err(invalid(format!(
                "{key} : `thickness` et `plane_state` ne servent qu'aux triangles : les retirer pour des tétraèdres"
            )));
        }
        if triangles && analysis == Analysis::Static && def.plane_state.is_none() {
            return Err(invalid(format!(
                "{key} : `plane_state` est obligatoire pour des triangles en statique : \"stress\" (plaque mince) ou \"strain\" (corps long)"
            )));
        }
        for &e in &set.elements {
            if let Some(first) = section_of[e].replace(i) {
                return Err(IoError::OverlappingSections {
                    element: e,
                    first,
                    second: i,
                });
            }
        }
    }
    if let Some(first) = section_of.iter().position(Option::is_none) {
        return Err(IoError::UncoveredElements {
            count: section_of.iter().filter(|s| s.is_none()).count(),
            first,
        });
    }

    mesh.cells
        .iter()
        .zip(&section_of)
        .map(|(nodes, section)| {
            let def = &file.sections[section.unwrap_or(0)];
            let m = materials[def.material.as_str()];
            let thickness = def.thickness.unwrap_or(1.0);
            let state = match def.plane_state {
                Some(PlaneStateDef::Strain) => PlaneState::Strain,
                _ => PlaneState::Stress,
            };
            Ok(match mesh.kind {
                CellKind::Tri3 => Element::Tri3(Tri3::new(array(nodes)?, m, thickness, state)?),
                CellKind::Tri6 => Element::Tri6(Tri6::new(array(nodes)?, m, thickness, state)?),
                CellKind::Tet4 => Element::Tet4(Tet4::new(array(nodes)?, m)),
                CellKind::Tet10 => Element::Tet10(Tet10::new(array(nodes)?, m)),
            })
        })
        .collect()
}

fn step(def: &StepDef, analysis: AnalysisKind, nodes: usize) -> Result<Run, IoError> {
    match def.kind {
        StepKind::Transient => {
            if analysis != AnalysisKind::Thermal {
                return Err(invalid(
                    "[step] kind = \"transient\" demande analysis = \"thermal\"".to_string(),
                ));
            }
            if def.amplitude.is_some() {
                return Err(invalid(
                    "[step] `amplitude` ne sert qu'à kind = \"quasi_static\"".to_string(),
                ));
            }
            let Some(t0) = def.initial_temperature else {
                return Err(invalid(
                    "[step] kind = \"transient\" demande `initial_temperature` (température de départ, uniforme)"
                        .to_string(),
                ));
            };
            let mut step = Step::transient(def.duration, def.increments)?;
            if let Some(theta) = def.theta {
                step = step.with_theta(theta)?;
            }
            Ok(Run::Transient {
                step,
                initial: vec![t0; nodes],
            })
        }
        StepKind::QuasiStatic => {
            if analysis != AnalysisKind::Static {
                return Err(invalid(
                    "[step] kind = \"quasi_static\" demande analysis = \"static\"".to_string(),
                ));
            }
            if def.theta.is_some() || def.initial_temperature.is_some() {
                return Err(invalid(
                    "[step] `theta` et `initial_temperature` ne servent qu'à kind = \"transient\""
                        .to_string(),
                ));
            }
            let step = Step::quasi_static(def.duration, def.increments)?;
            let amplitude = match &def.amplitude {
                Some(points) => Amplitude::new(points.iter().map(|p| (p[0], p[1])).collect())?,
                None => Amplitude::ramp(def.duration)?,
            };
            Ok(Run::QuasiStatic { step, amplitude })
        }
    }
}
