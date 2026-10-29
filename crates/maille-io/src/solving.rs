//! Résoudre un fichier de modèle : lit le `.toml` et le `.msh`, calcule, écrit les `.vtu` et le `.pvd`.
//! C'est le calcul de `maille solve`, que l'interface lance telle quelle (donc les mêmes fichiers).
//! Les chemins d'un modèle (`mesh`, `from`, `source`) sont relatifs à son dossier.

use std::error::Error;
use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use maille_core::model::Observer;
use maille_core::results::Field;

use crate::meshing;
use crate::model::{self, Loaded, MeshingDef, Physics, Run};
use crate::msh::{self, MshFile};
use crate::vtk::{self, Entry};

pub type Failure = Box<dyn Error>;

fn read(path: &Path) -> Result<String, Failure> {
    fs::read_to_string(path)
        .map_err(|e| format!("impossible de lire {} : {e}", path.display()).into())
}

fn write(path: &Path, text: &str) -> Result<(), Failure> {
    fs::write(path, text)
        .map_err(|e| format!("impossible d'écrire {} : {e}", path.display()).into())
}

/// Le maillage de la recette `[meshing]` ; la source est relative au dossier du modèle.
pub fn mesh_from_recipe(
    dir: &Path,
    def: &MeshingDef,
    log: &mut dyn FnMut(&str),
) -> Result<MshFile, Failure> {
    let path = dir.join(&def.source);
    let source = meshing::Source::of(&path).ok_or_else(|| {
        format!(
            "[meshing] source = \"{}\" : format non pris en charge, écrire un fichier .dxf, .stl ou .step",
            def.source
        )
    })?;
    let bytes =
        fs::read(&path).map_err(|e| format!("impossible de lire {} : {e}", path.display()))?;
    let meshed = meshing::run(source, &bytes, def)?;
    if let Some(unit) = &meshed.unit {
        log(&format!("{} : unité de longueur {unit}", path.display()));
    }
    Ok(meshed.file)
}

/// `observer` est prévenu pendant les calculs longs (thermique transitoire, quasi-statique) ; s'il
/// répond `Flow::Stop`, la résolution rend `MailleError::Cancelled` et rien n'est écrit pour ce
/// calcul. `log` reçoit une phrase à chaque étape (maillage fabriqué, calcul thermique…).
/// Résout le modèle et écrit un `.vtu` par Increment, plus le `.pvd` (toujours, même pour un calcul
/// linéaire : c'est lui que cite `[temperatures] from`). Un calcul couplé écrit d'abord le thermique
/// (`<préfixe>_thermal.pvd`), puis le mécanique dilaté par ses températures. Rend les fichiers
/// écrits, `.pvd` en dernier pour chaque calcul.
pub fn solve_file(
    model_path: &Path,
    output: Option<&Path>,
    observer: &mut Observer,
    log: &mut dyn FnMut(&str),
) -> Result<Vec<PathBuf>, Failure> {
    let dir = model_path.parent().unwrap_or(Path::new("."));
    let file = model::parse(&read(model_path)?)?;
    // le maillage manque : la recette [meshing] le fabrique (sinon `mesh` est lu tel quel)
    let mesh_path = dir.join(&file.mesh);
    if let (false, Some(def)) = (mesh_path.exists(), &file.meshing) {
        let made = mesh_from_recipe(dir, def, log)?;
        write(&mesh_path, &msh::write(&made)?)?;
        log(&format!(
            "{} fabriqué d'après [meshing]",
            mesh_path.display()
        ));
    }
    let mesh = msh::read(&read(&mesh_path)?)?;

    let stem = model_path.file_stem().ok_or("nom de modèle invalide")?;
    let prefix = output.map_or_else(|| dir.join(stem), Path::to_path_buf);
    let named = |suffix: &str| -> PathBuf {
        let mut name: OsString = prefix.as_os_str().to_owned();
        name.push(suffix);
        PathBuf::from(name)
    };

    let mut written = Vec::new();
    // couplage : les températures de l'Increment choisi, soit d'un calcul thermique déjà fait
    // (`from`), soit de la partie thermique de ce modèle, calculée d'abord
    let (file, temperatures) = match model::physics(&file)? {
        Physics::Coupled => {
            let (thermal, mechanical) = model::split(&file)?;
            let loaded = model::build(&thermal, &mesh, None)?;
            log("calcul thermique");
            let increments = run(&loaded, observer)?;
            written.extend(write_results(
                &named,
                "_thermal",
                &increments,
                &loaded,
                &mesh,
            )?);
            let def = mechanical
                .temperatures
                .as_ref()
                .ok_or("calcul couplé sans [temperatures]")?;
            let k = def.index(increments.len(), "le calcul thermique")?;
            (mechanical, Some(increments[k].1.temperatures.clone()))
        }
        _ => {
            let temperatures = match file.temperatures.as_ref() {
                None => None,
                Some(def) => {
                    let from = def.from.as_ref().ok_or("[temperatures] sans `from`")?;
                    let pvd_path = dir.join(from);
                    let entries = vtk::read_pvd(&read(&pvd_path)?)?;
                    let entry = def.pick(&entries)?;
                    let vtu_path = pvd_path
                        .parent()
                        .unwrap_or(Path::new("."))
                        .join(&entry.file);
                    let array = vtk::read_point_data(&read(&vtu_path)?, "temperature")?;
                    Some(array.values)
                }
            };
            (file, temperatures)
        }
    };

    let loaded = model::build(&file, &mesh, temperatures.as_deref())?;
    log("calcul");
    let increments = run(&loaded, observer)?;
    written.extend(write_results(&named, "", &increments, &loaded, &mesh)?);
    // contacts : la force sur chaque surface rigide à chaque incrément, pour la courbe force-enfoncement
    if increments.iter().any(|(_, f)| !f.rigid_forces.is_empty()) {
        let mut csv = String::from("time,contact,fx,fy,fz\n");
        for (time, field) in &increments {
            for (c, f) in field.rigid_forces.iter().enumerate() {
                csv.push_str(&format!("{time:e},{c},{:e},{:e},{:e}\n", f[0], f[1], f[2]));
            }
        }
        let path = named("_contacts.csv");
        write(&path, &csv)?;
        written.push(path);
    }
    Ok(written)
}

/// Résout le calcul demandé : un `(temps, Field)` par Increment.
fn run(loaded: &Loaded, observer: &mut Observer) -> Result<Vec<(f64, Field)>, Failure> {
    Ok(match &loaded.run {
        Run::Linear => vec![(0.0, loaded.model.solve()?)],
        Run::Transient { step, initial } => loaded
            .model
            .solve_step_with(step, initial, observer)?
            .into_iter()
            .map(|i| (i.time, i.field))
            .collect(),
        Run::QuasiStatic { step, amplitude } => loaded
            .model
            .solve_quasi_static_with(step, amplitude, observer)?
            .into_iter()
            .map(|i| (i.time, i.field))
            .collect(),
    })
}

/// Écrit un `.vtu` par Increment (`<préfixe><tag>.vtu` s'il n'y en a qu'un, calcul linéaire) et le
/// `.pvd` ; rend les fichiers écrits, `.pvd` en dernier.
fn write_results(
    named: &dyn Fn(&str) -> PathBuf,
    tag: &str,
    increments: &[(f64, Field)],
    loaded: &Loaded,
    mesh: &msh::MshFile,
) -> Result<Vec<PathBuf>, Failure> {
    let linear = loaded.run == Run::Linear;
    let mut written = Vec::new();
    let mut entries = Vec::new();
    for (k, (time, field)) in increments.iter().enumerate() {
        let path = named(&if linear {
            format!("{tag}.vtu")
        } else {
            format!("{tag}_{k:04}.vtu")
        });
        let arrays = vtk::result_arrays(&loaded.mesh, field)?;
        write(
            &path,
            &vtk::write_vtu(mesh, &loaded.units, loaded.large_strain, &arrays)?,
        )?;
        entries.push(Entry {
            time: *time,
            file: path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
        });
        written.push(path);
    }
    let pvd = named(&format!("{tag}.pvd"));
    write(&pvd, &vtk::write_pvd(&entries)?)?;
    written.push(pvd);
    Ok(written)
}
