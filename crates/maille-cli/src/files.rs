//! Commandes `solve` et `view` : lisent et écrivent les fichiers (`.msh`, `.toml`, `.vtu`, `.pvd`),
//! `maille-io` ne fait que du texte. Les chemins d'un modèle (`mesh`, `from`) sont relatifs à son
//! dossier.

use std::error::Error;
use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use maille_core::results::Field;
use maille_io::model::{self, Loaded, Physics, Run};
use maille_io::msh;
use maille_io::vtk::{self, Entry};
use maille_viewer::Shot;

pub(crate) type Failure = Box<dyn Error>;

pub(crate) fn read(path: &Path) -> Result<String, Failure> {
    fs::read_to_string(path)
        .map_err(|e| format!("impossible de lire {} : {e}", path.display()).into())
}

pub(crate) fn read_bytes(path: &Path) -> Result<Vec<u8>, Failure> {
    fs::read(path).map_err(|e| format!("impossible de lire {} : {e}", path.display()).into())
}

pub(crate) fn write(path: &Path, text: &str) -> Result<(), Failure> {
    fs::write(path, text)
        .map_err(|e| format!("impossible d'écrire {} : {e}", path.display()).into())
}

/// `maille solve <modele.toml> [-o <préfixe>]`.
pub fn solve(args: &[String]) -> Result<(), Failure> {
    let usage = "usage : maille solve <modele.toml> [-o <préfixe>]";
    let (mut model, mut output) = (None, None);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-o" | "--output" => {
                output = Some(PathBuf::from(it.next().ok_or(usage)?));
            }
            _ if model.is_none() => model = Some(PathBuf::from(arg)),
            _ => return Err(usage.into()),
        }
    }
    let model = model.ok_or(usage)?;
    let written = solve_file(&model, output.as_deref())?;
    println!("{} fichier(s) écrit(s) :", written.len());
    for path in &written {
        println!("  {}", path.display());
    }
    Ok(())
}

/// Résout le modèle et écrit un `.vtu` par Increment, plus le `.pvd` (toujours, même pour un calcul
/// linéaire : c'est lui que cite `[temperatures] from`). Un calcul couplé écrit d'abord le thermique
/// (`<préfixe>_thermal.pvd`), puis le mécanique dilaté par ses températures. Rend les fichiers
/// écrits, `.pvd` en dernier pour chaque calcul.
fn solve_file(model_path: &Path, output: Option<&Path>) -> Result<Vec<PathBuf>, Failure> {
    let dir = model_path.parent().unwrap_or(Path::new("."));
    let file = model::parse(&read(model_path)?)?;
    // le maillage manque : la recette [meshing] le fabrique (sinon `mesh` est lu tel quel)
    let mesh_path = dir.join(&file.mesh);
    if let (false, Some(def)) = (mesh_path.exists(), &file.meshing) {
        let made = crate::mesh::from_recipe(dir, def)?;
        write(&mesh_path, &msh::write(&made)?)?;
        println!("{} fabriqué d'après [meshing]", mesh_path.display());
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
            let increments = run(&loaded)?;
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
    let increments = run(&loaded)?;
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
fn run(loaded: &Loaded) -> Result<Vec<(f64, Field)>, Failure> {
    Ok(match &loaded.run {
        Run::Linear => vec![(0.0, loaded.model.solve()?)],
        Run::Transient { step, initial } => loaded
            .model
            .solve_step(step, initial)?
            .into_iter()
            .map(|i| (i.time, i.field))
            .collect(),
        Run::QuasiStatic { step, amplitude } => loaded
            .model
            .solve_quasi_static(step, amplitude)?
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

/// `maille view <resultats.pvd|.vtu> [indice] [--png image.png ...]`.
pub fn view(args: &[String]) -> Result<(), Failure> {
    let usage = "usage : maille view <resultats.pvd|resultats.vtu> [indice]\n        maille view <resultats.pvd|resultats.vtu> [indice] --png <image.png> [--size LxH] [--angles lacet,tangage] [--field k] [--deformed] [--cut x|y|z=position]";
    let (mut path, mut index, mut png) = (None, None, None);
    let mut shot = Shot::default();
    let mut it = args.iter();
    let number = |text: &str, what: &str| {
        text.parse::<f64>()
            .map_err(|_| format!("{what} : un nombre est attendu, trouvé « {text} »"))
    };
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--png" => png = Some(PathBuf::from(it.next().ok_or(usage)?)),
            "--size" => {
                let text = it.next().ok_or(usage)?;
                let (w, h) = text.split_once('x').ok_or_else(|| {
                    format!("--size attend largeurxhauteur, par exemple 900x650, trouvé « {text} »")
                })?;
                let pixels = |t: &str| {
                    t.parse::<u32>()
                        .ok()
                        .filter(|&n| (1..=8192).contains(&n))
                        .ok_or_else(|| {
                            format!("--size : « {t} » n'est pas une taille de 1 à 8192 pixels")
                        })
                };
                (shot.width, shot.height) = (pixels(w)?, pixels(h)?);
            }
            "--angles" => {
                let text = it.next().ok_or(usage)?;
                let (yaw, pitch) = text.split_once(',').ok_or_else(|| {
                    format!("--angles attend lacet,tangage en degrés, par exemple 30,20, trouvé « {text} »")
                })?;
                (shot.yaw, shot.pitch) = (number(yaw, "--angles")?, number(pitch, "--angles")?);
            }
            "--field" => {
                let text = it.next().ok_or(usage)?;
                shot.quantity = text
                    .parse::<usize>()
                    .ok()
                    .and_then(|k| k.checked_sub(1))
                    .ok_or_else(|| {
                        format!("--field attend le numéro de la touche (1, 2…), trouvé « {text} »")
                    })?;
            }
            "--deformed" => shot.deformed = true,
            "--cut" => {
                let text = it.next().ok_or(usage)?;
                let (axis, position) = text.split_once('=').ok_or_else(|| {
                    format!("--cut attend axe=position, par exemple x=0.5 (position de 0 à 1 dans la boîte), trouvé « {text} »")
                })?;
                let axis = ["x", "y", "z"]
                    .iter()
                    .position(|a| *a == axis)
                    .ok_or_else(|| format!("--cut : axe « {axis} » inconnu, choisir x, y ou z"))?;
                shot.cut = Some((axis, number(position, "--cut")?));
            }
            _ if path.is_none() => path = Some(arg.clone()),
            _ if index.is_none() => {
                index = Some(arg.parse::<usize>().map_err(|_| {
                    format!("indice « {arg} » : écrire un entier, depuis 0 ({usage})")
                })?)
            }
            _ => return Err(usage.into()),
        }
    }
    let path = path.ok_or(usage)?;
    let results =
        maille_gui::results::to_nodal(&maille_gui::results::open(Path::new(&path), index)?)?;
    match png {
        None => maille_viewer::show_nodal(&results),
        Some(file) => {
            fs::write(&file, maille_viewer::screenshot_nodal(&results, &shot)?)
                .map_err(|e| format!("impossible d'écrire {} : {e}", file.display()))?;
            println!("{}", file.display());
            Ok(())
        }
    }
}
