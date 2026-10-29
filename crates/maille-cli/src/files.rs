//! Commandes `solve` et `view` : lisent et écrivent les fichiers (`.msh`, `.toml`, `.vtu`, `.pvd`),
//! `maille-io` ne fait que du texte. Les chemins d'un modèle (`mesh`, `from`) sont relatifs à son
//! dossier.

use std::error::Error;
use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use maille_core::results::{Field, Stress};
use maille_io::model::{self, Run};
use maille_io::msh;
use maille_io::vtk::{self, Entry, VtuResults};
use maille_viewer::{NodalResults, Shot};

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
/// linéaire : c'est lui que cite `[temperatures]`). Rend les fichiers écrits, `.pvd` en dernier.
fn solve_file(model_path: &Path, output: Option<&Path>) -> Result<Vec<PathBuf>, Failure> {
    let dir = model_path.parent().unwrap_or(Path::new("."));
    let file = model::parse(&read(model_path)?)?;
    let mesh = msh::read(&read(&dir.join(&file.mesh))?)?;

    // couplage : les températures de l'Increment choisi d'un calcul thermique déjà fait
    let temperatures = match &file.temperatures {
        None => None,
        Some(def) => {
            let pvd_path = dir.join(&def.from);
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

    let loaded = model::build(&file, &mesh, temperatures.as_deref())?;
    let linear = loaded.run == Run::Linear;
    let increments: Vec<(f64, Field)> = match &loaded.run {
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
    };

    let stem = model_path.file_stem().ok_or("nom de modèle invalide")?;
    let prefix = output.map_or_else(|| dir.join(stem), Path::to_path_buf);
    let named = |suffix: &str| -> PathBuf {
        let mut name: OsString = prefix.as_os_str().to_owned();
        name.push(suffix);
        PathBuf::from(name)
    };
    let mut written = Vec::new();
    let mut entries = Vec::new();
    for (k, (time, field)) in increments.iter().enumerate() {
        let path = named(&if linear {
            ".vtu".to_string()
        } else {
            format!("_{k:04}.vtu")
        });
        let arrays = vtk::result_arrays(&loaded.mesh, field)?;
        write(
            &path,
            &vtk::write_vtu(&mesh, &loaded.units, loaded.large_strain, &arrays)?,
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
    let pvd = named(".pvd");
    write(&pvd, &vtk::write_pvd(&entries)?)?;
    written.push(pvd);
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
    let results = to_nodal(&open_results(Path::new(&path), index)?)?;
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

/// Le `.vtu` demandé : directement, ou l'Increment `index` d'un `.pvd` (le dernier par défaut).
fn open_results(path: &Path, index: Option<usize>) -> Result<VtuResults, Failure> {
    let vtu = match path.extension().and_then(|e| e.to_str()) {
        Some("pvd") => {
            let entries = vtk::read_pvd(&read(path)?)?;
            let k = index.unwrap_or(entries.len() - 1);
            let entry = entries.get(k).ok_or_else(|| {
                format!(
                    "incrément {k} : {} en contient {}, numérotés depuis 0",
                    path.display(),
                    entries.len()
                )
            })?;
            path.parent().unwrap_or(Path::new(".")).join(&entry.file)
        }
        Some("vtu") if index.is_none() => path.to_path_buf(),
        Some("vtu") => return Err("un .vtu n'a qu'un incrément : retirer l'indice".into()),
        _ => return Err(format!("{} : attendu un fichier .pvd ou .vtu", path.display()).into()),
    };
    Ok(vtk::read_results(&read(&vtu)?)?)
}

/// Résultats relus → ce que le viewer affiche (triangles et tétraèdres).
fn to_nodal(r: &VtuResults) -> Result<NodalResults, Failure> {
    let array = |name: &str, components: usize| {
        r.arrays
            .iter()
            .find(|a| a.name == name && a.components == components)
            .map(|a| a.values.as_slice())
    };
    let temperatures = array("temperature", 1).unwrap_or_default().to_vec();
    let stresses: Vec<Stress> = array("stress", 6)
        .map(|s| {
            s.chunks(6)
                .zip(&r.nodes)
                .map(|(c, p)| Stress {
                    x: p.x,
                    y: p.y,
                    z: p.z,
                    sigma_xx: c[0],
                    sigma_yy: c[1],
                    sigma_zz: c[2],
                    tau_xy: c[3],
                    tau_yz: c[4],
                    tau_xz: c[5],
                })
                .collect()
        })
        .unwrap_or_default();
    let displacements: Vec<[f64; 3]> = array("displacement", 3)
        .map(|d| d.chunks(3).map(|u| [u[0], u[1], u[2]]).collect())
        .unwrap_or_default();
    if temperatures.is_empty() && (stresses.is_empty() || displacements.is_empty()) {
        return Err(
            "ce fichier n'a ni « temperature », ni « displacement » et « stress » : rien à afficher"
                .into(),
        );
    }
    Ok(NodalResults {
        nodes: r.nodes.iter().map(|p| [p.x, p.y, p.z]).collect(),
        cells: r.cells.clone(),
        displacements,
        stresses,
        temperatures,
        plastic_strains: array("plastic_strain", 1).unwrap_or_default().to_vec(),
        large_strain: r.large_strain,
    })
}

#[cfg(test)]
mod tests {
    use maille_core::mesh::Node;
    use maille_io::msh::CellKind;
    use maille_io::vtk::PointArray;

    use super::*;

    fn results(kind: CellKind, arrays: Vec<PointArray>) -> VtuResults {
        let nodes = vec![
            Node::new_2d(0.0, 0.0),
            Node::new_2d(1.0, 0.0),
            Node::new_2d(0.0, 1.0),
        ];
        VtuResults {
            nodes,
            kind,
            cells: vec![vec![0, 1, 2]],
            arrays,
            units: None,
            large_strain: false,
        }
    }

    fn array(name: &str, components: usize) -> PointArray {
        PointArray {
            name: name.to_string(),
            components,
            values: (0..3 * components).map(|i| i as f64).collect(),
        }
    }

    #[test]
    fn mechanical_results_become_nodal_results() {
        let r = results(
            CellKind::Tri3,
            vec![array("displacement", 3), array("stress", 6)],
        );
        let n = to_nodal(&r).unwrap();
        assert_eq!(
            n.displacements,
            [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]
        );
        assert_eq!(n.stresses[1].sigma_xx, 6.0);
        assert_eq!((n.stresses[1].x, n.stresses[2].y), (1.0, 1.0));
        assert!(n.temperatures.is_empty() && n.plastic_strains.is_empty());
    }

    #[test]
    fn thermal_results_need_only_the_temperature() {
        let n = to_nodal(&results(CellKind::Tri3, vec![array("temperature", 1)])).unwrap();
        assert_eq!(n.temperatures, [0.0, 1.0, 2.0]);
    }

    #[test]
    fn tetrahedra_are_shown_and_empty_files_are_refused() {
        let solid = to_nodal(&results(CellKind::Tet4, vec![array("temperature", 1)])).unwrap();
        assert_eq!(solid.temperatures.len(), 3);
        let empty = to_nodal(&results(CellKind::Tri3, vec![array("stress", 6)]));
        assert!(empty.unwrap_err().to_string().contains("rien à afficher"));
    }
}
