//! Commandes `solve` et `view` : lisent et écrivent les fichiers (`.msh`, `.toml`, `.vtu`, `.pvd`),
//! `maille-io` ne fait que du texte. Les chemins d'un modèle (`mesh`, `from`) sont relatifs à son
//! dossier.

use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use maille_core::model::Flow;
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

/// Résout le modèle (voir `maille_io::solving`), en écrivant à l'écran ce que le calcul raconte.
fn solve_file(model_path: &Path, output: Option<&Path>) -> Result<Vec<PathBuf>, Failure> {
    maille_io::solving::solve_file(model_path, output, &mut |_| Flow::Continue, &mut |line| {
        println!("{line}")
    })
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
