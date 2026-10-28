//! Commande `mesh` : dessin DXF 2D, surface STL 3D ou pièce STEP → maillage Gmsh `.msh`. Le calque
//! d'une entité DXF, le nom d'un `solid` STL (ou `surface_N` à défaut), le nom d'une face STEP (ou
//! `face_N` à défaut) nomme le bord (voir `maille_io::dxf`, `maille_io::stl` et `maille_io::step`).

use std::path::{Path, PathBuf};

use maille_core::mesher::{
    tetrahedralize, tetrahedralize_brep, triangulate, MeshOptions, TetOptions,
};
use maille_io::msh::{self, MshFile, DOMAIN_SET};
use maille_io::{dxf, step, stl};

use crate::files::{read, read_bytes, write, Failure};

/// Angle vif par défaut du STL, en degrés.
const FEATURE_ANGLE: f64 = 30.0;

/// `maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]` ou
/// `maille mesh <piece.stl> [--size <h>] [--feature-angle <degrés>] [--quadratic] [-o <maillage.msh>]`
/// ou `maille mesh <piece.step> --size <h> [--quadratic] [-o <maillage.msh>]`.
pub fn mesh(args: &[String]) -> Result<(), Failure> {
    let usage = "usage : maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]\n        maille mesh <piece.stl> [--size <h>] [--feature-angle <degrés>] [--quadratic] [-o <maillage.msh>]\n        maille mesh <piece.step> --size <h> [--quadratic] [-o <maillage.msh>]";
    let (mut input, mut size, mut angle, mut output, mut quadratic) =
        (None, None, None, None, false);
    let mut it = args.iter();
    let number = |option: &str, unit: &str, text: Option<&String>| -> Result<f64, Failure> {
        let text = text.ok_or(usage)?;
        text.parse::<f64>()
            .map_err(|_| format!("{option} attend un nombre ({unit}), trouvé « {text} »").into())
    };
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--size" => {
                size = Some(number(
                    "--size",
                    "taille h des Elements, dans les unités du fichier",
                    it.next(),
                )?)
            }
            "--feature-angle" => {
                angle = Some(number(
                    "--feature-angle",
                    "angle vif, en degrés",
                    it.next(),
                )?)
            }
            "--quadratic" => quadratic = true,
            "-o" | "--output" => output = Some(PathBuf::from(it.next().ok_or(usage)?)),
            _ if input.is_none() => input = Some(PathBuf::from(arg)),
            _ => return Err(usage.into()),
        }
    }
    let input: PathBuf = input.ok_or(usage)?;
    let output = output.unwrap_or_else(|| input.with_extension("msh"));
    let extension = input
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    let (file, what) = match extension.as_deref() {
        Some("dxf") => {
            if angle.is_some() {
                return Err("--feature-angle ne sert qu'aux fichiers STL".into());
            }
            let size = size.ok_or(
                "--size manque : donner la taille h des Elements, dans les unités du dessin",
            )?;
            (from_dxf(&input, size, quadratic)?, "triangles")
        }
        Some("stl") => {
            let options = TetOptions {
                size,
                ..TetOptions::default()
            };
            let angle = angle.unwrap_or(FEATURE_ANGLE);
            (from_stl(&input, angle, &options, quadratic)?, "tétraèdres")
        }
        Some("step" | "stp") => {
            if angle.is_some() {
                return Err("--feature-angle ne sert qu'aux fichiers STL".into());
            }
            let options = TetOptions {
                size: Some(size.ok_or(
                    "--size manque : donner la taille h des Elements, dans l'unité de longueur du fichier STEP",
                )?),
                ..TetOptions::default()
            };
            (from_step(&input, &options, quadratic)?, "tétraèdres")
        }
        other => {
            return Err(format!(
                "extension {} non prise en charge : Maille maille un dessin .dxf (2D), une surface .stl ou une pièce .step (3D)",
                other.map_or("absente".to_string(), |e| format!("« .{e} »"))
            )
            .into())
        }
    };
    write(&output, &msh::write(&file)?)?;
    println!(
        "{} : {} Nodes, {} {what}, bords {}",
        output.display(),
        file.nodes.len(),
        file.cells.len(),
        file.sets
            .keys()
            .filter(|name| *name != DOMAIN_SET)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}

fn from_dxf(path: &Path, size: f64, quadratic: bool) -> Result<MshFile, Failure> {
    let geometry = dxf::read_geometry(&read(path)?, size)?;
    let triangulation = triangulate(&geometry, &MeshOptions::default())?;
    Ok(if quadratic {
        MshFile::from_triangulation(&triangulation.to_quadratic())?
    } else {
        MshFile::from_triangulation(&triangulation)?
    })
}

/// Sans `--size`, les tétraèdres ont la finesse des triangles du STL ; avec, la surface est
/// remaillée à la taille h puis l'intérieur raffiné (voir `TetOptions`).
/// Les faces sont maillées à la taille h sur leurs surfaces exactes, puis l'intérieur raffiné
/// (voir `tetrahedralize_brep`). L'unité de longueur du fichier est affichée : Maille ne convertit
/// rien, la taille h et le modèle sont dans cette unité.
fn from_step(path: &Path, options: &TetOptions, quadratic: bool) -> Result<MshFile, Failure> {
    let (brep, unit) = step::read_brep(&read(path)?)?;
    println!("{} : unité de longueur {unit}", path.display());
    let tetrahedra = tetrahedralize_brep(&brep, options)?;
    Ok(if quadratic {
        MshFile::from_tetrahedralization(&tetrahedra.to_quadratic())?
    } else {
        MshFile::from_tetrahedralization(&tetrahedra)?
    })
}

fn from_stl(
    path: &Path,
    feature_angle: f64,
    options: &TetOptions,
    quadratic: bool,
) -> Result<MshFile, Failure> {
    let surface = stl::read_surface(&read_bytes(path)?, feature_angle)?;
    let tetrahedra = tetrahedralize(&surface, options)?;
    Ok(if quadratic {
        MshFile::from_tetrahedralization(&tetrahedra.to_quadratic())?
    } else {
        MshFile::from_tetrahedralization(&tetrahedra)?
    })
}
