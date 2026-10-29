//! Commande `mesh` : dessin DXF 2D, surface STL 3D ou pièce STEP → maillage Gmsh `.msh`. Le calque
//! d'une entité DXF, le nom d'un `solid` STL (ou `surface_N` à défaut), le nom d'une face STEP (ou
//! `face_N` à défaut) nomme le bord (voir `maille_io::dxf`, `maille_io::stl` et `maille_io::step`).

use std::path::PathBuf;

use maille_io::meshing;
use maille_io::msh::{self, DOMAIN_SET};

use crate::files::{read_bytes, write, Failure};

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
    let source = meshing::Source::of(&input).ok_or_else(|| {
        format!(
            "extension {} non prise en charge : Maille maille un dessin .dxf (2D), une surface .stl ou une pièce .step (3D)",
            input
                .extension()
                .and_then(|e| e.to_str())
                .map_or("absente".to_string(), |e| format!("« .{} »", e.to_ascii_lowercase()))
        )
    })?;
    if angle.is_some() && source != meshing::Source::Stl {
        return Err("--feature-angle ne sert qu'aux fichiers STL".into());
    }
    if size.is_none() && source.needs_size() {
        return Err(match source {
            meshing::Source::Dxf => {
                "--size manque : donner la taille h des Elements, dans les unités du dessin"
            }
            _ => "--size manque : donner la taille h des Elements, dans l'unité de longueur du fichier STEP",
        }
        .into());
    }
    let recipe = meshing::Recipe {
        size,
        feature_angle: angle.unwrap_or(meshing::FEATURE_ANGLE),
        quadratic,
    };
    let meshed = meshing::mesh(source, &read_bytes(&input)?, &recipe)?;
    // l'unité de longueur du STEP est affichée : Maille ne convertit rien, h et le modèle sont dans cette unité
    if let Some(unit) = &meshed.unit {
        println!("{} : unité de longueur {unit}", input.display());
    }
    let file = meshed.file;
    let what = if file.kind.dim() == 2 {
        "triangles"
    } else {
        "tétraèdres"
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
