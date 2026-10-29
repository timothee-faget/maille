//! Commande `mesh` : dessin DXF 2D, surface STL 3D ou pièce STEP → maillage Gmsh `.msh`. Le calque
//! d'une entité DXF, le nom d'un `solid` STL (ou `surface_N` à défaut), le nom d'une face STEP (ou
//! `face_N` à défaut) nomme le bord (voir `maille_io::dxf`, `maille_io::stl` et `maille_io::step`).

use std::path::{Path, PathBuf};

use maille_io::meshing;
use maille_io::model::{self, MeshingDef};
use maille_io::msh::{self, MshFile, DOMAIN_SET};

use crate::files::{read, read_bytes, write, Failure};

/// `maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]` ou
/// `maille mesh <piece.stl> [--size <h>] [--feature-angle <degrés>] [--quadratic] [-o <maillage.msh>]`
/// ou `maille mesh <piece.step> --size <h> [--quadratic] [-o <maillage.msh>]`.
pub fn mesh(args: &[String]) -> Result<(), Failure> {
    let usage = "usage : maille mesh <modele.toml> [-o <maillage.msh>]\n        maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]\n        maille mesh <piece.stl> [--size <h>] [--feature-angle <degrés>] [--quadratic] [-o <maillage.msh>]\n        maille mesh <piece.step> --size <h> [--quadratic] [-o <maillage.msh>]";
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
    if input.extension().is_some_and(|e| e == "toml") {
        if size.is_some() || angle.is_some() || quadratic {
            return Err("--size, --feature-angle et --quadratic sont dans la section [meshing] du modèle : les y écrire".into());
        }
        return mesh_model(&input, output.as_deref());
    }
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

/// Maille depuis la recette `[meshing]` d'un modèle : le `.msh` de sa clé `mesh` (ou `-o`), avec les
/// Sets utilisateur.
fn mesh_model(model_path: &Path, output: Option<&Path>) -> Result<(), Failure> {
    let dir = model_path.parent().unwrap_or(Path::new("."));
    let file = model::parse(&read(model_path)?)?;
    let def = file.meshing.as_ref().ok_or_else(|| {
        format!(
            "{} n'a pas de section [meshing] : il n'y a pas de recette à suivre (source, size)",
            model_path.display()
        )
    })?;
    let output = output.map_or_else(|| dir.join(&file.mesh), Path::to_path_buf);
    let mesh = from_recipe(dir, def)?;
    write(&output, &msh::write(&mesh)?)?;
    println!(
        "{} : {} Nodes, {} Elements",
        output.display(),
        mesh.nodes.len(),
        mesh.cells.len()
    );
    Ok(())
}

/// Le maillage de la recette `[meshing]` ; la source est relative au dossier du modèle.
pub(crate) fn from_recipe(dir: &Path, def: &MeshingDef) -> Result<MshFile, Failure> {
    let path = dir.join(&def.source);
    let source = meshing::Source::of(&path).ok_or_else(|| {
        format!(
            "[meshing] source = \"{}\" : format non pris en charge, écrire un fichier .dxf, .stl ou .step",
            def.source
        )
    })?;
    let meshed = meshing::run(source, &read_bytes(&path)?, def)?;
    if let Some(unit) = &meshed.unit {
        println!("{} : unité de longueur {unit}", path.display());
    }
    Ok(meshed.file)
}
