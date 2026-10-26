//! Commande `mesh` : dessin DXF 2D → maillage Gmsh `.msh`. Le calque de chaque entité nomme le bord
//! (voir `maille_io::dxf`).

use std::path::PathBuf;

use maille_core::mesher::{triangulate, MeshOptions};
use maille_io::dxf;
use maille_io::msh::{self, MshFile, DOMAIN_SET};

use crate::files::{read, write, Failure};

/// `maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]`.
pub fn mesh(args: &[String]) -> Result<(), Failure> {
    let usage = "usage : maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]";
    let (mut drawing, mut size, mut output, mut quadratic) = (None, None, None, false);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--size" => {
                let text = it.next().ok_or(usage)?;
                size = Some(text.parse::<f64>().map_err(|_| {
                    format!("--size attend un nombre (taille h des Elements, dans les unités du dessin), trouvé « {text} »")
                })?);
            }
            "--quadratic" => quadratic = true,
            "-o" | "--output" => output = Some(PathBuf::from(it.next().ok_or(usage)?)),
            _ if drawing.is_none() => drawing = Some(PathBuf::from(arg)),
            _ => return Err(usage.into()),
        }
    }
    let drawing = drawing.ok_or(usage)?;
    let size =
        size.ok_or("--size manque : donner la taille h des Elements, dans les unités du dessin")?;
    let output = output.unwrap_or_else(|| drawing.with_extension("msh"));

    let geometry = dxf::read_geometry(&read(&drawing)?, size)?;
    let triangulation = triangulate(&geometry, &MeshOptions::default())?;
    let file = if quadratic {
        MshFile::from_triangulation(&triangulation.to_quadratic())?
    } else {
        MshFile::from_triangulation(&triangulation)?
    };
    write(&output, &msh::write(&file)?)?;
    println!(
        "{} : {} Nodes, {} triangles, bords {}",
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
