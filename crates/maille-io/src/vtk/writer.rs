//! Écrivain `.vtu` et `.pvd`. Les flottants s'écrivent avec leur représentation décimale exacte
//! (`{}` de Rust) : la relecture redonne les mêmes bits.

use std::fmt::Write;

use super::{vtk_error, Entry, PointArray};
use crate::msh::{CellKind, MshFile};
use crate::IoError;

fn cell_type(kind: CellKind) -> u8 {
    match kind {
        CellKind::Tri3 => 5,
        CellKind::Tri6 => 22,
        CellKind::Tet4 => 10,
        // sommets, puis milieux des arêtes (0,1) (1,2) (2,0) (0,3) (1,3) (2,3), comme Maille
        CellKind::Tet10 => 24,
    }
}

/// Un nom qui s'écrit tel quel dans un attribut XML.
fn plain(name: &str) -> bool {
    !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_alphanumeric() || matches!(c, '_' | '-' | '.' | ' ' | '/'))
}

fn push_array<T: std::fmt::Display>(
    out: &mut String,
    attributes: &str,
    components: usize,
    values: impl Iterator<Item = T>,
) {
    let _ = write!(out, "<DataArray {attributes} format=\"ascii\">");
    for (i, v) in values.enumerate() {
        let sep = if i % components == 0 { "\n" } else { " " };
        let _ = write!(out, "{sep}{v}");
    }
    let _ = writeln!(out, "\n</DataArray>");
}

/// Écrit le maillage et ses grandeurs aux Nodes. `units` est recopié en commentaire XML en tête ;
/// `large_strain` ajoute le commentaire `kinematics: large` (grandes déformations : les contraintes
/// sont de Cauchy et le viewer montre la déformée à l'échelle 1).
pub fn write_vtu(
    mesh: &MshFile,
    units: &str,
    large_strain: bool,
    arrays: &[PointArray],
) -> Result<String, IoError> {
    if units.contains("--") || units.contains('\n') || units.contains('\r') {
        return Err(vtk_error(
            "`units` ne peut pas contenir « -- » ni de saut de ligne".to_string(),
        ));
    }
    for a in arrays {
        if !plain(&a.name) {
            return Err(vtk_error(format!(
                "nom de grandeur « {} » : lettres, chiffres, espace, _ - . / seulement",
                a.name
            )));
        }
        if a.components == 0 || a.values.len() != a.components * mesh.nodes.len() {
            return Err(vtk_error(format!(
                "la grandeur « {} » a {} valeurs, il en faut {} par Node ({} Nodes)",
                a.name,
                a.values.len(),
                a.components,
                mesh.nodes.len()
            )));
        }
        if let Some(i) = a.values.iter().position(|v| !v.is_finite()) {
            return Err(vtk_error(format!(
                "la grandeur « {} » n'est pas finie au Node {} : le calcul a divergé ?",
                a.name,
                i / a.components
            )));
        }
    }

    let width = mesh.kind.node_count();
    let mut out = String::new();
    let _ = writeln!(out, "<?xml version=\"1.0\"?>");
    let _ = writeln!(out, "<!-- units: {units} -->");
    if large_strain {
        let _ = writeln!(out, "<!-- kinematics: large -->");
    }
    let _ = writeln!(
        out,
        "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">"
    );
    let _ = writeln!(out, "<UnstructuredGrid>");
    let _ = writeln!(
        out,
        "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
        mesh.nodes.len(),
        mesh.cells.len()
    );
    let _ = writeln!(out, "<PointData>");
    for a in arrays {
        push_array(
            &mut out,
            &format!(
                "type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\"",
                a.name, a.components
            ),
            a.components,
            a.values.iter(),
        );
    }
    let _ = writeln!(out, "</PointData>");
    let _ = writeln!(out, "<Points>");
    push_array(
        &mut out,
        "type=\"Float64\" NumberOfComponents=\"3\"",
        3,
        mesh.nodes.iter().flat_map(|p| [p.x, p.y, p.z]),
    );
    let _ = writeln!(out, "</Points>");
    let _ = writeln!(out, "<Cells>");
    push_array(
        &mut out,
        "type=\"Int64\" Name=\"connectivity\"",
        width,
        mesh.cells.iter().flatten(),
    );
    push_array(
        &mut out,
        "type=\"Int64\" Name=\"offsets\"",
        1,
        (1..=mesh.cells.len()).map(|k| k * width),
    );
    push_array(
        &mut out,
        "type=\"UInt8\" Name=\"types\"",
        1,
        mesh.cells.iter().map(|_| cell_type(mesh.kind)),
    );
    let _ = writeln!(out, "</Cells>");
    let _ = writeln!(out, "</Piece>\n</UnstructuredGrid>\n</VTKFile>");
    Ok(out)
}

/// Écrit le `.pvd` : les fichiers `.vtu` (noms relatifs au `.pvd`) et le temps de chacun.
pub fn write_pvd(entries: &[Entry]) -> Result<String, IoError> {
    if entries.is_empty() {
        return Err(vtk_error(
            "un .pvd sans fichier n'a pas de sens".to_string(),
        ));
    }
    let mut out = String::new();
    let _ = writeln!(out, "<?xml version=\"1.0\"?>");
    let _ = writeln!(
        out,
        "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>"
    );
    for e in entries {
        if e.file.contains(['"', '<', '>', '&']) || !e.time.is_finite() {
            return Err(vtk_error(format!(
                "entrée « {} » à t = {} : nom sans guillemet, < > ni &, temps fini",
                e.file, e.time
            )));
        }
        let _ = writeln!(
            out,
            "<DataSet timestep=\"{}\" part=\"0\" file=\"{}\"/>",
            e.time, e.file
        );
    }
    let _ = writeln!(out, "</Collection>\n</VTKFile>");
    Ok(out)
}
