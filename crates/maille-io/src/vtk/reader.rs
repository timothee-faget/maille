//! Lecteur minimal de `.vtu` et `.pvd` ASCII : ce que `writer` écrit, et ce que ParaView écrit en
//! ASCII. Un petit scanner de balises suffit, il ne comprend pas le XML en général.

use maille_core::mesh::Node;

use super::{vtk_error, Entry, PointArray, VtuResults};
use crate::msh::CellKind;
use crate::IoError;

/// Valeur de l'attribut `name` dans le texte d'une balise ouvrante.
fn attribute<'a>(tag: &'a str, name: &str) -> Option<&'a str> {
    let key = format!(" {name}=\"");
    let start = tag.find(&key)? + key.len();
    let len = tag[start..].find('"')?;
    Some(&tag[start..start + len])
}

/// Textes (de `<` à `>`, balise ouvrante) des balises `<name` du texte, avec ce qui les suit.
fn tags<'a>(text: &'a str, name: &str) -> Vec<(&'a str, &'a str)> {
    let open = format!("<{name}");
    let mut found = Vec::new();
    let mut rest = text;
    while let Some(i) = rest.find(&open) {
        rest = &rest[i..];
        let Some(end) = rest.find('>') else { break };
        found.push((&rest[..end], &rest[end + 1..]));
        rest = &rest[end + 1..];
    }
    found
}

/// Liste des `.vtu` d'un `.pvd` avec leur temps, dans l'ordre du fichier.
pub fn read_pvd(text: &str) -> Result<Vec<Entry>, IoError> {
    let mut entries = Vec::new();
    for (tag, _) in tags(text, "DataSet") {
        let file = attribute(tag, "file")
            .ok_or_else(|| vtk_error("un <DataSet> du .pvd n'a pas d'attribut file".to_string()))?;
        let time = attribute(tag, "timestep").unwrap_or("0");
        let time = time
            .parse()
            .map_err(|_| vtk_error(format!("timestep=\"{time}\" n'est pas un nombre valable")))?;
        entries.push(Entry {
            time,
            file: file.to_string(),
        });
    }
    if entries.is_empty() {
        return Err(vtk_error(
            "le .pvd ne contient aucun <DataSet> : ce n'est pas un fichier de résultats de Maille ?"
                .to_string(),
        ));
    }
    Ok(entries)
}

/// Partie du texte entre `<name` et `</name>`.
fn section<'a>(text: &'a str, name: &str) -> Result<&'a str, IoError> {
    let start = text.find(&format!("<{name}")).ok_or_else(|| {
        vtk_error(format!(
            "le .vtu n'a pas de <{name}> : ce n'est pas un fichier de résultats de Maille ?"
        ))
    })?;
    let rest = &text[start..];
    Ok(&rest[..rest.find(&format!("</{name}>")).unwrap_or(rest.len())])
}

/// Valeurs d'un `<DataArray>` ASCII (`tag` : sa balise ouvrante, `body` : le texte qui la suit).
fn parse_array(tag: &str, body: &str, name: &str) -> Result<PointArray, IoError> {
    if attribute(tag, "format") != Some("ascii") {
        return Err(vtk_error(format!(
            "la grandeur « {name} » n'est pas en ASCII : Maille lit les .vtu ASCII (dans ParaView, enregistrer en ASCII)"
        )));
    }
    let components = match attribute(tag, "NumberOfComponents") {
        None => 1,
        Some(c) => c.parse().map_err(|_| {
            vtk_error(format!(
                "NumberOfComponents=\"{c}\" n'est pas un nombre valable"
            ))
        })?,
    };
    let body = if tag.ends_with('/') {
        ""
    } else {
        &body[..body.find("</DataArray>").unwrap_or(body.len())]
    };
    let values = body
        .split_whitespace()
        .map(|v| {
            v.parse::<f64>().map_err(|_| {
                vtk_error(format!(
                    "« {v} » n'est pas un nombre valable dans « {name} »"
                ))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if components == 0 || !values.len().is_multiple_of(components) {
        return Err(vtk_error(format!(
            "« {name} » : {} valeurs pour {components} composantes",
            values.len()
        )));
    }
    Ok(PointArray {
        name: name.to_string(),
        components,
        values,
    })
}

/// Grandeur aux Nodes `name` d'un `.vtu` ASCII. Erreur, avec les noms présents, si elle manque.
pub fn read_point_data(text: &str, name: &str) -> Result<PointArray, IoError> {
    let section = section(text, "PointData")?;
    let mut names = Vec::new();
    for (tag, body) in tags(section, "DataArray") {
        let Some(found) = attribute(tag, "Name") else {
            continue;
        };
        names.push(found);
        if found == name {
            return parse_array(tag, body, name);
        }
    }
    Err(vtk_error(format!(
        "pas de grandeur « {name} » aux Nodes ; présentes : {}",
        if names.is_empty() {
            "aucune".to_string()
        } else {
            names
                .iter()
                .map(|n| format!("« {n} »"))
                .collect::<Vec<_>>()
                .join(", ")
        }
    )))
}

fn cell_kind(code: usize) -> Option<CellKind> {
    match code {
        5 => Some(CellKind::Tri3),
        22 => Some(CellKind::Tri6),
        10 => Some(CellKind::Tet4),
        24 => Some(CellKind::Tet10),
        _ => None,
    }
}

/// Un `.vtu` en entier : maillage et toutes les grandeurs aux Nodes.
pub fn read_results(text: &str) -> Result<VtuResults, IoError> {
    let points = section(text, "Points")?;
    let (tag, body) = *tags(points, "DataArray")
        .first()
        .ok_or_else(|| vtk_error("<Points> sans coordonnées".to_string()))?;
    let coordinates = parse_array(tag, body, "Points")?;
    if coordinates.components != 3 {
        return Err(vtk_error(
            "les points doivent avoir 3 coordonnées (NumberOfComponents=\"3\")".to_string(),
        ));
    }
    let nodes: Vec<Node> = coordinates
        .values
        .chunks(3)
        .map(|p| Node::new(p[0], p[1], p[2]))
        .collect();

    let cells_section = section(text, "Cells")?;
    let cell_array = |name: &str| -> Result<Vec<usize>, IoError> {
        let (tag, body) = tags(cells_section, "DataArray")
            .into_iter()
            .find(|(tag, _)| attribute(tag, "Name") == Some(name))
            .ok_or_else(|| vtk_error(format!("<Cells> sans « {name} »")))?;
        // des entiers : exacts en f64 jusqu'à 2^53
        Ok(parse_array(tag, body, name)?
            .values
            .into_iter()
            .map(|v| v as usize)
            .collect())
    };
    let (connectivity, offsets, types) = (
        cell_array("connectivity")?,
        cell_array("offsets")?,
        cell_array("types")?,
    );
    let Some(&first) = types.first() else {
        return Err(vtk_error("le .vtu ne contient aucune cellule".to_string()));
    };
    let kind = cell_kind(first).ok_or_else(|| {
        vtk_error(format!(
            "type de cellule VTK {first} non pris en charge : Maille lit les triangles (5, 22) et les tétraèdres (10, 24)"
        ))
    })?;
    let width = kind.node_count();
    let regular = types.iter().all(|&t| t == first)
        && offsets
            .iter()
            .enumerate()
            .all(|(k, &o)| o == (k + 1) * width)
        && offsets.len() == types.len()
        && connectivity.len() == types.len() * width;
    if !regular {
        return Err(vtk_error(
            "les cellules ne sont pas toutes du même type, ou leurs offsets sont irréguliers"
                .to_string(),
        ));
    }
    if connectivity.iter().any(|&i| i >= nodes.len()) {
        return Err(vtk_error(
            "la connectivité cite un point qui n'existe pas".to_string(),
        ));
    }

    let mut arrays = Vec::new();
    if let Ok(point_data) = section(text, "PointData") {
        for (tag, body) in tags(point_data, "DataArray") {
            if let Some(name) = attribute(tag, "Name") {
                arrays.push(parse_array(tag, body, name)?);
            }
        }
    }
    let units = text.find("<!-- units: ").and_then(|i| {
        let rest = &text[i + "<!-- units: ".len()..];
        rest.find(" -->").map(|j| rest[..j].to_string())
    });
    Ok(VtuResults {
        nodes,
        kind,
        cells: connectivity.chunks(width).map(<[usize]>::to_vec).collect(),
        arrays,
        units,
    })
}
