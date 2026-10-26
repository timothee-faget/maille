//! Lecteur minimal de `.vtu` et `.pvd` ASCII : ce que `writer` écrit, et ce que ParaView écrit en
//! ASCII. Un petit scanner de balises suffit, il ne comprend pas le XML en général.

use super::{vtk_error, Entry, PointArray};
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

/// Grandeur aux Nodes `name` d'un `.vtu` ASCII. Erreur, avec les noms présents, si elle manque.
pub fn read_point_data(text: &str, name: &str) -> Result<PointArray, IoError> {
    let start = text.find("<PointData").ok_or_else(|| {
        vtk_error("le .vtu n'a pas de <PointData> : aucune grandeur aux Nodes".to_string())
    })?;
    let section = &text[start..];
    let section = &section[..section.find("</PointData>").unwrap_or(section.len())];

    let mut names = Vec::new();
    for (tag, body) in tags(section, "DataArray") {
        let Some(found) = attribute(tag, "Name") else {
            continue;
        };
        names.push(found);
        if found != name {
            continue;
        }
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
        return Ok(PointArray {
            name: name.to_string(),
            components,
            values,
        });
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
