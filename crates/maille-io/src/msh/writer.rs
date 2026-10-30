//! Écrivain `.msh` 4.1 ASCII.
//!
//! Chaque Set devient un groupe physique. Les entités Gmsh se regroupent par appartenance : les
//! Elements qui sont dans les mêmes Sets forment une entité, ce qui permet à un Element d'être
//! dans plusieurs Sets. Un Set de Nodes s'écrit en points, un Set de faces en éléments de bord
//! (ses Nodes se déduisent des faces), un Set d'arêtes (3D) en lignes (ses Nodes se déduisent des
//! arêtes), un Set d'Elements en éléments de la dimension du maillage.

use std::collections::BTreeMap;
use std::fmt::Write;

use super::numbering::{edge_count, edge_positions, face_count, face_positions, swap_tet10_mids};
use super::{CellKind, MshFile, Set};
use crate::IoError;

/// Entité Gmsh à écrire : groupes physiques, type d'élément et Nodes (indices depuis 0, ordre Gmsh).
struct Entity {
    physical: Vec<usize>,
    gmsh_type: u32,
    elements: Vec<Vec<usize>>,
}

fn invalid(problem: String) -> IoError {
    IoError::Invalid { problem }
}

/// Écrit le fichier. `Invalid` si le `MshFile` est incohérent (indice hors limites, Set mélangeant
/// Elements et faces, nom avec des guillemets…).
pub fn write(file: &MshFile) -> Result<String, IoError> {
    validate(file)?;
    let dim = file.kind.dim();

    // groupes physiques numérotés dans l'ordre des noms ; (dimension, nom) pour $PhysicalNames
    let physical: Vec<(&String, &Set)> = file.sets.iter().collect();
    let tag_dim = |set: &Set| {
        if !set.elements.is_empty() {
            dim
        } else if !set.faces.is_empty() {
            dim - 1
        } else if !set.edges.is_empty() {
            1
        } else {
            0
        }
    };

    let mut cell_sets = vec![Vec::new(); file.cells.len()];
    let mut face_sets: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    let mut edge_sets: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    let mut point_sets: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, (_, set)) in physical.iter().enumerate() {
        let tag = i + 1;
        if !set.elements.is_empty() {
            for &e in &set.elements {
                cell_sets[e].push(tag);
            }
        } else if !set.faces.is_empty() {
            for &face in &set.faces {
                face_sets.entry(face).or_default().push(tag);
            }
        } else if !set.edges.is_empty() {
            for &edge in &set.edges {
                edge_sets.entry(edge).or_default().push(tag);
            }
        } else {
            for &n in &set.nodes {
                point_sets.entry(n).or_default().push(tag);
            }
        }
    }

    let gmsh_cell = |cell: &Vec<usize>| {
        let mut nodes = cell.clone();
        if file.kind == CellKind::Tet10 {
            swap_tet10_mids(&mut nodes);
        }
        nodes
    };
    // les Elements gardent leur ordre : une entité par suite d'Elements consécutifs de mêmes Sets
    let mut cell_entities: Vec<Entity> = Vec::new();
    for (cell, tags) in file.cells.iter().zip(&mut cell_sets) {
        tags.dedup();
        match cell_entities.last_mut() {
            Some(last) if last.physical == *tags => last.elements.push(gmsh_cell(cell)),
            _ => cell_entities.push(Entity {
                physical: tags.clone(),
                gmsh_type: file.kind.gmsh_type(),
                elements: vec![gmsh_cell(cell)],
            }),
        }
    }
    let mut face_groups: BTreeMap<Vec<usize>, Vec<Vec<usize>>> = BTreeMap::new();
    for ((e, f), mut tags) in face_sets {
        tags.dedup();
        let cell = &file.cells[e];
        let nodes = face_positions(file.kind, f)
            .map(|pos| pos.iter().map(|&p| cell[p]).collect())
            .unwrap_or_default();
        face_groups.entry(tags).or_default().push(nodes);
    }

    let mut edge_groups: BTreeMap<Vec<usize>, Vec<Vec<usize>>> = BTreeMap::new();
    for ((e, k), mut tags) in edge_sets {
        tags.dedup();
        let cell = &file.cells[e];
        let nodes = edge_positions(file.kind, k)
            .map(|pos| pos.iter().map(|&p| cell[p]).collect())
            .unwrap_or_default();
        edge_groups.entry(tags).or_default().push(nodes);
    }

    let face_entities: Vec<Entity> = face_groups
        .into_iter()
        .map(|(physical, elements)| Entity {
            physical,
            gmsh_type: file.kind.face_gmsh_type(),
            elements,
        })
        .collect();
    let edge_entities: Vec<Entity> = edge_groups
        .into_iter()
        .map(|(physical, elements)| Entity {
            physical,
            gmsh_type: file.kind.edge_gmsh_type(),
            elements,
        })
        .collect();
    let point_entities: Vec<Entity> = point_sets
        .into_iter()
        .map(|(n, mut physical)| {
            physical.dedup();
            Entity {
                physical,
                gmsh_type: 15,
                elements: vec![vec![n]],
            }
        })
        .collect();

    let mut out = String::new();
    let _ = writeln!(out, "$MeshFormat\n4.1 0 8\n$EndMeshFormat");

    let _ = writeln!(out, "$PhysicalNames\n{}", physical.len());
    for (i, (name, set)) in physical.iter().enumerate() {
        let _ = writeln!(out, "{} {} \"{name}\"", tag_dim(set), i + 1);
    }
    let _ = writeln!(out, "$EndPhysicalNames");

    // entités par dimension : points (0), lignes (1, arêtes d'un maillage 3D), faces (dim - 1),
    // cellules (dim)
    let mut by_dim: [&[Entity]; 4] = [&[]; 4];
    by_dim[0] = &point_entities;
    by_dim[1] = &edge_entities;
    by_dim[dim - 1] = &face_entities;
    by_dim[dim] = &cell_entities;
    let _ = writeln!(
        out,
        "$Entities\n{} {} {} {}",
        by_dim[0].len(),
        by_dim[1].len(),
        by_dim[2].len(),
        by_dim[3].len()
    );
    for (d, list) in by_dim.iter().enumerate() {
        for (k, entity) in list.iter().enumerate() {
            let (lo, hi) = bounding_box(file, entity);
            let tags: String = entity.physical.iter().map(|t| format!(" {t}")).collect();
            if d == 0 {
                let _ = writeln!(
                    out,
                    "{} {} {} {} {}{tags}",
                    k + 1,
                    lo[0],
                    lo[1],
                    lo[2],
                    entity.physical.len()
                );
            } else {
                let _ = writeln!(
                    out,
                    "{} {} {} {} {} {} {} {}{tags} 0",
                    k + 1,
                    lo[0],
                    lo[1],
                    lo[2],
                    hi[0],
                    hi[1],
                    hi[2],
                    entity.physical.len()
                );
            }
        }
    }
    let _ = writeln!(out, "$EndEntities");

    // tous les Nodes dans un seul bloc, rattaché à la première entité de cellules
    let n = file.nodes.len();
    let _ = writeln!(out, "$Nodes\n1 {n} 1 {n}\n{dim} 1 0 {n}");
    for i in 1..=n {
        let _ = writeln!(out, "{i}");
    }
    for p in &file.nodes {
        let _ = writeln!(out, "{} {} {}", p.x, p.y, p.z);
    }
    let _ = writeln!(out, "$EndNodes");

    let blocks: Vec<(usize, usize, &Entity)> = by_dim
        .iter()
        .enumerate()
        .flat_map(|(d, list)| list.iter().enumerate().map(move |(k, e)| (d, k + 1, e)))
        .collect();
    let total: usize = blocks.iter().map(|b| b.2.elements.len()).sum();
    let _ = writeln!(out, "$Elements\n{} {total} 1 {total}", blocks.len());
    let mut tag = 1;
    for (d, k, entity) in blocks {
        let _ = writeln!(
            out,
            "{d} {k} {} {}",
            entity.gmsh_type,
            entity.elements.len()
        );
        for nodes in &entity.elements {
            let list: String = nodes.iter().map(|&i| format!(" {}", i + 1)).collect();
            let _ = writeln!(out, "{tag}{list}");
            tag += 1;
        }
    }
    let _ = writeln!(out, "$EndElements");
    Ok(out)
}

fn bounding_box(file: &MshFile, entity: &Entity) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &i in entity.elements.iter().flatten() {
        let p = file.nodes[i];
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[k] = lo[k].min(v);
            hi[k] = hi[k].max(v);
        }
    }
    (lo, hi)
}

fn validate(file: &MshFile) -> Result<(), IoError> {
    if file.nodes.is_empty() || file.cells.is_empty() {
        return Err(invalid(
            "le maillage n'a aucun Node ou aucun Element".to_string(),
        ));
    }
    for (e, cell) in file.cells.iter().enumerate() {
        if cell.len() != file.kind.node_count() {
            return Err(invalid(format!(
                "l'Element {e} a {} Nodes, {:?} en demande {}",
                cell.len(),
                file.kind,
                file.kind.node_count()
            )));
        }
        if let Some(&n) = cell.iter().find(|&&n| n >= file.nodes.len()) {
            return Err(invalid(format!(
                "l'Element {e} utilise le Node {n}, qui n'existe pas"
            )));
        }
    }
    for (name, set) in &file.sets {
        if name.contains(['"', '\n', '\r']) {
            return Err(invalid(format!(
                "le nom de Set « {name} » ne peut pas contenir de guillemet ni de saut de ligne"
            )));
        }
        if !set.elements.is_empty()
            && !(set.faces.is_empty() && set.edges.is_empty() && set.nodes.is_empty())
        {
            return Err(invalid(format!(
                "le Set « {name} » mélange des Elements avec des faces, des arêtes ou des Nodes : un Set a un seul rôle"
            )));
        }
        if !set.edges.is_empty() && !set.faces.is_empty() {
            return Err(invalid(format!(
                "le Set « {name} » mélange des faces et des arêtes : un Set a un seul rôle"
            )));
        }
        if !set.edges.is_empty() && file.kind.dim() != 3 {
            return Err(invalid(format!(
                "le Set « {name} » a des arêtes, qui n'existent en Set que dans un maillage 3D (en 2D, une arête de bord est une face)"
            )));
        }
        if set.elements.is_empty()
            && set.faces.is_empty()
            && set.edges.is_empty()
            && set.nodes.is_empty()
        {
            return Err(invalid(format!("le Set « {name} » est vide")));
        }
        if let Some(&e) = set.elements.iter().find(|&&e| e >= file.cells.len()) {
            return Err(invalid(format!(
                "le Set « {name} » cite l'Element {e}, qui n'existe pas"
            )));
        }
        for &(e, f) in &set.faces {
            if e >= file.cells.len() || f >= face_count(file.kind) {
                return Err(invalid(format!(
                    "le Set « {name} » cite la face {f} de l'Element {e}, qui n'existe pas"
                )));
            }
        }
        for &(e, k) in &set.edges {
            if e >= file.cells.len() || k >= edge_count(file.kind) {
                return Err(invalid(format!(
                    "le Set « {name} » cite l'arête {k} de l'Element {e}, qui n'existe pas"
                )));
            }
        }
        if let Some(&n) = set.nodes.iter().find(|&&n| n >= file.nodes.len()) {
            return Err(invalid(format!(
                "le Set « {name} » cite le Node {n}, qui n'existe pas"
            )));
        }
    }
    Ok(())
}
