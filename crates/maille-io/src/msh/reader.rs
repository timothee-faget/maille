//! Lecteur `.msh` 4.1 ASCII : `$MeshFormat`, `$PhysicalNames`, `$Entities`, `$Nodes`, `$Elements`.
//! Les autres sections sont ignorées.

use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::str::FromStr;

use maille_core::mesh::Node;

use super::numbering::{face_count, face_key, face_positions, orient, swap_tet10_mids};
use super::{CellKind, MshFile, Set};
use crate::IoError;

/// Lignes non vides d'une section, avec leur numéro dans le fichier (depuis 1).
type Body<'a> = Vec<(usize, &'a str)>;

/// Entité Gmsh : (dimension, numéro).
type Entity = (usize, i64);

/// Lit un `.msh` 4.1 ASCII. Les erreurs donnent la ligne et ce qu'il faut corriger.
pub fn read(text: &str) -> Result<MshFile, IoError> {
    let sections = split_sections(text)?;
    let section = |name: &'static str| sections.get(name).ok_or(IoError::MissingSection { name });

    check_format(section("MeshFormat")?)?;
    let names = match sections.get("PhysicalNames") {
        Some(body) => parse_names(body)?,
        None => HashMap::new(),
    };
    let entities = match sections.get("Entities") {
        Some(body) => parse_entities(body)?,
        None => HashMap::new(),
    };
    let (points, node_index) = parse_nodes(section("Nodes")?)?;
    let blocks = parse_elements(section("Elements")?, &node_index)?;
    assemble(points, &blocks, &names, &entities)
}

fn syntax(line: usize, message: String) -> IoError {
    IoError::Syntax { line, message }
}

/// Découpe le fichier en sections `$Nom` … `$EndNom`.
fn split_sections(text: &str) -> Result<HashMap<&str, Body<'_>>, IoError> {
    let mut sections = HashMap::new();
    let mut current: Option<(&str, Body)> = None;
    let mut last = 0;
    for (i, raw) in text.lines().enumerate() {
        let (n, line) = (i + 1, raw.trim());
        if line.is_empty() {
            continue;
        }
        last = n;
        match (&mut current, line.strip_prefix('$')) {
            (None, Some(name)) if !name.starts_with("End") => current = Some((name, Vec::new())),
            (None, _) => {}
            (Some((name, _)), Some(end)) if end.strip_prefix("End") == Some(name) => {
                if let Some((name, body)) = current.take() {
                    sections.insert(name, body);
                }
            }
            (Some((name, _)), Some(other)) => {
                return Err(syntax(
                    n,
                    format!("« ${other} » dans la section ${name} : il manque $End{name}"),
                ));
            }
            (Some((_, body)), _) => body.push((n, line)),
        }
    }
    match current {
        Some((name, _)) => Err(syntax(
            last,
            format!("la section ${name} n'est pas fermée : il manque $End{name}"),
        )),
        None => Ok(sections),
    }
}

/// Parcours des lignes d'une section.
struct Cursor<'a> {
    name: &'static str,
    lines: &'a [(usize, &'a str)],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(name: &'static str, body: &'a Body<'a>) -> Self {
        Self {
            name,
            lines: body,
            pos: 0,
        }
    }

    fn next_line(&mut self) -> Result<(usize, &'a str), IoError> {
        let Some(&line) = self.lines.get(self.pos) else {
            let at = self.lines.last().map_or(0, |l| l.0);
            return Err(syntax(
                at,
                format!("la section ${} se termine trop tôt", self.name),
            ));
        };
        self.pos += 1;
        Ok(line)
    }

    fn next(&mut self) -> Result<(usize, Vec<&'a str>), IoError> {
        let (n, line) = self.next_line()?;
        Ok((n, line.split_whitespace().collect()))
    }
}

fn num<T: FromStr>(line: usize, toks: &[&str], i: usize, what: &str) -> Result<T, IoError> {
    let tok = toks
        .get(i)
        .ok_or_else(|| syntax(line, format!("il manque {what}")))?;
    tok.parse().map_err(|_| {
        syntax(
            line,
            format!("{what} : « {tok} » n'est pas un nombre valable"),
        )
    })
}

fn check_format(body: &Body) -> Result<(), IoError> {
    let mut c = Cursor::new("MeshFormat", body);
    let (line, toks) = c.next()?;
    let version = toks.first().copied().unwrap_or("");
    if version != "4.1" {
        return Err(IoError::UnsupportedVersion {
            line,
            found: version.to_string(),
        });
    }
    if toks.get(1) != Some(&"0") {
        return Err(IoError::BinaryFormat { line });
    }
    Ok(())
}

fn parse_names(body: &Body) -> Result<HashMap<Entity, String>, IoError> {
    let mut c = Cursor::new("PhysicalNames", body);
    let (line, toks) = c.next()?;
    let count: usize = num(line, &toks, 0, "le nombre de groupes physiques")?;
    let mut names = HashMap::new();
    for _ in 0..count {
        let (line, text) = c.next_line()?;
        let (Some(open), Some(close)) = (text.find('"'), text.rfind('"')) else {
            return Err(syntax(
                line,
                "le nom du groupe doit être entre guillemets".to_string(),
            ));
        };
        if close <= open {
            return Err(syntax(
                line,
                "le nom du groupe doit être entre guillemets".to_string(),
            ));
        }
        let head: Vec<&str> = text[..open].split_whitespace().collect();
        let dim: usize = num(line, &head, 0, "la dimension du groupe")?;
        let tag: i64 = num(line, &head, 1, "le numéro du groupe")?;
        names.insert((dim, tag), text[open + 1..close].to_string());
    }
    Ok(names)
}

/// Numéros de groupes physiques de chaque entité.
fn parse_entities(body: &Body) -> Result<HashMap<Entity, Vec<i64>>, IoError> {
    let mut c = Cursor::new("Entities", body);
    let (line, toks) = c.next()?;
    let mut counts = [0usize; 4];
    for (dim, n) in counts.iter_mut().enumerate() {
        *n = num(line, &toks, dim, "le nombre d'entités")?;
    }
    let mut entities = HashMap::new();
    for (dim, &n) in counts.iter().enumerate() {
        for _ in 0..n {
            let (line, toks) = c.next()?;
            let tag: i64 = num(line, &toks, 0, "le numéro de l'entité")?;
            // un point a X Y Z, les autres une boîte min/max
            let at = if dim == 0 { 4 } else { 7 };
            let count: usize = num(
                line,
                &toks,
                at,
                "le nombre de groupes physiques de l'entité",
            )?;
            let phys = (0..count)
                .map(|k| num(line, &toks, at + 1 + k, "un numéro de groupe physique"))
                .collect::<Result<Vec<i64>, _>>()?;
            entities.insert((dim, tag), phys);
        }
    }
    Ok(entities)
}

/// Nodes dans l'ordre du fichier, et indice de chaque numéro Gmsh.
fn parse_nodes(body: &Body) -> Result<(Vec<Node>, HashMap<i64, usize>), IoError> {
    let mut c = Cursor::new("Nodes", body);
    let (line, toks) = c.next()?;
    let blocks: usize = num(line, &toks, 0, "le nombre de blocs de Nodes")?;
    let total: usize = num(line, &toks, 1, "le nombre de Nodes")?;
    let mut nodes = Vec::new();
    let mut index = HashMap::new();
    for _ in 0..blocks {
        let (line, toks) = c.next()?;
        let count: usize = num(line, &toks, 3, "le nombre de Nodes du bloc")?;
        let mut tags = Vec::new();
        for _ in 0..count {
            let (line, toks) = c.next()?;
            tags.push((line, num::<i64>(line, &toks, 0, "le numéro du Node")?));
        }
        // 3 coordonnées ; les paramètres d'une entité paramétrique, s'il y en a, sont ignorés
        for (tag_line, tag) in tags {
            let (line, toks) = c.next()?;
            let node = Node::new(
                num(line, &toks, 0, "x")?,
                num(line, &toks, 1, "y")?,
                num(line, &toks, 2, "z")?,
            );
            if !(node.x.is_finite() && node.y.is_finite() && node.z.is_finite()) {
                return Err(syntax(
                    line,
                    format!("une coordonnée du Node {tag} n'est pas finie (NaN ou infinie)"),
                ));
            }
            if index.insert(tag, nodes.len()).is_some() {
                return Err(IoError::DuplicateNodeTag {
                    line: tag_line,
                    tag,
                });
            }
            nodes.push(node);
        }
    }
    if nodes.len() != total {
        return Err(syntax(
            line,
            format!("{total} Nodes annoncés, {} trouvés", nodes.len()),
        ));
    }
    Ok((nodes, index))
}

/// Dimension et nombre de Nodes d'un type d'élément Gmsh géré (point, ligne, triangle, tétraèdre).
fn element_info(gmsh_type: u32) -> Option<(usize, usize)> {
    match gmsh_type {
        15 => Some((0, 1)),
        1 => Some((1, 2)),
        8 => Some((1, 3)),
        2 => Some((2, 3)),
        9 => Some((2, 6)),
        4 => Some((3, 4)),
        11 => Some((3, 10)),
        _ => None,
    }
}

/// Bloc d'éléments d'une entité.
struct Block {
    entity: Entity,
    gmsh_type: u32,
    /// Dimension de l'élément.
    dim: usize,
    /// (ligne, Nodes en indices) de chaque élément, dans la numérotation de Gmsh.
    elements: Vec<(usize, Vec<usize>)>,
}

fn parse_elements(body: &Body, node_index: &HashMap<i64, usize>) -> Result<Vec<Block>, IoError> {
    let mut c = Cursor::new("Elements", body);
    let (line, toks) = c.next()?;
    let blocks: usize = num(line, &toks, 0, "le nombre de blocs d'éléments")?;
    let mut out = Vec::new();
    for _ in 0..blocks {
        let (line, toks) = c.next()?;
        let entity_dim: usize = num(line, &toks, 0, "la dimension de l'entité")?;
        let entity_tag: i64 = num(line, &toks, 1, "le numéro de l'entité")?;
        let gmsh_type: u32 = num(line, &toks, 2, "le type d'élément")?;
        let count: usize = num(line, &toks, 3, "le nombre d'éléments du bloc")?;
        let Some((dim, width)) = element_info(gmsh_type) else {
            return Err(IoError::UnsupportedElement { line, gmsh_type });
        };
        let mut elements = Vec::new();
        for _ in 0..count {
            let (line, toks) = c.next()?;
            if toks.len() != width + 1 {
                return Err(syntax(
                    line,
                    format!(
                        "l'élément de type {gmsh_type} doit avoir {width} Nodes, il en a {}",
                        toks.len().saturating_sub(1)
                    ),
                ));
            }
            let mut nodes = Vec::with_capacity(width);
            for k in 1..=width {
                let tag: i64 = num(line, &toks, k, "le numéro d'un Node")?;
                let i = *node_index
                    .get(&tag)
                    .ok_or(IoError::UnknownNode { line, tag })?;
                nodes.push(i);
            }
            elements.push((line, nodes));
        }
        out.push(Block {
            entity: (entity_dim, entity_tag),
            gmsh_type,
            dim,
            elements,
        });
    }
    Ok(out)
}

/// Éléments → cellules, Sets.
fn assemble(
    points: Vec<Node>,
    blocks: &[Block],
    names: &HashMap<Entity, String>,
    entities: &HashMap<Entity, Vec<i64>>,
) -> Result<MshFile, IoError> {
    let has = |types: [u32; 2]| blocks.iter().any(|b| types.contains(&b.gmsh_type));
    let dim = if has([4, 11]) {
        3
    } else if has([2, 9]) {
        2
    } else {
        return Err(IoError::NoCells);
    };
    if dim == 2 {
        check_planar(&points)?;
    }

    let mut kind: Option<CellKind> = None;
    let mut cells: Vec<Vec<usize>> = Vec::new();
    let mut entity_cells: HashMap<Entity, Vec<usize>> = HashMap::new();
    for b in blocks.iter().filter(|b| b.dim == dim) {
        let Some(k) = CellKind::from_gmsh_type(b.gmsh_type) else {
            continue;
        };
        match kind {
            Some(first) if first != k => {
                return Err(IoError::MixedCellKinds { first, found: k });
            }
            _ => kind = Some(k),
        }
        for (_, nodes) in &b.elements {
            let mut nodes = nodes.clone();
            if k == CellKind::Tet10 {
                swap_tet10_mids(&mut nodes);
            }
            orient(k, &mut nodes, &points);
            entity_cells.entry(b.entity).or_default().push(cells.len());
            cells.push(nodes);
        }
    }
    let kind = kind.ok_or(IoError::NoCells)?;

    // la première face rencontrée gagne : un bord interne est rattaché à l'Element de plus petit indice
    let mut faces: HashMap<[usize; 3], (usize, usize)> = HashMap::new();
    for (e, cell) in cells.iter().enumerate() {
        for f in 0..face_count(kind) {
            if let Some(pos) = face_positions(kind, f) {
                let vertices: Vec<usize> = pos[..dim].iter().map(|&p| cell[p]).collect();
                faces.entry(face_key(&vertices)).or_insert((e, f));
            }
        }
    }

    let mut sets: BTreeMap<String, Set> = BTreeMap::new();
    let mut owner: HashMap<String, Entity> = HashMap::new();
    for b in blocks {
        let Some(physical) = entities.get(&b.entity) else {
            continue;
        };
        for &tag in physical {
            let group = (b.entity.0, tag);
            let name = names
                .get(&group)
                .cloned()
                .unwrap_or_else(|| tag.to_string());
            match owner.entry(name.clone()) {
                Entry::Occupied(o) if *o.get() != group => {
                    return Err(IoError::DuplicateSetName { name });
                }
                Entry::Occupied(_) => {}
                Entry::Vacant(v) => {
                    v.insert(group);
                }
            }
            let set = sets.entry(name.clone()).or_default();
            if b.dim == dim {
                set.elements
                    .extend(entity_cells.get(&b.entity).into_iter().flatten());
            } else if b.dim + 1 == dim {
                for (line, nodes) in &b.elements {
                    let found = faces.get(&face_key(&nodes[..dim])).copied();
                    let (e, f) = found.ok_or_else(|| IoError::BoundaryFaceNotFound {
                        line: *line,
                        set: name.clone(),
                    })?;
                    set.faces.push((e, f));
                    if let Some(pos) = face_positions(kind, f) {
                        set.nodes.extend(pos.iter().map(|&p| cells[e][p]));
                    }
                }
            } else {
                set.nodes
                    .extend(b.elements.iter().flat_map(|(_, n)| n.iter().copied()));
            }
        }
    }
    for set in sets.values_mut() {
        set.nodes.sort_unstable();
        set.nodes.dedup();
        set.faces.sort_unstable();
        set.faces.dedup();
        set.elements.sort_unstable();
        set.elements.dedup();
    }
    Ok(MshFile {
        nodes: points,
        kind,
        cells,
        sets,
    })
}

/// Un maillage de triangles doit être dans le plan z = 0 (à la précision de sa taille près).
fn check_planar(points: &[Node]) -> Result<(), IoError> {
    let span = points
        .iter()
        .flat_map(|p| [p.x.abs(), p.y.abs()])
        .fold(0.0_f64, f64::max);
    let tolerance = 1e-9 * span.max(f64::MIN_POSITIVE);
    match points.iter().find(|p| p.z.abs() > tolerance) {
        Some(p) => Err(IoError::NonPlanar { z: p.z }),
        None => Ok(()),
    }
}
