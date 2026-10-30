//! Désigner des morceaux du Mesh lui-même, sans passer par la géométrie : un cadre tracé à l'écran
//! (Nodes ou Elements), ou une région de la peau qu'on étend d'une face à l'autre tant que l'angle
//! entre voisines reste petit. Sans fenêtre : comme `sets`, on projette avec la caméra du viewer.
//! Le résultat est un `Set` du `.msh` (Nodes seuls, Elements seuls, ou faces avec leurs Nodes).

use std::collections::{BTreeSet, HashMap, VecDeque};

use maille_io::msh::{edge_positions, face_positions, MshFile, Set};

use crate::sets::{segment_distance, triangle_depth, EDGE_TOLERANCE_PX};

type Point = [f64; 3];

/// Cadre en pixels de la texture (origine en haut à gauche), coins dans n'importe quel ordre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    min: [f64; 2],
    max: [f64; 2],
}

impl Rect {
    pub fn from_corners(a: [f64; 2], b: [f64; 2]) -> Self {
        Self {
            min: [a[0].min(b[0]), a[1].min(b[1])],
            max: [a[0].max(b[0]), a[1].max(b[1])],
        }
    }

    fn contains(&self, p: Point) -> bool {
        (self.min[0]..=self.max[0]).contains(&p[0]) && (self.min[1]..=self.max[1]).contains(&p[1])
    }
}

fn position(file: &MshFile, node: usize) -> Point {
    let p = &file.nodes[node];
    [p.x, p.y, p.z]
}

/// Les sommets d'un Element (3 pour un triangle, 4 pour un tétraèdre), sans les milieux.
fn corners<'a>(file: &MshFile, cell: &'a [usize]) -> &'a [usize] {
    &cell[..file.kind.dim() + 1]
}

/// Les Nodes dont la projection tombe dans le cadre. En 3D, le cadre prend aussi ceux qui sont
/// cachés derrière la pièce.
pub fn nodes_in_rect(file: &MshFile, project: &dyn Fn(Point) -> Point, rect: Rect) -> Vec<usize> {
    (0..file.nodes.len())
        .filter(|&n| rect.contains(project(position(file, n))))
        .collect()
}

/// Les Elements dont le centre (moyenne des sommets) se projette dans le cadre.
pub fn elements_in_rect(
    file: &MshFile,
    project: &dyn Fn(Point) -> Point,
    rect: Rect,
) -> Vec<usize> {
    file.cells
        .iter()
        .enumerate()
        .filter(|(_, cell)| rect.contains(project(centroid(file, cell))))
        .map(|(e, _)| e)
        .collect()
}

fn centroid(file: &MshFile, cell: &[usize]) -> Point {
    let c = corners(file, cell);
    let mut sum = [0.0; 3];
    for &n in c {
        let p = position(file, n);
        (0..3).for_each(|k| sum[k] += p[k]);
    }
    sum.map(|v| v / c.len() as f64)
}

/// Un Set de Nodes seuls.
pub fn node_set(nodes: &[usize]) -> Set {
    let mut nodes = nodes.to_vec();
    nodes.sort_unstable();
    nodes.dedup();
    Set {
        nodes,
        ..Set::default()
    }
}

/// Un Set d'Elements seuls (le `.msh` ne mélange pas les rôles).
pub fn element_set(elements: &[usize]) -> Set {
    let mut elements = elements.to_vec();
    elements.sort_unstable();
    elements.dedup();
    Set {
        elements,
        ..Set::default()
    }
}

/// Une face de la peau : celle d'un seul Element.
#[derive(Debug, Clone, PartialEq)]
pub struct SkinFace {
    pub element: usize,
    pub face: usize,
    /// Ses sommets : 2 en 2D, 3 en 3D.
    pub corners: Vec<usize>,
}

/// La peau du maillage : les faces qui n'appartiennent qu'à un Element.
pub fn skin(file: &MshFile) -> Vec<SkinFace> {
    let dim = file.kind.dim();
    let mut seen: HashMap<Vec<usize>, Vec<SkinFace>> = HashMap::new();
    for (element, cell) in file.cells.iter().enumerate() {
        for face in 0.. {
            let Some(positions) = face_positions(file.kind, face) else {
                break;
            };
            let corners: Vec<usize> = positions[..dim].iter().map(|&p| cell[p]).collect();
            let mut key = corners.clone();
            key.sort_unstable();
            seen.entry(key).or_default().push(SkinFace {
                element,
                face,
                corners,
            });
        }
    }
    let mut skin: Vec<SkinFace> = seen
        .into_values()
        .filter(|faces| faces.len() == 1)
        .flatten()
        .collect();
    skin.sort_by_key(|f| (f.element, f.face));
    skin
}

/// La face de la peau sous le pixel `cursor` : la plus proche de l'observateur en 3D, le segment le
/// plus proche (à moins de `EDGE_TOLERANCE_PX`) en 2D.
pub fn pick_skin(
    file: &MshFile,
    skin: &[SkinFace],
    project: &dyn Fn(Point) -> Point,
    cursor: [f64; 2],
) -> Option<usize> {
    let at = |f: &SkinFace, k: usize| project(position(file, f.corners[k]));
    let scored = skin.iter().enumerate().filter_map(|(i, f)| {
        if file.kind.dim() == 2 {
            let d = segment_distance(cursor, at(f, 0), at(f, 1));
            (d <= EDGE_TOLERANCE_PX).then_some((d, i))
        } else {
            triangle_depth(cursor, [at(f, 0), at(f, 1), at(f, 2)]).map(|depth| (depth, i))
        }
    });
    scored.min_by(|a, b| a.0.total_cmp(&b.0)).map(|(_, i)| i)
}

fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn unit(v: Point) -> Point {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n > 0.0 {
        v.map(|c| c / n)
    } else {
        v
    }
}

/// Normale unitaire sortante d'une face de la peau (en 2D, la normale au segment, dans le plan).
fn outward_normal(file: &MshFile, face: &SkinFace) -> Point {
    let p: Vec<Point> = face.corners.iter().map(|&n| position(file, n)).collect();
    let raw = if p.len() == 2 {
        let t = sub(p[1], p[0]);
        [t[1], -t[0], 0.0]
    } else {
        let (a, b) = (sub(p[1], p[0]), sub(p[2], p[0]));
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let mid = {
        let mut m = [0.0; 3];
        for q in &p {
            (0..3).for_each(|k| m[k] += q[k] / p.len() as f64);
        }
        m
    };
    let inside = sub(mid, centroid(file, &file.cells[face.element]));
    let n = unit(raw);
    if n[0] * inside[0] + n[1] * inside[1] + n[2] * inside[2] < 0.0 {
        n.map(|c| -c)
    } else {
        n
    }
}

/// Le Set de faces obtenu en partant de la face `seed` de la peau et en passant aux voisines (qui
/// partagent une arête en 3D, un sommet en 2D) tant que leur normale ne tourne pas de plus de
/// `max_angle` degrés par rapport à la face d'où l'on vient : un plan ou une surface lisse se
/// sélectionne d'un clic, un coin vif arrête la région. Ses Nodes (milieux compris) sont dans le Set.
pub fn angle_region(file: &MshFile, skin: &[SkinFace], seed: usize, max_angle: f64) -> Set {
    let min_cos = max_angle.to_radians().cos();
    let normals: Vec<Point> = skin.iter().map(|f| outward_normal(file, f)).collect();
    // les voisines : on indexe les faces par leurs arêtes (3D) ou leurs sommets (2D)
    let mut by_part: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (i, f) in skin.iter().enumerate() {
        for part in parts(&f.corners) {
            by_part.entry(part).or_default().push(i);
        }
    }
    let mut region = BTreeSet::from([seed]);
    let mut queue = VecDeque::from([seed]);
    while let Some(i) = queue.pop_front() {
        for part in parts(&skin[i].corners) {
            for &j in &by_part[&part] {
                let cos = (0..3).map(|k| normals[i][k] * normals[j][k]).sum::<f64>();
                if cos >= min_cos && region.insert(j) {
                    queue.push_back(j);
                }
            }
        }
    }
    let mut faces: Vec<(usize, usize)> = region
        .iter()
        .map(|&i| (skin[i].element, skin[i].face))
        .collect();
    faces.sort_unstable();
    let mut nodes = Vec::new();
    for &(element, face) in &faces {
        if let Some(positions) = face_positions(file.kind, face) {
            nodes.extend(positions.iter().map(|&p| file.cells[element][p]));
        }
    }
    nodes.sort_unstable();
    nodes.dedup();
    Set {
        nodes,
        faces,
        ..Set::default()
    }
}

/// Ce que deux faces voisines ont en commun : chaque arête d'un triangle, chaque sommet d'un segment.
fn parts(corners: &[usize]) -> Vec<Vec<usize>> {
    let sorted = |mut v: Vec<usize>| {
        v.sort_unstable();
        v
    };
    if corners.len() == 2 {
        corners.iter().map(|&c| vec![c]).collect()
    } else {
        (0..3)
            .map(|k| sorted(vec![corners[k], corners[(k + 1) % 3]]))
            .collect()
    }
}

/// Les triangles à surligner pour une sélection de Nodes : ceux (en 2D, les Elements ; en 3D, les
/// faces de la peau) qui en touchent au moins un. Un cadre posé sur une ligne de Nodes se voit ainsi.
pub fn highlight_nodes(file: &MshFile, skin: &[SkinFace], nodes: &[usize]) -> Vec<[Point; 3]> {
    let chosen: BTreeSet<usize> = nodes.iter().copied().collect();
    let tri = |c: &[usize]| [0, 1, 2].map(|k| position(file, c[k]));
    if file.kind.dim() == 2 {
        file.cells
            .iter()
            .filter(|c| corners(file, c).iter().any(|n| chosen.contains(n)))
            .map(|c| tri(c))
            .collect()
    } else {
        skin.iter()
            .filter(|f| f.corners.iter().any(|n| chosen.contains(n)))
            .map(|f| tri(&f.corners))
            .collect()
    }
}

/// Les triangles à surligner pour une sélection d'Elements : eux-mêmes en 2D, leurs faces sur la
/// peau en 3D (l'intérieur ne se voit pas).
pub fn highlight_elements(
    file: &MshFile,
    skin: &[SkinFace],
    elements: &[usize],
) -> Vec<[Point; 3]> {
    let chosen: BTreeSet<usize> = elements.iter().copied().collect();
    let tri = |c: &[usize]| [0, 1, 2].map(|k| position(file, c[k]));
    if file.kind.dim() == 2 {
        chosen.iter().map(|&e| tri(&file.cells[e])).collect()
    } else {
        skin.iter()
            .filter(|f| chosen.contains(&f.element))
            .map(|f| tri(&f.corners))
            .collect()
    }
}

/// Une arête vive de la peau d'un maillage 3D : celle que se partagent deux faces dont les normales
/// font plus que l'angle demandé (le pli d'un cube, le bord d'un trou).
#[derive(Debug, Clone, PartialEq)]
pub struct SharpEdge {
    /// Ses deux sommets, `nodes[0] < nodes[1]`.
    pub nodes: [usize; 2],
}

/// Les arêtes vives de la peau : entre deux faces dont l'angle des normales dépasse `max_angle`
/// degrés. Vide en 2D, où un bord est déjà une face.
pub fn sharp_edges(file: &MshFile, skin: &[SkinFace], max_angle: f64) -> Vec<SharpEdge> {
    if file.kind.dim() != 3 {
        return Vec::new();
    }
    let min_cos = max_angle.to_radians().cos();
    let normals: Vec<Point> = skin.iter().map(|f| outward_normal(file, f)).collect();
    let mut by_edge: HashMap<[usize; 2], Vec<usize>> = HashMap::new();
    for (i, f) in skin.iter().enumerate() {
        for part in parts(&f.corners) {
            by_edge.entry([part[0], part[1]]).or_default().push(i);
        }
    }
    let mut edges: Vec<SharpEdge> = by_edge
        .into_iter()
        .filter(|(_, faces)| {
            faces.len() == 2
                && (0..3)
                    .map(|k| normals[faces[0]][k] * normals[faces[1]][k])
                    .sum::<f64>()
                    < min_cos
        })
        .map(|(nodes, _)| SharpEdge { nodes })
        .collect();
    edges.sort_by_key(|e| e.nodes);
    edges
}

/// Point de `a`–`b` le plus proche de `p` à l'écran : sa distance et sa profondeur.
fn nearest_on_segment(p: [f64; 2], a: Point, b: Point) -> (f64, [f64; 2], f64) {
    let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
    let len2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if len2 > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let at = [a[0] + t * ab[0], a[1] + t * ab[1]];
    let depth = a[2] + t * (b[2] - a[2]);
    ((p[0] - at[0]).hypot(p[1] - at[1]), at, depth)
}

/// L'arête vive sous le pixel `cursor` : la plus proche (à moins de `EDGE_TOLERANCE_PX`) parmi celles
/// qu'aucune face ne cache.
pub fn pick_edge(
    file: &MshFile,
    skin: &[SkinFace],
    edges: &[SharpEdge],
    project: &dyn Fn(Point) -> Point,
    cursor: [f64; 2],
) -> Option<usize> {
    let tri = |f: &SkinFace| [0, 1, 2].map(|k| project(position(file, f.corners[k])));
    let mut best: Option<(f64, usize)> = None;
    for (i, e) in edges.iter().enumerate() {
        let (a, b) = (
            project(position(file, e.nodes[0])),
            project(position(file, e.nodes[1])),
        );
        let (distance, at, depth) = nearest_on_segment(cursor, a, b);
        if distance > EDGE_TOLERANCE_PX || best.is_some_and(|(d, _)| d <= distance) {
            continue;
        }
        // une face plus proche de l'observateur au même pixel cache l'arête
        let front = skin
            .iter()
            .filter_map(|f| triangle_depth(at, tri(f)))
            .fold(f64::INFINITY, f64::min);
        if depth <= front + 1e-6 + 1e-4 * front.abs() {
            best = Some((distance, i));
        }
    }
    best.map(|(_, i)| i)
}

/// Le Set d'arêtes qui prolonge l'arête vive `seed` : on avance de sommet en sommet tant que le
/// sommet touche exactement deux arêtes vives, et on s'arrête à un coin (trois arêtes ou plus) ou à un
/// bout. Un pli de cube donne ses segments, le bord d'un trou toute sa boucle.
pub fn edge_chain(file: &MshFile, edges: &[SharpEdge], seed: usize) -> Set {
    let mut at_node: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, e) in edges.iter().enumerate() {
        for n in e.nodes {
            at_node.entry(n).or_default().push(i);
        }
    }
    let mut chain = BTreeSet::from([seed]);
    let mut queue = VecDeque::from([seed]);
    while let Some(i) = queue.pop_front() {
        for n in edges[i].nodes {
            if let [a, b] = at_node[&n][..] {
                let next = if a == i { b } else { a };
                if chain.insert(next) {
                    queue.push_back(next);
                }
            }
        }
    }
    let pairs: Vec<[usize; 2]> = chain.iter().map(|&i| edges[i].nodes).collect();
    edge_set(file, &pairs)
}

/// Le Set d'arêtes de ces paires de sommets, chacune rattachée à l'Element de plus petit indice qui la
/// porte (comme le fait la lecture d'un `.msh`), avec les Nodes des arêtes, milieux compris.
pub fn edge_set(file: &MshFile, pairs: &[[usize; 2]]) -> Set {
    let wanted: BTreeSet<[usize; 2]> = pairs
        .iter()
        .map(|p| [p[0].min(p[1]), p[0].max(p[1])])
        .collect();
    let mut found: HashMap<[usize; 2], (usize, usize)> = HashMap::new();
    let mut nodes = Vec::new();
    for (element, cell) in file.cells.iter().enumerate() {
        for k in 0.. {
            let Some(positions) = edge_positions(file.kind, k) else {
                break;
            };
            let (a, b) = (cell[positions[0]], cell[positions[1]]);
            let key = [a.min(b), a.max(b)];
            if wanted.contains(&key) && !found.contains_key(&key) {
                found.insert(key, (element, k));
                nodes.extend(positions.iter().map(|&p| cell[p]));
            }
        }
    }
    let mut edges: Vec<(usize, usize)> = found.into_values().collect();
    edges.sort_unstable();
    nodes.sort_unstable();
    nodes.dedup();
    Set {
        nodes,
        edges,
        ..Set::default()
    }
}

/// Les triangles à surligner pour des arêtes : les faces de la peau qui les touchent, un ruban de
/// chaque côté (un trait d'un pixel ne se verrait pas).
pub fn highlight_edges(
    file: &MshFile,
    skin: &[SkinFace],
    edges: &[(usize, usize)],
) -> Vec<[Point; 3]> {
    let mut at_node: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, f) in skin.iter().enumerate() {
        for &n in &f.corners {
            at_node.entry(n).or_default().push(i);
        }
    }
    let mut touched = BTreeSet::new();
    for &(element, k) in edges {
        let Some(positions) = edge_positions(file.kind, k) else {
            continue;
        };
        let (a, b) = (
            file.cells[element][positions[0]],
            file.cells[element][positions[1]],
        );
        for &i in at_node.get(&a).into_iter().flatten() {
            if skin[i].corners.contains(&b) {
                touched.insert(i);
            }
        }
    }
    touched
        .into_iter()
        .map(|i| [0, 1, 2].map(|k| position(file, skin[i].corners[k])))
        .collect()
}

/// Les triangles à surligner pour un Set du `.msh`, quel que soit son rôle (Elements, faces, Nodes).
pub fn highlight_set(file: &MshFile, skin: &[SkinFace], set: &Set) -> Vec<[Point; 3]> {
    if !set.elements.is_empty() {
        highlight_elements(file, skin, &set.elements)
    } else if !set.faces.is_empty() {
        let faces: BTreeSet<(usize, usize)> = set.faces.iter().copied().collect();
        let tri = |c: &[usize]| [0, 1, 2].map(|k| position(file, c[k]));
        if file.kind.dim() == 2 {
            let elements: BTreeSet<usize> = faces.iter().map(|f| f.0).collect();
            elements.iter().map(|&e| tri(&file.cells[e])).collect()
        } else {
            skin.iter()
                .filter(|f| faces.contains(&(f.element, f.face)))
                .map(|f| tri(&f.corners))
                .collect()
        }
    } else if !set.edges.is_empty() {
        highlight_edges(file, skin, &set.edges)
    } else {
        highlight_nodes(file, skin, &set.nodes)
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::Path;

    use maille_io::msh::read;

    use super::*;
    use crate::mesh_space::{import, mesh, Imported};

    fn example(name: &str) -> MshFile {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(name);
        read(&fs::read_to_string(path).unwrap()).unwrap()
    }

    /// Vue orthographique de face : 100 pixels par unité, y vers le haut, profondeur −z.
    fn front(p: Point) -> Point {
        [p[0] * 100.0, 300.0 - p[1] * 100.0, -p[2]]
    }

    fn at(x: f64, y: f64) -> [f64; 2] {
        [x * 100.0, 300.0 - y * 100.0]
    }

    #[test]
    fn a_box_picks_the_nodes_and_the_elements_inside() {
        let file = example("plaque.msh"); // 4 × 2
        let strip = Rect::from_corners(at(-0.01, 2.1), at(0.01, -0.1));
        let nodes = nodes_in_rect(&file, &front, strip);
        let left = &file.sets["gauche"];
        assert!(!nodes.is_empty());
        assert!(nodes.iter().all(|&n| file.nodes[n].x.abs() < 1e-9));
        // tous les Nodes du bord gauche, coins compris, et eux seuls (les Nodes milieux aussi)
        assert_eq!(node_set(&nodes).nodes, left.nodes);

        let half = Rect::from_corners(at(-1.0, 3.0), at(2.0, -1.0));
        let elements = elements_in_rect(&file, &front, half);
        let expected: Vec<usize> = (0..file.cells.len())
            .filter(|&e| {
                let c = &file.cells[e][..3];
                c.iter().map(|&n| file.nodes[n].x).sum::<f64>() / 3.0 <= 2.0
            })
            .collect();
        assert_eq!(elements, expected);
        assert!(!elements.is_empty() && elements.len() < file.cells.len());
        // le Set d'Elements s'écrit et se relit tel quel
        let mut with = file.clone();
        with.sets.insert("moitie".into(), element_set(&elements));
        let again = read(&maille_io::msh::write(&with).unwrap()).unwrap();
        assert_eq!(again.sets["moitie"].elements, elements);
    }

    #[test]
    fn in_2d_the_angle_region_stops_at_the_corners() {
        let file = example("plaque.msh");
        let skin = skin(&file);
        // 2 × (4 + 2) arêtes de bord au moins : chaque Set de bord est fait de faces de peau
        let total: usize = ["bas", "droite", "haut", "gauche"]
            .iter()
            .map(|n| file.sets[*n].faces.len())
            .sum();
        assert_eq!(skin.len(), total);
        let i = pick_skin(&file, &skin, &front, at(2.0, 0.04)).expect("le bord bas");
        let region = angle_region(&file, &skin, i, 30.0);
        assert_eq!(region.faces, file.sets["bas"].faces);
        assert_eq!(region.nodes, file.sets["bas"].nodes);
        // à 100° les coins à 90° se franchissent : tout le contour
        assert_eq!(angle_region(&file, &skin, i, 100.0).faces.len(), total);
        // au milieu de la plaque : aucun bord à portée
        assert_eq!(pick_skin(&file, &skin, &front, at(2.0, 1.0)), None);
    }

    #[test]
    fn in_3d_the_angle_region_follows_a_smooth_surface_and_stops_at_the_edges() {
        // un huitième de sphère creuse, maillé depuis le STEP (le .msh n'est pas versionné) : la
        // surface extérieure est lisse, ses bords sont à 90° des plans de symétrie
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/sphere.step");
        let Imported::Cad(cad) = import(&path).unwrap() else {
            panic!("sphere.step est une pièce");
        };
        let file = mesh(&cad, 5.0, false, &Default::default()).unwrap().file;
        let skin = skin(&file);
        let outer = &file.sets["exterieur"];
        let seed = skin
            .iter()
            .position(|f| outer.faces.contains(&(f.element, f.face)))
            .unwrap();
        let region = angle_region(&file, &skin, seed, 30.0);
        assert_eq!(region.faces, outer.faces);
        assert_eq!(region.nodes, outer.nodes);
        // le Set de faces s'écrit et se relit
        let mut with = file.clone();
        with.sets.insert("region".into(), region.clone());
        let again = read(&maille_io::msh::write(&with).unwrap()).unwrap();
        assert_eq!(again.sets["region"].faces, region.faces);
    }

    #[test]
    fn in_3d_a_flat_face_is_a_plane_and_the_nearest_face_wins() {
        let file = example("cube.msh");
        let skin = skin(&file);
        // vue de face : au milieu du cube on voit la face z = 1 (la plus proche)
        let i = pick_skin(&file, &skin, &front, [50.0, 250.0]).expect("une face");
        assert!(skin[i].corners.iter().all(|&n| file.nodes[n].z > 0.99));
        let region = angle_region(&file, &skin, i, 30.0);
        let plane = skin
            .iter()
            .filter(|f| f.corners.iter().all(|&n| file.nodes[n].z > 0.99))
            .count();
        assert_eq!(region.faces.len(), plane);
        assert!(region.faces.len() < skin.len());
    }

    /// Les coordonnées du cube unité qui valent 0 ou 1 pour tous les Nodes de l'ensemble : une arête
    /// du cube est une droite où deux de ces coordonnées sont fixées.
    fn fixed_coordinates(file: &MshFile, nodes: &[usize]) -> usize {
        (0..3)
            .filter(|&c| {
                let v = |n: usize| [file.nodes[n].x, file.nodes[n].y, file.nodes[n].z][c];
                let first = v(nodes[0]);
                (first.abs() < 1e-9 || (first - 1.0).abs() < 1e-9)
                    && nodes.iter().all(|&n| (v(n) - first).abs() < 1e-9)
            })
            .count()
    }

    #[test]
    fn in_3d_the_sharp_edges_of_a_cube_follow_its_creases() {
        let file = example("cube.msh");
        let skin = skin(&file);
        let sharp = sharp_edges(&file, &skin, 30.0);
        assert!(!sharp.is_empty());
        // aucune diagonale de face : tous les segments suivent un pli du cube
        for e in &sharp {
            assert_eq!(fixed_coordinates(&file, &e.nodes), 2, "{e:?}");
        }
        // à 100° le pli à 90° n'est plus vif ; en 2D il n'y a pas d'arête
        assert!(sharp_edges(&file, &skin, 100.0).is_empty());
        let flat = example("plaque.msh");
        assert!(sharp_edges(&flat, &super::skin(&flat), 30.0).is_empty());
    }

    #[test]
    fn a_click_picks_the_visible_sharp_edge_and_the_chain_stops_at_the_corners() {
        let file = example("cube.msh");
        let skin = skin(&file);
        let sharp = sharp_edges(&file, &skin, 30.0);
        // vue de face : l'arête y = 0, z = 1 (devant) cache celle de y = 0, z = 0 au même pixel
        let i = pick_edge(&file, &skin, &sharp, &front, [50.0, 300.0]).expect("une arête");
        assert!(sharp[i].nodes.iter().all(|&n| {
            let p = file.nodes[n];
            p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9
        }));
        // loin de toute arête vive : rien
        assert_eq!(pick_edge(&file, &skin, &sharp, &front, [50.0, 250.0]), None);

        let chain = edge_chain(&file, &sharp, i);
        // tout le pli, rien que lui : autant d'arêtes que de segments du maillage sur cette droite
        let on_line = |n: usize| {
            let p = file.nodes[n];
            p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9
        };
        let expected = sharp
            .iter()
            .filter(|e| e.nodes.iter().all(|&n| on_line(n)))
            .count();
        assert!(expected >= 1);
        assert_eq!(chain.edges.len(), expected);
        assert!(chain.nodes.iter().all(|&n| on_line(n)));
        assert!(chain.faces.is_empty() && chain.elements.is_empty());
        // chaque arête (élément, numéro) retombe sur les Nodes du Set
        for &(e, k) in &chain.edges {
            let positions = edge_positions(file.kind, k).unwrap();
            assert!(positions
                .iter()
                .all(|&p| chain.nodes.contains(&file.cells[e][p])));
        }
    }

    #[test]
    fn an_edge_set_is_written_read_and_highlighted() {
        let file = example("cube.msh");
        let skin = skin(&file);
        let sharp = sharp_edges(&file, &skin, 30.0);
        let i = pick_edge(&file, &skin, &sharp, &front, [50.0, 300.0]).unwrap();
        let chain = edge_chain(&file, &sharp, i);
        let mut with = file.clone();
        with.sets.insert("pli".into(), chain.clone());
        let again = read(&maille_io::msh::write(&with).unwrap()).unwrap();
        assert_eq!(again.sets["pli"], chain);
        // le ruban : des faces de la peau qui touchent l'arête, et pas la peau entière
        let ribbon = highlight_set(&file, &skin, &chain);
        assert!(!ribbon.is_empty() && ribbon.len() < skin.len());
    }

    #[test]
    fn the_highlight_covers_what_is_selected() {
        let file = example("plaque.msh");
        let skin = skin(&file);
        assert_eq!(highlight_elements(&file, &skin, &[0, 1]).len(), 2);
        let nodes = &file.sets["gauche"].nodes;
        let touched = highlight_nodes(&file, &skin, nodes);
        assert!(!touched.is_empty() && touched.len() < file.cells.len());
        let cube = example("cube.msh");
        let cskin = super::skin(&cube);
        let faces = highlight_elements(&cube, &cskin, &(0..cube.cells.len()).collect::<Vec<_>>());
        assert_eq!(faces.len(), cskin.len());
    }

    /// Le parcours sans géométrie : un `.msh` dont on a retiré les Sets, des Sets tirés au cadre et par
    /// angle, écrits dans le `.msh` puis relus, donnent la solution exacte de `examples/traction.toml`.
    #[test]
    fn sets_picked_on_a_bare_mesh_give_the_exact_traction_solution() {
        use maille_io::model::{self, DofDef, ValueDef};

        use crate::mesh_space::Built;
        use crate::model_space::Draft;

        let mut bare = example("plaque.msh");
        bare.sets
            .retain(|name, _| name == maille_io::msh::DOMAIN_SET);
        let mut built = Built::new(bare, None);
        let skin = skin(&built.file);

        let left = Rect::from_corners(at(-0.01, 2.1), at(0.01, -0.1));
        let right = Rect::from_corners(at(3.99, 2.1), at(4.01, -0.1));
        let nodes = |r| node_set(&nodes_in_rect(&built.file, &front, r));
        let (gauche, droite) = (nodes(left), nodes(right));
        let bottom = pick_skin(&built.file, &skin, &front, at(2.0, 0.04)).unwrap();
        let bas = angle_region(&built.file, &skin, bottom, 30.0);
        built.add_mesh_set("gauche", gauche).unwrap();
        built.add_mesh_set("droite", droite).unwrap();
        built.add_mesh_set("bas", bas).unwrap();

        // écrit puis relu, comme le fait `maille solve`
        let file = read(&maille_io::msh::write(&built.file).unwrap()).unwrap();
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.thickness = 2.0;
        let bc = |set: &str, dof, value| ValueDef {
            set: set.into(),
            dof,
            value,
        };
        d.bcs = vec![
            bc("gauche", DofDef::Ux, 0.0),
            bc("droite", DofDef::Ux, 0.004),
            bc("bas", DofDef::Uy, 0.0),
        ];
        let text = d.to_toml("plaque_sets.msh", &file, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &file, None).unwrap();
        let field = loaded.model.solve().unwrap();
        for (i, p) in file.nodes.iter().enumerate() {
            let (ux, uy) = (field.displacements[2 * i], field.displacements[2 * i + 1]);
            assert!((ux - 0.001 * p.x).abs() < 1e-9, "ux du Node {i}");
            assert!((uy + 0.3 * 0.001 * p.y).abs() < 1e-9, "uy du Node {i}");
        }
    }

    /// Le parcours en 3D : des arêtes tirées sur un cube, écrites dans le `.msh` puis relues, portent une
    /// force linéique dont la résultante est q × la longueur des arêtes (4 pour le contour de x = 1).
    #[test]
    fn edges_picked_on_a_cube_carry_a_line_load_with_the_exact_resultant() {
        use maille_io::model::{self, DofDef, EdgeLoadDef, ValueDef};

        use crate::mesh_space::Built;
        use crate::model_space::{diagnose, Draft};

        let mut built = Built::new(example("cube.msh"), None);
        let skin = skin(&built.file);
        let sharp = sharp_edges(&built.file, &skin, 30.0);
        // les quatre plis de la face x = 1, chacun pris en entier par une chaîne
        let on_end = |e: &SharpEdge| e.nodes.iter().all(|&n| built.file.nodes[n].x > 0.99);
        let mut rim = Set::default();
        for i in (0..sharp.len()).filter(|&i| on_end(&sharp[i])) {
            let chain = edge_chain(&built.file, &sharp, i);
            rim.edges.extend(chain.edges);
            rim.nodes.extend(chain.nodes);
        }
        rim.edges.sort_unstable();
        rim.edges.dedup();
        rim.nodes.sort_unstable();
        rim.nodes.dedup();
        assert!(rim.edges.len() >= 4);
        built.add_mesh_set("contour", rim).unwrap();
        assert_eq!(built.edge_sets(), ["contour"]);

        let file = read(&maille_io::msh::write(&built.file).unwrap()).unwrap();
        let q = 25.0;
        let mut d = Draft::default();
        d.material = Some(d.steel());
        let bc = |set: &str, dof| ValueDef {
            set: set.into(),
            dof,
            value: 0.0,
        };
        d.bcs = vec![
            bc("x0", DofDef::Ux),
            bc("y0", DofDef::Uy),
            bc("z0", DofDef::Uz),
        ];
        // sur un Set de Nodes, le diagnostic le dit ; sur le contour, il se tait
        d.edge_loads = vec![EdgeLoadDef {
            set: "x1".into(),
            value: vec![q, 0.0, 0.0],
        }];
        assert!(diagnose(&d, &file, &[])
            .iter()
            .any(|x| x.message.contains("arêtes")));
        d.edge_loads[0].set = "contour".into();
        assert!(!diagnose(&d, &file, &[])
            .iter()
            .any(|x| x.message.contains("arêtes")));

        let text = d.to_toml("cube_contour.msh", &file, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &file, None).unwrap();
        let field = loaded.model.solve().unwrap();
        let reaction: f64 = file.sets["x0"]
            .nodes
            .iter()
            .map(|&n| field.reactions[3 * n])
            .sum();
        // une arête diagonale de la face x = 1 ajouterait sa longueur : le compte est exact
        assert!((reaction + 4.0 * q).abs() < 1e-9, "réaction {reaction}");
    }

    /// En 2D, le bord d'une géométrie est déjà un Set de faces : la force linéique s'y pose, et q = σ · t
    /// (l'épaisseur n'est pas multipliée) donne la traction exacte.
    #[test]
    fn a_boundary_carries_a_line_load_in_2d() {
        use maille_io::model::{self, DofDef, EdgeLoadDef, ValueDef};

        use crate::mesh_space::Built;
        use crate::model_space::{diagnose, Draft};

        let file = example("plaque.msh");
        let built = Built::new(file.clone(), None);
        assert!(built.edge_sets().contains(&"droite"));
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.thickness = 2.0;
        let young = d.steel().young;
        let bc = |set: &str, dof| ValueDef {
            set: set.into(),
            dof,
            value: 0.0,
        };
        d.bcs = vec![bc("gauche", DofDef::Ux), bc("bas", DofDef::Uy)];
        d.edge_loads = vec![EdgeLoadDef {
            set: "droite".into(),
            value: vec![200.0, 0.0],
        }];
        assert!(diagnose(&d, &file, &[])
            .iter()
            .all(|x| !x.message.contains("force linéique")));
        let text = d.to_toml("plaque_force.msh", &file, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &file, None).unwrap();
        let field = loaded.model.solve().unwrap();
        // σ = 200 / 2 = 100 MPa en contraintes planes : ux = σ x / E
        for (i, p) in file.nodes.iter().enumerate() {
            let ux = field.displacements[2 * i];
            assert!((ux - 100.0 * p.x / young).abs() < 1e-9, "ux du Node {i}");
        }
    }
}
