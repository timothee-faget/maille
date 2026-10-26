//! Des facettes brutes à la `Surface` : soudure des sommets, contrôles, orientation, patches.

use std::collections::HashMap;

use super::reader::Facets;
use super::{cross, dot, error, norm, sub, Patch, Surface};
use crate::IoError;

/// Sommets confondus s'ils sont plus proches que cette fraction de la diagonale de la pièce.
const WELD: f64 = 1e-6;
/// Facette dégénérée : aire sous cette fraction du carré de la diagonale.
const DEGENERATE: f64 = 1e-12;

pub(super) fn build(facets: Facets, feature_angle: f64) -> Result<Surface, IoError> {
    let (vertices, triangles, diagonal) = weld(&facets.corners);
    for (i, t) in triangles.iter().enumerate() {
        let [a, b, c] = t.map(|v| vertices[v]);
        let area = 0.5 * norm(cross(sub(b, a), sub(c, a)));
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] || area <= DEGENERATE * diagonal * diagonal
        {
            return Err(error(format!(
                "facette {} dégénérée (sommets confondus ou alignés) : la supprimer ou réparer le maillage dans l'outil de CAO",
                i + 1
            )));
        }
    }
    let edges = check_edges(&triangles, &vertices)?;
    check_one_piece(&triangles, vertices.len())?;

    let mut surface = Surface {
        vertices,
        triangles,
        patches: Vec::new(),
    };
    let volume = surface.volume();
    if volume.abs() <= DEGENERATE * diagonal.powi(3) {
        return Err(error(
            "la surface n'enferme aucun volume (volume nul) : vérifier le fichier".to_string(),
        ));
    }
    if volume < 0.0 {
        for t in &mut surface.triangles {
            t.swap(1, 2);
        }
    }
    surface.patches = if facets.solids.len() >= 2 {
        patches_by_solid(&facets)
    } else {
        patches_by_angle(&surface, &edges, feature_angle)
    };
    Ok(surface)
}

/// Sommets soudés (grille de côté la tolérance, voisines comprises), triangles, diagonale.
fn weld(corners: &[[[f64; 3]; 3]]) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, f64) {
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in corners.iter().flatten() {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let diagonal = norm(sub(hi, lo));
    let tol = (WELD * diagonal).max(f64::MIN_POSITIVE);
    let cell = |p: &[f64; 3]| p.map(|x| (x / tol).floor() as i64);
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut grid: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let mut triangles = Vec::with_capacity(corners.len());
    for tri in corners {
        let mut t = [0; 3];
        for (slot, p) in t.iter_mut().zip(tri) {
            let c = cell(p);
            let found = (-1..=1)
                .flat_map(|i| (-1..=1).flat_map(move |j| (-1..=1).map(move |k| (i, j, k))))
                .filter_map(|(i, j, k)| grid.get(&[c[0] + i, c[1] + j, c[2] + k]))
                .flatten()
                .copied()
                .find(|&v| norm(sub(vertices[v], *p)) <= tol);
            *slot = found.unwrap_or_else(|| {
                vertices.push(*p);
                grid.entry(c).or_default().push(vertices.len() - 1);
                vertices.len() - 1
            });
        }
        triangles.push(t);
    }
    (vertices, triangles, diagonal)
}

/// Arête (petit sommet, grand sommet) vers les triangles qui la bordent, avec `true` si le
/// triangle la parcourt du petit au grand sommet.
type Edges = HashMap<(usize, usize), Vec<(usize, bool)>>;

fn check_edges(triangles: &[[usize; 3]], vertices: &[[f64; 3]]) -> Result<Edges, IoError> {
    let mut edges = Edges::new();
    for (i, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            edges
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((i, a < b));
        }
    }
    let at = |e: &(usize, usize)| {
        let p = vertices[e.0];
        format!("({}, {}, {})", p[0], p[1], p[2])
    };
    // les arêtes sont parcourues dans l'ordre des sommets pour un message qui ne change pas d'un lancement à l'autre
    let mut keys: Vec<_> = edges.keys().copied().collect();
    keys.sort_unstable();
    if let Some(e) = keys.iter().find(|e| edges[e].len() > 2) {
        return Err(error(format!(
            "arête partagée par {} triangles près de {} : la surface n'est pas une variété (parois qui se touchent ou se recoupent)",
            edges[e].len(),
            at(e)
        )));
    }
    let open: Vec<_> = keys.iter().filter(|e| edges[*e].len() == 1).collect();
    if let Some(e) = open.first() {
        return Err(error(format!(
            "surface ouverte : {} arête(s) n'ont qu'un triangle, dont une près de {} ; boucher les trous dans l'outil de CAO",
            open.len(),
            at(e)
        )));
    }
    if let Some(e) = keys.iter().find(|e| edges[*e][0].1 == edges[*e][1].1) {
        return Err(error(format!(
            "orientation incohérente des triangles près de {} : certaines normales sortent, d'autres rentrent ; recalculer les normales dans l'outil de CAO",
            at(e)
        )));
    }
    Ok(edges)
}

fn check_one_piece(triangles: &[[usize; 3]], count: usize) -> Result<(), IoError> {
    let mut parent: Vec<usize> = (0..count).collect();
    fn root(parent: &mut [usize], mut v: usize) -> usize {
        while parent[v] != v {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        v
    }
    for t in triangles {
        for k in 1..3 {
            let (a, b) = (root(&mut parent, t[0]), root(&mut parent, t[k]));
            parent[a] = b;
        }
    }
    let pieces = (0..count).filter(|&v| root(&mut parent, v) == v).count();
    if pieces > 1 {
        return Err(error(format!(
            "la surface compte {pieces} pièces séparées (ou une cavité) : Maille maille une seule pièce pleine, séparer le fichier"
        )));
    }
    Ok(())
}

fn patches_by_solid(facets: &Facets) -> Vec<Patch> {
    let mut patches: Vec<Patch> = Vec::new();
    for (i, &s) in facets.solid.iter().enumerate() {
        let name = match facets.solids[s].as_str() {
            "" => format!("surface_{}", s + 1),
            name => name.to_string(),
        };
        match patches.iter_mut().find(|p| p.name == name) {
            Some(p) => p.triangles.push(i),
            None => patches.push(Patch {
                name,
                triangles: vec![i],
            }),
        }
    }
    patches
}

fn patches_by_angle(surface: &Surface, edges: &Edges, feature_angle: f64) -> Vec<Patch> {
    let normals: Vec<[f64; 3]> = surface
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| surface.vertices[v]);
            let n = cross(sub(b, a), sub(c, a));
            n.map(|x| x / norm(n))
        })
        .collect();
    let smooth = feature_angle.to_radians().cos();
    let mut neighbours = vec![Vec::new(); surface.triangles.len()];
    for pair in edges.values() {
        let (a, b) = (pair[0].0, pair[1].0);
        if dot(normals[a], normals[b]) >= smooth {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
    }
    let mut patch_of = vec![usize::MAX; surface.triangles.len()];
    let mut patches: Vec<Patch> = Vec::new();
    for start in 0..surface.triangles.len() {
        if patch_of[start] != usize::MAX {
            continue;
        }
        let id = patches.len();
        let mut members = Vec::new();
        let mut stack = vec![start];
        patch_of[start] = id;
        while let Some(t) = stack.pop() {
            members.push(t);
            for &n in &neighbours[t] {
                if patch_of[n] == usize::MAX {
                    patch_of[n] = id;
                    stack.push(n);
                }
            }
        }
        members.sort_unstable();
        patches.push(Patch {
            name: format!("surface_{}", id + 1),
            triangles: members,
        });
    }
    patches
}
