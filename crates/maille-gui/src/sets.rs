//! Désigner des entités de la géométrie à l'écran, pour en faire des Sets. Sans fenêtre : on projette
//! les faces de bord du maillage avec la caméra du viewer et on cherche celle qui est sous le curseur.
//! En 3D, c'est la face la plus proche de l'observateur ; en 2D, le bord le plus proche du curseur.

use maille_io::msh::face_positions;

use crate::mesh_space::Built;

/// En 2D, à moins de tant de pixels d'un bord, le curseur le désigne.
pub const EDGE_TOLERANCE_PX: f64 = 10.0;

/// Une face de bord désignable.
#[derive(Debug, Clone, PartialEq)]
pub struct PickFace {
    /// Rang de l'entité qui la porte, dans `Built::entities`.
    pub entity: usize,
    /// Element dont c'est une face.
    pub element: usize,
    /// Ses sommets : 2 en 2D (un segment), 3 en 3D (un triangle).
    pub corners: Vec<usize>,
}

/// Toutes les faces de bord que portent les entités du maillage.
pub fn pickable(built: &Built) -> Vec<PickFace> {
    let kind = built.file.kind;
    let corners = kind.dim();
    let mut faces = Vec::new();
    for (entity, name) in built.entities.iter().enumerate() {
        let Some(set) = built.file.sets.get(name) else {
            continue;
        };
        for &(element, face) in &set.faces {
            let Some(positions) = face_positions(kind, face) else {
                continue;
            };
            let cell = &built.file.cells[element];
            faces.push(PickFace {
                entity,
                element,
                corners: positions[..corners].iter().map(|&p| cell[p]).collect(),
            });
        }
    }
    faces
}

type Point = [f64; 3];

fn position(built: &Built, node: usize) -> Point {
    let p = &built.file.nodes[node];
    [p.x, p.y, p.z]
}

/// Distance du point `p` au segment `a`–`b`, dans le plan de l'écran.
pub(crate) fn segment_distance(p: [f64; 2], a: Point, b: Point) -> f64 {
    let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
    let len2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if len2 > 0.0 {
        ((ap[0] * ab[0] + ap[1] * ab[1]) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (ap[0] - t * ab[0]).hypot(ap[1] - t * ab[1])
}

/// Profondeur du triangle `t` (pixels, profondeur) au pixel `p`, s'il le couvre.
pub(crate) fn triangle_depth(p: [f64; 2], t: [Point; 3]) -> Option<f64> {
    let d = (t[1][1] - t[2][1]) * (t[0][0] - t[2][0]) + (t[2][0] - t[1][0]) * (t[0][1] - t[2][1]);
    if d.abs() < 1e-12 {
        return None;
    }
    let a = ((t[1][1] - t[2][1]) * (p[0] - t[2][0]) + (t[2][0] - t[1][0]) * (p[1] - t[2][1])) / d;
    let b = ((t[2][1] - t[0][1]) * (p[0] - t[2][0]) + (t[0][0] - t[2][0]) * (p[1] - t[2][1])) / d;
    let c = 1.0 - a - b;
    let eps = -1e-9;
    (a >= eps && b >= eps && c >= eps).then(|| a * t[0][2] + b * t[1][2] + c * t[2][2])
}

/// L'entité sous le pixel `cursor` (origine en haut à gauche de la texture). `project` donne, pour un
/// point du maillage, ses pixels (x, y) et sa profondeur (0 : devant), comme `ResultsView::project`.
pub fn pick<'a>(
    built: &'a Built,
    faces: &[PickFace],
    project: &dyn Fn(Point) -> Point,
    cursor: [f64; 2],
) -> Option<&'a str> {
    let projected = |f: &PickFace, k: usize| project(position(built, f.corners[k]));
    let best = if built.file.kind.dim() == 2 {
        faces
            .iter()
            .map(|f| {
                (
                    segment_distance(cursor, projected(f, 0), projected(f, 1)),
                    f,
                )
            })
            .filter(|(d, _)| *d <= EDGE_TOLERANCE_PX)
            .min_by(|a, b| a.0.total_cmp(&b.0))
    } else {
        faces
            .iter()
            .filter_map(|f| {
                let t = [projected(f, 0), projected(f, 1), projected(f, 2)];
                triangle_depth(cursor, t).map(|depth| (depth, f))
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
    };
    best.map(|(_, f)| built.entities[f.entity].as_str())
}

/// Les triangles à surligner pour ces entités : en 3D leurs faces, en 2D les triangles qui touchent
/// leurs bords (un bord seul est trop fin pour se voir).
pub fn highlight(built: &Built, entities: &[String]) -> Vec<[Point; 3]> {
    let mut triangles = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for name in entities {
        let Some(set) = built.file.sets.get(name) else {
            continue;
        };
        for &(element, face) in &set.faces {
            let cell = &built.file.cells[element];
            if built.file.kind.dim() == 2 {
                if seen.insert(element) {
                    triangles.push([0, 1, 2].map(|k| position(built, cell[k])));
                }
            } else if let Some(positions) = face_positions(built.file.kind, face) {
                triangles.push([0, 1, 2].map(|k| position(built, cell[positions[k]])));
            }
        }
    }
    triangles
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::path::{Path, PathBuf};

    use super::*;
    use crate::mesh_space::{import, mesh, Imported};

    fn built(name: &str) -> Built {
        let path: PathBuf = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(name);
        let Imported::Cad(cad) = import(&path).unwrap() else {
            panic!("{name}");
        };
        mesh(&cad, cad.suggested, false, &BTreeMap::new()).unwrap()
    }

    /// Vue orthographique de face : 100 pixels par unité, y vers le haut, profondeur −z.
    fn front(p: Point) -> Point {
        [p[0] * 100.0, 300.0 - p[1] * 100.0, -p[2]]
    }

    #[test]
    fn in_2d_the_nearest_edge_is_picked_within_a_few_pixels() {
        let b = built("plaque.dxf");
        let faces = pickable(&b);
        assert!(!faces.is_empty() && faces.iter().all(|f| f.corners.len() == 2));
        // plaque 4 × 2 : bord gauche en x = 0, droit en x = 4, bas en y = 0, haut en y = 2
        let at = |x: f64, y: f64| [x * 100.0, 300.0 - y * 100.0];
        let pick_at = |c| pick(&b, &faces, &front, c);
        assert_eq!(pick_at(at(0.05, 1.0)), Some("gauche"));
        assert_eq!(pick_at(at(3.95, 1.0)), Some("droite"));
        assert_eq!(pick_at(at(2.0, 0.04)), Some("bas"));
        assert_eq!(pick_at(at(2.0, 1.96)), Some("haut"));
        // au milieu de la plaque : aucun bord à portée
        assert_eq!(pick_at(at(2.0, 1.0)), None);
    }

    #[test]
    fn in_3d_the_face_nearest_the_viewer_is_picked() {
        let b = built("cube.stl");
        let faces = pickable(&b);
        assert!(!faces.is_empty() && faces.iter().all(|f| f.corners.len() == 3));
        // cube unité vu de face : au milieu, la face z = 1 cache celle de derrière
        let name = pick(&b, &faces, &front, [50.0, 250.0]).expect("une face");
        let set = &b.file.sets[name];
        assert!(!set.faces.is_empty());
        for f in faces.iter().filter(|f| b.entities[f.entity] == name) {
            for &n in &f.corners {
                assert!(
                    (b.file.nodes[n].z - 1.0).abs() < 1e-9,
                    "{name} : z = {}",
                    b.file.nodes[n].z
                );
            }
        }
        // vu de derrière, c'est l'autre face
        let back = |p: Point| [p[0] * 100.0, 300.0 - p[1] * 100.0, p[2]];
        let other = pick(&b, &faces, &back, [50.0, 250.0]).expect("une face");
        assert_ne!(other, name);
        // hors du cube : rien
        assert_eq!(pick(&b, &faces, &front, [500.0, 250.0]), None);
    }

    #[test]
    fn the_highlight_covers_the_faces_of_the_entities() {
        let b3 = built("cube.stl");
        let name = b3.entities[0].clone();
        let t = highlight(&b3, std::slice::from_ref(&name));
        assert_eq!(t.len(), b3.file.sets[&name].faces.len());
        // en 2D : les triangles voisins du bord, chacun une seule fois
        let b2 = built("plaque.dxf");
        let t = highlight(&b2, &["gauche".to_string()]);
        let mut elements: Vec<usize> = b2.file.sets["gauche"].faces.iter().map(|f| f.0).collect();
        elements.sort_unstable();
        elements.dedup();
        assert_eq!(t.len(), elements.len());
        assert!(highlight(&b2, &["fantome".to_string()]).is_empty());
    }
}
