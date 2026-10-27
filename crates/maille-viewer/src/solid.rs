//! Peau d'un Mesh de tétraèdres : les faces portées par un seul tétraèdre, ce qu'on voit de
//! l'extérieur. Tet4 et Tet10 (les nombres de Nodes d'une cellule : 4 ou 10).

use std::collections::{HashMap, HashSet};

use maille_core::shape::{TET10_EDGES, TET_FACES};

/// Triangles à remplir (indices de Nodes, normale sortante) et segments du contour des faces. Les
/// faces d'un Tet10 sont découpées en 4 triangles, mais leur contour reste celui de la face (avec
/// ses milieux) : les arêtes internes de ces sous-triangles ne sont pas dessinées.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Skin {
    pub triangles: Vec<[u32; 3]>,
    pub edges: Vec<[u32; 2]>,
}

/// Peau des cellules `cells`. `None` si une cellule n'a ni 4 ni 10 Nodes.
pub fn skin(cells: &[Vec<usize>]) -> Option<Skin> {
    if cells.iter().any(|c| c.len() != 4 && c.len() != 10) {
        return None;
    }
    // chaque face, par ses trois sommets triés : combien de tétraèdres la portent
    let mut count: HashMap<[usize; 3], usize> = HashMap::new();
    for cell in cells {
        for face in TET_FACES {
            *count.entry(sorted(face.map(|k| cell[k]))).or_default() += 1;
        }
    }
    let mut triangles = Vec::new();
    let mut edges = HashSet::new();
    for cell in cells {
        for face in TET_FACES {
            if count[&sorted(face.map(|k| cell[k]))] != 1 {
                continue;
            }
            // TET_FACES regarde vers l'intérieur d'un tétraèdre positif : on retourne la face
            let local = [face[0], face[2], face[1]];
            let corner = local.map(|k| cell[k] as u32);
            if cell.len() == 4 {
                triangles.push(corner);
                edges.extend((0..3).map(|k| pair(corner[k], corner[(k + 1) % 3])));
                continue;
            }
            let mid = |a: usize, b: usize| {
                let k = TET10_EDGES.iter().position(|&e| e == (a, b) || e == (b, a));
                cell[4 + k.unwrap_or(0)] as u32
            };
            let [p, q, r] = corner;
            let (mpq, mqr, mrp) = (
                mid(local[0], local[1]),
                mid(local[1], local[2]),
                mid(local[2], local[0]),
            );
            triangles.extend([[p, mpq, mrp], [mpq, q, mqr], [mrp, mqr, r], [mpq, mqr, mrp]]);
            edges.extend(
                [(p, mpq), (mpq, q), (q, mqr), (mqr, r), (r, mrp), (mrp, p)]
                    .map(|(a, b)| pair(a, b)),
            );
        }
    }
    let mut edges: Vec<[u32; 2]> = edges.into_iter().collect();
    edges.sort_unstable();
    Some(Skin { triangles, edges })
}

fn sorted(mut v: [usize; 3]) -> [usize; 3] {
    v.sort_unstable();
    v
}

fn pair(a: u32, b: u32) -> [u32; 2] {
    [a.min(b), a.max(b)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testmesh::{cube, cube_quadratic};
    use crate::vec3::{cross, dot, norm, scale, sub, V3};

    fn area_and_orientation(positions: &[V3], skin: &Skin) -> f64 {
        let mut area = 0.0;
        for t in &skin.triangles {
            let [a, b, c] = t.map(|i| positions[i as usize]);
            let n = cross(sub(b, a), sub(c, a));
            // le cube est centré à l'origine : la normale sortante regarde du centre vers le triangle
            let centroid = scale(
                [a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]],
                1.0 / 3.0,
            );
            assert!(
                dot(n, centroid) > 0.0,
                "triangle {t:?} tourné vers l'intérieur"
            );
            area += 0.5 * norm(n);
        }
        area
    }

    #[test]
    fn the_skin_of_a_tet4_cube_is_its_six_faces() {
        let (positions, cells) = cube();
        let skin = skin(&cells).unwrap();
        // aire du cube de côté 2 : 6 · 4
        assert!((area_and_orientation(&positions, &skin) - 24.0).abs() < 1e-12);
        // chaque triangle de peau est sur une face du cube
        for t in &skin.triangles {
            assert!(
                (0..3).any(|k| t.iter().all(|&i| positions[i as usize][k].abs() == 1.0
                    && positions[i as usize][k] == positions[t[0] as usize][k]))
            );
        }
        // chaque arête de peau est portée par deux triangles de la peau (surface fermée)
        let mut count: HashMap<[u32; 2], usize> = HashMap::new();
        for t in &skin.triangles {
            for k in 0..3 {
                *count.entry(pair(t[k], t[(k + 1) % 3])).or_default() += 1;
            }
        }
        assert!(count.values().all(|&c| c == 2));
        assert_eq!(count.len(), skin.edges.len());
    }

    #[test]
    fn the_skin_of_a_tet10_cube_has_four_triangles_per_face_and_the_outline_of_the_face() {
        let (positions, cells) = cube();
        let linear = skin(&cells).unwrap();
        let (positions10, cells10) = cube_quadratic();
        let quadratic = skin(&cells10).unwrap();
        assert!((area_and_orientation(&positions10, &quadratic) - 24.0).abs() < 1e-12);
        assert_eq!(quadratic.triangles.len(), 4 * linear.triangles.len());
        // chaque arête de peau Tet4 devient deux segments, corner-milieu et milieu-corner
        assert_eq!(quadratic.edges.len(), 2 * linear.edges.len());
        // les milieux sont bien au milieu : les segments ont la moitié de la longueur des arêtes
        let longest = |p: &[V3], e: &[[u32; 2]]| {
            e.iter()
                .map(|e| norm(sub(p[e[0] as usize], p[e[1] as usize])))
                .fold(0.0, f64::max)
        };
        assert!(
            (longest(&positions10, &quadratic.edges) - longest(&positions, &linear.edges) / 2.0)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn other_cells_are_refused() {
        assert!(skin(&[vec![0, 1, 2]]).is_none());
        assert!(skin(&[]).unwrap().triangles.is_empty());
    }
}
