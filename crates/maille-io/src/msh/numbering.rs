//! Numérotation des Nodes : Gmsh et Maille, orientation, faces.

use maille_core::mesh::Node;
use maille_core::shape::{TET10_EDGES, TET_FACES, TRI6_EDGES};

use super::CellKind;

/// Passe un Tet10 de la numérotation de Gmsh à celle de Maille, ou l'inverse (la permutation est
/// son propre inverse). Gmsh place les milieux 7 = (3,0), 8 = (3,2), 9 = (3,1) ; Maille 7 = (0,3),
/// 8 = (1,3), 9 = (2,3). Le Tri6 est identique (3 = 0-1, 4 = 1-2, 5 = 2-0).
pub fn swap_tet10_mids(nodes: &mut [usize]) {
    if nodes.len() == 10 {
        nodes.swap(8, 9);
    }
}

/// Permute deux sommets (et les milieux en conséquence) si l'Element est retourné : triangle
/// horaire ou tétraèdre de volume négatif. `nodes` est en numérotation de Maille. Un Element plat
/// n'est pas touché : `Mesh::new` le refuse.
pub fn orient(kind: CellKind, nodes: &mut [usize], points: &[Node]) {
    let p = |i: usize| points[nodes[i]];
    let reversed = match kind {
        CellKind::Tri3 | CellKind::Tri6 => {
            (p(1).x - p(0).x) * (p(2).y - p(0).y) - (p(2).x - p(0).x) * (p(1).y - p(0).y) < 0.0
        }
        CellKind::Tet4 | CellKind::Tet10 => {
            let e = |k: usize| [p(k).x - p(0).x, p(k).y - p(0).y, p(k).z - p(0).z];
            let (a, b, c) = (e(1), e(2), e(3));
            a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
                < 0.0
        }
    };
    if !reversed {
        return;
    }
    nodes.swap(0, 1);
    match kind {
        CellKind::Tri6 => nodes.swap(4, 5),
        CellKind::Tet10 => {
            nodes.swap(5, 6);
            nodes.swap(7, 8);
        }
        CellKind::Tri3 | CellKind::Tet4 => {}
    }
}

/// Nombre de faces d'un Element : 3 arêtes pour un triangle, 4 faces pour un tétraèdre.
pub fn face_count(kind: CellKind) -> usize {
    match kind {
        CellKind::Tri3 | CellKind::Tri6 => 3,
        CellKind::Tet4 | CellKind::Tet10 => 4,
    }
}

/// Positions, dans les Nodes de l'Element, de ceux de la face `face` : sommets (2 en 2D, 3 en 3D)
/// puis milieux. Même numérotation que `Element::face_nodes`, et même ordre que les éléments de
/// bord de Gmsh (ligne à 3 Nodes, triangle à 6 Nodes).
pub fn face_positions(kind: CellKind, face: usize) -> Option<Vec<usize>> {
    let tet_mid = |a: usize, b: usize| {
        let k = TET10_EDGES.iter().position(|&e| e == (a, b) || e == (b, a));
        4 + k.unwrap_or(0)
    };
    match kind {
        CellKind::Tri3 => TRI6_EDGES.get(face).map(|&(a, b)| vec![a, b]),
        CellKind::Tri6 => TRI6_EDGES.get(face).map(|&(a, b)| vec![a, b, 3 + face]),
        CellKind::Tet4 => TET_FACES.get(face).map(|f| f.to_vec()),
        CellKind::Tet10 => TET_FACES
            .get(face)
            .map(|&[a, b, c]| vec![a, b, c, tet_mid(a, b), tet_mid(b, c), tet_mid(c, a)]),
    }
}

/// Clé d'une face indépendante de l'ordre de ses sommets (2 en 2D, 3 en 3D, complétée par MAX).
pub fn face_key(vertices: &[usize]) -> [usize; 3] {
    let mut key = [usize::MAX; 3];
    key[..vertices.len()].copy_from_slice(vertices);
    key[..vertices.len()].sort_unstable();
    key
}

#[cfg(test)]
mod tests {
    use maille_core::material::Material;
    use maille_core::mesh::{Element, PlaneState, Tet10, Tet4, Tri3, Tri6};

    use super::*;

    fn element(kind: CellKind) -> Element {
        let m = Material::with_poisson(1.0, 0.3).expect("matériau valide");
        let ids: Vec<usize> = (0..kind.node_count()).collect();
        match kind {
            CellKind::Tri3 => Element::Tri3(
                Tri3::new([0, 1, 2], m, 1.0, PlaneState::Stress).expect("Tri3 valide"),
            ),
            CellKind::Tri6 => Element::Tri6(
                Tri6::new(ids.try_into().unwrap(), m, 1.0, PlaneState::Stress)
                    .expect("Tri6 valide"),
            ),
            CellKind::Tet4 => Element::Tet4(Tet4::new([0, 1, 2, 3], m)),
            CellKind::Tet10 => Element::Tet10(Tet10::new(ids.try_into().unwrap(), m)),
        }
    }

    #[test]
    fn faces_as_in_core() {
        for kind in [
            CellKind::Tri3,
            CellKind::Tri6,
            CellKind::Tet4,
            CellKind::Tet10,
        ] {
            let e = element(kind);
            for f in 0..face_count(kind) {
                assert_eq!(
                    face_positions(kind, f),
                    e.face_nodes(f),
                    "{kind:?} face {f}"
                );
            }
            assert_eq!(face_positions(kind, face_count(kind)), None);
        }
    }

    #[test]
    fn reversed_triangle_is_flipped() {
        let pts = [
            Node::new_2d(0.0, 0.0),
            Node::new_2d(1.0, 0.0),
            Node::new_2d(0.0, 1.0),
        ];
        let mut ccw = vec![0, 1, 2];
        orient(CellKind::Tri3, &mut ccw, &pts);
        assert_eq!(ccw, [0, 1, 2]);
        let mut cw = vec![0, 2, 1];
        orient(CellKind::Tri3, &mut cw, &pts);
        assert_eq!(cw, [2, 0, 1]);
    }
}
