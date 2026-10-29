//! Qualité des Elements d'un maillage : 1 pour un triangle équilatéral ou un tétraèdre régulier, vers
//! 0 pour un Element écrasé. Un Element écrasé donne un calcul peu précis, donc c'est ce que
//! l'espace Maillage montre. Seuls les sommets comptent (les Nodes milieux d'un Tri6 ou d'un Tet10 ne
//! changent pas la forme).

use maille_core::mesh::Node;
use maille_io::msh::MshFile;

/// En dessous, un Element est signalé comme mauvais.
pub const POOR: f64 = 0.3;

#[derive(Debug, Clone, PartialEq)]
pub struct Quality {
    /// Pire qualité des Elements qui touchent le Node : ce que le viewer colorie.
    pub per_node: Vec<f64>,
    pub elements: usize,
    pub min: f64,
    pub mean: f64,
    /// Nombre d'Elements sous `POOR`.
    pub poor: usize,
}

fn sub(a: &Node, b: &Node) -> [f64; 3] {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 4√3 A / (somme des carrés des côtés) : 1 pour le triangle équilatéral.
fn triangle(p: [&Node; 3]) -> f64 {
    let (a, b) = (sub(p[1], p[0]), sub(p[2], p[0]));
    let area = 0.5 * dot(cross(a, b), cross(a, b)).sqrt();
    let sides = dot(a, a) + dot(b, b) + dot(sub(p[2], p[1]), sub(p[2], p[1]));
    if sides > 0.0 {
        4.0 * 3f64.sqrt() * area / sides
    } else {
        0.0
    }
}

/// 12 (3V)^(2/3) / (somme des carrés des 6 arêtes) : 1 pour le tétraèdre régulier.
fn tetrahedron(p: [&Node; 4]) -> f64 {
    let (a, b, c) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
    let volume = dot(a, cross(b, c)).abs() / 6.0;
    let edges: f64 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        .iter()
        .map(|&(i, j)| {
            let d = sub(p[j], p[i]);
            dot(d, d)
        })
        .sum();
    if edges > 0.0 {
        12.0 * (3.0 * volume).powf(2.0 / 3.0) / edges
    } else {
        0.0
    }
}

/// Qualité de chaque Element et des Nodes du maillage.
pub fn of(file: &MshFile) -> Quality {
    let node = |i: usize| &file.nodes[i];
    let by_element: Vec<f64> = file
        .cells
        .iter()
        .map(|c| {
            if file.kind.dim() == 2 {
                triangle([node(c[0]), node(c[1]), node(c[2])])
            } else {
                tetrahedron([node(c[0]), node(c[1]), node(c[2]), node(c[3])])
            }
        })
        .collect();
    let mut per_node = vec![1.0_f64; file.nodes.len()];
    for (cell, q) in file.cells.iter().zip(&by_element) {
        for &i in cell {
            per_node[i] = per_node[i].min(*q);
        }
    }
    let elements = by_element.len();
    Quality {
        per_node,
        elements,
        min: by_element.iter().copied().fold(1.0, f64::min),
        mean: by_element.iter().sum::<f64>() / elements.max(1) as f64,
        poor: by_element.iter().filter(|&&q| q < POOR).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_io::msh::CellKind;
    use std::collections::BTreeMap;

    fn file(kind: CellKind, nodes: Vec<Node>, cells: Vec<Vec<usize>>) -> MshFile {
        MshFile {
            nodes,
            kind,
            cells,
            sets: BTreeMap::new(),
        }
    }

    fn close(got: f64, want: f64) {
        assert!((got - want).abs() < 1e-12, "{got} contre {want}");
    }

    #[test]
    fn an_equilateral_triangle_is_perfect() {
        let h = 3f64.sqrt() / 2.0;
        let f = file(
            CellKind::Tri3,
            vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.5, h),
            ],
            vec![vec![0, 1, 2]],
        );
        let q = of(&f);
        close(q.min, 1.0);
        close(q.mean, 1.0);
        assert_eq!((q.elements, q.poor), (1, 0));
        assert!(q.per_node.iter().all(|&v| (v - 1.0).abs() < 1e-12));
    }

    #[test]
    fn a_flat_triangle_is_poor() {
        // presque aplati : de l'aire, mais pas de forme
        let f = file(
            CellKind::Tri3,
            vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.5, 0.01),
            ],
            vec![vec![0, 1, 2]],
        );
        let q = of(&f);
        assert!(q.min < 0.05 && q.poor == 1, "{q:?}");
    }

    #[test]
    fn a_right_isoceles_triangle_has_the_known_quality() {
        // côtés 1, 1, √2 : 4√3 (1/2) / 4 = √3 / 2
        let f = file(
            CellKind::Tri3,
            vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.0, 1.0),
            ],
            vec![vec![0, 1, 2]],
        );
        close(of(&f).min, 3f64.sqrt() / 2.0);
    }

    #[test]
    fn a_regular_tetrahedron_is_perfect_and_a_corner_one_is_not() {
        // sommets d'un cube : tétraèdre régulier d'arête √2
        let regular = file(
            CellKind::Tet4,
            vec![
                Node::new(0.0, 0.0, 0.0),
                Node::new(1.0, 1.0, 0.0),
                Node::new(1.0, 0.0, 1.0),
                Node::new(0.0, 1.0, 1.0),
            ],
            vec![vec![0, 1, 2, 3]],
        );
        close(of(&regular).min, 1.0);
        // coin d'un cube : arêtes 1, 1, 1, √2, √2, √2 ; V = 1/6
        let corner = file(
            CellKind::Tet4,
            vec![
                Node::new(0.0, 0.0, 0.0),
                Node::new(1.0, 0.0, 0.0),
                Node::new(0.0, 1.0, 0.0),
                Node::new(0.0, 0.0, 1.0),
            ],
            vec![vec![0, 1, 2, 3]],
        );
        close(of(&corner).min, 12.0 * 0.5f64.powf(2.0 / 3.0) / 9.0);
    }

    #[test]
    fn a_node_takes_the_worst_quality_of_the_elements_around_it() {
        let h = 3f64.sqrt() / 2.0;
        let f = file(
            CellKind::Tri3,
            vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.5, h),
                Node::new_2d(1.5, 0.01),
            ],
            vec![vec![0, 1, 2], vec![1, 3, 2]],
        );
        let q = of(&f);
        close(q.per_node[0], 1.0);
        assert!(q.per_node[1] < 1.0);
        close(q.per_node[1], q.per_node[3]);
        assert!(q.mean > q.min);
    }
}
