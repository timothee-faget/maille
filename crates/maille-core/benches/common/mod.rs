//! Maillages 3D structurés pour les bancs : copie de `tests/common/mod.rs` (les tests et les bancs
//! ne partagent pas de module), réduite aux tétraèdres.
#![allow(dead_code)]

use std::collections::HashMap;

use maille_core::mesh::Node;
use maille_core::shape::TET10_EDGES;

/// Pavé lx × ly × lz de nx × ny × nz cubes, chacun coupé en 6 tétraèdres (découpage de Kuhn : les
/// 6 chemins de (0,0,0) à (1,1,1) le long des arêtes, conforme entre cubes voisins), orientés comme
/// `Tet4`. Le Node (i, j, k) a pour indice i + (nx + 1) (j + (ny + 1) k).
pub fn box_grid(n: [usize; 3], size: [f64; 3]) -> (Vec<Node>, Vec<[usize; 4]>) {
    let idx = |p: [usize; 3]| p[0] + (n[0] + 1) * (p[1] + (n[1] + 1) * p[2]);
    let mut nodes = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                let at = |c: usize, v: usize| size[c] * v as f64 / n[c] as f64;
                nodes.push(Node::new(at(0, i), at(1, j), at(2, k)));
            }
        }
    }
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut tets = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                for axes in permutations {
                    let mut p = [i, j, k];
                    let mut tet = [idx(p); 4];
                    for (v, &axis) in tet.iter_mut().skip(1).zip(&axes) {
                        p[axis] += 1;
                        *v = idx(p);
                    }
                    // une permutation impaire des axes retourne le tétraèdre
                    if signed_volume(&nodes, &tet) < 0.0 {
                        tet.swap(1, 2);
                    }
                    tets.push(tet);
                }
            }
        }
    }
    (nodes, tets)
}

fn signed_volume(nodes: &[Node], t: &[usize]) -> f64 {
    let e = |k: usize| {
        let (a, b) = (nodes[t[0]], nodes[t[k]]);
        [b.x - a.x, b.y - a.y, b.z - a.z]
    };
    let (a, b, c) = (e(1), e(2), e(3));
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// Maillage Tet10 obtenu en ajoutant un Node au milieu de chaque arête (partagé entre voisins).
pub struct QuadraticTet {
    pub nodes: Vec<Node>,
    /// Connectivité Tet10 : sommets puis milieux d'arêtes 0-1, 1-2, 2-0, 0-3, 1-3, 2-3.
    pub elements: Vec<[usize; 10]>,
    mids: HashMap<(usize, usize), usize>,
}

impl QuadraticTet {
    /// Indice du Node milieu de l'arête (a, b), dans un sens ou dans l'autre.
    pub fn mid(&self, a: usize, b: usize) -> usize {
        self.mids[&(a.min(b), a.max(b))]
    }
}

pub fn to_quadratic_tet(nodes: &[Node], tets: &[[usize; 4]]) -> QuadraticTet {
    let mut q = QuadraticTet {
        nodes: nodes.to_vec(),
        elements: Vec::new(),
        mids: HashMap::new(),
    };
    for t in tets {
        let mut element = [0; 10];
        element[..4].copy_from_slice(t);
        for (k, (a, b)) in TET10_EDGES.into_iter().enumerate() {
            let key = (t[a].min(t[b]), t[a].max(t[b]));
            element[4 + k] = *q.mids.entry(key).or_insert_with(|| {
                let (pa, pb) = (q.nodes[t[a]], q.nodes[t[b]]);
                q.nodes.push(Node::new(
                    (pa.x + pb.x) / 2.0,
                    (pa.y + pb.y) / 2.0,
                    (pa.z + pb.z) / 2.0,
                ));
                q.nodes.len() - 1
            });
        }
        q.elements.push(element);
    }
    q
}
