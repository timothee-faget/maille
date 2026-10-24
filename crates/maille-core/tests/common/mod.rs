//! Helpers de maillage pour les tests : ce ne sont pas des mailleurs (le mailleur est le jalon M4).
#![allow(dead_code)]

use std::collections::HashMap;

use maille_core::mesh::Node;

/// Grille w × h de nx × ny cases, chacune coupée en 2 triangles antihoraires.
/// Le Node (i, j) a pour indice i + j (nx + 1).
pub fn rect_grid(nx: usize, ny: usize, w: f64, h: f64) -> (Vec<Node>, Vec<[usize; 3]>) {
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push(Node {
                x: w * i as f64 / nx as f64,
                y: h * j as f64 / ny as f64,
            });
        }
    }
    let mut tris = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
            tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
        }
    }
    (nodes, tris)
}

/// Maillage quadratique obtenu en ajoutant un Node au milieu de chaque arête (partagé entre voisins).
pub struct Quadratic {
    pub nodes: Vec<Node>,
    /// Connectivité Tri6 : sommets puis milieux d'arêtes 0-1, 1-2, 2-0.
    pub elements: Vec<[usize; 6]>,
    mids: HashMap<(usize, usize), usize>,
}

impl Quadratic {
    /// Indice du Node milieu de l'arête (a, b), dans un sens ou dans l'autre.
    pub fn mid(&self, a: usize, b: usize) -> usize {
        self.mids[&(a.min(b), a.max(b))]
    }
}

pub fn to_quadratic(nodes: &[Node], tris: &[[usize; 3]]) -> Quadratic {
    let mut q = Quadratic {
        nodes: nodes.to_vec(),
        elements: Vec::new(),
        mids: HashMap::new(),
    };
    for t in tris {
        let mut mid = |a: usize, b: usize| -> usize {
            let key = (a.min(b), a.max(b));
            if let Some(&m) = q.mids.get(&key) {
                return m;
            }
            let (pa, pb) = (q.nodes[a], q.nodes[b]);
            q.nodes.push(Node {
                x: (pa.x + pb.x) / 2.0,
                y: (pa.y + pb.y) / 2.0,
            });
            q.mids.insert(key, q.nodes.len() - 1);
            q.nodes.len() - 1
        };
        let (m01, m12, m20) = (mid(t[0], t[1]), mid(t[1], t[2]), mid(t[2], t[0]));
        q.elements.push([t[0], t[1], t[2], m01, m12, m20]);
    }
    q
}
