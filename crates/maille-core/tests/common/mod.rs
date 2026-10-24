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
            nodes.push(Node::new_2d(
                w * i as f64 / nx as f64,
                h * j as f64 / ny as f64,
            ));
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
            q.nodes
                .push(Node::new_2d((pa.x + pb.x) / 2.0, (pa.y + pb.y) / 2.0));
            q.mids.insert(key, q.nodes.len() - 1);
            q.nodes.len() - 1
        };
        let (m01, m12, m20) = (mid(t[0], t[1]), mid(t[1], t[2]), mid(t[2], t[0]));
        q.elements.push([t[0], t[1], t[2], m01, m12, m20]);
    }
    q
}

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
        for (k, (a, b)) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
            .into_iter()
            .enumerate()
        {
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
