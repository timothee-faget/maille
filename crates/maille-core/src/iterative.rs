//! Solveur itératif : gradients conjugués préconditionnés par un Cholesky incomplet IC(0), sur une
//! matrice symétrique creuse en CSR. La mémoire croît en O(n) (le Cholesky en ligne de ciel la
//! fait croître en n^(5/3) en 3D), ce qui rend les gros modèles 3D possibles.

use crate::error::MailleError;
use crate::solve::{dot, Adjacency};

/// Résidu relatif ‖b − A x‖ / ‖b‖ visé.
const TOLERANCE: f64 = 1e-12;
/// Nombre maximal d'itérations avant `IterativeNotConverged`.
const MAX_ITERATIONS: usize = 10_000;
/// Décalages relatifs de la diagonale essayés quand IC(0) rencontre un pivot négatif ou nul
/// (une matrice d'élasticité n'est pas une M-matrice : IC(0) peut échouer sans décalage).
const SHIFTS: [f64; 6] = [0.0, 1e-3, 1e-2, 1e-1, 1.0, 10.0];

/// Graphe en CSR : les voisins du sommet v sont `targets[offsets[v]..offsets[v + 1]]`, dans
/// aucun ordre particulier.
#[derive(Debug, Clone, PartialEq)]
pub struct Graph {
    offsets: Vec<usize>,
    targets: Vec<usize>,
}

impl Adjacency for Graph {
    fn len(&self) -> usize {
        self.offsets.len() - 1
    }

    fn neighbours(&self, v: usize) -> &[usize] {
        &self.targets[self.offsets[v]..self.offsets[v + 1]]
    }
}

impl Graph {
    /// Graphe sur `n` sommets où deux sommets sont voisins s'ils sont dans un même élément de
    /// `elements` (liste de sommets distincts, tous < n). Pas de liste de voisins avec doublons :
    /// chaque rangée est dédupliquée au fil de l'eau, la mémoire est celle du résultat.
    pub fn from_elements(n: usize, elements: &[Vec<usize>]) -> Self {
        // sommet -> éléments qui le contiennent
        let mut incident = vec![0; n + 1];
        for element in elements {
            for &v in element {
                incident[v + 1] += 1;
            }
        }
        for v in 0..n {
            incident[v + 1] += incident[v];
        }
        let mut cursor = incident.clone();
        let mut of_vertex = vec![0; incident[n]];
        for (e, element) in elements.iter().enumerate() {
            for &v in element {
                of_vertex[cursor[v]] = e;
                cursor[v] += 1;
            }
        }

        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        let mut targets = Vec::new();
        // mark[w] == v : w est déjà voisin de v (v lui-même compris)
        let mut mark = vec![usize::MAX; n];
        for v in 0..n {
            mark[v] = v;
            for &e in &of_vertex[incident[v]..incident[v + 1]] {
                for &w in &elements[e] {
                    if mark[w] != v {
                        mark[w] = v;
                        targets.push(w);
                    }
                }
            }
            offsets.push(targets.len());
        }
        Self { offsets, targets }
    }
}

/// Triangle inférieur d'une matrice symétrique en CSR : les colonnes de chaque rangée sont
/// croissantes, la diagonale en dernier.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    offsets: Vec<usize>,
    cols: Vec<usize>,
    values: Vec<f64>,
}

impl Csr {
    /// Matrice nulle de profil `graph` dans la numérotation `rank` (`rank[ancien] = nouveau`,
    /// `order` est son inverse).
    pub fn zeros(graph: &Graph, order: &[usize], rank: &[usize]) -> Self {
        let mut offsets = Vec::with_capacity(order.len() + 1);
        offsets.push(0);
        let mut cols = Vec::new();
        for (r, &v) in order.iter().enumerate() {
            let start = cols.len();
            cols.extend(
                graph
                    .neighbours(v)
                    .iter()
                    .map(|&w| rank[w])
                    .filter(|&c| c < r),
            );
            cols[start..].sort_unstable();
            cols.push(r);
            offsets.push(cols.len());
        }
        let values = vec![0.0; cols.len()];
        Self {
            offsets,
            cols,
            values,
        }
    }

    pub fn n(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Ajoute `v` au terme (i, j) et à son symétrique, qui doit être dans le profil.
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        let (i, j) = if j > i { (j, i) } else { (i, j) };
        let row = self.offsets[i]..self.offsets[i + 1];
        let k = self.cols[row.clone()]
            .binary_search(&j)
            .expect("terme hors du profil");
        self.values[row.start + k] += v;
    }

    /// y = A x.
    fn multiply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for i in 0..self.n() {
            let row = self.offsets[i]..self.offsets[i + 1] - 1;
            let mut s = self.values[row.end] * x[i];
            for (&j, &a) in self.cols[row.clone()].iter().zip(&self.values[row]) {
                s += a * x[j];
                y[j] += a * x[i];
            }
            y[i] += s;
        }
    }

    /// Facteur L de Cholesky incomplet IC(0) de A avec la diagonale multipliée par 1 + `shift` :
    /// L a le profil de A. `None` si un pivot n'est pas strictement positif.
    fn incomplete_cholesky(&self, shift: f64) -> Option<Vec<f64>> {
        let mut l = vec![0.0; self.values.len()];
        // w[c] : l_ic de la rangée i en cours, 0 hors des colonnes déjà calculées
        let mut w = vec![0.0; self.n()];
        for i in 0..self.n() {
            let (start, diag) = (self.offsets[i], self.offsets[i + 1] - 1);
            for p in start..diag {
                let j = self.cols[p];
                let (js, jd) = (self.offsets[j], self.offsets[j + 1] - 1);
                let s: f64 = (js..jd).map(|q| l[q] * w[self.cols[q]]).sum();
                l[p] = (self.values[p] - s) / l[jd];
                w[j] = l[p];
            }
            let squares: f64 = l[start..diag].iter().map(|x| x * x).sum();
            let pivot = self.values[diag] * (1.0 + shift) - squares;
            if pivot.is_nan() || pivot <= 1e-12 * self.values[diag] {
                return None;
            }
            l[diag] = pivot.sqrt();
            for p in start..diag {
                w[self.cols[p]] = 0.0;
            }
        }
        Some(l)
    }

    /// Résout L Lᵀ z = r (`l` : valeurs de L dans le profil de `self`).
    fn precondition(&self, l: &[f64], r: &[f64]) -> Vec<f64> {
        let mut z = r.to_vec();
        for i in 0..self.n() {
            let (start, diag) = (self.offsets[i], self.offsets[i + 1] - 1);
            let s: f64 = (start..diag).map(|p| l[p] * z[self.cols[p]]).sum();
            z[i] = (z[i] - s) / l[diag];
        }
        for i in (0..self.n()).rev() {
            let (start, diag) = (self.offsets[i], self.offsets[i + 1] - 1);
            z[i] /= l[diag];
            for p in start..diag {
                z[self.cols[p]] -= l[p] * z[i];
            }
        }
        z
    }
}

/// Gradients conjugués préconditionnés par IC(0) : prêt à résoudre A x = b pour autant de
/// seconds membres qu'on veut.
#[derive(Debug, Clone, PartialEq)]
pub struct Pcg {
    a: Csr,
    l: Vec<f64>,
}

impl Pcg {
    /// Renvoie `Singular` si IC(0) échoue même avec le plus grand décalage de la diagonale.
    pub fn new(a: Csr) -> Result<Self, MailleError> {
        let l = SHIFTS
            .iter()
            .find_map(|&shift| a.incomplete_cholesky(shift))
            .ok_or(MailleError::Singular)?;
        Ok(Self { a, l })
    }

    /// Résout A x = b jusqu'à un résidu relatif de 1e-12.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>, MailleError> {
        let norm = |v: &[f64]| dot(v, v).sqrt();
        let b_norm = norm(b);
        let mut x = vec![0.0; b.len()];
        if b_norm == 0.0 {
            return Ok(x);
        }
        let mut r = b.to_vec();
        let mut z = self.a.precondition(&self.l, &r);
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let mut ap = vec![0.0; b.len()];
        for _ in 0..MAX_ITERATIONS {
            self.a.multiply(&p, &mut ap);
            let curvature = dot(&p, &ap);
            if curvature.is_nan() || curvature <= 0.0 {
                return Err(MailleError::Singular);
            }
            let alpha = rz / curvature;
            for k in 0..x.len() {
                x[k] += alpha * p[k];
                r[k] -= alpha * ap[k];
            }
            if norm(&r) <= TOLERANCE * b_norm {
                return Ok(x);
            }
            z = self.a.precondition(&self.l, &r);
            let rz_next = dot(&r, &z);
            let beta = rz_next / rz;
            rz = rz_next;
            for k in 0..p.len() {
                p[k] = z[k] + beta * p[k];
            }
        }
        Err(MailleError::IterativeNotConverged {
            iterations: MAX_ITERATIONS,
            residual: norm(&r) / b_norm,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solve::reverse_cuthill_mckee;

    /// Laplacien (4 + degré, −1) d'une grille n × n, numéroté par RCM, et sa solution choisie.
    fn grid_system(n: usize) -> (Csr, Vec<f64>, Vec<f64>) {
        let id = |i: usize, j: usize| i + n * j;
        let mut elements = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                elements.push(vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]);
            }
        }
        let graph = Graph::from_elements(n * n, &elements);
        let order = reverse_cuthill_mckee(&graph);
        let mut rank = vec![0; n * n];
        for (r, &v) in order.iter().enumerate() {
            rank[v] = r;
        }
        let mut a = Csr::zeros(&graph, &order, &rank);
        for v in 0..n * n {
            let neighbours = graph.neighbours(v);
            a.add(rank[v], rank[v], 1.0 + neighbours.len() as f64);
            for &w in neighbours.iter().filter(|&&w| w < v) {
                a.add(rank[v], rank[w], -1.0);
            }
        }
        let x: Vec<f64> = (0..n * n).map(|k| (k as f64 * 0.37).sin()).collect();
        let mut b = vec![0.0; n * n];
        a.multiply(&x, &mut b);
        (a, x, b)
    }

    #[test]
    fn graph_lists_each_neighbour_once() {
        // deux triangles (0 1 2) et (1 2 3) : 1 et 2 sont dans les deux éléments
        let graph = Graph::from_elements(4, &[vec![0, 1, 2], vec![1, 2, 3]]);
        let mut n1 = graph.neighbours(1).to_vec();
        n1.sort_unstable();
        assert_eq!(n1, vec![0, 2, 3]);
        assert_eq!(graph.neighbours(0).len(), 2);
    }

    #[test]
    fn rcm_orders_a_csr_graph_like_the_vec_graph() {
        let graph = Graph::from_elements(5, &[vec![0, 1], vec![1, 2], vec![2, 3]]);
        let lists: Vec<Vec<usize>> = (0..5).map(|v| graph.neighbours(v).to_vec()).collect();
        assert_eq!(reverse_cuthill_mckee(&graph), reverse_cuthill_mckee(&lists));
    }

    #[test]
    fn pcg_recovers_the_chosen_solution() {
        let (a, x, b) = grid_system(30);
        let got = Pcg::new(a).unwrap().solve(&b).unwrap();
        for (g, e) in got.iter().zip(&x) {
            assert!((g - e).abs() < 1e-8, "{g} contre {e}");
        }
    }

    #[test]
    fn pcg_solves_several_right_hand_sides_and_a_zero_one() {
        let (a, _, b) = grid_system(10);
        let pcg = Pcg::new(a.clone()).unwrap();
        assert!(pcg
            .solve(&vec![0.0; 100])
            .unwrap()
            .iter()
            .all(|&v| v == 0.0));
        let x = pcg.solve(&b).unwrap();
        let mut ax = vec![0.0; 100];
        a.multiply(&x, &mut ax);
        assert!(ax.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-8));
    }

    #[test]
    fn incomplete_cholesky_of_a_tridiagonal_matrix_is_exact() {
        // sans remplissage possible, IC(0) est le Cholesky : L Lᵀ = A
        let graph = Graph::from_elements(4, &[vec![0, 1], vec![1, 2], vec![2, 3]]);
        let id: Vec<usize> = (0..4).collect();
        let mut a = Csr::zeros(&graph, &id, &id);
        for i in 0..4 {
            a.add(i, i, 2.0);
            if i > 0 {
                a.add(i, i - 1, -1.0);
            }
        }
        let l = a.incomplete_cholesky(0.0).unwrap();
        let z = a.precondition(&l, &[1.0, 0.0, 0.0, 1.0]);
        let mut back = vec![0.0; 4];
        a.multiply(&z, &mut back);
        assert!(back
            .iter()
            .zip([1.0, 0.0, 0.0, 1.0])
            .all(|(p, q)| (p - q).abs() < 1e-13));
    }

    #[test]
    fn singular_matrix_is_reported() {
        // −Laplacien d'une chaîne sans appui : singulière (le vecteur constant est dans le noyau), avec
        // un effort net non nul : pas de solution
        let graph = Graph::from_elements(3, &[vec![0, 1], vec![1, 2]]);
        let id: Vec<usize> = (0..3).collect();
        let mut a = Csr::zeros(&graph, &id, &id);
        a.add(0, 0, 1.0);
        a.add(1, 1, 2.0);
        a.add(2, 2, 1.0);
        a.add(1, 0, -1.0);
        a.add(2, 1, -1.0);
        let result = Pcg::new(a).and_then(|p| p.solve(&[1.0, 1.0, 1.0]));
        assert!(result.is_err(), "{result:?}");
    }
}
