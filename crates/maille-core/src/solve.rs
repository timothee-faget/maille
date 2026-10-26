//! Solveur creux direct : stockage en ligne de ciel (skyline), renumérotation de Cuthill-McKee
//! inverse pour réduire le profil, puis Cholesky en place.

use crate::error::MailleError;

/// Matrice symétrique stockée en ligne de ciel : pour chaque ligne i, seule la partie
/// triangulaire inférieure entre la colonne `first[i]` et la diagonale est gardée. Le remplissage
/// du Cholesky reste dans ce profil, d'où l'intérêt de le réduire (voir `reverse_cuthill_mckee`).
#[derive(Debug, Clone, PartialEq)]
pub struct SkylineMatrix {
    first: Vec<usize>,
    /// La ligne i occupe `data[start[i]..start[i + 1]]`, de la colonne `first[i]` à i.
    start: Vec<usize>,
    data: Vec<f64>,
}

impl SkylineMatrix {
    /// Matrice nulle de profil `first` (`first[i] ≤ i`, colonne du premier terme non nul de la ligne i).
    pub fn zeros(first: Vec<usize>) -> Self {
        let mut start = Vec::with_capacity(first.len() + 1);
        start.push(0);
        for (i, &f) in first.iter().enumerate() {
            start.push(start[i] + i + 1 - f);
        }
        let len = start[first.len()];
        Self {
            first,
            start,
            data: vec![0.0; len],
        }
    }

    pub fn n(&self) -> usize {
        self.first.len()
    }

    /// Nombre de termes stockés : la mémoire du Cholesky, que la renumérotation cherche à réduire.
    pub fn stored(&self) -> usize {
        self.data.len()
    }

    fn row(&self, i: usize) -> &[f64] {
        &self.data[self.start[i]..self.start[i + 1]]
    }

    /// Terme (i, j), nul hors du profil.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let (i, j) = if j > i { (j, i) } else { (i, j) };
        if j < self.first[i] {
            0.0
        } else {
            self.row(i)[j - self.first[i]]
        }
    }

    /// Ajoute `v` au terme (i, j) et à son symétrique, qui doit être dans le profil.
    pub fn add(&mut self, i: usize, j: usize, v: f64) {
        let (i, j) = if j > i { (j, i) } else { (i, j) };
        self.data[self.start[i] + j - self.first[i]] += v;
    }

    /// Résout A x = b par Cholesky, voir `factorize`.
    pub fn cholesky_solve(self, b: &[f64]) -> Result<Vec<f64>, MailleError> {
        Ok(self.factorize()?.solve(b))
    }

    /// Factorisation de Cholesky A = L Lᵀ, L écrasant A dans le même profil : une fois faite, elle
    /// résout A x = b pour autant de seconds membres b qu'on veut.
    /// Renvoie `Singular` si un pivot est négligeable devant la plus grande diagonale.
    pub fn factorize(mut self) -> Result<Cholesky, MailleError> {
        let n = self.n();
        let max_diag = (0..n).map(|i| self.get(i, i)).fold(0.0_f64, f64::max);
        let tol = 1e-12 * max_diag;

        for i in 0..n {
            let fi = self.first[i];
            for j in fi..=i {
                // Σ_k L_ik L_jk sur les colonnes communes aux deux profils
                let k0 = fi.max(self.first[j]);
                let (ri, rj) = (self.start[i] - fi, self.start[j] - self.first[j]);
                let s = dot(&self.data[ri + k0..ri + j], &self.data[rj + k0..rj + j]);
                let a = self.data[ri + j] - s;
                if i == j {
                    if a.is_nan() || a <= tol {
                        return Err(MailleError::Singular);
                    }
                    self.data[ri + j] = a.sqrt();
                } else {
                    self.data[ri + j] = a / self.data[rj + j];
                }
            }
        }
        Ok(Cholesky(self))
    }
}

/// Produit scalaire de deux tranches de même longueur. Quatre sommes partielles indépendantes
/// : le Cholesky passe l'essentiel de son temps ici, et une seule somme force chaque addition à
/// attendre la précédente.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let ((a4, a_tail), (b4, b_tail)) = (a.as_chunks::<8>(), b.as_chunks::<8>());
    let tail: f64 = a_tail.iter().zip(b_tail).map(|(x, y)| x * y).sum();
    let mut s = [0.0; 8];
    for (x, y) in a4.iter().zip(b4) {
        for k in 0..8 {
            s[k] += x[k] * y[k];
        }
    }
    ((s[0] + s[1]) + (s[2] + s[3])) + ((s[4] + s[5]) + (s[6] + s[7])) + tail
}

/// Facteur L de Cholesky, dans le profil de la matrice factorisée.
#[derive(Debug, Clone, PartialEq)]
pub struct Cholesky(SkylineMatrix);

impl Cholesky {
    /// Résout L Lᵀ x = b.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let l = &self.0;
        // L y = b
        let mut x = b.to_vec();
        for i in 0..l.n() {
            let row = l.row(i);
            let (fi, diag) = (l.first[i], row[i - l.first[i]]);
            let s: f64 = (fi..i).map(|k| row[k - fi] * x[k]).sum();
            x[i] = (x[i] - s) / diag;
        }
        // Lᵀ x = y, par colonne : la ligne i de L est la colonne i de Lᵀ
        for i in (0..l.n()).rev() {
            let row = l.row(i);
            let fi = l.first[i];
            x[i] /= row[i - fi];
            for k in fi..i {
                x[k] -= row[k - fi] * x[i];
            }
        }
        x
    }
}

/// Renumérotation de Cuthill-McKee inverse du graphe `adjacency` (voisins de chaque sommet) :
/// `order[nouveau] = ancien`. Parcours en largeur depuis un sommet pseudo-périphérique, voisins
/// par degré croissant, puis ordre inversé. Chaque composante connexe est traitée à part.
pub fn reverse_cuthill_mckee(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let degree = |v: usize| adjacency[v].len();
    let mut visited = vec![false; n];
    let mut order = Vec::with_capacity(n);
    while let Some(seed) = (0..n).filter(|&v| !visited[v]).min_by_key(|&v| degree(v)) {
        // sommet pseudo-périphérique (George et Liu) : tant que le sommet de plus petit degré
        // du dernier niveau est plus excentré, on repart de lui
        let mut start = seed;
        let mut levels = bfs_levels(adjacency, start, &visited);
        loop {
            let last = levels.last().map_or(&[][..], |l| l.as_slice());
            let Some(&far) = last.iter().min_by_key(|&&v| degree(v)) else {
                break;
            };
            let far_levels = bfs_levels(adjacency, far, &visited);
            if far_levels.len() <= levels.len() {
                break;
            }
            (start, levels) = (far, far_levels);
        }

        let first = order.len();
        visited[start] = true;
        order.push(start);
        let mut head = first;
        while head < order.len() {
            let v = order[head];
            head += 1;
            let mut next: Vec<usize> = adjacency[v]
                .iter()
                .copied()
                .filter(|&w| !visited[w])
                .collect();
            next.sort_by_key(|&w| degree(w));
            for w in next {
                if !std::mem::replace(&mut visited[w], true) {
                    order.push(w);
                }
            }
        }
    }
    order.reverse();
    order
}

/// Niveaux du parcours en largeur depuis `start`, sans entrer dans les sommets `excluded`.
fn bfs_levels(adjacency: &[Vec<usize>], start: usize, excluded: &[bool]) -> Vec<Vec<usize>> {
    let mut seen = excluded.to_vec();
    seen[start] = true;
    let mut levels = vec![vec![start]];
    loop {
        let mut next = Vec::new();
        for &v in levels.last().map_or(&[][..], |l| l.as_slice()) {
            for &w in &adjacency[v] {
                if !std::mem::replace(&mut seen[w], true) {
                    next.push(w);
                }
            }
        }
        if next.is_empty() {
            return levels;
        }
        levels.push(next);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Profil plein (dense) n × n rempli avec `a`.
    fn dense(a: &[&[f64]]) -> SkylineMatrix {
        let mut m = SkylineMatrix::zeros(vec![0; a.len()]);
        for (i, row) in a.iter().enumerate() {
            for (j, &v) in row.iter().enumerate().take(i + 1) {
                m.add(i, j, v);
            }
        }
        m
    }

    #[test]
    fn cholesky_solves_a_known_spd_system() {
        // A x = b avec x = (1, −2, 3, 0.5) choisi d'avance
        let a: [&[f64]; 4] = [
            &[4.0, 1.0, 0.0, 0.5],
            &[1.0, 5.0, 2.0, 0.0],
            &[0.0, 2.0, 6.0, 1.0],
            &[0.5, 0.0, 1.0, 3.0],
        ];
        let x = [1.0, -2.0, 3.0, 0.5];
        let b: Vec<f64> = a
            .iter()
            .map(|row| row.iter().zip(x).map(|(a, x)| a * x).sum())
            .collect();
        let got = dense(&a).cholesky_solve(&b).unwrap();
        for (g, e) in got.iter().zip(x) {
            assert!((g - e).abs() < 1e-13, "{got:?}");
        }
    }

    #[test]
    fn variable_profile_gives_the_same_solution_as_the_full_one() {
        // tridiagonale (−1, 2, −1) + un terme lointain (4, 1) : profil [0, 0, 1, 2, 1]
        let n = 5;
        let mut full = SkylineMatrix::zeros(vec![0; n]);
        let mut sky = SkylineMatrix::zeros(vec![0, 0, 1, 2, 1]);
        for m in [&mut full, &mut sky] {
            for i in 0..n {
                m.add(i, i, 2.0);
                if i > 0 {
                    m.add(i, i - 1, -1.0);
                }
            }
            m.add(1, 4, 0.5);
        }
        assert_eq!(sky.stored(), 1 + 2 + 2 + 2 + 4);
        assert!((sky.get(1, 4) - 0.5).abs() < 1e-15 && sky.get(0, 3).abs() < 1e-15);
        let b = [1.0, 0.0, 2.0, -1.0, 3.0];
        let (x_full, x_sky) = (
            full.cholesky_solve(&b).unwrap(),
            sky.cholesky_solve(&b).unwrap(),
        );
        for (a, b) in x_full.iter().zip(&x_sky) {
            assert!((a - b).abs() < 1e-13);
        }
    }

    #[test]
    fn one_factorization_solves_several_right_hand_sides() {
        let a: [&[f64]; 3] = [&[4.0, 1.0, 0.0], &[1.0, 3.0, 1.0], &[0.0, 1.0, 2.0]];
        let factor = dense(&a).factorize().unwrap();
        for b in [[1.0, 0.0, 0.0], [0.5, -2.0, 3.0]] {
            let (x, direct) = (factor.solve(&b), dense(&a).cholesky_solve(&b).unwrap());
            for (row, bi) in a.iter().zip(b) {
                let ax: f64 = row.iter().zip(&x).map(|(a, x)| a * x).sum();
                assert!((ax - bi).abs() < 1e-14);
            }
            assert!(x.iter().zip(&direct).all(|(a, b)| (a - b).abs() < 1e-15));
        }
    }

    #[test]
    fn dot_matches_the_plain_sum_for_every_length_around_the_chunk_size() {
        for n in 0..40 {
            let a: Vec<f64> = (0..n).map(|k| 0.5 + k as f64).collect();
            let b: Vec<f64> = (0..n).map(|k| 2.0 - 0.25 * k as f64).collect();
            let plain: f64 = a.iter().zip(&b).map(|(x, y)| x * y).sum();
            assert!(
                (dot(&a, &b) - plain).abs() <= 1e-12 * (1.0 + plain.abs()),
                "n = {n}"
            );
        }
    }

    #[test]
    fn singular_matrix_is_reported() {
        let a: [&[f64]; 2] = [&[1.0, 1.0], &[1.0, 1.0]];
        assert_eq!(
            dense(&a).cholesky_solve(&[1.0, 1.0]),
            Err(MailleError::Singular)
        );
    }

    /// Profil (nombre de termes stockés) de la matrice du graphe numérotée selon `order`.
    fn profile(adjacency: &[Vec<usize>], order: &[usize]) -> usize {
        let mut rank = vec![0; order.len()];
        for (new, &old) in order.iter().enumerate() {
            rank[old] = new;
        }
        (0..adjacency.len())
            .map(|v| {
                let lowest = adjacency[v]
                    .iter()
                    .map(|&w| rank[w])
                    .fold(rank[v], usize::min);
                rank[v] - lowest + 1
            })
            .sum()
    }

    #[test]
    fn rcm_is_a_permutation_that_shrinks_the_profile_of_a_shuffled_grid() {
        // grille 6 × 6 (voisins à 4), numérotée en mélangeant les lignes et les colonnes
        let side = 6;
        let shuffle = [3, 0, 5, 1, 4, 2];
        let id = |i: usize, j: usize| shuffle[i] * side + shuffle[j];
        let mut adjacency = vec![Vec::new(); side * side];
        for i in 0..side {
            for j in 0..side {
                if i + 1 < side {
                    adjacency[id(i, j)].push(id(i + 1, j));
                    adjacency[id(i + 1, j)].push(id(i, j));
                }
                if j + 1 < side {
                    adjacency[id(i, j)].push(id(i, j + 1));
                    adjacency[id(i, j + 1)].push(id(i, j));
                }
            }
        }
        let order = reverse_cuthill_mckee(&adjacency);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..side * side).collect::<Vec<_>>());

        let identity: Vec<usize> = (0..side * side).collect();
        let (before, after) = (profile(&adjacency, &identity), profile(&adjacency, &order));
        // numérotation ligne par ligne d'une grille propre : demi-largeur de bande = side
        assert!(after <= side * side * (side + 1), "après RCM : {after}");
        assert!(after < before, "avant {before}, après {after}");
    }

    #[test]
    fn rcm_handles_disconnected_graphs_and_isolated_vertices() {
        // deux arêtes 0-2 et 1-3, sommet 4 isolé
        let adjacency = vec![vec![2], vec![3], vec![0], vec![1], vec![]];
        let mut order = reverse_cuthill_mckee(&adjacency);
        order.sort_unstable();
        assert_eq!(order, vec![0, 1, 2, 3, 4]);
    }
}
