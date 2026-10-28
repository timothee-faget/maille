//! Raffinement de l'intérieur d'un maillage 3D dont la surface est figée : tant qu'un tétraèdre du
//! domaine a un rapport rayon circonscrit / plus petite arête au-delà de `MAX_RATIO` ou dépasse la
//! taille locale h(x), on insère le centre de sa sphère circonscrite. Un centre hors du domaine, ou
//! qui empiète sur une face de la surface, est refusé : la surface ne bouge pas, et le tétraèdre
//! reste tel quel. La qualité est donc visée, pas garantie près de la surface ; les slivers
//! (tétraèdres plats de bon rapport rayon / arête) ne sont pas traités.

use std::collections::VecDeque;

use super::delaunay3::Delaunay3;
use super::predicates::distance3;
use super::surface::{cross, dot, sub, Surface};
use super::MAX_NODES;
use crate::error::MailleError;

/// Rapport rayon circonscrit / plus petite arête visé (2 : le seuil de Shewchuk).
pub(super) const MAX_RATIO: f64 = 2.0;

/// Rayon circonscrit du tétraèdre régulier de côté 1 : un tétraèdre de taille h a un rayon d'au
/// plus h √6 / 4.
const REGULAR_RADIUS: f64 = 0.612_372_435_695_794_5;

/// Raffine les tétraèdres du domaine de `del` (déjà classés, surface en murs) ; les points insérés
/// sont ajoutés à `points` (même numérotation que `del`).
pub(super) fn refine(
    del: &mut Delaunay3,
    points: &mut Vec<[f64; 3]>,
    size: &SizeField3,
) -> Result<(), MailleError> {
    let mut queue: VecDeque<usize> = (0..del.tet_slots()).collect();
    while let Some(t) = queue.pop_front() {
        let Some(v) = del.inside_tet(t) else {
            continue;
        };
        let p = v.map(|i| points[i]);
        let (center, radius) = circumsphere(p);
        let shortest = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .map(|(i, j)| distance3(p[i], p[j]))
            .fold(f64::INFINITY, f64::min);
        let centroid = std::array::from_fn(|k| p.iter().map(|q| q[k]).sum::<f64>() / 4.0);
        let limit = size.at(centroid) * REGULAR_RADIUS * (1.0 + 1e-9);
        if radius <= MAX_RATIO * shortest && radius <= limit {
            continue;
        }
        if points.len() >= MAX_NODES {
            return Err(MailleError::TooManyNodes { max: MAX_NODES });
        }
        if del.insert_inside(center).is_some() {
            points.push(center);
            queue.extend(del.new_tets());
        }
    }
    Ok(())
}

/// Centre et rayon de la sphère circonscrite d'un tétraèdre.
pub(super) fn circumsphere([a, b, c, d]: [[f64; 3]; 4]) -> ([f64; 3], f64) {
    let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
    let (x, y, z) = (cross(v, w), cross(w, u), cross(u, v));
    let (uu, vv, ww) = (dot(u, u), dot(v, v), dot(w, w));
    let det = 2.0 * dot(u, x);
    let o: [f64; 3] = std::array::from_fn(|k| (uu * x[k] + vv * y[k] + ww * z[k]) / det);
    (std::array::from_fn(|k| a[k] + o[k]), dot(o, o).sqrt())
}

/// Taille locale h(x) = min(h_max, minᵢ(hᵢ + g |x − pᵢ|)) sur les sommets pᵢ de la surface, hᵢ la
/// longueur moyenne de leurs arêtes (au plus h_max) et g la gradation : comme le `SizeField` du 2D,
/// h vaut la finesse de la surface la plus proche et croît au plus de g par unité de longueur. Les
/// sommets sont rangés dans une grille, parcourue par couronnes jusqu'à ce qu'aucune case plus
/// lointaine ne puisse faire mieux.
pub(super) struct SizeField3 {
    max: f64,
    grading: f64,
    smallest: f64,
    lo: [f64; 3],
    cell: f64,
    n: [usize; 3],
    cells: Vec<Vec<([f64; 3], f64)>>,
}

/// Cases au plus par direction de la grille.
const MAX_CELLS: usize = 64;

impl SizeField3 {
    pub fn new(surface: &Surface, max: f64, grading: f64) -> Self {
        let mut sum = vec![(0.0, 0usize); surface.vertices.len()];
        for t in &surface.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let l = distance3(surface.vertices[a], surface.vertices[b]);
                for v in [a, b] {
                    sum[v].0 += l;
                    sum[v].1 += 1;
                }
            }
        }
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for p in &surface.vertices {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let extent = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        // au plus MAX_CELLS cases par direction, sans tronquer la dernière
        let cell = max
            .max(extent / (MAX_CELLS - 1) as f64)
            .max(f64::MIN_POSITIVE);
        let n =
            std::array::from_fn(|k| (((hi[k] - lo[k]) / cell).floor() as usize + 1).min(MAX_CELLS));
        let mut field = Self {
            max,
            grading,
            smallest: max,
            lo,
            cell,
            n,
            cells: vec![Vec::new(); n[0] * n[1] * n[2]],
        };
        for (p, &(total, count)) in surface.vertices.iter().zip(&sum) {
            let h = if count == 0 {
                max
            } else {
                (total / count as f64).min(max)
            };
            field.smallest = field.smallest.min(h);
            let i = field.index(field.cell_of(*p));
            field.cells[i].push((*p, h));
        }
        field
    }

    fn cell_of(&self, p: [f64; 3]) -> [usize; 3] {
        std::array::from_fn(|k| {
            (((p[k] - self.lo[k]) / self.cell).floor().max(0.0) as usize).min(self.n[k] - 1)
        })
    }

    fn index(&self, c: [usize; 3]) -> usize {
        (c[0] * self.n[1] + c[1]) * self.n[2] + c[2]
    }

    pub fn at(&self, x: [f64; 3]) -> f64 {
        let c = self.cell_of(x);
        let mut best = self.max;
        let rings = self.n.iter().copied().max().unwrap_or(1);
        for r in 0..rings {
            // tout sommet de la couronne r est à plus de (r − 1) cases de x
            let reach = self.smallest + self.grading * (r.saturating_sub(1) as f64) * self.cell;
            if reach >= best {
                break;
            }
            let range = |k: usize| c[k].saturating_sub(r)..=(c[k] + r).min(self.n[k] - 1);
            for i in range(0) {
                for j in range(1) {
                    for l in range(2) {
                        let ring = [i, j, l]
                            .iter()
                            .zip(&c)
                            .map(|(&a, &b)| a.abs_diff(b))
                            .max()
                            .unwrap_or(0);
                        if ring != r {
                            continue;
                        }
                        for &(p, h) in &self.cells[self.index([i, j, l])] {
                            best = best.min(h + self.grading * distance3(x, p));
                        }
                    }
                }
            }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesher::testsurface::*;

    #[test]
    fn sphere_circonscrite_equidistante() {
        let p = [
            [0.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [0.3, 1.5, 0.2],
            [0.4, 0.2, 1.7],
        ];
        let (c, r) = circumsphere(p);
        for q in p {
            assert!((distance3(c, q) - r).abs() < 1e-12);
        }
        let side = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3f64.sqrt() / 2.0, 0.0],
        ];
        let apex = [0.5, 3f64.sqrt() / 6.0, (2.0f64 / 3.0).sqrt()];
        let (_, r) = circumsphere([side[0], side[1], side[2], apex]);
        assert!((r - REGULAR_RADIUS).abs() < 1e-12);
    }

    #[test]
    fn taille_locale_comme_la_recherche_exhaustive() {
        let s = surface(&cube(4));
        let (h, g) = (0.8, 0.3);
        let field = SizeField3::new(&s, h, g);
        // référence : tous les sommets, arêtes de 0,25 (grille 4 × 4) et diagonales
        let mut rng = 1u64;
        for _ in 0..200 {
            let x: [f64; 3] = std::array::from_fn(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                (rng >> 11) as f64 / (1u64 << 53) as f64
            });
            let want = field
                .cells
                .iter()
                .flatten()
                .map(|&(p, hp)| hp + g * distance3(x, p))
                .fold(h, f64::min);
            assert!((field.at(x) - want).abs() < 1e-12);
        }
        // au bord, la finesse de la surface (arêtes de 0,25 et diagonales de 0,354) ; au centre,
        // plus grand
        assert!(field.at([0.0, 0.0, 0.0]) < 0.354);
        assert!(field.at([0.5, 0.5, 0.5]) > field.at([0.0, 0.0, 0.0]));
    }
}
