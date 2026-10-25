//! Bord respecté, puis raffinement de Ruppert : tant qu'un triangle du domaine a un angle trop
//! petit ou dépasse la taille locale h(x), on insère le centre de son cercle circonscrit. Si ce
//! centre empiète sur un morceau de bord (tombe dans son cercle diamétral), on coupe ce morceau en
//! son milieu à la place. Un morceau de bord dont le cercle diamétral contient le sommet opposé de
//! son triangle est coupé lui aussi.

use std::collections::{HashMap, VecDeque};

use super::delaunay::Delaunay;
use super::geometry::{Boundary, Piece, Segment};
use super::{MeshOptions, MAX_NODES};
use crate::error::MailleError;

pub(crate) struct Mesher<'a> {
    pub d: Delaunay,
    /// Morceaux de bord actuels (un morceau coupé garde son indice pour sa première moitié).
    pub segments: Vec<Segment>,
    pieces: &'a [Piece],
    /// Morceau de bord porté par chaque arête contrainte (indices de points, le plus petit d'abord).
    by_edge: HashMap<(usize, usize), usize>,
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

impl<'a> Mesher<'a> {
    /// Triangulation des points du bord, chaque morceau de bord présent et contraint, domaine délimité.
    pub fn new(boundary: &'a Boundary) -> Result<Self, MailleError> {
        let (mut min, mut max) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in &boundary.points {
            for j in 0..2 {
                min[j] = min[j].min(p[j]);
                max[j] = max[j].max(p[j]);
            }
        }
        let mut d = Delaunay::new(min, max);
        for &p in &boundary.points {
            d.insert(p)?; // indice dans la triangulation = indice dans boundary.points
        }
        let mut m = Self {
            d,
            segments: boundary.segments.clone(),
            pieces: &boundary.pieces,
            by_edge: HashMap::new(),
        };
        m.recover(
            (0..m.segments.len()).collect(),
            &mut Vec::new(),
            &mut Vec::new(),
        )?;
        m.carve()?;
        Ok(m)
    }

    /// Délimite le domaine à partir des morceaux de bord actuels.
    pub fn carve(&mut self) -> Result<(), MailleError> {
        self.d.carve(self.segments.iter().map(|s| (s.a, s.b)))
    }

    /// Raffine jusqu'à ce que chaque triangle du domaine ait ses angles d'au moins
    /// `options.min_angle` et un rayon circonscrit d'au plus h/√3 (celui du triangle équilatéral de
    /// côté h), h étant pris au centre de gravité.
    pub fn refine(&mut self, options: &MeshOptions, size: &SizeField) -> Result<(), MailleError> {
        let sin_min = options.min_angle.to_radians().sin();
        let mut tris: VecDeque<usize> = (0..self.d.triangle_slots()).collect();
        let mut segs: Vec<usize> = (0..self.segments.len()).collect();
        loop {
            let mut created = Vec::new();
            if let Some(s) = segs.pop() {
                if self.encroached(s) {
                    self.split(s, &mut created)?;
                }
            } else {
                let Some(t) = tris.pop_front() else {
                    return Ok(());
                };
                let Some(v) = self.d.triangle(t) else {
                    continue;
                };
                let p = v.map(|i| self.d.point(i));
                if !bad(p, sin_min, size) {
                    continue;
                }
                let c = circumcenter(p);
                let hit: Vec<usize> = (0..self.segments.len())
                    .filter(|&s| self.live(s) && self.diametral(s, c) < 1.0 - 1e-12)
                    .collect();
                if !hit.is_empty() {
                    for s in hit {
                        self.split(s, &mut created)?;
                    }
                } else if !self.d.inside_at(c)? {
                    // centre hors du domaine sans empiéter sur aucun morceau : seulement possible à
                    // l'arrondi près, on coupe le morceau dont le cercle diamétral est le plus proche
                    let nearest = (0..self.segments.len())
                        .filter(|&s| self.live(s))
                        .min_by(|&r, &s| self.diametral(r, c).total_cmp(&self.diametral(s, c)))
                        .ok_or(MailleError::MeshingFailed)?;
                    self.split(nearest, &mut created)?;
                } else {
                    self.insert(c, &mut created)?;
                }
                tris.push_back(t); // revu s'il a survécu
            }
            // tout morceau de bord bordé par un nouveau triangle peut être devenu empiété
            for &t in &created {
                if let Some(v) = self.d.triangle(t) {
                    for k in 0..3 {
                        if let Some(&s) = self.by_edge.get(&key(v[k], v[(k + 1) % 3])) {
                            segs.push(s);
                        }
                    }
                }
            }
            tris.extend(created);
        }
    }

    /// Insère p (`MeshingFailed` au-delà de `MAX_NODES`) et note les triangles créés.
    fn insert(&mut self, p: [f64; 2], created: &mut Vec<usize>) -> Result<usize, MailleError> {
        if self.d.point_count() >= MAX_NODES {
            return Err(MailleError::MeshingFailed);
        }
        let i = self.d.insert(p)?;
        created.extend(self.d.new_triangles());
        Ok(i)
    }

    /// Coupe le morceau de bord présent `s` en son milieu, respecte ses deux moitiés, puis corrige
    /// l'appartenance au domaine des triangles qui les bordent.
    fn split(&mut self, s: usize, created: &mut Vec<usize>) -> Result<(), MailleError> {
        let seg = self.segments[s];
        self.by_edge.remove(&key(seg.a, seg.b));
        self.d.set_constrained(seg.a, seg.b, false);
        let second = self.halve(s, created)?;
        let mut done = Vec::new();
        self.recover(vec![s, second], &mut done, created)?;
        for s in done {
            let seg = self.segments[s];
            self.d.mark_sides(seg.a, seg.b, created);
        }
        Ok(())
    }

    /// Respecte les morceaux `pending` : un morceau présent devient contraint (et noté dans `done`),
    /// un morceau absent est coupé en deux et ses moitiés sont revues.
    fn recover(
        &mut self,
        mut pending: Vec<usize>,
        done: &mut Vec<usize>,
        created: &mut Vec<usize>,
    ) -> Result<(), MailleError> {
        while let Some(s) = pending.pop() {
            let seg = self.segments[s];
            if self.d.has_edge(seg.a, seg.b) {
                self.d.set_constrained(seg.a, seg.b, true);
                self.by_edge.insert(key(seg.a, seg.b), s);
                done.push(s);
            } else {
                let second = self.halve(s, created)?;
                pending.extend([s, second]);
            }
        }
        Ok(())
    }

    /// Insère le milieu du morceau `s` sur sa Curve : `s` devient la première moitié, la seconde est
    /// ajoutée à la fin et son indice renvoyé.
    fn halve(&mut self, s: usize, created: &mut Vec<usize>) -> Result<usize, MailleError> {
        let seg = self.segments[s];
        let t = (seg.t[0] + seg.t[1]) / 2.0;
        let m = self.insert(self.pieces[seg.piece].point(t), created)?;
        self.segments[s].b = m;
        self.segments[s].t[1] = t;
        self.segments.push(Segment {
            a: m,
            t: [t, seg.t[1]],
            ..seg
        });
        Ok(self.segments.len() - 1)
    }

    /// Vrai si le morceau `s` est présent et contraint.
    fn live(&self, s: usize) -> bool {
        let seg = self.segments[s];
        self.by_edge.get(&key(seg.a, seg.b)) == Some(&s)
    }

    /// Distance de p au milieu du morceau `s`, rapportée au rayon de son cercle diamétral : moins
    /// de 1 si p y est.
    fn diametral(&self, s: usize, p: [f64; 2]) -> f64 {
        let seg = self.segments[s];
        let (a, b) = (self.d.point(seg.a), self.d.point(seg.b));
        let mid = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
        distance(p, mid) / (distance(a, b) / 2.0)
    }

    /// Vrai si le sommet opposé au morceau `s`, dans son triangle du domaine, est dans son cercle
    /// diamétral (angle obtus en ce sommet).
    fn encroached(&self, s: usize) -> bool {
        if !self.live(s) {
            return false;
        }
        let seg = self.segments[s];
        let Some(v) = self.d.left_vertex(seg.a, seg.b) else {
            return false;
        };
        let (a, b, c) = (self.d.point(seg.a), self.d.point(seg.b), self.d.point(v));
        let (u, w) = ([a[0] - c[0], a[1] - c[1]], [b[0] - c[0], b[1] - c[1]]);
        u[0] * w[0] + u[1] * w[1] < -1e-12 * u[0].hypot(u[1]) * w[0].hypot(w[1])
    }
}

/// Taille locale h(x) = min(h_max, minᵢ(hᵢ + g |x − pᵢ|)) sur les points pᵢ du bord découpé, hᵢ
/// la taille de leur Curve et g la gradation : h vaut la taille du bord le plus proche et croît
/// au plus de g par unité de longueur en s'en éloignant.
pub(crate) struct SizeField {
    points: Vec<([f64; 2], f64)>,
    max: f64,
    grading: f64,
}

impl SizeField {
    pub fn new(boundary: &Boundary, grading: f64) -> Self {
        let points: Vec<([f64; 2], f64)> = boundary
            .segments
            .iter()
            .map(|s| (boundary.points[s.a], boundary.pieces[s.piece].size))
            .collect();
        let max = points.iter().map(|&(_, h)| h).fold(0.0, f64::max);
        Self {
            points,
            max,
            grading,
        }
    }

    pub fn at(&self, x: [f64; 2]) -> f64 {
        self.points
            .iter()
            .map(|&(p, h)| h + self.grading * distance(x, p))
            .fold(self.max, f64::min)
    }
}

/// Vrai si le triangle p a un angle de sinus inférieur à `sin_min` ou un rayon circonscrit
/// au-delà de h/√3. Le plus petit angle fait face au plus petit côté : sin θ = l_min / 2R.
fn bad(p: [[f64; 2]; 3], sin_min: f64, size: &SizeField) -> bool {
    let l = [0, 1, 2].map(|k| distance(p[k], p[(k + 1) % 3]));
    let area2 =
        (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
    let radius = l[0] * l[1] * l[2] / (2.0 * area2);
    let lmin = l.iter().copied().fold(f64::INFINITY, f64::min);
    let centroid = [0, 1].map(|j| (p[0][j] + p[1][j] + p[2][j]) / 3.0);
    lmin / (2.0 * radius) < sin_min || radius > (1.0 + 1e-9) * size.at(centroid) / 3f64.sqrt()
}

fn circumcenter(p: [[f64; 2]; 3]) -> [f64; 2] {
    let (bx, by) = (p[1][0] - p[0][0], p[1][1] - p[0][1]);
    let (cx, cy) = (p[2][0] - p[0][0], p[2][1] - p[0][1]);
    let d = 2.0 * (bx * cy - by * cx);
    let (b2, c2) = (bx * bx + by * by, cx * cx + cy * cy);
    [
        p[0][0] + (cy * b2 - by * c2) / d,
        p[0][1] + (bx * c2 - cx * b2) / d,
    ]
}

fn distance(p: [f64; 2], q: [f64; 2]) -> f64 {
    (p[0] - q[0]).hypot(p[1] - q[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circumcenter_is_equidistant() {
        let p = [[0.0, 0.0], [2.0, 0.0], [0.5, 1.5]];
        let c = circumcenter(p);
        let r = distance(c, p[0]);
        assert!((distance(c, p[1]) - r).abs() < 1e-12 && (distance(c, p[2]) - r).abs() < 1e-12);
    }

    #[test]
    fn equilateral_triangle_of_side_h_is_just_good() {
        let h = 0.1;
        let p = [[0.0, 0.0], [h, 0.0], [h / 2.0, h * 3f64.sqrt() / 2.0]];
        let size = SizeField {
            points: vec![([0.0, 0.0], h)],
            max: h,
            grading: 0.3,
        };
        let sin25 = 25f64.to_radians().sin();
        // h(x) = h partout (h_max = h)
        assert!(!bad(p, sin25, &size));
        let bigger = p.map(|q| [q[0] * 1.01, q[1] * 1.01]);
        assert!(bad(bigger, sin25, &size));
        // triangle plat : angle de ~11°
        assert!(bad(
            [[0.0, 0.0], [h, 0.0], [h / 2.0, h / 10.0]],
            sin25,
            &size
        ));
    }
}
