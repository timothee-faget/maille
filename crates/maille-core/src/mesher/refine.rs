//! Bord respecté, puis raffinement de Ruppert : tant qu'un triangle du domaine a un angle trop
//! petit ou dépasse la taille locale h(x), on insère le centre de son cercle circonscrit. Si ce
//! centre empiète sur un morceau de bord (tombe dans son cercle diamétral), on coupe ce morceau en
//! son milieu à la place. Un morceau de bord dont le cercle diamétral contient le sommet opposé de
//! son triangle est coupé lui aussi, sauf si ce sommet est sur un autre bord parti du même coin
//! aigu.

use std::collections::{HashMap, HashSet, VecDeque};

use super::delaunay::Delaunay;
use super::geometry::{Boundary, Piece, Segment};
use super::predicates::distance;
use super::{MeshOptions, MAX_NODES};
use crate::error::MailleError;

pub(crate) struct Mesher<'a> {
    pub d: Delaunay,
    /// Morceaux de bord actuels (un morceau coupé garde son indice pour sa première moitié).
    pub segments: Vec<Segment>,
    pieces: &'a [Piece],
    /// Morceau de bord porté par chaque arête contrainte (indices de points, le plus petit d'abord).
    by_edge: HashMap<(usize, usize), usize>,
    /// Points du bord où le domaine fait un angle aigu (facettes du mailleur 3D ; vide en 2D, où
    /// `MIN_CORNER` les refuse) : un triangle qui en a un n'est raffiné que pour sa taille, pas
    /// pour ses angles, sinon Ruppert ne terminerait pas.
    pub acute: HashSet<usize>,
    /// Groupe de chaque Piece pour les coins aigus : les Pieces d'un même bord d'entrée (une arête
    /// découpée en cordes) partent ensemble du coin. Vide : chaque Piece est son propre groupe.
    pub groups: Vec<usize>,
    /// Coins aigus de chaque groupe qui en touche un, et (groupe, coin) des points posés sur ses
    /// Pieces : tenus à jour seulement s'il y a des coins aigus.
    piece_apex: HashMap<usize, Vec<usize>>,
    apex_at: HashMap<usize, Vec<(usize, usize)>>,
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
            acute: HashSet::new(),
            groups: Vec::new(),
            piece_apex: HashMap::new(),
            apex_at: HashMap::new(),
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
        for seg in &self.segments {
            let ends = [(seg.t[0] == 0.0, seg.a), (seg.t[1] == 1.0, seg.b)];
            for (end, v) in ends {
                if end && self.acute.contains(&v) {
                    let g = self.group(seg.piece);
                    self.piece_apex.entry(g).or_default().push(v);
                }
            }
        }
        for seg in &self.segments {
            let g = self.group(seg.piece);
            for &apex in self.piece_apex.get(&g).into_iter().flatten() {
                for v in [seg.a, seg.b] {
                    self.apex_at.entry(v).or_default().push((g, apex));
                }
            }
        }
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
                let sin = if v.iter().any(|i| self.acute.contains(i)) || self.in_wedge(p, v) {
                    0.0
                } else {
                    sin_min
                };
                if !bad(p, sin, size) {
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

    /// Insère p (`TooManyNodes` au-delà de `MAX_NODES`) et note les triangles créés.
    fn insert(&mut self, p: [f64; 2], created: &mut Vec<usize>) -> Result<usize, MailleError> {
        if self.d.point_count() >= MAX_NODES {
            return Err(MailleError::TooManyNodes { max: MAX_NODES });
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
        let g = self.group(seg.piece);
        for &apex in self.piece_apex.get(&g).into_iter().flatten() {
            self.apex_at.entry(m).or_default().push((g, apex));
        }
        self.segments[s].b = m;
        self.segments[s].t[1] = t;
        self.segments.push(Segment {
            a: m,
            t: [t, seg.t[1]],
            ..seg
        });
        Ok(self.segments.len() - 1)
    }

    fn group(&self, piece: usize) -> usize {
        self.groups.get(piece).copied().unwrap_or(piece)
    }

    /// Vrai si le plus petit côté du triangle relie deux groupes qui partent du même coin aigu : son
    /// petit angle vient de celui du coin, Ruppert ne peut pas l'améliorer (règle de Shewchuk).
    fn in_wedge(&self, p: [[f64; 2]; 3], v: [usize; 3]) -> bool {
        if self.apex_at.is_empty() {
            return false;
        }
        let k = (0..3)
            .min_by(|&i, &j| {
                distance(p[i], p[(i + 1) % 3]).total_cmp(&distance(p[j], p[(j + 1) % 3]))
            })
            .unwrap_or(0);
        let (Some(u), Some(w)) = (self.apex_at.get(&v[k]), self.apex_at.get(&v[(k + 1) % 3]))
        else {
            return false;
        };
        u.iter()
            .any(|&(pu, au)| w.iter().any(|&(pw, aw)| pu != pw && au == aw))
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
        // sommet posé sur un autre groupe parti du même coin aigu : couper ne ferait que
        // rapprocher les deux bords du coin, sans fin (règle de Shewchuk)
        let g = self.group(seg.piece);
        if let (Some(apexes), Some(at)) = (self.piece_apex.get(&g), self.apex_at.get(&v)) {
            if at
                .iter()
                .any(|&(other, apex)| other != g && apexes.contains(&apex))
            {
                return false;
            }
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesher::remesh::acute_corners;

    /// Deux pointes : la droite x = 1 et la droite y = 0 sont tangentes au cercle de centre (0, 1)
    /// et de rayon 1, en (1, 1) et en (0, 0). Chaque bord est découpé en cordes, un groupe par bord :
    /// le raffinement termine et couvre le polygone.
    #[test]
    fn pointes_entre_bords_tangents() {
        let n = 8;
        let mut ring = Vec::new();
        let mut groups = Vec::new();
        for i in 0..n {
            ring.push(([i as f64 / n as f64, 0.0], 0.2));
            groups.push(0);
        }
        for i in 0..n {
            ring.push(([1.0, i as f64 / n as f64], 0.2));
            groups.push(1);
        }
        for i in 0..n {
            let a = -std::f64::consts::FRAC_PI_2 * i as f64 / n as f64;
            ring.push(([a.cos(), 1.0 + a.sin()], 0.2));
            groups.push(2);
        }
        let rings = [ring];
        let boundary = Boundary::polygons(&rings, "pointes");
        let mut m = Mesher::new(&boundary).unwrap();
        m.acute = acute_corners(&rings);
        assert_eq!(m.acute.len(), 2);
        m.groups = groups;
        m.refine(&MeshOptions::default(), &SizeField::new(&boundary, 0.3))
            .unwrap();
        m.carve().unwrap();
        // aire du polygone de cordes
        let r = &rings[0];
        let polygon: f64 = (0..r.len())
            .map(|i| {
                let (a, b) = (r[i].0, r[(i + 1) % r.len()].0);
                a[0] * b[1] - a[1] * b[0]
            })
            .sum::<f64>()
            / 2.0;
        let area: f64 =
            m.d.inside_triangles()
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|v| m.d.point(v));
                    ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0
                })
                .sum();
        assert!((area - polygon).abs() < 1e-12, "{area} {polygon}");
    }

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
