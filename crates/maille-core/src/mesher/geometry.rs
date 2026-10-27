//! Géométrie d'entrée du mailleur : un contour extérieur et des trous, chacun fermé et fait de
//! Curves (segments ou arcs de cercle). Vérification, orientation, puis découpage des Curves en
//! morceaux de longueur au plus h.

use std::f64::consts::PI;
use std::fmt;

use super::predicates::{distance, orient2d};
use super::MAX_NODES;
use crate::error::MailleError;

/// Angle intérieur minimal entre deux Curves qui se suivent, en degrés : en dessous, Ruppert peut ne
/// pas terminer.
pub const MIN_CORNER: f64 = 60.0;

/// Morceau de bord qui va de `start` au `start` de la Curve suivante du même contour (la dernière
/// revient à la première).
#[derive(Debug, Clone, PartialEq)]
pub struct Curve {
    pub start: [f64; 2],
    pub kind: CurveKind,
    /// Taille cible h des triangles le long de cette Curve, strictement positive.
    pub size: f64,
    /// Étiquette du bord, pour y poser des BC ; plusieurs Curves peuvent la partager.
    pub name: String,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveKind {
    Line,
    /// Arc de cercle de centre `center`, parcouru dans le sens horaire si `clockwise`. Le point
    /// d'arrivée doit être à la même distance du centre que `start` ; s'il est confondu avec
    /// `start`, l'arc est un cercle complet.
    Arc {
        center: [f64; 2],
        clockwise: bool,
    },
}

impl Curve {
    pub fn line(start: [f64; 2], size: f64, name: &str) -> Self {
        Self {
            start,
            kind: CurveKind::Line,
            size,
            name: name.to_string(),
        }
    }

    pub fn arc(start: [f64; 2], center: [f64; 2], clockwise: bool, size: f64, name: &str) -> Self {
        Self {
            start,
            kind: CurveKind::Arc { center, clockwise },
            size,
            name: name.to_string(),
        }
    }
}

/// Domaine à mailler. Le sens de parcours des contours est libre : le mailleur l'oriente lui-même.
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    pub outer: Vec<Curve>,
    pub holes: Vec<Vec<Curve>>,
}

/// Repère d'une Curve dans la Geometry, pour les messages d'erreur.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CurveId {
    /// `None` pour le contour extérieur, sinon l'indice du trou.
    pub hole: Option<usize>,
    /// Rang de la Curve dans son contour, tel que l'utilisateur l'a donné.
    pub index: usize,
    pub name: String,
}

impl fmt::Display for CurveId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "la Curve {} « {} » ", self.index, self.name)?;
        match self.hole {
            None => write!(f, "du contour extérieur"),
            Some(h) => write!(f, "du trou {h}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Path {
    Line,
    /// Angle de départ `theta0`, balayage `sweep` signé (positif dans le sens antihoraire).
    Arc {
        center: [f64; 2],
        radius: f64,
        theta0: f64,
        sweep: f64,
    },
}

/// Curve vérifiée, paramétrée par t ∈ [0, 1] de `start` à `end`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Piece {
    pub id: CurveId,
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub size: f64,
    path: Path,
}

impl Piece {
    fn new(id: CurveId, curve: &Curve, end: [f64; 2]) -> Result<Self, MailleError> {
        if !(curve.size.is_finite() && curve.size > 0.0) {
            return Err(MailleError::InvalidCurveSize {
                curve: id,
                value: curve.size,
            });
        }
        let (start, size) = (curve.start, curve.size);
        let scale = start[0].hypot(start[1]).max(end[0].hypot(end[1]));
        let path = match curve.kind {
            CurveKind::Line => {
                if distance(start, end) <= 1e-12 * scale {
                    return Err(MailleError::ZeroLengthCurve { curve: id });
                }
                Path::Line
            }
            CurveKind::Arc { center, clockwise } => {
                finite(center, &id, "centre")?;
                let radius = distance(start, center);
                if radius <= 1e-12 * scale || (distance(end, center) - radius).abs() > 1e-9 * radius
                {
                    return Err(MailleError::ArcRadiusMismatch { curve: id });
                }
                let angle = |p: [f64; 2]| (p[1] - center[1]).atan2(p[0] - center[0]);
                let theta0 = angle(start);
                let full = distance(start, end) <= 1e-9 * radius;
                // balayage antihoraire dans ]0, 2π[ hors cercle complet
                let ccw = (angle(end) - theta0).rem_euclid(2.0 * PI);
                let sweep = match (full, clockwise) {
                    (true, false) => 2.0 * PI,
                    (true, true) => -2.0 * PI,
                    (false, false) => ccw,
                    (false, true) => ccw - 2.0 * PI,
                };
                Path::Arc {
                    center,
                    radius,
                    theta0,
                    sweep,
                }
            }
        };
        Ok(Self {
            id,
            start,
            end,
            size,
            path,
        })
    }

    /// Point de paramètre t. Les extrémités sont exactes (pas recalculées sur le cercle).
    pub fn point(&self, t: f64) -> [f64; 2] {
        if t <= 0.0 {
            return self.start;
        }
        if t >= 1.0 {
            return self.end;
        }
        match self.path {
            Path::Line => {
                std::array::from_fn(|j| self.start[j] + t * (self.end[j] - self.start[j]))
            }
            Path::Arc {
                center,
                radius,
                theta0,
                sweep,
            } => {
                let a = theta0 + t * sweep;
                [center[0] + radius * a.cos(), center[1] + radius * a.sin()]
            }
        }
    }

    /// Direction de parcours (non normée) au paramètre t.
    fn tangent(&self, t: f64) -> [f64; 2] {
        match self.path {
            Path::Line => [self.end[0] - self.start[0], self.end[1] - self.start[1]],
            Path::Arc { theta0, sweep, .. } => {
                let a = theta0 + t * sweep;
                [-sweep * a.sin(), sweep * a.cos()]
            }
        }
    }

    fn length(&self) -> f64 {
        match self.path {
            Path::Line => distance(self.start, self.end),
            Path::Arc { radius, sweep, .. } => radius * sweep.abs(),
        }
    }

    /// Nombre de morceaux de longueur au plus h ; au moins 2 pour un arc et 3 pour un cercle
    /// complet, pour que tout contour découpé soit un vrai polygone.
    fn divisions(&self) -> Result<usize, MailleError> {
        let n = (self.length() / self.size).ceil();
        if n > MAX_NODES as f64 {
            return Err(MailleError::TooManyNodes { max: MAX_NODES });
        }
        let min = match self.path {
            Path::Line => 1,
            Path::Arc { sweep, .. } if sweep.abs() > 1.5 * PI => 3,
            Path::Arc { .. } => 2,
        };
        Ok((n as usize).max(min))
    }

    fn reversed(&self) -> Self {
        let path = match self.path {
            Path::Line => Path::Line,
            Path::Arc {
                center,
                radius,
                theta0,
                sweep,
            } => Path::Arc {
                center,
                radius,
                theta0: theta0 + sweep,
                sweep: -sweep,
            },
        };
        Self {
            id: self.id.clone(),
            start: self.end,
            end: self.start,
            size: self.size,
            path,
        }
    }
}

/// Morceau de bord entre les points `a` et `b` de `Boundary::points`, portion [t₀, t₁] de la Piece
/// `piece`. Le domaine est à gauche de a → b.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Segment {
    pub a: usize,
    pub b: usize,
    pub piece: usize,
    pub t: [f64; 2],
}

/// Contours découpés : contour extérieur dans le sens antihoraire, trous dans le sens horaire.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Boundary {
    pub points: Vec<[f64; 2]>,
    pub segments: Vec<Segment>,
    pub pieces: Vec<Piece>,
}

impl Geometry {
    /// Vérifie la géométrie, oriente les contours et découpe chaque Curve en morceaux de
    /// longueur au plus h.
    pub(crate) fn discretize(&self) -> Result<Boundary, MailleError> {
        let loops = std::iter::once((None, &self.outer))
            .chain(self.holes.iter().enumerate().map(|(h, c)| (Some(h), c)));
        let mut b = Boundary::default();
        let mut ranges = Vec::new(); // points de chaque contour : b.points[range]
        for (hole, curves) in loops {
            let mut pieces = loop_pieces(hole, curves)?;
            let ring: Vec<[f64; 2]> = sample(&pieces)?.into_iter().map(|s| s.p).collect();
            if (signed_area(&ring) > 0.0) != hole.is_none() {
                pieces = pieces.iter().rev().map(Piece::reversed).collect();
            }
            let first = b.points.len();
            // les paramètres t changent avec le sens : on redécoupe après une éventuelle inversion
            let samples = sample(&pieces)?;
            if first + samples.len() > MAX_NODES {
                return Err(MailleError::TooManyNodes { max: MAX_NODES });
            }
            for (k, s) in samples.iter().enumerate() {
                // le dernier morceau d'une Piece s'arrête au début de la suivante, en t = 1
                let t1 = match samples.get(k + 1) {
                    Some(next) if next.piece == s.piece => next.t,
                    _ => 1.0,
                };
                b.points.push(s.p);
                b.segments.push(Segment {
                    a: first + k,
                    b: first + (k + 1) % samples.len(),
                    piece: b.pieces.len() + s.piece,
                    t: [s.t, t1],
                });
            }
            b.pieces.extend(pieces.iter().cloned());
            ranges.push((first..b.points.len(), pieces));
        }

        check_intersections(&b)?;
        let outer = &b.points[ranges[0].0.clone()];
        for (h, (range, _)) in ranges.iter().enumerate().skip(1) {
            let p = b.points[range.start];
            let in_hole = ranges
                .iter()
                .enumerate()
                .skip(1)
                .any(|(o, (r, _))| o != h && inside(p, &b.points[r.clone()]));
            if !inside(p, outer) || in_hole {
                return Err(MailleError::HoleOutsideBoundary { hole: h - 1 });
            }
        }
        for (_, pieces) in &ranges {
            check_corners(pieces)?;
        }
        Ok(b)
    }
}

impl Boundary {
    /// Contours polygonaux déjà découpés et orientés (domaine à gauche), pour le remaillage d'une
    /// facette en 3D : chaque côté devient une Piece droite, sans redécoupage ni contrôle des angles.
    /// Chaque sommet porte la taille h du côté qui en part.
    pub(crate) fn polygons(loops: &[Vec<([f64; 2], f64)>], name: &str) -> Self {
        let mut b = Self::default();
        for (l, ring) in loops.iter().enumerate() {
            let first = b.points.len();
            for (index, &(start, size)) in ring.iter().enumerate() {
                let k = (index + 1) % ring.len();
                b.points.push(start);
                b.segments.push(Segment {
                    a: first + index,
                    b: first + k,
                    piece: b.pieces.len(),
                    t: [0.0, 1.0],
                });
                b.pieces.push(Piece {
                    id: CurveId {
                        hole: l.checked_sub(1),
                        index,
                        name: name.to_string(),
                    },
                    start,
                    end: ring[k].0,
                    size,
                    path: Path::Line,
                });
            }
        }
        b
    }
}

/// Pieces d'un contour, chaque Curve finissant au début de la suivante.
fn loop_pieces(hole: Option<usize>, curves: &[Curve]) -> Result<Vec<Piece>, MailleError> {
    if curves.is_empty() {
        return Err(MailleError::EmptyBoundary { hole });
    }
    curves
        .iter()
        .enumerate()
        .map(|(index, c)| {
            let id = CurveId {
                hole,
                index,
                name: c.name.clone(),
            };
            // l'arrivée est le départ de la Curve suivante : vérifiée avec elle (et la dernière
            // revient à la première, vérifiée d'abord)
            finite(c.start, &id, "départ")?;
            Piece::new(id, c, curves[(index + 1) % curves.len()].start)
        })
        .collect()
}

/// Point de découpage : position `p`, de paramètre `t` sur la Piece `piece` de son contour.
struct Sample {
    piece: usize,
    t: f64,
    p: [f64; 2],
}

/// Points de découpage d'un contour (t = 1 exclu : c'est le premier point de la Piece suivante).
fn sample(pieces: &[Piece]) -> Result<Vec<Sample>, MailleError> {
    let mut out = Vec::new();
    for (piece, p) in pieces.iter().enumerate() {
        let n = p.divisions()?;
        for k in 0..n {
            let t = k as f64 / n as f64;
            out.push(Sample {
                piece,
                t,
                p: p.point(t),
            });
        }
    }
    Ok(out)
}

/// Deux morceaux de bord qui ne se suivent pas ne doivent ni se couper ni se toucher.
fn check_intersections(b: &Boundary) -> Result<(), MailleError> {
    let s = &b.segments;
    for i in 0..s.len() {
        for j in i + 1..s.len() {
            let shared = [s[i].a, s[i].b]
                .iter()
                .any(|v| *v == s[j].a || *v == s[j].b);
            let p = |v: usize| b.points[v];
            if !shared && segments_touch(p(s[i].a), p(s[i].b), p(s[j].a), p(s[j].b)) {
                return Err(MailleError::SelfIntersectingBoundary {
                    first: b.pieces[s[i].piece].id.clone(),
                    second: b.pieces[s[j].piece].id.clone(),
                });
            }
        }
    }
    Ok(())
}

/// Angle intérieur (côté domaine, à gauche du parcours) à chaque raccord entre deux Curves.
fn check_corners(pieces: &[Piece]) -> Result<(), MailleError> {
    for (i, before) in pieces.iter().enumerate() {
        let after = &pieces[(i + 1) % pieces.len()];
        let (u, v) = (before.tangent(1.0), after.tangent(0.0));
        let turn = (u[0] * v[1] - u[1] * v[0]).atan2(u[0] * v[0] + u[1] * v[1]);
        // un demi-tour (turn = ±π) est un retour en arrière, d'angle 0 : le signe du zéro du
        // produit vectoriel le ferait sinon passer pour un angle de 360°
        let angle = (PI - turn).to_degrees();
        let angle = if angle > 360.0 - 1e-6 { 0.0 } else { angle };
        if angle < MIN_CORNER {
            return Err(MailleError::SharpCorner {
                first: before.id.clone(),
                second: after.id.clone(),
                angle,
            });
        }
    }
    Ok(())
}

/// `NonFiniteCurve` si une coordonnée du point `what` de la Curve `id` n'est pas finie.
fn finite(p: [f64; 2], id: &CurveId, what: &'static str) -> Result<(), MailleError> {
    if p.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(MailleError::NonFiniteCurve {
            curve: id.clone(),
            what,
        })
    }
}

/// Aire signée d'un polygone (positive dans le sens antihoraire).
fn signed_area(ring: &[[f64; 2]]) -> f64 {
    (0..ring.len())
        .map(|i| {
            let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
            p[0] * q[1] - q[0] * p[1]
        })
        .sum::<f64>()
        / 2.0
}

/// Vrai si p est strictement dans le polygone `ring` (lancer de rayon vers les x croissants).
fn inside(p: [f64; 2], ring: &[[f64; 2]]) -> bool {
    let mut odd = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if x > p[0] {
                odd = !odd;
            }
        }
    }
    odd
}

/// Vrai si les segments pq et rs se coupent ou se touchent (extrémité posée sur l'autre segment,
/// recouvrement de segments alignés).
fn segments_touch(p: [f64; 2], q: [f64; 2], r: [f64; 2], s: [f64; 2]) -> bool {
    let (o1, o2) = (orient2d(p, q, r), orient2d(p, q, s));
    let (o3, o4) = (orient2d(r, s, p), orient2d(r, s, q));
    if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
        return true;
    }
    // point aligné avec un segment et dans sa boîte englobante : posé dessus. Le test == 0.0 est
    // exact : `orient2d` renvoie 0 tout rond sous l'erreur d'arrondi possible
    let on = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        (0..2).all(|j| a[j].min(b[j]) <= c[j] && c[j] <= a[j].max(b[j]))
    };
    (o1 == 0.0 && on(p, q, r))
        || (o2 == 0.0 && on(p, q, s))
        || (o3 == 0.0 && on(r, s, p))
        || (o4 == 0.0 && on(r, s, q))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(size: f64, name: &str) -> Vec<Curve> {
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
            .iter()
            .map(|&p| Curve::line(p, size, name))
            .collect()
    }

    /// Carré [lo, hi]², parcouru dans le sens antihoraire ou horaire.
    fn small_square(lo: f64, hi: f64, ccw: bool) -> Vec<Curve> {
        let mut c = vec![[lo, lo], [hi, lo], [hi, hi], [lo, hi]];
        if !ccw {
            c.reverse();
        }
        c.into_iter().map(|p| Curve::line(p, 1.0, "trou")).collect()
    }

    fn id(hole: Option<usize>, index: usize, name: &str) -> CurveId {
        CurveId {
            hole,
            index,
            name: name.to_string(),
        }
    }

    fn ring(b: &Boundary, range: std::ops::Range<usize>) -> Vec<[f64; 2]> {
        b.points[range].to_vec()
    }

    /// Quart de couronne de Kirsch : a = 1, R = 4, comme le cas du CLI.
    fn kirsch(size_hole: f64, size_out: f64) -> Geometry {
        Geometry {
            outer: vec![
                Curve::line([1.0, 0.0], size_out, "sym_y"),
                Curve::arc([4.0, 0.0], [0.0, 0.0], false, size_out, "exterieur"),
                Curve::line([0.0, 4.0], size_out, "sym_x"),
                Curve::arc([0.0, 1.0], [0.0, 0.0], true, size_hole, "trou"),
            ],
            holes: vec![],
        }
    }

    #[test]
    fn square_is_cut_into_pieces_of_size_h() {
        let g = Geometry {
            outer: square(0.25, "bord"),
            holes: vec![],
        };
        let b = g.discretize().unwrap();
        assert_eq!(b.points.len(), 16);
        assert_eq!(b.segments.len(), 16);
        // chaîne fermée, morceaux de longueur 0,25, contour dans le sens antihoraire
        for (k, s) in b.segments.iter().enumerate() {
            assert_eq!(s.b, b.segments[(k + 1) % 16].a);
            assert!((distance(b.points[s.a], b.points[s.b]) - 0.25).abs() < 1e-12);
        }
        assert!((signed_area(&b.points) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn loops_are_oriented() {
        // contour donné dans le sens horaire, trou dans le sens antihoraire : tous deux retournés
        let mut outer = square(0.5, "bord");
        outer.reverse();
        let g = Geometry {
            outer,
            holes: vec![small_square(0.4, 0.6, true)],
        };
        let b = g.discretize().unwrap();
        assert!((signed_area(&ring(&b, 0..8)) - 1.0).abs() < 1e-12);
        assert!((signed_area(&ring(&b, 8..12)) + 0.04).abs() < 1e-12);
        // le domaine est à gauche de chaque morceau : son milieu décalé à gauche est dans le carré
        // et hors du trou
        for s in &b.segments {
            let (p, q) = (b.points[s.a], b.points[s.b]);
            let m = [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
            let left = [m[0] - 1e-3 * (q[1] - p[1]), m[1] + 1e-3 * (q[0] - p[0])];
            assert!(inside(left, &ring(&b, 0..8)) && !inside(left, &ring(&b, 8..12)));
        }
        // un retournement garde le nom et le rang de la Curve donnés par l'utilisateur
        assert_eq!(b.pieces[0].id, id(None, 3, "bord"));
    }

    #[test]
    fn arcs_are_cut_on_the_circle() {
        let b = kirsch(0.1, 0.5).discretize().unwrap();
        for s in &b.segments {
            let piece = &b.pieces[s.piece];
            let radius = match piece.id.name.as_str() {
                "trou" => 1.0,
                "exterieur" => 4.0,
                _ => continue,
            };
            for p in [
                b.points[s.a],
                b.points[s.b],
                piece.point((s.t[0] + s.t[1]) / 2.0),
            ] {
                assert!((p[0].hypot(p[1]) - radius).abs() < 1e-12);
            }
            assert!(distance(b.points[s.a], b.points[s.b]) <= piece.size);
        }
        // trou : quart de cercle de longueur π/2 en morceaux de 0,1 au plus, soit 16
        let hole = b
            .segments
            .iter()
            .filter(|s| b.pieces[s.piece].id.name == "trou");
        assert_eq!(hole.count(), 16);
        // le paramètre suit le parcours : t₁ du dernier morceau d'une Piece vaut 1
        let s = b.segments.iter().rfind(|s| s.piece == 1).unwrap();
        assert!((s.t[1] - 1.0).abs() < 1e-15);
        assert!(distance(b.points[s.b], [0.0, 4.0]) < 1e-15);
    }

    #[test]
    fn full_circle_hole() {
        for clockwise in [false, true] {
            full_circle(clockwise);
        }
    }

    fn full_circle(clockwise: bool) {
        let g = Geometry {
            outer: square(0.5, "bord"),
            holes: vec![vec![Curve::arc(
                [0.7, 0.5],
                [0.5, 0.5],
                clockwise,
                1.0,
                "trou",
            )]],
        };
        let b = g.discretize().unwrap();
        // h plus grand que le cercle : 3 morceaux au moins, tournés dans le sens horaire
        let hole = ring(&b, 8..b.points.len());
        assert_eq!(hole.len(), 3);
        assert!(signed_area(&hole) < 0.0);
        for p in &hole {
            assert!((distance(*p, [0.5, 0.5]) - 0.2).abs() < 1e-12);
        }
        // 3 points distincts : triangle équilatéral inscrit, d'aire 3√3/4 r²
        assert!((signed_area(&hole) + 3f64.sqrt() * 0.75 * 0.04).abs() < 1e-12);
    }

    #[test]
    fn reflex_corner_is_accepted() {
        // L : angle intérieur de 270° en (1, 1)
        let l = [[0., 0.], [2., 0.], [2., 1.], [1., 1.], [1., 2.], [0., 2.]];
        let g = Geometry {
            outer: l.iter().map(|&p| Curve::line(p, 1.0, "bord")).collect(),
            holes: vec![],
        };
        assert!(g.discretize().is_ok());
    }

    fn err(g: Geometry) -> MailleError {
        g.discretize().unwrap_err()
    }

    #[test]
    fn invalid_geometries() {
        let outer = |curves: Vec<Curve>| Geometry {
            outer: curves,
            holes: vec![],
        };
        assert_eq!(
            err(outer(vec![])),
            MailleError::EmptyBoundary { hole: None }
        );
        assert_eq!(
            err(outer(square(0.0, "bord"))),
            MailleError::InvalidCurveSize {
                curve: id(None, 0, "bord"),
                value: 0.0
            }
        );
        let mut nan = square(1.0, "bord");
        nan[2].start[0] = f64::NAN;
        assert_eq!(
            err(outer(nan)),
            MailleError::NonFiniteCurve {
                curve: id(None, 2, "bord"),
                what: "départ"
            }
        );
        let bad_center = vec![
            Curve::line([0.0, 0.0], 1.0, "a"),
            Curve::arc([1.0, 0.0], [f64::INFINITY, 0.0], false, 1.0, "arc"),
            Curve::line([1.0, 1.0], 1.0, "b"),
        ];
        assert_eq!(
            err(outer(bad_center)),
            MailleError::NonFiniteCurve {
                curve: id(None, 1, "arc"),
                what: "centre"
            }
        );

        let mut doubled = square(1.0, "bord");
        doubled.insert(2, Curve::line([1.0, 0.0], 1.0, "double"));
        assert_eq!(
            err(outer(doubled)),
            MailleError::ZeroLengthCurve {
                curve: id(None, 1, "bord")
            }
        );

        // arc de centre (0, 0) qui part du rayon 1 et arrive au rayon 2
        let bad_arc = vec![
            Curve::line([0.0, 0.0], 1.0, "a"),
            Curve::line([1.0, 0.0], 1.0, "b"),
            Curve::arc([1.0, 1.0], [0.0, 0.0], false, 1.0, "arc"),
        ];
        assert!(matches!(
            err(outer(vec![
                bad_arc[0].clone(),
                bad_arc[2].clone(),
                Curve::line([0.0, 2.0], 1.0, "c")
            ])),
            MailleError::ArcRadiusMismatch { .. }
        ));
        let mut wrong_end = bad_arc;
        wrong_end.push(Curve::line([0.0, 2.0], 1.0, "c"));
        assert_eq!(
            err(outer(wrong_end)),
            MailleError::ArcRadiusMismatch {
                curve: id(None, 2, "arc")
            }
        );

        // nœud papillon : les côtés 1 et 3 se croisent en (½, ½)
        let bowtie = [[0., 0.], [1., 1.], [1., 0.], [0., 1.]];
        assert!(matches!(
            err(outer(
                bowtie.iter().map(|&p| Curve::line(p, 1.0, "b")).collect()
            )),
            MailleError::SelfIntersectingBoundary { .. }
        ));

        // trou qui touche le bord, hors du contour, dans un autre trou
        let with_holes = |holes| Geometry {
            outer: square(1.0, "bord"),
            holes,
        };
        assert!(matches!(
            err(with_holes(vec![small_square(0.0, 0.5, false)])),
            MailleError::SelfIntersectingBoundary { .. }
        ));
        assert_eq!(
            err(with_holes(vec![small_square(2.0, 3.0, false)])),
            MailleError::HoleOutsideBoundary { hole: 0 }
        );
        assert_eq!(
            err(with_holes(vec![
                small_square(0.2, 0.8, false),
                small_square(0.4, 0.6, false)
            ])),
            MailleError::HoleOutsideBoundary { hole: 1 }
        );
    }

    /// Demi-tour au raccord de deux Curves (fissure de largeur nulle). Le produit vectoriel des
    /// tangentes vaut −0 pour (0, 1) puis (0, −1) et +0 pour (1, 0) puis (−1, 0) : les deux
    /// raccords doivent donner un angle de 0°. (Sur un contour complet, la vérification des
    /// croisements refuse déjà ce cas : on teste ici `check_corners` seul.)
    #[test]
    fn folding_back_is_a_zero_angle() {
        for (a, b, c) in [
            ([1.0, 0.0], [1.0, 1.0], [1.0, 0.5]),
            ([0.0, 0.0], [1.0, 0.0], [0.5, 0.0]),
        ] {
            let curves = [a, b, c].map(|p| Curve::line(p, 1.0, "bord"));
            let pieces = loop_pieces(None, &curves).unwrap();
            match check_corners(&pieces[..2]) {
                Err(MailleError::SharpCorner { angle, .. }) => assert!(angle.abs() < 1e-9),
                other => panic!("{other:?}"),
            }
        }
    }

    #[test]
    fn sharp_corner_is_refused() {
        // triangle isocèle d'angles 30°, 75°, 75°
        let a = 30f64.to_radians();
        let t = [[0.0, 0.0], [1.0, 0.0], [a.cos(), a.sin()]];
        let g = Geometry {
            outer: t.iter().map(|&p| Curve::line(p, 1.0, "bord")).collect(),
            holes: vec![],
        };
        match err(g) {
            MailleError::SharpCorner {
                first,
                second,
                angle,
            } => {
                assert_eq!((first.index, second.index), (2, 0));
                assert!((angle - 30.0).abs() < 1e-9);
            }
            e => panic!("{e:?}"),
        }
        // arc tangent à une droite : raccord à 180°, accepté (demi-disque arrondi)
        let g = Geometry {
            outer: vec![
                Curve::line([0.0, 0.0], 1.0, "bas"),
                Curve::line([1.0, 0.0], 1.0, "droite"),
                Curve::arc([1.0, 1.0], [0.5, 1.0], false, 0.1, "haut"),
                Curve::line([0.0, 1.0], 1.0, "gauche"),
            ],
            holes: vec![],
        };
        assert!(g.discretize().is_ok());
    }
}
