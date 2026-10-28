//! Pièce B-rep : le bord exact d'un solide, lu d'un STEP (`maille_io::step::read_brep`). Des
//! sommets, des arêtes portées par des courbes (`EdgeCurve`), et des faces bordées par des contours
//! d'arêtes orientées, chacune sur sa surface (`FaceSurface`). C'est l'entrée du mailleur des faces
//! (`tetrahedralize_brep`). Le nom d'une face devient celui d'un Patch, donc du Set où poser les BC.

use crate::error::MailleError;

/// Pièce B-rep d'un seul solide.
#[derive(Debug, Clone, PartialEq)]
pub struct Brep {
    pub vertices: Vec<[f64; 3]>,
    pub edges: Vec<BrepEdge>,
    pub faces: Vec<BrepFace>,
}

/// Arête de `start` à `end` (indices dans `Brep::vertices`, égaux pour une courbe fermée), le
/// long de `curve`, parcourue de `start` à `end` quand t va de 0 à 1.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepEdge {
    pub start: usize,
    pub end: usize,
    pub curve: EdgeCurve,
}

/// Face de la pièce (ce n'est pas la face d'un Element).
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFace {
    /// Nom de la face dans le fichier, ou `face_1`, `face_2`… s'il est vide.
    pub name: String,
    /// Surface qui porte la face.
    pub surface: FaceSurface,
    /// Contours fermés : indices dans `Brep::edges`, avec `true` pour une arête parcourue de `end`
    /// à `start`. Le premier est le contour extérieur. Vue de dehors, la face est à gauche de ses
    /// contours : ce sont eux qui l'orientent, pas la surface.
    pub loops: Vec<Vec<(usize, bool)>>,
}

/// Surface d'une face (ce n'est pas la `Surface` triangulée d'un STL). Les directions sont
/// unitaires ; son sens ne compte pas, la face est orientée par ses contours.
#[derive(Debug, Clone, PartialEq)]
pub enum FaceSurface {
    /// Plan passant par `origin`, de normale `normal`.
    Plane { origin: [f64; 3], normal: [f64; 3] },
    /// Cylindre de rayon `radius` autour de l'axe qui passe par `origin` dans la direction `axis`.
    Cylinder {
        origin: [f64; 3],
        axis: [f64; 3],
        radius: f64,
    },
    /// Cône d'axe (`origin`, `axis`) : rayon `radius` (positif ou nul) dans le plan de `origin`,
    /// qui croît de tan(`half_angle`) par unité de longueur le long de `axis` ; demi-angle au
    /// sommet dans ]0, π/2[, en radians.
    Cone {
        origin: [f64; 3],
        axis: [f64; 3],
        radius: f64,
        half_angle: f64,
    },
}

/// Courbe d'une arête, paramétrée par t ∈ [0, 1] du début à la fin de l'arête.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeCurve {
    /// Segment de `from` à `to`.
    Line { from: [f64; 3], to: [f64; 3] },
    /// Arc d'ellipse, de cercle si les deux demi-axes sont égaux : centre + a cos θ u + b sin θ v,
    /// θ allant de `start_angle` à `start_angle + sweep` (signé ; ±2π pour une courbe fermée).
    /// `u` et `v` sont unitaires et orthogonaux.
    Arc {
        center: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        radii: [f64; 2],
        start_angle: f64,
        sweep: f64,
    },
    /// Morceau de B-spline, de `range[0]` à `range[1]` dans son paramètre (décroissant si l'arête
    /// la parcourt à rebours).
    BSpline { curve: BSpline, range: [f64; 2] },
}

/// Morceaux d'une longueur calculée par cordes (extrapolation de Richardson entre N et 2N cordes).
const LENGTH_CHORDS: usize = 2048;

impl EdgeCurve {
    /// Point de paramètre t ∈ [0, 1].
    pub fn point(&self, t: f64) -> [f64; 3] {
        match self {
            Self::Line { from, to } => std::array::from_fn(|k| from[k] + t * (to[k] - from[k])),
            Self::Arc {
                center,
                u,
                v,
                radii,
                start_angle,
                sweep,
            } => {
                let a = start_angle + t * sweep;
                let (c, s) = (radii[0] * a.cos(), radii[1] * a.sin());
                std::array::from_fn(|k| center[k] + c * u[k] + s * v[k])
            }
            Self::BSpline { curve, range } => curve.at(range[0] + t * (range[1] - range[0])),
        }
    }

    /// Longueur : exacte pour un segment et un arc de cercle, par cordes extrapolées sinon (écart
    /// relatif de l'ordre de 1e-12 sur une courbe lisse).
    pub fn length(&self) -> f64 {
        match self {
            Self::Line { from, to } => distance(*from, *to),
            Self::Arc { radii, sweep, .. } if radii[0] == radii[1] => radii[0] * sweep.abs(),
            _ => {
                let (n, m) = (self.chords(LENGTH_CHORDS), self.chords(2 * LENGTH_CHORDS));
                (4.0 * m - n) / 3.0
            }
        }
    }

    fn chords(&self, n: usize) -> f64 {
        (0..n)
            .map(|i| {
                distance(
                    self.point(i as f64 / n as f64),
                    self.point((i + 1) as f64 / n as f64),
                )
            })
            .sum()
    }
}

/// B-spline de degré p, rationnelle si elle a des poids : C(u) = Σ wᵢ Nᵢ,ₚ(u) Pᵢ / Σ wᵢ Nᵢ,ₚ(u),
/// pour u dans [nœud p, nœud n] (n points de contrôle, n + p + 1 nœuds, répétés).
#[derive(Debug, Clone, PartialEq)]
pub struct BSpline {
    degree: usize,
    points: Vec<[f64; 3]>,
    weights: Vec<f64>,
    knots: Vec<f64>,
}

impl BSpline {
    /// `weights` vide pour une B-spline non rationnelle, sinon un poids strictement positif par
    /// point. `InvalidBrep` si les tailles ne vont pas ensemble, si les nœuds décroissent ou si une
    /// valeur n'est pas finie.
    pub fn new(
        degree: usize,
        points: Vec<[f64; 3]>,
        weights: Vec<f64>,
        knots: Vec<f64>,
    ) -> Result<Self, MailleError> {
        let invalid = |problem: String| MailleError::InvalidBrep { problem };
        let n = points.len();
        if degree == 0 || n <= degree {
            return Err(invalid(format!(
                "B-spline de degré {degree} à {n} points de contrôle : il en faut au moins degré + 1"
            )));
        }
        if knots.len() != n + degree + 1 {
            return Err(invalid(format!(
                "B-spline à {n} points de degré {degree} : {} nœuds au lieu de {}",
                knots.len(),
                n + degree + 1
            )));
        }
        if !(weights.is_empty() || weights.len() == n) {
            return Err(invalid(format!(
                "B-spline rationnelle : {} poids pour {n} points",
                weights.len()
            )));
        }
        if points
            .iter()
            .flatten()
            .chain(&knots)
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "B-spline : coordonnée ou nœud non fini".to_string(),
            ));
        }
        if weights.iter().any(|w| !(w.is_finite() && *w > 0.0)) {
            return Err(invalid(
                "B-spline rationnelle : un poids n'est pas strictement positif".to_string(),
            ));
        }
        if knots.windows(2).any(|w| w[1] < w[0]) || knots[degree] >= knots[n] {
            return Err(invalid(
                "B-spline : nœuds décroissants ou intervalle de paramètre vide".to_string(),
            ));
        }
        Ok(Self {
            degree,
            points,
            weights,
            knots,
        })
    }

    /// Intervalle du paramètre.
    pub fn domain(&self) -> [f64; 2] {
        [self.knots[self.degree], self.knots[self.points.len()]]
    }

    /// Point de paramètre u (ramené dans le domaine), par l'algorithme de de Boor en coordonnées
    /// homogènes.
    pub fn at(&self, u: f64) -> [f64; 3] {
        let (p, n) = (self.degree, self.points.len());
        let [lo, hi] = self.domain();
        let u = u.clamp(lo, hi);
        // intervalle k : nœud k ≤ u < nœud k + 1, le dernier non vide en u = hi
        let mut k = p;
        while k + 1 < n && self.knots[k + 1] <= u {
            k += 1;
        }
        let mut d: Vec<[f64; 4]> = (0..=p)
            .map(|j| {
                let i = k - p + j;
                let w = self.weights.get(i).copied().unwrap_or(1.0);
                let q = self.points[i];
                [w * q[0], w * q[1], w * q[2], w]
            })
            .collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = k - p + j;
                let span = self.knots[i + p + 1 - r] - self.knots[i];
                let a = if span > 0.0 {
                    (u - self.knots[i]) / span
                } else {
                    0.0
                };
                d[j] = std::array::from_fn(|c| (1.0 - a) * d[j - 1][c] + a * d[j][c]);
            }
        }
        let h = d[p];
        [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
    }

    /// Paramètre du point de la courbe le plus proche de `q`, et sa distance : échantillonnage
    /// puis affinage par section dorée autour du meilleur échantillon.
    pub fn parameter_of(&self, q: [f64; 3]) -> (f64, f64) {
        let [lo, hi] = self.domain();
        let samples = 64 * self.points.len();
        let at = |i: usize| lo + (hi - lo) * i as f64 / samples as f64;
        let best = (0..=samples)
            .min_by(|&i, &j| distance(self.at(at(i)), q).total_cmp(&distance(self.at(at(j)), q)))
            .unwrap_or(0);
        let (mut a, mut b) = (at(best.saturating_sub(1)), at((best + 1).min(samples)));
        let g = (5f64.sqrt() - 1.0) / 2.0;
        for _ in 0..100 {
            let (x, y) = (b - g * (b - a), a + g * (b - a));
            if distance(self.at(x), q) < distance(self.at(y), q) {
                b = y;
            } else {
                a = x;
            }
        }
        let u = (a + b) / 2.0;
        (u, distance(self.at(u), q))
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, PI};

    #[test]
    fn segment_et_arc_de_cercle() {
        let line = EdgeCurve::Line {
            from: [1.0, 2.0, 3.0],
            to: [4.0, 6.0, 3.0],
        };
        assert_eq!(line.length(), 5.0);
        assert_eq!(line.point(0.5), [2.5, 4.0, 3.0]);
        // trois quarts de cercle de rayon 2 autour de z, parcourus à rebours
        let arc = EdgeCurve::Arc {
            center: [0.0; 3],
            u: [1.0, 0.0, 0.0],
            v: [0.0, 1.0, 0.0],
            radii: [2.0, 2.0],
            start_angle: PI / 2.0,
            sweep: -1.5 * PI,
        };
        assert!((arc.length() - 3.0 * PI).abs() < 1e-12);
        let end = arc.point(1.0);
        assert!((end[0] + 2.0).abs() < 1e-12 && end[1].abs() < 1e-12);
    }

    #[test]
    fn perimetre_d_ellipse() {
        // demi-axes 2 et 1 : périmètre 9,688448220547675… (intégrale elliptique complète)
        let e = EdgeCurve::Arc {
            center: [0.0; 3],
            u: [1.0, 0.0, 0.0],
            v: [0.0, 0.0, 1.0],
            radii: [2.0, 1.0],
            start_angle: 0.0,
            sweep: 2.0 * PI,
        };
        assert!((e.length() - 9.688_448_220_547_675).abs() < 1e-9);
    }

    /// Quart de cercle exact : B-spline rationnelle quadratique, poids (1, √2/2, 1).
    fn quarter_circle(r: f64) -> BSpline {
        BSpline::new(
            2,
            vec![[r, 0.0, 0.0], [r, r, 0.0], [0.0, r, 0.0]],
            vec![1.0, FRAC_1_SQRT_2, 1.0],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        )
        .unwrap()
    }

    #[test]
    fn b_spline_rationnelle_exacte() {
        let c = quarter_circle(3.0);
        for i in 0..=20 {
            let p = c.at(i as f64 / 20.0);
            assert!((p[0].hypot(p[1]) - 3.0).abs() < 1e-14);
        }
        let edge = EdgeCurve::BSpline {
            curve: c.clone(),
            range: c.domain(),
        };
        assert!((edge.length() - 1.5 * PI).abs() < 1e-10);
        // la même parcourue à rebours, et un morceau
        let back = EdgeCurve::BSpline {
            curve: c.clone(),
            range: [1.0, 0.0],
        };
        assert_eq!(back.point(0.0), c.at(1.0));
        // sommet sur la courbe (le cas d'une arête rognée) : paramètre précis
        let (u, d) = c.parameter_of([3.0 * FRAC_1_SQRT_2, 3.0 * FRAC_1_SQRT_2, 0.0]);
        assert!((u - 0.5).abs() < 1e-12 && d < 1e-12, "{u} {d}");
        // point hors de la courbe : la distance est plate au minimum, u connu à √ε près
        let (u, d) = c.parameter_of([3.0 * FRAC_1_SQRT_2, 3.0 * FRAC_1_SQRT_2, 0.5]);
        assert!((u - 0.5).abs() < 1e-7 && (d - 0.5).abs() < 1e-12, "{u} {d}");
    }

    #[test]
    fn b_spline_cubique_a_noeuds_interieurs() {
        // cubique non rationnelle dont les points sont alignés à pas égal : la droite, parcourue
        // à vitesse constante (nœuds uniformes à l'intérieur)
        let points: Vec<[f64; 3]> = (0..6).map(|i| [i as f64, 2.0 * i as f64, 0.0]).collect();
        let c = BSpline::new(
            3,
            points,
            vec![],
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0],
        )
        .unwrap();
        assert_eq!(c.domain(), [0.0, 3.0]);
        let edge = EdgeCurve::BSpline {
            curve: c.clone(),
            range: c.domain(),
        };
        assert!((edge.length() - 5.0 * 5f64.sqrt()).abs() < 1e-10);
        let end = c.at(3.0);
        assert!((end[0] - 5.0).abs() < 1e-14 && (end[1] - 10.0).abs() < 1e-14);
    }

    #[test]
    fn b_splines_invalides_refusees() {
        let p = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let k = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let bad = [
            BSpline::new(2, p[..2].to_vec(), vec![], k.clone()),
            BSpline::new(2, p.clone(), vec![], k[..5].to_vec()),
            BSpline::new(2, p.clone(), vec![1.0], k.clone()),
            BSpline::new(2, p.clone(), vec![1.0, -1.0, 1.0], k.clone()),
            BSpline::new(2, p.clone(), vec![], vec![0.0, 0.0, 1.0, 0.5, 1.0, 1.0]),
            BSpline::new(2, p, vec![], vec![0.0; 6]),
        ];
        for b in bad {
            assert!(matches!(b, Err(MailleError::InvalidBrep { .. })), "{b:?}");
        }
    }
}
