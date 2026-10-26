//! Chaînage des pièces d'un dessin en contours fermés, puis `Geometry`.

use super::error;
use super::reader::{sweep, Arc, Piece};
use crate::IoError;
use maille_core::mesher::{Curve, Geometry};

/// Nom du bord pour une entité sur le calque « 0 » ou sans calque.
const DEFAULT_NAME: &str = "bord";

fn dist(p: [f64; 2], q: [f64; 2]) -> f64 {
    (p[0] - q[0]).hypot(p[1] - q[1])
}

fn reverse(p: &Piece) -> Piece {
    Piece {
        layer: p.layer.clone(),
        a: p.b,
        b: p.a,
        arc: p.arc.map(|a| Arc { ccw: !a.ccw, ..a }),
    }
}

/// Tolérance de raccord : un millionième de la diagonale du dessin.
fn tolerance(pieces: &[Piece]) -> f64 {
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for p in pieces {
        let mut points = vec![p.a, p.b];
        if let Some(arc) = p.arc {
            let r = dist(arc.center, p.a);
            points.push([arc.center[0] - r, arc.center[1] - r]);
            points.push([arc.center[0] + r, arc.center[1] + r]);
        }
        for q in points {
            for k in 0..2 {
                lo[k] = lo[k].min(q[k]);
                hi[k] = hi[k].max(q[k]);
            }
        }
    }
    1e-6 * dist(lo, hi).max(f64::MIN_POSITIVE)
}

/// Pièces d'un contour fermé, orientées bout à bout.
fn chain(pieces: &[Piece], tol: f64) -> Result<Vec<Vec<Piece>>, IoError> {
    let mut used = vec![false; pieces.len()];
    let mut loops = Vec::new();
    for first in 0..pieces.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut contour = vec![pieces[first].clone()];
        loop {
            let (start, end) = (contour[0].a, contour[contour.len() - 1].b);
            if dist(start, end) < tol {
                break;
            }
            let candidates: Vec<usize> = (0..pieces.len())
                .filter(|&k| {
                    !used[k] && (dist(end, pieces[k].a) < tol || dist(end, pieces[k].b) < tol)
                })
                .collect();
            match candidates[..] {
                [k] => {
                    used[k] = true;
                    let next = if dist(end, pieces[k].a) < tol { pieces[k].clone() } else { reverse(&pieces[k]) };
                    contour.push(next);
                }
                [] => {
                    return Err(error(format!(
                        "contour ouvert en ({}, {}) : aucune entité ne continue le tracé à moins de {tol:.1e}. Joindre les extrémités dans le logiciel de dessin",
                        end[0], end[1]
                    )))
                }
                _ => {
                    return Err(error(format!(
                        "{} entités se rejoignent en ({}, {}) : un contour ne peut pas se ramifier ; séparer les contours",
                        candidates.len() + 1,
                        end[0],
                        end[1]
                    )))
                }
            }
        }
        loops.push(contour);
    }
    Ok(loops)
}

/// Aire signée (positive si antihoraire) du contour, calculée exactement avec les arcs.
fn area(contour: &[Piece]) -> f64 {
    let twice: f64 = contour
        .iter()
        .map(|p| match p.arc {
            None => p.a[0] * p.b[1] - p.b[0] * p.a[1],
            Some(Arc { center: c, ccw }) => {
                let r = dist(c, p.a);
                let t1 = (p.a[1] - c[1]).atan2(p.a[0] - c[0]);
                let t2 = (p.b[1] - c[1]).atan2(p.b[0] - c[0]);
                r * r * sweep(t1, t2, ccw) + r * c[0] * (t2.sin() - t1.sin())
                    - r * c[1] * (t2.cos() - t1.cos())
            }
        })
        .sum();
    twice / 2.0
}

fn curves(contour: &[Piece], size: f64) -> Vec<Curve> {
    let n = contour.len();
    (0..n)
        .map(|i| {
            let (prev, p) = (&contour[(i + n - 1) % n], &contour[i]);
            // le raccord est calculé sur l'arc quand il y en a un : le mailleur exige que
            // l'arrivée d'un arc soit sur son cercle
            let start = if prev.arc.is_some() { prev.b } else { p.a };
            let name = match p.layer.as_str() {
                "" | "0" => DEFAULT_NAME,
                layer => layer,
            };
            match p.arc {
                None => Curve::line(start, size, name),
                Some(arc) => Curve::arc(start, arc.center, !arc.ccw, size, name),
            }
        })
        .collect()
}

pub(super) fn build(pieces: Vec<Piece>, size: f64) -> Result<Geometry, IoError> {
    let tol = tolerance(&pieces);
    let mut loops = chain(&pieces, tol)?;
    let outer = (0..loops.len())
        .max_by(|&i, &j| area(&loops[i]).abs().total_cmp(&area(&loops[j]).abs()))
        .ok_or_else(|| error("aucun contour fermé dans le dessin".to_string()))?;
    let outer = curves(&loops.swap_remove(outer), size);
    let holes = loops.iter().map(|l| curves(l, size)).collect();
    Ok(Geometry { outer, holes })
}
