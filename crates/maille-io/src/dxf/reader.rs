//! Lecture du DXF ASCII : paires (code, valeur), section ENTITIES, une `Piece` par segment ou arc.

use std::f64::consts::TAU;

use super::error;
use crate::IoError;

/// Morceau de dessin de `a` à `b` : un segment, ou un arc de cercle (`a == b` pour un cercle
/// complet) parcouru dans le sens antihoraire si `ccw`.
#[derive(Debug, Clone)]
pub(super) struct Piece {
    pub layer: String,
    pub a: [f64; 2],
    pub b: [f64; 2],
    pub arc: Option<Arc>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Arc {
    pub center: [f64; 2],
    pub ccw: bool,
}

struct Group<'a> {
    line: usize,
    code: i32,
    value: &'a str,
}

fn groups(text: &str) -> Result<Vec<Group<'_>>, IoError> {
    if text.starts_with("AutoCAD Binary DXF") {
        return Err(error(
            "fichier DXF binaire : Maille lit le DXF en ASCII, ré-exporter sans l'option « binaire »"
                .to_string(),
        ));
    }
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::new();
    for (i, pair) in lines.chunks(2).enumerate() {
        let line = 2 * i + 1;
        let code = pair[0].trim().parse::<i32>().map_err(|_| {
            error(format!(
                "ligne {line} : code de groupe attendu, trouvé « {} » : ce n'est pas un DXF ASCII",
                pair[0].trim()
            ))
        })?;
        let Some(value) = pair.get(1) else {
            return Err(error(format!(
                "ligne {line} : fichier tronqué après le code {code}"
            )));
        };
        if code != 999 {
            out.push(Group {
                line,
                code,
                value: value.trim(),
            });
        }
    }
    Ok(out)
}

fn number(g: &Group) -> Result<f64, IoError> {
    match g.value.parse::<f64>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(error(format!(
            "ligne {} : nombre attendu pour le code {}, trouvé « {} »",
            g.line, g.code, g.value
        ))),
    }
}

/// Première valeur du code `code` d'une entité, obligatoire.
fn required(entity: &[Group], kind: &str, code: i32) -> Result<f64, IoError> {
    match entity.iter().find(|g| g.code == code) {
        Some(g) => number(g),
        None => Err(error(format!(
            "ligne {} : {kind} sans valeur de code {code} (entité incomplète)",
            entity[0].line
        ))),
    }
}

pub(super) fn read_pieces(text: &str) -> Result<Vec<Piece>, IoError> {
    let all = groups(text)?;
    let is = |g: &Group, code: i32, value: &str| g.code == code && g.value == value;
    let start = (0..all.len().saturating_sub(1))
        .find(|&i| is(&all[i], 0, "SECTION") && is(&all[i + 1], 2, "ENTITIES"))
        .ok_or_else(|| {
            error("section ENTITIES absente : ce fichier ne contient pas de dessin".to_string())
        })?;
    let body = &all[start + 2..];
    let end = body
        .iter()
        .position(|g| is(g, 0, "ENDSEC"))
        .unwrap_or(body.len());
    let body = &body[..end];

    let mut pieces = Vec::new();
    let mut from = 0;
    while from < body.len() {
        let to = body[from + 1..]
            .iter()
            .position(|g| g.code == 0)
            .map_or(body.len(), |k| from + 1 + k);
        let entity = &body[from..to];
        if entity[0].code != 0 {
            return Err(error(format!(
                "ligne {} : entité sans type (code 0 attendu)",
                entity[0].line
            )));
        }
        entity_pieces(entity, &mut pieces)?;
        from = to;
    }
    if pieces.is_empty() {
        return Err(error(
            "aucune entité dessinable : dessiner des LINE, ARC, CIRCLE ou LWPOLYLINE".to_string(),
        ));
    }
    Ok(pieces)
}

fn entity_pieces(entity: &[Group], out: &mut Vec<Piece>) -> Result<(), IoError> {
    let kind = entity[0].value;
    let layer = entity
        .iter()
        .find(|g| g.code == 8)
        .map_or("", |g| g.value)
        .to_string();
    let point = |cx: i32, cy: i32| -> Result<[f64; 2], IoError> {
        Ok([required(entity, kind, cx)?, required(entity, kind, cy)?])
    };
    match kind {
        "LINE" => out.push(Piece { layer, a: point(10, 20)?, b: point(11, 21)?, arc: None }),
        "CIRCLE" => {
            let (c, r) = (point(10, 20)?, required(entity, kind, 40)?);
            let a = [c[0] + r, c[1]];
            out.push(Piece { layer, a, b: a, arc: Some(Arc { center: c, ccw: true }) });
        }
        "ARC" => {
            let (c, r) = (point(10, 20)?, required(entity, kind, 40)?);
            let (s, e) = (required(entity, kind, 50)?.to_radians(), required(entity, kind, 51)?.to_radians());
            let at = |t: f64| [c[0] + r * t.cos(), c[1] + r * t.sin()];
            out.push(Piece { layer, a: at(s), b: at(e), arc: Some(Arc { center: c, ccw: true }) });
        }
        "LWPOLYLINE" => polyline_pieces(entity, &layer, out)?,
        other => {
            return Err(error(format!(
                "ligne {} : entité {other} (calque « {layer} ») non prise en charge ; Maille lit LINE, ARC, CIRCLE et LWPOLYLINE : supprimer cette entité ou la remplacer par l'un de ces types",
                entity[0].line
            )))
        }
    }
    Ok(())
}

fn polyline_pieces(entity: &[Group], layer: &str, out: &mut Vec<Piece>) -> Result<(), IoError> {
    let mut closed = false;
    // (point, arrondi du côté qui part de ce sommet)
    let mut vertices: Vec<([f64; 2], f64)> = Vec::new();
    for g in entity {
        match g.code {
            70 => closed = number(g)? as i64 & 1 == 1,
            10 => vertices.push(([number(g)?, 0.0], 0.0)),
            20 | 42 => match vertices.last_mut() {
                Some(v) if g.code == 20 => v.0[1] = number(g)?,
                Some(v) => v.1 = number(g)?,
                None => {
                    return Err(error(format!(
                        "ligne {} : LWPOLYLINE : coordonnée avant le premier sommet",
                        g.line
                    )))
                }
            },
            _ => {}
        }
    }
    let n = vertices.len();
    if n < 2 {
        return Err(error(format!(
            "ligne {} : LWPOLYLINE avec moins de 2 sommets",
            entity[0].line
        )));
    }
    let sides = if closed { n } else { n - 1 };
    for i in 0..sides {
        let ((a, bulge), (b, _)) = (vertices[i], vertices[(i + 1) % n]);
        out.push(Piece {
            layer: layer.to_string(),
            a,
            b,
            arc: bulge_arc(a, b, bulge),
        });
    }
    Ok(())
}

/// Arc d'arrondi `bulge` = tan(θ/4) de `a` à `b` : positif, il tourne dans le sens antihoraire.
fn bulge_arc(a: [f64; 2], b: [f64; 2], bulge: f64) -> Option<Arc> {
    if bulge == 0.0 {
        return None;
    }
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let chord = dx.hypot(dy);
    // le centre est sur la médiatrice, à gauche de la corde pour un arc antihoraire
    let offset = chord * (1.0 - bulge * bulge) / (4.0 * bulge);
    let center = [
        (a[0] + b[0]) / 2.0 - dy / chord * offset,
        (a[1] + b[1]) / 2.0 + dx / chord * offset,
    ];
    Some(Arc {
        center,
        ccw: bulge > 0.0,
    })
}

/// Angle de rotation de `from` à `to` dans le sens demandé, dans ]0, 2π] (2π si les angles
/// coïncident : cercle complet).
pub(super) fn sweep(from: f64, to: f64, ccw: bool) -> f64 {
    let d = if ccw { to - from } else { from - to };
    let d = d.rem_euclid(TAU);
    let d = if d < 1e-12 { TAU } else { d };
    if ccw {
        d
    } else {
        -d
    }
}
