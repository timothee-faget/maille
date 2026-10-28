//! Du contenu Part 21 à la `Brep` : le solide unique (MANIFOLD_SOLID_BREP), sa coque fermée, ses
//! faces (ADVANCED_FACE), leurs contours (EDGE_LOOP d'ORIENTED_EDGE), les arêtes (EDGE_CURVE) et
//! leurs courbes (LINE, CIRCLE, ELLIPSE, B_SPLINE_CURVE_WITH_KNOTS, rationnelle ou non), puis
//! l'unité de longueur du fichier.

use std::collections::HashMap;
use std::f64::consts::TAU;

use maille_core::mesher::{BSpline, Brep, BrepEdge, BrepFace, EdgeCurve};

use super::error;
use super::part21::{Instance, Part21, Record, Value};
use crate::IoError;

/// Écart toléré entre un sommet et la courbe de son arête, rapporté à la diagonale de la pièce.
const ON_CURVE: f64 = 1e-6;

pub(super) fn build(p: &Part21) -> Result<(Brep, String), IoError> {
    let solids: Vec<(&u64, &Instance)> = p
        .instances
        .iter()
        .filter(|(_, i)| {
            i.record("MANIFOLD_SOLID_BREP").is_some() || i.record("BREP_WITH_VOIDS").is_some()
        })
        .collect();
    let (&id, solid) = match solids[..] {
        [one] => one,
        [] => {
            return Err(error(
                "aucun solide (MANIFOLD_SOLID_BREP) : le fichier décrit une surface ou un fil de fer, qui ne se maille pas en volume ; exporter la pièce en solide".to_string(),
            ))
        }
        _ => {
            let lines: Vec<String> = solids.iter().map(|(_, i)| i.line.to_string()).collect();
            return Err(error(format!(
                "{} solides (lignes {}) : Maille maille une seule pièce ; exporter chaque pièce dans son fichier",
                solids.len(),
                lines.join(", ")
            )));
        }
    };
    if solid.record("BREP_WITH_VOIDS").is_some() {
        return Err(error(format!(
            "ligne {} : solide à cavités internes (BREP_WITH_VOIDS), pas encore pris en charge",
            solid.line
        )));
    }
    let mut r = Reader {
        p,
        vertices: Vec::new(),
        vertex_of: HashMap::new(),
        edges: Vec::new(),
        edge_of: HashMap::new(),
        tolerance: 0.0,
    };
    r.tolerance = ON_CURVE * r.diagonal();
    let record = r.only(id, solid.line, &["MANIFOLD_SOLID_BREP"])?;
    let shell_id = r.reference(&record.params, 1, solid.line, "coque du solide")?;
    let shell = r.only(shell_id, solid.line, &["CLOSED_SHELL", "OPEN_SHELL"])?;
    let shell_line = p.get(shell_id, solid.line)?.line;
    if shell.name == "OPEN_SHELL" {
        return Err(error(format!(
            "ligne {shell_line} : coque ouverte (OPEN_SHELL) : la pièce n'est pas un volume fermé"
        )));
    }
    let mut faces = Vec::new();
    for (k, face_id) in r
        .references(&shell.params, 1, shell_line, "faces de la coque")?
        .into_iter()
        .enumerate()
    {
        faces.push(r.face(face_id, shell_line, k)?);
    }
    let brep = Brep {
        vertices: r.vertices,
        edges: r.edges,
        faces,
    };
    check_closed(&brep, &r.edge_of)?;
    Ok((brep, length_unit(p)))
}

struct Reader<'a> {
    p: &'a Part21,
    vertices: Vec<[f64; 3]>,
    vertex_of: HashMap<u64, usize>,
    edges: Vec<BrepEdge>,
    /// Arête de la `Brep` de chaque EDGE_CURVE, et son numéro dans le fichier.
    edge_of: HashMap<u64, usize>,
    tolerance: f64,
}

impl<'a> Reader<'a> {
    /// Diagonale de la boîte des points cartésiens 3D du fichier.
    fn diagonal(&self) -> f64 {
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for i in self.p.instances.values() {
            if let Some(Record { params, .. }) = i.record("CARTESIAN_POINT") {
                if let Some(Value::List(c)) = params.get(1) {
                    for (k, v) in c.iter().take(3).enumerate() {
                        if let Some(x) = number(v) {
                            lo[k] = lo[k].min(x);
                            hi[k] = hi[k].max(x);
                        }
                    }
                }
            }
        }
        let d = (0..3)
            .map(|k| (hi[k] - lo[k]).max(0.0).powi(2))
            .sum::<f64>()
            .sqrt();
        if d.is_finite() && d > 0.0 {
            d
        } else {
            1.0
        }
    }

    /// Enregistrement de l'instance simple `id`, dont le nom doit être l'un de `names`.
    fn only(&self, id: u64, from: usize, names: &[&str]) -> Result<&'a Record, IoError> {
        let i = self.p.get(id, from)?;
        i.records
            .iter()
            .find(|r| names.contains(&r.name.as_str()))
            .ok_or_else(|| {
                let found: Vec<&str> = i.records.iter().map(|r| r.name.as_str()).collect();
                error(format!(
                    "ligne {} : #{id} est {} ; attendu : {}",
                    i.line,
                    found.join(" + "),
                    names.join(" ou ")
                ))
            })
    }

    fn reference(
        &self,
        params: &[Value],
        k: usize,
        line: usize,
        what: &str,
    ) -> Result<u64, IoError> {
        match params.get(k) {
            Some(Value::Ref(id)) => Ok(*id),
            _ => Err(error(format!(
                "ligne {line} : {what} : référence « #n » attendue"
            ))),
        }
    }

    fn references(
        &self,
        params: &[Value],
        k: usize,
        line: usize,
        what: &str,
    ) -> Result<Vec<u64>, IoError> {
        match params.get(k) {
            Some(Value::List(l)) => l
                .iter()
                .map(|v| match v {
                    Value::Ref(id) => Ok(*id),
                    _ => Err(error(format!(
                        "ligne {line} : {what} : liste de « #n » attendue"
                    ))),
                })
                .collect(),
            _ => Err(error(format!("ligne {line} : {what} : liste attendue"))),
        }
    }

    fn face(&mut self, id: u64, from: usize, k: usize) -> Result<BrepFace, IoError> {
        let line = self.p.get(id, from)?.line;
        let face = self.only(id, from, &["ADVANCED_FACE", "FACE_SURFACE"])?;
        let name = match face.params.first() {
            Some(Value::String(s)) if !s.trim().is_empty() => s.trim().to_string(),
            _ => format!("face_{}", k + 1),
        };
        let mut outer = Vec::new();
        let mut inner = Vec::new();
        for bound_id in self.references(&face.params, 1, line, "contours de la face")? {
            let bound_line = self.p.get(bound_id, line)?.line;
            let bound = self.only(bound_id, line, &["FACE_OUTER_BOUND", "FACE_BOUND"])?;
            let loop_id = self.reference(&bound.params, 1, bound_line, "contour")?;
            let same = boolean(bound.params.get(2), bound_line, "sens du contour")?;
            let lp = self.only(
                loop_id,
                bound_line,
                &["EDGE_LOOP", "VERTEX_LOOP", "POLY_LOOP"],
            )?;
            let loop_line = self.p.get(loop_id, bound_line)?.line;
            match lp.name.as_str() {
                // contour réduit à un sommet (pôle d'une sphère) : pas d'arête
                "VERTEX_LOOP" => continue,
                "POLY_LOOP" => {
                    return Err(error(format!(
                    "ligne {loop_line} : contour polygonal (POLY_LOOP) pas encore pris en charge"
                )))
                }
                _ => {}
            }
            let mut edges = Vec::new();
            for oriented in self.references(&lp.params, 1, loop_line, "arêtes du contour")? {
                let o_line = self.p.get(oriented, loop_line)?.line;
                let o = self.only(oriented, loop_line, &["ORIENTED_EDGE"])?;
                let edge_id = self.reference(&o.params, 3, o_line, "arête orientée")?;
                let forward = boolean(o.params.get(4), o_line, "sens de l'arête")?;
                edges.push((self.edge(edge_id, o_line)?, !forward));
            }
            if !same {
                edges.reverse();
                for e in &mut edges {
                    e.1 = !e.1;
                }
            }
            if bound.name == "FACE_OUTER_BOUND" {
                outer.push(edges);
            } else {
                inner.push(edges);
            }
        }
        outer.extend(inner);
        Ok(BrepFace { name, loops: outer })
    }

    fn vertex(&mut self, id: u64, from: usize) -> Result<usize, IoError> {
        if let Some(&v) = self.vertex_of.get(&id) {
            return Ok(v);
        }
        let line = self.p.get(id, from)?.line;
        let vp = self.only(id, from, &["VERTEX_POINT"])?;
        let point = self.reference(&vp.params, 1, line, "point du sommet")?;
        let p = self.point(point, line)?;
        self.vertices.push(p);
        self.vertex_of.insert(id, self.vertices.len() - 1);
        Ok(self.vertices.len() - 1)
    }

    fn point(&self, id: u64, from: usize) -> Result<[f64; 3], IoError> {
        let line = self.p.get(id, from)?.line;
        let r = self.only(id, from, &["CARTESIAN_POINT"])?;
        self.triple(r.params.get(1), line, "coordonnées d'un point")
    }

    fn direction(&self, id: u64, from: usize) -> Result<[f64; 3], IoError> {
        let line = self.p.get(id, from)?.line;
        let r = self.only(id, from, &["DIRECTION"])?;
        let d = self.triple(r.params.get(1), line, "direction")?;
        let n = norm(d);
        if n == 0.0 {
            return Err(error(format!("ligne {line} : direction nulle")));
        }
        Ok(d.map(|x| x / n))
    }

    fn triple(&self, v: Option<&Value>, line: usize, what: &str) -> Result<[f64; 3], IoError> {
        let bad = || {
            error(format!(
                "ligne {line} : {what} : trois nombres attendus (pièce 3D)"
            ))
        };
        match v {
            Some(Value::List(l)) if l.len() == 3 => {
                let c: Vec<f64> = l
                    .iter()
                    .map(number)
                    .collect::<Option<_>>()
                    .ok_or_else(bad)?;
                Ok([c[0], c[1], c[2]])
            }
            _ => Err(bad()),
        }
    }

    /// Repère (origine, u, v, axe) d'un AXIS2_PLACEMENT_3D ; axe z et direction x par défaut.
    fn placement(&self, id: u64, from: usize) -> Result<[[f64; 3]; 4], IoError> {
        let line = self.p.get(id, from)?.line;
        let r = self.only(id, from, &["AXIS2_PLACEMENT_3D"])?;
        let origin = self.point(
            self.reference(&r.params, 1, line, "origine du repère")?,
            line,
        )?;
        let axis = match r.params.get(2) {
            Some(Value::Ref(d)) => self.direction(*d, line)?,
            _ => [0.0, 0.0, 1.0],
        };
        let reference = match r.params.get(3) {
            Some(Value::Ref(d)) => self.direction(*d, line)?,
            _ if axis[0].abs() < 0.9 => [1.0, 0.0, 0.0],
            _ => [0.0, 1.0, 0.0],
        };
        let d = dot(reference, axis);
        let u: [f64; 3] = std::array::from_fn(|k| reference[k] - d * axis[k]);
        let n = norm(u);
        if n < 1e-12 {
            return Err(error(format!(
                "ligne {line} : direction de référence parallèle à l'axe"
            )));
        }
        let u = u.map(|x| x / n);
        Ok([origin, u, cross(axis, u), axis])
    }

    fn edge(&mut self, id: u64, from: usize) -> Result<usize, IoError> {
        if let Some(&e) = self.edge_of.get(&id) {
            return Ok(e);
        }
        let line = self.p.get(id, from)?.line;
        let ec = self.only(id, from, &["EDGE_CURVE"])?;
        let start = self.vertex(
            self.reference(&ec.params, 1, line, "premier sommet de l'arête")?,
            line,
        )?;
        let end = self.vertex(
            self.reference(&ec.params, 2, line, "second sommet de l'arête")?,
            line,
        )?;
        let curve_id = self.reference(&ec.params, 3, line, "courbe de l'arête")?;
        let same = boolean(ec.params.get(4), line, "sens de la courbe")?;
        let curve = self.curve(curve_id, line, start, end, same)?;
        self.edges.push(BrepEdge { start, end, curve });
        self.edge_of.insert(id, self.edges.len() - 1);
        Ok(self.edges.len() - 1)
    }

    /// Courbe de l'arête de `start` à `end` (`same` : dans le sens de la courbe).
    fn curve(
        &self,
        id: u64,
        from: usize,
        start: usize,
        end: usize,
        same: bool,
    ) -> Result<EdgeCurve, IoError> {
        let i = self.p.get(id, from)?;
        let line = i.line;
        let (a, b) = (self.vertices[start], self.vertices[end]);
        let names: Vec<&str> = i.records.iter().map(|r| r.name.as_str()).collect();
        if let Some(r) = i.record("SURFACE_CURVE").or(i.record("SEAM_CURVE")) {
            // courbe 3D d'une courbe tracée sur une surface
            let curve = self.reference(&r.params, 1, line, "courbe 3D")?;
            return self.curve(curve, line, start, end, same);
        }
        if i.record("LINE").is_some() {
            if norm(sub(b, a)) <= self.tolerance {
                return Err(error(format!("ligne {line} : segment de longueur nulle")));
            }
            return Ok(EdgeCurve::Line { from: a, to: b });
        }
        if let Some(r) = i.record("CIRCLE").or(i.record("ELLIPSE")) {
            let [center, u, v, _] =
                self.placement(self.reference(&r.params, 1, line, "repère")?, line)?;
            let radii = if r.name == "CIRCLE" {
                let radius = positive(r.params.get(2), line, "rayon")?;
                [radius, radius]
            } else {
                [
                    positive(r.params.get(2), line, "premier demi-axe")?,
                    positive(r.params.get(3), line, "second demi-axe")?,
                ]
            };
            let angle = |p: [f64; 3]| -> Result<f64, IoError> {
                let d = sub(p, center);
                let (x, y) = (dot(d, u), dot(d, v));
                let t = (y / radii[1]).atan2(x / radii[0]);
                let on = [0, 1, 2]
                    .map(|k| center[k] + radii[0] * t.cos() * u[k] + radii[1] * t.sin() * v[k]);
                if norm(sub(on, p)) > self.tolerance {
                    return Err(error(format!(
                        "ligne {line} : un sommet de l'arête n'est pas sur son {} (écart {:e})",
                        if r.name == "CIRCLE" {
                            "cercle"
                        } else {
                            "ellipse"
                        },
                        norm(sub(on, p))
                    )));
                }
                Ok(t)
            };
            let (t0, t1) = (angle(a)?, angle(b)?);
            let closed = start == end;
            let sweep = match (same, closed) {
                (true, true) => TAU,
                (false, true) => -TAU,
                (true, false) => (t1 - t0).rem_euclid(TAU),
                (false, false) => -(t0 - t1).rem_euclid(TAU),
            };
            return Ok(EdgeCurve::Arc {
                center,
                u,
                v,
                radii,
                start_angle: t0,
                sweep,
            });
        }
        if i.record("B_SPLINE_CURVE_WITH_KNOTS").is_some() {
            let curve = self.bspline(i)?;
            let [lo, hi] = curve.domain();
            let param = |p: [f64; 3]| -> Result<f64, IoError> {
                let (t, d) = curve.parameter_of(p);
                if d > self.tolerance {
                    return Err(error(format!(
                        "ligne {line} : un sommet de l'arête n'est pas sur sa B-spline (écart {d:e})"
                    )));
                }
                Ok(t)
            };
            let range = if start == end {
                if same {
                    [lo, hi]
                } else {
                    [hi, lo]
                }
            } else {
                let (t0, t1) = (param(a)?, param(b)?);
                if (t1 > t0) != same {
                    return Err(error(format!(
                        "ligne {line} : arête à cheval sur la jointure d'une B-spline fermée, pas encore prise en charge"
                    )));
                }
                [t0, t1]
            };
            return Ok(EdgeCurve::BSpline { curve, range });
        }
        Err(error(format!(
            "ligne {line} : courbe {} pas encore prise en charge (Maille lit LINE, CIRCLE, ELLIPSE et B_SPLINE_CURVE_WITH_KNOTS)",
            names.join(" + ")
        )))
    }

    /// B-spline simple (B_SPLINE_CURVE_WITH_KNOTS) ou complexe (B_SPLINE_CURVE +
    /// B_SPLINE_CURVE_WITH_KNOTS + RATIONAL_B_SPLINE_CURVE…).
    fn bspline(&self, i: &Instance) -> Result<BSpline, IoError> {
        let line = i.line;
        let knots_record = i
            .record("B_SPLINE_CURVE_WITH_KNOTS")
            .ok_or_else(|| error(format!("ligne {line} : B-spline sans nœuds")))?;
        // instance simple : nom, degré, points, forme, fermée, auto-intersection, multiplicités,
        // nœuds ; complexe : chaque partie n'a que ses propres attributs
        let (curve_params, knot_params): (&[Value], &[Value]) = match i.record("B_SPLINE_CURVE") {
            Some(c) => (&c.params[..], &knots_record.params[..]),
            None => {
                let p = &knots_record.params;
                if p.len() < 8 {
                    return Err(error(format!(
                        "ligne {line} : B_SPLINE_CURVE_WITH_KNOTS incomplète"
                    )));
                }
                (&p[1..6], &p[6..])
            }
        };
        let degree = match curve_params.first() {
            Some(Value::Integer(d)) if *d > 0 => *d as usize,
            _ => {
                return Err(error(format!(
                    "ligne {line} : degré de la B-spline attendu"
                )))
            }
        };
        let points = self
            .references(curve_params, 1, line, "points de contrôle")?
            .into_iter()
            .map(|p| self.point(p, line))
            .collect::<Result<Vec<_>, _>>()?;
        let integers = |v: Option<&Value>, what: &str| -> Result<Vec<usize>, IoError> {
            match v {
                Some(Value::List(l)) => l
                    .iter()
                    .map(|x| match x {
                        Value::Integer(m) if *m > 0 => Ok(*m as usize),
                        _ => Err(error(format!(
                            "ligne {line} : {what} : entiers positifs attendus"
                        ))),
                    })
                    .collect(),
                _ => Err(error(format!("ligne {line} : {what} : liste attendue"))),
            }
        };
        let reals = |v: Option<&Value>, what: &str| -> Result<Vec<f64>, IoError> {
            match v {
                Some(Value::List(l)) => l
                    .iter()
                    .map(|x| {
                        number(x).ok_or_else(|| {
                            error(format!("ligne {line} : {what} : nombres attendus"))
                        })
                    })
                    .collect(),
                _ => Err(error(format!("ligne {line} : {what} : liste attendue"))),
            }
        };
        let multiplicities = integers(knot_params.first(), "multiplicités des nœuds")?;
        let values = reals(knot_params.get(1), "nœuds")?;
        if multiplicities.len() != values.len() {
            return Err(error(format!(
                "ligne {line} : autant de multiplicités que de nœuds attendu"
            )));
        }
        let knots = values
            .iter()
            .zip(&multiplicities)
            .flat_map(|(&u, &m)| std::iter::repeat_n(u, m))
            .collect();
        let weights = match i.record("RATIONAL_B_SPLINE_CURVE") {
            Some(r) => reals(r.params.first(), "poids")?,
            None => Vec::new(),
        };
        BSpline::new(degree, points, weights, knots)
            .map_err(|e| error(format!("ligne {line} : {e}")))
    }
}

/// Chaque arête doit être parcourue deux fois, une fois dans chaque sens : sinon la coque n'est
/// pas fermée, ou ses faces ne sont pas orientées de façon cohérente.
fn check_closed(brep: &Brep, edge_of: &HashMap<u64, usize>) -> Result<(), IoError> {
    let mut uses = vec![[0usize; 2]; brep.edges.len()];
    for f in &brep.faces {
        for &(e, reversed) in f.loops.iter().flatten() {
            uses[e][reversed as usize] += 1;
        }
    }
    if let Some(e) = uses.iter().position(|u| *u != [1, 1]) {
        let id = edge_of
            .iter()
            .find(|&(_, &v)| v == e)
            .map_or(0, |(&id, _)| id);
        return Err(error(format!(
            "l'arête #{id} est parcourue {} fois dans un sens et {} dans l'autre : la coque n'est pas fermée ou mal orientée",
            uses[e][0], uses[e][1]
        )));
    }
    Ok(())
}

/// Unité de longueur du fichier : `mm`, `m`, `inch`… ou `inconnue`.
fn length_unit(p: &Part21) -> String {
    for i in p.instances.values() {
        if i.record("LENGTH_UNIT").is_none() {
            continue;
        }
        if let Some(si) = i.record("SI_UNIT") {
            let prefix = match si.params.first() {
                Some(Value::Enum(e)) => match e.as_str() {
                    "MILLI" => "m",
                    "CENTI" => "c",
                    "DECI" => "d",
                    "KILO" => "k",
                    "MICRO" => "µ",
                    "NANO" => "n",
                    _ => "?",
                },
                _ => "",
            };
            return format!("{prefix}m");
        }
        if let Some(c) = i.record("CONVERSION_BASED_UNIT") {
            if let Some(Value::String(name)) = c.params.first() {
                // une unité convertie peut porter le nom d'une unité SI (Onshape écrit « METRE »)
                let name = name.trim().to_lowercase();
                return match name.as_str() {
                    "metre" | "meter" => "m".to_string(),
                    "millimetre" | "millimeter" => "mm".to_string(),
                    "centimetre" | "centimeter" => "cm".to_string(),
                    _ => name,
                };
            }
        }
    }
    "inconnue".to_string()
}

fn number(v: &Value) -> Option<f64> {
    match v {
        Value::Real(x) => Some(*x),
        Value::Integer(i) => Some(*i as f64),
        _ => None,
    }
}

fn positive(v: Option<&Value>, line: usize, what: &str) -> Result<f64, IoError> {
    match v.and_then(number) {
        Some(x) if x.is_finite() && x > 0.0 => Ok(x),
        _ => Err(error(format!(
            "ligne {line} : {what} strictement positif attendu"
        ))),
    }
}

fn boolean(v: Option<&Value>, line: usize, what: &str) -> Result<bool, IoError> {
    match v {
        Some(Value::Enum(e)) if e == "T" => Ok(true),
        Some(Value::Enum(e)) if e == "F" => Ok(false),
        _ => Err(error(format!("ligne {line} : {what} : .T. ou .F. attendu"))),
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
