//! Maillage des faces d'une `Brep` à la taille h, première étape de `tetrahedralize_brep`. Chaque
//! arête est découpée une seule fois, en 3D sur sa courbe : en morceaux d'au plus h, et d'au plus un
//! seizième de tour là où elle tourne. Chaque face est mise à plat, maillée par le mailleur 2D avec
//! ses contours, puis remise sur sa surface : plan, cylindre et cône se déroulent sans déformer les
//! longueurs ; la sphère (projection stéréographique) et le tore (coordonnées isothermes) gardent
//! les angles, donc la qualité des triangles, et la taille h y est divisée par l'échelle locale de
//! la carte. Sur une face courbe, la taille est en plus bornée par un seizième de tour
//! (flèche des cordes d'environ 2 % du rayon). Comme pour les facettes (`remesh`), un morceau de
//! bord coupé par le mailleur 2D rejoint l'arête, en un point de sa courbe, et les faces qui la
//! bordent sont remaillées avec lui : la surface reste fermée.

use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::{PI, TAU};

use super::brep::{Brep, BrepEdge, EdgeCurve, FaceSurface};
use super::geometry::Boundary;
use super::refine::{Mesher, SizeField};
use super::remesh::{acute_corners, basis};
use super::surface::{cross, dot, norm, sub, Patch, Surface};
use super::{MeshOptions, MAX_NODES};
use crate::error::MailleError;
use crate::material::positive;

/// Morceaux par tour au moins, sur une arête qui tourne et sur une face courbe.
const SIDES: f64 = 16.0;
/// Maillages de face au plus, rapportés au nombre de faces : au-delà, les coupes ne se
/// stabilisent pas d'une face à l'autre.
const MAX_REMESH: usize = 20;
/// Échantillons d'une arête pour la découper.
const SAMPLES: usize = 256;
/// Écart toléré entre un sommet et la surface de sa face, rapporté à la diagonale de la pièce.
const ON_SURFACE: f64 = 1e-6;
/// Points des arêtes au plus, rapportés à ceux du premier découpage : au-delà, les coupes
/// s'emballent d'une face à l'autre (mesuré : 1,2 fois sur les pièces qui se maillent, croissance
/// exponentielle autour d'une pointe entre faces).
const MAX_GROWTH: usize = 4;
const RUNAWAY: &str = "la pièce ne se laisse pas mailler à cette taille h : les coupes s'emballent entre des faces qui se rejoignent tangentes (pointe, par exemple le chanfrein d'un écrou contre un de ses pans) ; essayer une autre taille h";
const STUCK: &str = "la pièce ne se laisse pas mailler à la taille h : les coupes ne se stabilisent pas entre deux faces (coin ou dièdre de quelques degrés) ; arrondir ou épaissir ce coin";

/// Surface fermée des faces de `brep` maillées à la taille `size` avec la gradation `grading` :
/// un Patch par nom de face, triangles de rayon circonscrit au plus h/√3 dans le plan déroulé,
/// sommets sur les surfaces des faces. `InvalidProperty` si h ou la gradation n'est pas
/// strictement positive, `InvalidBrep` si une face ne se laisse pas mailler.
pub(crate) fn triangulate(brep: &Brep, size: f64, grading: f64) -> Result<Surface, MailleError> {
    let size = positive("taille h", size)?;
    let grading = positive("gradation", grading)?;
    let invalid = |problem: String| MailleError::InvalidBrep { problem };
    let tolerance = ON_SURFACE * diagonal(&brep.vertices);
    check_indices(brep)?;
    let brep = &prepare_closed(brep)?;

    // carte de chaque face, sommets sur sa surface ; taille de la face, bornée par un seizième de
    // tour sur son plus petit rayon de courbure (le long de ses arêtes, le cône rétrécit)
    let mut charts = brep
        .faces
        .iter()
        .map(|f| Chart::new(&f.surface, &f.name))
        .collect::<Result<Vec<_>, _>>()?;
    let mut face_sizes = vec![size; brep.faces.len()];
    for ((face, chart), h) in brep.faces.iter().zip(&charts).zip(&mut face_sizes) {
        match chart.kind {
            Kind::Sphere { radius, .. } => *h = h.min(radius * TAU / SIDES),
            Kind::Torus { major, minor, .. } => *h = h.min(minor.min(major - minor) * TAU / SIDES),
            _ => {}
        }
        for &(e, _) in face.loops.iter().flatten() {
            let edge = &brep.edges[e];
            for v in [edge.start, edge.end] {
                let p = brep.vertices[v];
                let gap = chart.gap(p);
                if gap > tolerance {
                    return Err(invalid(format!(
                        "face « {} » : un sommet est à {gap:e} de sa surface",
                        face.name
                    )));
                }
            }
            if let Kind::Round { radius, sin, .. } = chart.kind {
                let rho = if sin == 0.0 {
                    radius
                } else {
                    (0..=16)
                        .map(|i| chart.radial(edge.curve.point(i as f64 / 16.0)))
                        .fold(f64::INFINITY, f64::min)
                };
                if rho <= tolerance {
                    return Err(invalid(format!(
                        "face « {} » : cône jusqu'à son sommet, pas encore pris en charge",
                        face.name
                    )));
                }
                *h = h.min(rho * TAU / SIDES);
            }
        }
    }
    let brep = &add_seams(brep, &charts)?;
    // une face sphérique est vue depuis l'opposé de son milieu
    for (face, chart) in brep.faces.iter().zip(&mut charts) {
        let ring: Vec<[f64; 3]> = face.loops[0]
            .iter()
            .flat_map(|&(e, reversed)| {
                (0..16).map(move |i| {
                    let t = i as f64 / 16.0;
                    brep.edges[e]
                        .curve
                        .point(if reversed { 1.0 - t } else { t })
                })
            })
            .collect();
        chart.aim(&ring, &face.name)?;
    }

    // taille d'une arête : la plus petite de ses faces
    let mut sizes = vec![size; brep.edges.len()];
    let mut faces_of = vec![Vec::new(); brep.edges.len()];
    for (f, face) in brep.faces.iter().enumerate() {
        for &(e, _) in face.loops.iter().flatten() {
            sizes[e] = sizes[e].min(face_sizes[f]);
            faces_of[e].push(f);
        }
    }

    // points de chaque arête, du début à la fin, avec leur paramètre ; taille de ses morceaux
    if brep.vertices.len() > MAX_NODES {
        return Err(MailleError::TooManyNodes { max: MAX_NODES });
    }
    let mut points = brep.vertices.clone();
    let mut chains: Vec<Vec<(f64, usize)>> = Vec::with_capacity(brep.edges.len());
    for (e, edge) in brep.edges.iter().enumerate() {
        let (ts, length) = divide(&edge.curve, sizes[e], edge.start == edge.end)?;
        let n = ts.len() - 1;
        sizes[e] = sizes[e].min(length / n as f64);
        let mut chain = Vec::with_capacity(n + 1);
        for (i, &t) in ts.iter().enumerate() {
            let v = match i {
                0 => edge.start,
                i if i == n => edge.end,
                _ => add(&mut points, edge.curve.point(t))?,
            };
            chain.push((t, v));
        }
        chains.push(chain);
    }

    let initial: usize = chains.iter().map(Vec::len).sum();
    let mut meshes: Vec<Vec<[usize; 3]>> = vec![Vec::new(); brep.faces.len()];
    let mut pending: BTreeSet<usize> = (0..brep.faces.len()).collect();
    let mut budget = MAX_REMESH * brep.faces.len();
    while let Some(f) = pending.pop_first() {
        budget = budget
            .checked_sub(1)
            .ok_or_else(|| invalid(STUCK.to_string()))?;
        let face = &brep.faces[f];
        if chains.iter().map(Vec::len).sum::<usize>() > MAX_GROWTH * initial {
            return Err(invalid(RUNAWAY.to_string()));
        }
        let (triangles, origins) = mesh_face(
            brep,
            f,
            charts[f],
            (&chains, &sizes, &points),
            face_sizes[f],
            grading,
        )
            .map_err(|e| match e {
                MailleError::MeshingFailed => invalid(format!(
                    "la face « {} » ne se laisse pas mailler à la taille h (contour qui se recoupe une fois déroulé, ou coin très aigu)",
                    face.name
                )),
                e => e,
            })?;
        let mut splits: BTreeMap<(usize, usize), Vec<(f64, usize)>> = BTreeMap::new();
        // une couture est bordée deux fois par sa face : une coupe vue des deux côtés est un point
        let mut seen: HashMap<(usize, usize, u64), usize> = HashMap::new();
        let mut global = Vec::with_capacity(origins.len());
        for origin in origins {
            global.push(match origin {
                Origin::Known(i) => i,
                Origin::Interior(p) => add(&mut points, p)?,
                Origin::Split { edge, k, t } => match seen.entry((edge, k, t.to_bits())) {
                    Entry::Occupied(o) => *o.get(),
                    Entry::Vacant(slot) => {
                        let c = &chains[edge];
                        let s = c[k].0 + t * (c[k + 1].0 - c[k].0);
                        let i = add(&mut points, brep.edges[edge].curve.point(s))?;
                        splits.entry((edge, k)).or_default().push((s, i));
                        *slot.insert(i)
                    }
                },
            });
        }
        // du dernier morceau au premier, pour que les rangs k restent justes. La face est remaillée
        // aussi si le point n'est pas sur la corde où le mailleur 2D l'a posé (arête courbe) ou si
        // l'arête est sa couture : l'autre côté de la couture n'a pas encore le point
        for ((edge, k), mut new) in splits.into_iter().rev() {
            new.sort_by(|a, b| a.0.total_cmp(&b.0));
            chains[edge].splice(k + 1..k + 1, new);
            let curved = !matches!(brep.edges[edge].curve, EdgeCurve::Line { .. });
            let seam = faces_of[edge].iter().filter(|&&g| g == f).count() > 1;
            pending.extend(faces_of[edge].iter().filter(|&&g| g != f || curved || seam));
        }
        meshes[f] = triangles
            .into_iter()
            .map(|t| t.map(|v| global[v]))
            .collect();
    }

    // surface finale : un Patch par nom de face, sans les points des maillages remplacés
    let mut out = Surface {
        vertices: Vec::new(),
        triangles: Vec::new(),
        patches: Vec::new(),
    };
    let mut number = vec![usize::MAX; points.len()];
    for (face, triangles) in brep.faces.iter().zip(&meshes) {
        let patch = match out.patches.iter().position(|p| p.name == face.name) {
            Some(p) => p,
            None => {
                out.patches.push(Patch {
                    name: face.name.clone(),
                    triangles: Vec::new(),
                });
                out.patches.len() - 1
            }
        };
        for t in triangles {
            let t = t.map(|v| {
                if number[v] == usize::MAX {
                    number[v] = out.vertices.len();
                    out.vertices.push(points[v]);
                }
                number[v]
            });
            out.patches[patch].triangles.push(out.triangles.len());
            out.triangles.push(t);
        }
    }
    Ok(out)
}

/// Indices des arêtes et des sommets dans les bornes.
fn check_indices(brep: &Brep) -> Result<(), MailleError> {
    let invalid = |problem: String| MailleError::InvalidBrep { problem };
    for face in &brep.faces {
        for &(e, _) in face.loops.iter().flatten() {
            if e >= brep.edges.len() {
                return Err(invalid(format!(
                    "face « {} » : arête {e} absente (la pièce en a {})",
                    face.name,
                    brep.edges.len()
                )));
            }
        }
    }
    for (e, edge) in brep.edges.iter().enumerate() {
        for v in [edge.start, edge.end] {
            if v >= brep.vertices.len() {
                return Err(invalid(format!(
                    "arête {e} : sommet {v} absent (la pièce en a {})",
                    brep.vertices.len()
                )));
            }
        }
    }
    Ok(())
}

/// Faces fermées prêtes pour leur carte. Sur une sphère, la couture (arête parcourue dans les deux
/// sens par la face) n'est pas un bord : elle est retirée avec les contours qu'elle vide, et une
/// sphère entière, sans plus aucun contour, est coupée en deux hémisphères par un équateur, sous
/// le même nom. Un tore entier sans arête (OCCT l'écrit avec un seul contour réduit à un sommet)
/// reçoit ses deux coutures, un parallèle et un méridien, orientées pour que la matière soit dans
/// le tube. Une autre face sans contour est refusée.
fn prepare_closed(brep: &Brep) -> Result<Brep, MailleError> {
    let mut out = brep.clone();
    let mut halves = Vec::new();
    for f in 0..out.faces.len() {
        let face = &mut out.faces[f];
        if let FaceSurface::Sphere { .. } = face.surface {
            let uses: Vec<(usize, bool)> = face.loops.iter().flatten().copied().collect();
            let seam = |e: usize| uses.contains(&(e, false)) && uses.contains(&(e, true));
            for lp in &mut face.loops {
                lp.retain(|&(e, _)| !seam(e));
            }
            face.loops.retain(|lp| !lp.is_empty());
        }
        if !out.faces[f].loops.is_empty() {
            continue;
        }
        match out.faces[f].surface {
            FaceSurface::Sphere { center, radius, .. } => {
                out.vertices
                    .push([center[0] + radius, center[1], center[2]]);
                let v = out.vertices.len() - 1;
                let equator = full_circle(
                    &mut out,
                    v,
                    center,
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    radius,
                );
                // vu de dehors, l'hémisphère nord est à gauche de l'équateur antihoraire
                out.faces[f].loops = vec![vec![(equator, false)]];
                let mut south = out.faces[f].clone();
                south.loops = vec![vec![(equator, true)]];
                halves.push(south);
            }
            FaceSurface::Torus {
                center,
                axis,
                major,
                minor,
            } => {
                let l = norm(axis);
                let axis = axis.map(|x| x / l);
                let (u, w) = basis(axis);
                let at = |r: f64| std::array::from_fn(|k| center[k] + r * u[k]);
                out.vertices.push(at(major + minor));
                let v = out.vertices.len() - 1;
                let parallel = full_circle(&mut out, v, center, [u, w], major + minor);
                let meridian = full_circle(&mut out, v, at(major), [u, axis], minor);
                // (θ, φ) croissants : sens direct autour de la normale qui sort du tube
                out.faces[f].loops = vec![vec![
                    (parallel, false),
                    (meridian, false),
                    (parallel, true),
                    (meridian, true),
                ]];
            }
            _ => {
                return Err(MailleError::InvalidBrep {
                    problem: format!("face « {} » sans contour", out.faces[f].name),
                })
            }
        }
    }
    out.faces.extend(halves);
    Ok(out)
}

/// Cercle complet de centre `center` dans le plan (u, v), du sommet `v` au même sommet.
fn full_circle(
    brep: &mut Brep,
    vertex: usize,
    center: [f64; 3],
    [u, v]: [[f64; 3]; 2],
    r: f64,
) -> usize {
    brep.edges.push(BrepEdge {
        start: vertex,
        end: vertex,
        curve: EdgeCurve::Arc {
            center,
            u,
            v,
            radii: [r, r],
            start_angle: 0.0,
            sweep: TAU,
        },
    });
    brep.edges.len() - 1
}

/// Couture d'une face de révolution fermée sans couture (bordée de deux contours qui font le tour
/// de l'axe, deux cercles par exemple, comme l'écrivent certains modeleurs) : une arête le long
/// d'une génératrice (un méridien sur un tore), et les deux contours fondus en un seul qui la
/// parcourt dans les deux sens, comme la couture d'OCCT. La génératrice est prise au milieu du plus
/// grand intervalle d'angle sans sommet, loin des détails des contours ; chaque contour y est coupé
/// (la face voisine suit) et ne doit la croiser qu'une fois.
fn add_seams(brep: &Brep, charts: &[Chart]) -> Result<Brep, MailleError> {
    let mut out = brep.clone();
    for (f, chart) in charts.iter().enumerate() {
        if chart.periods()[0].is_none() {
            continue;
        }
        let name = out.faces[f].name.clone();
        let refuse = |why: &str| MailleError::InvalidBrep {
            problem: format!(
                "face « {name} » : face fermée sans couture {why}, pas encore prise en charge"
            ),
        };
        let turns: Vec<f64> = out.faces[f]
            .loops
            .iter()
            .map(|lp| winding(&out, lp, chart))
            .collect();
        let around: Vec<usize> = (0..turns.len()).filter(|&l| turns[l].abs() > PI).collect();
        let (a, b) = match around[..] {
            [] => continue,
            [a, b] if turns[a] * turns[b] < 0.0 => (a, b),
            _ => return Err(refuse("(un contour fait le tour de l'axe)")),
        };
        // angle de la couture : milieu du plus grand intervalle entre les sommets des deux contours
        let mut angles: Vec<f64> = [a, b]
            .iter()
            .flat_map(|&l| out.faces[f].loops[l].iter())
            .map(|entry| chart.coords(out.vertices[start(&out, entry)])[0])
            .collect();
        angles.sort_by(f64::total_cmp);
        let theta = (0..angles.len())
            .map(|k| {
                let next = angles.get(k + 1).copied().unwrap_or(angles[0] + TAU);
                (next - angles[k], (angles[k] + next) / 2.0)
            })
            .max_by(|x, y| x.0.total_cmp(&y.0))
            .map_or(0.0, |(_, mid)| mid);
        let mut ends = [0; 2];
        for (end, l) in ends.iter_mut().zip([a, b]) {
            let lp = out.faces[f].loops[l].clone();
            if crossings(&out, &lp, chart, theta) != 1 {
                return Err(refuse(
                    "dont un contour recroise la génératrice de la couture",
                ));
            }
            let (e, t) = lp
                .iter()
                .find_map(|&(e, _)| angle_on(&out.edges[e].curve, chart, theta).map(|t| (e, t)))
                .ok_or_else(|| {
                    refuse("dont un contour ne croise pas la génératrice de la couture")
                })?;
            split_edge(&mut out, e, t);
            *end = out.vertices.len() - 1;
        }
        let [va, vb] = ends;
        let face = &out.faces[f];
        let position = |l: usize, v: usize| {
            face.loops[l]
                .iter()
                .position(|entry| start(&out, entry) == v)
                .unwrap_or(0)
        };
        let mut merged = face.loops[a].clone();
        merged.rotate_left(position(a, va));
        merged.push((out.edges.len(), false));
        let mut other = face.loops[b].clone();
        other.rotate_left(position(b, vb));
        merged.extend(other);
        merged.push((out.edges.len(), true));
        let mut loops = vec![merged];
        loops.extend(
            (0..face.loops.len())
                .filter(|&l| l != a && l != b)
                .map(|l| face.loops[l].clone()),
        );
        out.faces[f].loops = loops;
        let curve = chart.generator(out.vertices[va], out.vertices[vb]);
        out.edges.push(BrepEdge {
            start: va,
            end: vb,
            curve,
        });
    }
    Ok(out)
}

/// Combien de fois un contour croise la génératrice d'angle `theta`.
fn crossings(brep: &Brep, lp: &[(usize, bool)], chart: &Chart, theta: f64) -> usize {
    let mut gaps = Vec::new();
    for &(e, reversed) in lp {
        for i in 0..SAMPLES {
            let t = i as f64 / SAMPLES as f64;
            let t = if reversed { 1.0 - t } else { t };
            gaps.push(wrap(chart.coords(brep.edges[e].curve.point(t))[0] - theta));
        }
    }
    (0..gaps.len())
        .filter(|&k| {
            let (g, h) = (gaps[k], gaps[(k + 1) % gaps.len()]);
            (g < 0.0) != (h < 0.0) && (g - h).abs() < PI
        })
        .count()
}

/// Sommet de départ d'une arête parcourue dans un contour.
fn start(brep: &Brep, &(e, reversed): &(usize, bool)) -> usize {
    let edge = &brep.edges[e];
    if reversed {
        edge.end
    } else {
        edge.start
    }
}

/// Paramètre t ∈ ]0, 1[ où la courbe passe à l'angle `theta` autour de l'axe de la carte :
/// échantillons, puis dichotomie sur l'écart d'angle.
fn angle_on(curve: &EdgeCurve, chart: &Chart, theta: f64) -> Option<f64> {
    let gap = |t: f64| wrap(chart.coords(curve.point(t))[0] - theta);
    let (mut lo, mut hi) = (0..SAMPLES)
        .map(|i| i as f64 / SAMPLES as f64)
        .zip((1..=SAMPLES).map(|i| i as f64 / SAMPLES as f64))
        .find(|&(a, b)| {
            let (ga, gb) = (gap(a), gap(b));
            ga * gb <= 0.0 && (ga - gb).abs() < PI
        })?;
    for _ in 0..60 {
        let m = (lo + hi) / 2.0;
        if gap(lo) * gap(m) <= 0.0 {
            hi = m;
        } else {
            lo = m;
        }
    }
    let t = (lo + hi) / 2.0;
    (t > 1e-9 && t < 1.0 - 1e-9).then_some(t)
}

/// Coupe l'arête `e` en son paramètre t : un sommet de plus, l'arête garde son début, une nouvelle
/// arête prend la suite, et chaque contour qui passait par `e` passe par les deux.
fn split_edge(brep: &mut Brep, e: usize, t: f64) {
    let edge = brep.edges[e].clone();
    let q = edge.curve.point(t);
    brep.vertices.push(q);
    let v = brep.vertices.len() - 1;
    let (first, second) = match edge.curve {
        EdgeCurve::Line { from, to } => (
            EdgeCurve::Line { from, to: q },
            EdgeCurve::Line { from: q, to },
        ),
        EdgeCurve::Arc {
            center,
            u,
            v,
            radii,
            start_angle,
            sweep,
        } => (
            EdgeCurve::Arc {
                center,
                u,
                v,
                radii,
                start_angle,
                sweep: sweep * t,
            },
            EdgeCurve::Arc {
                center,
                u,
                v,
                radii,
                start_angle: start_angle + sweep * t,
                sweep: sweep * (1.0 - t),
            },
        ),
        EdgeCurve::BSpline { curve, range } => {
            let m = range[0] + t * (range[1] - range[0]);
            (
                EdgeCurve::BSpline {
                    curve: curve.clone(),
                    range: [range[0], m],
                },
                EdgeCurve::BSpline {
                    curve,
                    range: [m, range[1]],
                },
            )
        }
    };
    brep.edges[e] = BrepEdge {
        start: edge.start,
        end: v,
        curve: first,
    };
    brep.edges.push(BrepEdge {
        start: v,
        end: edge.end,
        curve: second,
    });
    let next = brep.edges.len() - 1;
    for face in &mut brep.faces {
        for lp in &mut face.loops {
            *lp = lp
                .iter()
                .flat_map(|&(k, reversed)| match (k == e, reversed) {
                    (false, _) => vec![(k, reversed)],
                    (true, false) => vec![(e, false), (next, false)],
                    (true, true) => vec![(next, true), (e, true)],
                })
                .collect();
        }
    }
}

/// De combien un contour tourne autour de l'axe d'un cylindre ou d'un cône : 0 ou ±2π.
fn winding(brep: &Brep, lp: &[(usize, bool)], chart: &Chart) -> f64 {
    let mut angles = Vec::new();
    for &(e, reversed) in lp {
        for i in 0..32 {
            let t = i as f64 / 32.0;
            let t = if reversed { 1.0 - t } else { t };
            angles.push(chart.coords(brep.edges[e].curve.point(t))[0]);
        }
    }
    (0..angles.len())
        .map(|k| wrap(angles[(k + 1) % angles.len()] - angles[k]))
        .sum()
}

fn add(points: &mut Vec<[f64; 3]>, p: [f64; 3]) -> Result<usize, MailleError> {
    if points.len() >= MAX_NODES {
        return Err(MailleError::TooManyNodes { max: MAX_NODES });
    }
    points.push(p);
    Ok(points.len() - 1)
}

/// Paramètres des points d'une arête (0 et 1 compris) et sa longueur : morceaux de même « coût »,
/// le coût d'un bout de courbe étant sa longueur rapportée à h, ou son virage rapporté à un
/// seizième de tour s'il est plus grand. Une courbe fermée a au moins trois morceaux.
fn divide(curve: &EdgeCurve, size: f64, closed: bool) -> Result<(Vec<f64>, f64), MailleError> {
    let p: Vec<[f64; 3]> = (0..=SAMPLES)
        .map(|i| curve.point(i as f64 / SAMPLES as f64))
        .collect();
    let chords: Vec<[f64; 3]> = (0..SAMPLES).map(|i| sub(p[i + 1], p[i])).collect();
    let turn = |a: [f64; 3], b: [f64; 3]| {
        if norm(a) == 0.0 || norm(b) == 0.0 {
            0.0
        } else {
            norm(cross(a, b)).atan2(dot(a, b))
        }
    };
    let mut cost = vec![0.0; SAMPLES + 1];
    for i in 0..SAMPLES {
        let bend = chords.get(i + 1).map_or(0.0, |&next| turn(chords[i], next));
        cost[i + 1] = cost[i] + (norm(chords[i]) / size).max(bend * SIDES / TAU);
    }
    let length: f64 = chords.iter().map(|&c| norm(c)).sum();
    let total = cost[SAMPLES];
    let n = total.ceil().max(if closed { 3.0 } else { 1.0 });
    if !n.is_finite() || n > MAX_NODES as f64 {
        return Err(MailleError::TooManyNodes { max: MAX_NODES });
    }
    let n = n as usize;
    let mut ts = vec![0.0];
    let mut i = 0;
    for j in 1..n {
        let target = total * j as f64 / n as f64;
        while i + 1 < SAMPLES && cost[i + 1] < target {
            i += 1;
        }
        let span = cost[i + 1] - cost[i];
        let frac = if span > 0.0 {
            (target - cost[i]) / span
        } else {
            0.0
        };
        ts.push((i as f64 + frac) / SAMPLES as f64);
    }
    ts.push(1.0);
    Ok((ts, length))
}

/// D'où vient un point d'une face maillée : d'une arête, d'une coupe du morceau `k` de l'arête
/// `edge` (paramètre t dans le sens de l'arête), ou de l'intérieur.
enum Origin {
    Known(usize),
    Split { edge: usize, k: usize, t: f64 },
    Interior([f64; 3]),
}

/// Points de chaque arête (paramètre, indice), taille de leurs morceaux, coordonnées des points.
type EdgePoints<'a> = (&'a [Vec<(f64, usize)>], &'a [f64], &'a [[f64; 3]]);

/// Maille la face `f` déroulée à partir des points actuels de ses arêtes ; renvoie les triangles
/// (indices locaux) et l'origine de chaque point local.
fn mesh_face(
    brep: &Brep,
    f: usize,
    mut chart: Chart,
    (chains, sizes, points): EdgePoints,
    face_size: f64,
    grading: f64,
) -> Result<(Vec<[usize; 3]>, Vec<Origin>), MailleError> {
    let face = &brep.faces[f];
    // morceaux de bord dans le sens de la face : (arête, rang du morceau, à rebours)
    // et groupe de chaque morceau : son arête dans le contour (une couture y passe deux fois)
    let mut pieces = Vec::new();
    let mut groups = Vec::new();
    let mut known = Vec::new();
    let mut rings = Vec::new();
    for lp in &face.loops {
        let mut ring = Vec::new();
        for &(e, reversed) in lp {
            let m = chains[e].len() - 1;
            let group = groups.last().map_or(0, |g| g + 1);
            for j in 0..m {
                let k = if reversed { m - 1 - j } else { j };
                let start = chains[e][if reversed { k + 1 } else { k }].1;
                ring.push((chart.coords(points[start]), sizes[e]));
                pieces.push((e, k, reversed));
                groups.push(group);
                known.push(start);
            }
        }
        rings.push(ring);
    }
    chart.unroll(&mut rings, &face.name)?;
    // tailles du plan : celles de la surface divisées par l'échelle de la carte
    for p in rings.iter_mut().flatten() {
        p.1 /= chart.scale(p.0);
    }
    let boundary = Boundary::polygons(&rings, &format!("face « {} »", face.name));
    let mut m = Mesher::new(&boundary)?;
    m.acute = acute_corners(&rings);
    m.groups = groups;
    let options = MeshOptions {
        grading,
        ..MeshOptions::default()
    };
    let field = SizeField::new(&boundary, grading);
    let field = match chart.kind {
        Kind::Sphere { .. } | Kind::Torus { .. } => {
            field.limited(move |x| face_size / chart.scale(x))
        }
        _ => field,
    };
    m.refine(&options, &field)?;
    m.carve()?;

    let split: HashMap<usize, (usize, f64)> = m
        .segments
        .iter()
        .filter(|s| s.t[0] > 0.0)
        .map(|s| (s.a, (s.piece, s.t[0])))
        .collect();
    let origins = (0..m.d.point_count())
        .map(|i| {
            if i < known.len() {
                Origin::Known(known[i])
            } else if let Some(&(piece, t)) = split.get(&i) {
                let (edge, k, reversed) = pieces[piece];
                let t = if reversed { 1.0 - t } else { t };
                Origin::Split { edge, k, t }
            } else {
                Origin::Interior(chart.lift(m.d.point(i)))
            }
        })
        .collect();
    Ok((m.d.inside_triangles(), origins))
}

/// Carte d'une face : de sa surface au plan du mailleur 2D. Plan, cylindre et cône se déroulent
/// sans déformer les longueurs ; la sphère (projection stéréographique) et le tore (coordonnées
/// isothermes) gardent les angles, et leurs longueurs sont multipliées par `scale` en revenant sur
/// la surface.
#[derive(Debug, Clone, Copy)]
struct Chart {
    kind: Kind,
    /// x changé en −x, pour que le contour extérieur tourne dans le sens direct.
    mirror: bool,
    /// Angle déroulé (cône) amené en 0 dans le plan.
    shift: f64,
}

#[derive(Debug, Clone, Copy)]
enum Kind {
    /// Coordonnées (x, y) dans le repère (e1, e2) du plan.
    Plane {
        origin: [f64; 3],
        e1: [f64; 3],
        e2: [f64; 3],
        normal: [f64; 3],
    },
    /// Cylindre (sin = 0) ou cône : point o + ρ(v) (cos θ u + sin θ w) + v axe, ρ(v) = r + v tan α,
    /// coordonnées (θ, v). Le cylindre se déroule en (r θ, v), le cône en coordonnées polaires
    /// (ρ / sin α, θ sin α) autour de son sommet.
    Round {
        origin: [f64; 3],
        u: [f64; 3],
        w: [f64; 3],
        axis: [f64; 3],
        radius: f64,
        tan: f64,
        sin: f64,
    },
    /// Sphère projetée depuis le point c − R n sur le plan de son centre (repère e1, e2) :
    /// coordonnées R (q·e1, q·e2) / (1 + q·n), q = (p − c) / R. `n` vise le milieu de la face.
    Sphere {
        center: [f64; 3],
        radius: f64,
        convex: bool,
        n: [f64; 3],
        e1: [f64; 3],
        e2: [f64; 3],
    },
    /// Tore : point c + (R + r cos φ) (cos θ u + sin θ w) + r sin φ axe, coordonnées (θ, ψ) avec
    /// dψ = r dφ / (R + r cos φ), où la métrique est (R + r cos φ)² (dθ² + dψ²) ; mis à plat en
    /// (R θ, R ψ). `k` = √(R² − r²).
    Torus {
        center: [f64; 3],
        u: [f64; 3],
        w: [f64; 3],
        axis: [f64; 3],
        major: f64,
        minor: f64,
        k: f64,
    },
}

impl Chart {
    fn new(surface: &FaceSurface, name: &str) -> Result<Self, MailleError> {
        let invalid = |what: &str| MailleError::InvalidBrep {
            problem: format!("face « {name} » : {what}"),
        };
        let unit = |v: [f64; 3], what: &str| {
            let n = norm(v);
            if v.iter().all(|x| x.is_finite()) && n > 0.0 {
                Ok(v.map(|x| x / n))
            } else {
                Err(invalid(what))
            }
        };
        let point = |p: [f64; 3]| {
            if p.iter().all(|x| x.is_finite()) {
                Ok(p)
            } else {
                Err(invalid("origine non finie"))
            }
        };
        let kind = match *surface {
            FaceSurface::Plane { origin, normal } => {
                let normal = unit(normal, "normale du plan nulle ou non finie")?;
                let (e1, e2) = basis(normal);
                Kind::Plane {
                    origin: point(origin)?,
                    e1,
                    e2,
                    normal,
                }
            }
            FaceSurface::Cylinder {
                origin,
                axis,
                radius,
            } => {
                if !(radius.is_finite() && radius > 0.0) {
                    return Err(invalid("rayon du cylindre non strictement positif"));
                }
                let axis = unit(axis, "axe du cylindre nul ou non fini")?;
                let (u, w) = basis(axis);
                Kind::Round {
                    origin: point(origin)?,
                    u,
                    w,
                    axis,
                    radius,
                    tan: 0.0,
                    sin: 0.0,
                }
            }
            FaceSurface::Cone {
                origin,
                axis,
                radius,
                half_angle,
            } => {
                if !(radius.is_finite() && radius >= 0.0) {
                    return Err(invalid("rayon du cône négatif ou non fini"));
                }
                if !(half_angle > 0.0 && half_angle < PI / 2.0) {
                    return Err(invalid("demi-angle du cône hors de ]0, π/2["));
                }
                let axis = unit(axis, "axe du cône nul ou non fini")?;
                let (u, w) = basis(axis);
                Kind::Round {
                    origin: point(origin)?,
                    u,
                    w,
                    axis,
                    radius,
                    tan: half_angle.tan(),
                    sin: half_angle.sin(),
                }
            }
            FaceSurface::Sphere {
                center,
                radius,
                convex,
            } => {
                if !(radius.is_finite() && radius > 0.0) {
                    return Err(invalid("rayon de la sphère non strictement positif"));
                }
                let n = [0.0, 0.0, 1.0];
                let (e1, e2) = basis(n);
                Kind::Sphere {
                    center: point(center)?,
                    radius,
                    convex,
                    n,
                    e1,
                    e2,
                }
            }
            FaceSurface::Torus {
                center,
                axis,
                major,
                minor,
            } => {
                if !(minor.is_finite() && minor > 0.0 && major.is_finite()) {
                    return Err(invalid("petit rayon du tore non strictement positif"));
                }
                if minor >= major {
                    return Err(invalid("tore croisé (petit rayon au moins égal au grand : la surface se recoupe sur l'axe), pas pris en charge"));
                }
                let axis = unit(axis, "axe du tore nul ou non fini")?;
                let (u, w) = basis(axis);
                Kind::Torus {
                    center: point(center)?,
                    u,
                    w,
                    axis,
                    major,
                    minor,
                    k: (major * major - minor * minor).sqrt(),
                }
            }
        };
        Ok(Self {
            kind,
            mirror: false,
            shift: 0.0,
        })
    }

    /// Périodes des deux coordonnées de `coords` : θ (2π) d'un cylindre, d'un cône ou d'un tore,
    /// ψ du tore (2π r / k).
    fn periods(&self) -> [Option<f64>; 2] {
        match self.kind {
            Kind::Plane { .. } | Kind::Sphere { .. } => [None, None],
            Kind::Round { .. } => [Some(TAU), None],
            Kind::Torus { minor, k, .. } => [Some(TAU), Some(TAU * minor / k)],
        }
    }

    /// Vise le milieu d'une face sphérique : `n` suit l'aire orientée du contour extérieur `ring`
    /// (points de la sphère). Vue de dehors, la face est à gauche du contour : l'aire pointe vers
    /// elle si la face est bombée (`convex`), à l'opposé pour une cavité.
    fn aim(&mut self, ring: &[[f64; 3]], name: &str) -> Result<(), MailleError> {
        let Kind::Sphere {
            center,
            ref mut n,
            ref mut e1,
            ref mut e2,
            radius,
            convex,
        } = self.kind
        else {
            return Ok(());
        };
        let mut area = [0.0; 3];
        for i in 0..ring.len() {
            let c = cross(
                sub(ring[i], center),
                sub(ring[(i + 1) % ring.len()], center),
            );
            area = std::array::from_fn(|k| area[k] + c[k]);
        }
        let l = norm(area);
        let refuse = || {
            MailleError::InvalidBrep {
            problem: format!("face « {name} » : face sphérique qui couvre presque toute la sphère, pas encore prise en charge"),
        }
        };
        if l.is_nan() || l <= 1e-9 * radius * radius {
            return Err(refuse());
        }
        let side = if convex { 1.0 } else { -1.0 };
        *n = area.map(|x| side * x / l);
        (*e1, *e2) = basis(*n);
        // le point de projection ne doit pas toucher la face
        if ring
            .iter()
            .any(|&p| 1.0 + dot(sub(p, center), *n) / radius < 1e-3)
        {
            return Err(refuse());
        }
        Ok(())
    }

    /// Génératrice de la surface de `p` à `q`, deux points au même angle θ : un segment sur un
    /// cylindre ou un cône, le plus court arc de méridien sur un tore.
    fn generator(&self, p: [f64; 3], q: [f64; 3]) -> EdgeCurve {
        let Kind::Torus {
            center,
            axis,
            major,
            minor,
            ..
        } = self.kind
        else {
            return EdgeCurve::Line { from: p, to: q };
        };
        // repère du méridien : direction radiale de p, axe du tore
        let d = sub(p, center);
        let z = dot(d, axis);
        let radial: [f64; 3] = std::array::from_fn(|k| d[k] - z * axis[k]);
        let l = norm(radial);
        let radial = radial.map(|x| x / l);
        let middle: [f64; 3] = std::array::from_fn(|k| center[k] + major * radial[k]);
        let angle = |x: [f64; 3]| {
            let d = sub(x, middle);
            dot(d, axis).atan2(dot(d, radial))
        };
        let (a, b) = (angle(p), angle(q));
        EdgeCurve::Arc {
            center: middle,
            u: radial,
            v: axis,
            radii: [minor, minor],
            start_angle: a,
            sweep: wrap(b - a),
        }
    }

    /// Coordonnées (x, y) d'un plan, (θ, v) d'un cylindre ou d'un cône, θ dans ]−π, π] ;
    /// projection d'une sphère ; (θ, ψ) d'un tore, ψ dans ]−π r / k, π r / k].
    fn coords(&self, p: [f64; 3]) -> [f64; 2] {
        match self.kind {
            Kind::Plane { origin, e1, e2, .. } => {
                let d = sub(p, origin);
                [dot(d, e1), dot(d, e2)]
            }
            Kind::Round {
                origin, u, w, axis, ..
            } => {
                let d = sub(p, origin);
                [dot(d, w).atan2(dot(d, u)), dot(d, axis)]
            }
            Kind::Sphere {
                center,
                radius,
                n,
                e1,
                e2,
                ..
            } => {
                let q = sub(p, center).map(|x| x / radius);
                let d = 1.0 + dot(q, n);
                [radius * dot(q, e1) / d, radius * dot(q, e2) / d]
            }
            Kind::Torus {
                center,
                u,
                w,
                axis,
                major,
                minor,
                k,
            } => {
                let d = sub(p, center);
                let (x, y) = (dot(d, u), dot(d, w));
                let phi = dot(d, axis).atan2(x.hypot(y) - major);
                let psi = 2.0 * minor / k
                    * (((major - minor) / (major + minor)).sqrt() * (phi / 2.0).tan()).atan();
                [y.atan2(x), psi]
            }
        }
    }

    /// Distance de p à l'axe d'un cylindre, d'un cône ou d'un tore.
    fn radial(&self, p: [f64; 3]) -> f64 {
        match self.kind {
            Kind::Plane { .. } | Kind::Sphere { .. } => 0.0,
            Kind::Round { origin, axis, .. }
            | Kind::Torus {
                center: origin,
                axis,
                ..
            } => {
                let d = sub(p, origin);
                let a = dot(d, axis);
                norm(std::array::from_fn(|k| d[k] - a * axis[k]))
            }
        }
    }

    /// Distance de p à la surface (à l'ordre un pour le cône).
    fn gap(&self, p: [f64; 3]) -> f64 {
        match self.kind {
            Kind::Plane { origin, normal, .. } => dot(sub(p, origin), normal).abs(),
            Kind::Round {
                origin,
                axis,
                radius,
                tan,
                ..
            } => {
                let v = dot(sub(p, origin), axis);
                (self.radial(p) - (radius + v * tan)).abs() / (1.0 + tan * tan).sqrt()
            }
            Kind::Sphere { center, radius, .. } => (norm(sub(p, center)) - radius).abs(),
            Kind::Torus {
                center,
                axis,
                major,
                minor,
                ..
            } => {
                let z = dot(sub(p, center), axis);
                ((self.radial(p) - major).hypot(z) - minor).abs()
            }
        }
    }

    /// Longueur sur la surface d'une unité de longueur du plan, au point x du plan.
    fn scale(&self, [x, y]: [f64; 2]) -> f64 {
        match self.kind {
            Kind::Plane { .. } | Kind::Round { .. } => 1.0,
            Kind::Sphere { radius, .. } => 2.0 / (1.0 + (x * x + y * y) / (radius * radius)),
            Kind::Torus { major, minor, .. } => {
                (major + minor * self.tube_angle(y / major).cos()) / major
            }
        }
    }

    /// Angle φ autour du tube d'un tore, de sa coordonnée ψ (à 2π près).
    fn tube_angle(&self, psi: f64) -> f64 {
        let Kind::Torus {
            major, minor, k, ..
        } = self.kind
        else {
            return 0.0;
        };
        let period = TAU * minor / k;
        let psi = psi - period * (psi / period).round();
        2.0 * (((major + minor) / (major - minor)).sqrt() * (k * psi / (2.0 * minor)).tan()).atan()
    }

    /// Déroule les contours, donnés par `coords` : chaque coordonnée périodique (angle) est rendue
    /// continue le long de chaque contour (une couture apparaît deux fois, à θ et à θ + 2π), les
    /// trous sont ramenés dans l'intervalle du contour extérieur, puis mise à plat et miroir si le
    /// contour extérieur tourne en sens horaire. Refuse un contour qui fait le tour d'une période :
    /// la face n'a pas de couture.
    fn unroll(
        &mut self,
        rings: &mut [Vec<([f64; 2], f64)>],
        name: &str,
    ) -> Result<(), MailleError> {
        let mut range = [[0.0; 2]; 2];
        for (c, period) in self.periods().into_iter().enumerate() {
            let Some(period) = period else {
                continue;
            };
            for (l, ring) in rings.iter_mut().enumerate() {
                for i in 1..ring.len() {
                    ring[i].0[c] =
                        ring[i - 1].0[c] + wrap_by(ring[i].0[c] - ring[i - 1].0[c], period);
                }
                let (first, last) = (ring[0].0[c], ring[ring.len() - 1].0[c]);
                if (last + wrap_by(first - last, period) - first).abs() > period / 2.0 {
                    let around = if c == 0 {
                        "de l'axe"
                    } else {
                        "du tube du tore"
                    };
                    return Err(MailleError::InvalidBrep {
                        problem: format!("face « {name} » : un contour fait le tour {around} (face fermée sans couture), pas encore pris en charge"),
                    });
                }
                let values = ring.iter().map(|p| p.0[c]);
                if l == 0 {
                    range[c] = [
                        values.clone().fold(f64::INFINITY, f64::min),
                        values.fold(f64::NEG_INFINITY, f64::max),
                    ];
                } else {
                    let shift = period * ((range[c][0] - first) / period).ceil();
                    for p in ring.iter_mut() {
                        p.0[c] += shift;
                    }
                }
            }
        }
        if let Kind::Round { sin, .. } = self.kind {
            self.shift = (range[0][0] + range[0][1]) / 2.0 * sin;
        }
        for p in rings.iter_mut().flatten() {
            p.0 = self.flat(p.0);
        }
        let outer = &rings[0];
        let area: f64 = (0..outer.len())
            .map(|i| {
                let (a, b) = (outer[i].0, outer[(i + 1) % outer.len()].0);
                a[0] * b[1] - a[1] * b[0]
            })
            .sum();
        if area < 0.0 {
            self.mirror = true;
            for p in rings.iter_mut().flatten() {
                p.0[0] = -p.0[0];
            }
        }
        Ok(())
    }

    /// Point du plan déroulé (sans le miroir) des coordonnées de `coords`, angles rendus continus.
    fn flat(&self, [a, b]: [f64; 2]) -> [f64; 2] {
        match self.kind {
            Kind::Plane { .. } | Kind::Sphere { .. } => [a, b],
            Kind::Round {
                radius, sin: 0.0, ..
            } => [radius * a, b],
            Kind::Round {
                radius, tan, sin, ..
            } => {
                let s = (radius + b * tan) / sin;
                let phi = a * sin - self.shift;
                [s * phi.cos(), s * phi.sin()]
            }
            Kind::Torus { major, .. } => [major * a, major * b],
        }
    }

    /// Point de la surface d'un point du plan déroulé.
    fn lift(&self, [x, y]: [f64; 2]) -> [f64; 3] {
        let x = if self.mirror { -x } else { x };
        match self.kind {
            Kind::Plane { origin, e1, e2, .. } => {
                std::array::from_fn(|k| origin[k] + x * e1[k] + y * e2[k])
            }
            Kind::Round {
                origin,
                u,
                w,
                axis,
                radius,
                tan,
                sin,
            } => {
                let (theta, rho, v) = if sin == 0.0 {
                    (x / radius, radius, y)
                } else {
                    let rho = x.hypot(y) * sin;
                    ((y.atan2(x) + self.shift) / sin, rho, (rho - radius) / tan)
                };
                let (c, s) = (rho * theta.cos(), rho * theta.sin());
                std::array::from_fn(|k| origin[k] + c * u[k] + s * w[k] + v * axis[k])
            }
            Kind::Sphere {
                center,
                radius,
                n,
                e1,
                e2,
                ..
            } => {
                let (a, b) = (x / radius, y / radius);
                let s = a * a + b * b;
                std::array::from_fn(|k| {
                    center[k]
                        + radius * (2.0 * a * e1[k] + 2.0 * b * e2[k] + (1.0 - s) * n[k])
                            / (1.0 + s)
                })
            }
            Kind::Torus {
                center,
                u,
                w,
                axis,
                major,
                minor,
                ..
            } => {
                let theta = x / major;
                let phi = self.tube_angle(y / major);
                let rho = major + minor * phi.cos();
                let (c, s, z) = (rho * theta.cos(), rho * theta.sin(), minor * phi.sin());
                std::array::from_fn(|k| center[k] + c * u[k] + s * w[k] + z * axis[k])
            }
        }
    }
}

/// Angle ramené dans ]−π, π].
fn wrap(a: f64) -> f64 {
    wrap_by(a, TAU)
}

/// Écart ramené dans ]−période / 2, période / 2].
fn wrap_by(a: f64, period: f64) -> f64 {
    let a = a.rem_euclid(period);
    if a > period / 2.0 {
        a - period
    } else {
        a
    }
}

/// Diagonale de la boîte des points (1 s'ils sont confondus).
fn diagonal(points: &[[f64; 3]]) -> f64 {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let d = norm(sub(hi, lo));
    if d.is_finite() && d > 0.0 {
        d
    } else {
        1.0
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::mesher::{BSpline, BrepFace};

    fn line(brep: &mut Brep, a: usize, b: usize) -> usize {
        let (from, to) = (brep.vertices[a], brep.vertices[b]);
        brep.edges.push(BrepEdge {
            start: a,
            end: b,
            curve: EdgeCurve::Line { from, to },
        });
        brep.edges.len() - 1
    }

    /// Cercle fermé de rayon r autour de l'axe z, à la hauteur z, depuis son sommet en (r, 0, z),
    /// antihoraire vu de dessus.
    fn circle(brep: &mut Brep, v: usize, r: f64, z: f64) -> usize {
        brep.edges.push(BrepEdge {
            start: v,
            end: v,
            curve: EdgeCurve::Arc {
                center: [0.0, 0.0, z],
                u: [1.0, 0.0, 0.0],
                v: [0.0, 1.0, 0.0],
                radii: [r, r],
                start_angle: 0.0,
                sweep: TAU,
            },
        });
        brep.edges.len() - 1
    }

    fn face(brep: &mut Brep, name: &str, surface: FaceSurface, loops: Vec<Vec<(usize, bool)>>) {
        brep.faces.push(BrepFace {
            name: name.to_string(),
            surface,
            loops,
        });
    }

    fn plane(z: f64, up: bool) -> FaceSurface {
        FaceSurface::Plane {
            origin: [0.0, 0.0, z],
            normal: [0.0, 0.0, if up { 1.0 } else { -1.0 }],
        }
    }

    /// Pavé [lo, hi], faces `x0`, `x1`, `y0`, `y1`, `z0`, `z1` ; arêtes des faces z0 et z1 dans
    /// `extra` en plus de leur contour (des trous).
    fn block(lo: [f64; 3], hi: [f64; 3], extra: [Vec<Vec<(usize, bool)>>; 2], brep: &mut Brep) {
        let first = brep.vertices.len();
        for k in 0..8 {
            brep.vertices.push(std::array::from_fn(|j| {
                if (k >> j) & 1 == 1 {
                    hi[j]
                } else {
                    lo[j]
                }
            }));
        }
        let mut edges = HashMap::new();
        for i in 0..8usize {
            for bit in [1, 2, 4] {
                if i & bit == 0 {
                    edges.insert((i, i | bit), line(brep, first + i, first + (i | bit)));
                }
            }
        }
        // cycles antihoraires vus de dehors
        let cycles = [
            ("x0", [0, 4, 6, 2], 0, false),
            ("x1", [1, 3, 7, 5], 0, true),
            ("y0", [0, 1, 5, 4], 1, false),
            ("y1", [2, 6, 7, 3], 1, true),
            ("z0", [0, 2, 3, 1], 2, false),
            ("z1", [4, 5, 7, 6], 2, true),
        ];
        let [mut holes0, mut holes1] = extra;
        for (name, cycle, axis, up) in cycles {
            let outer: Vec<(usize, bool)> = (0..4)
                .map(|k| {
                    let (p, q) = (cycle[k], cycle[(k + 1) % 4]);
                    (edges[&(p.min(q), p.max(q))], p > q)
                })
                .collect();
            let mut normal = [0.0; 3];
            normal[axis] = if up { 1.0 } else { -1.0 };
            let origin = if up { hi } else { lo };
            let mut loops = vec![outer];
            if name == "z0" {
                loops.append(&mut holes0);
            } else if name == "z1" {
                loops.append(&mut holes1);
            }
            face(brep, name, FaceSurface::Plane { origin, normal }, loops);
        }
    }

    fn empty() -> Brep {
        Brep {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
        }
    }

    /// Cylindre de rayon r de z = 0 à z = h, ou tronc de cône du rayon r à `top` : faces `cote`
    /// (couture comme chez OCCT), `fond` et `haut`.
    pub(in crate::mesher) fn round(r: f64, top: f64, h: f64) -> Brep {
        let mut b = empty();
        b.vertices = vec![[r, 0.0, 0.0], [top, 0.0, h]];
        let bottom = circle(&mut b, 0, r, 0.0);
        let upper = circle(&mut b, 1, top, h);
        let seam = line(&mut b, 0, 1);
        let surface = if r == top {
            FaceSurface::Cylinder {
                origin: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: r,
            }
        } else {
            FaceSurface::Cone {
                origin: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: r,
                half_angle: ((top - r) / h).atan(),
            }
        };
        face(
            &mut b,
            "cote",
            surface,
            vec![vec![
                (bottom, false),
                (seam, false),
                (upper, true),
                (seam, true),
            ]],
        );
        face(
            &mut b,
            "fond",
            plane(0.0, false),
            vec![vec![(bottom, true)]],
        );
        face(&mut b, "haut", plane(h, true), vec![vec![(upper, false)]]);
        b
    }

    /// Pavé [−2, 2]² × [0, 2] percé d'un trou de rayon 1 selon z : faces du pavé et `trou`.
    fn holed_block() -> Brep {
        let mut b = empty();
        b.vertices = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 2.0]];
        let bottom = circle(&mut b, 0, 1.0, 0.0);
        let upper = circle(&mut b, 1, 1.0, 2.0);
        let seam = line(&mut b, 0, 1);
        // vu de dehors, le trou tourne en sens horaire autour de la face
        block(
            [-2.0, -2.0, 0.0],
            [2.0, 2.0, 2.0],
            [vec![vec![(bottom, false)]], vec![vec![(upper, true)]]],
            &mut b,
        );
        face(
            &mut b,
            "trou",
            FaceSurface::Cylinder {
                origin: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: 1.0,
            },
            vec![vec![
                (seam, false),
                (upper, false),
                (seam, true),
                (bottom, true),
            ]],
        );
        b
    }

    fn triangle_area(s: &Surface, t: [usize; 3]) -> f64 {
        let [a, b, c] = t.map(|v| s.vertices[v]);
        0.5 * norm(cross(sub(b, a), sub(c, a)))
    }

    /// Surface maillée : fermée (volume de signe positif), sommets de chaque Patch sur la surface
    /// de sa face, rayon circonscrit au plus h/√3.
    fn check(brep: &Brep, h: f64) -> Surface {
        check_within(brep, h, 1.0 + 1e-9)
    }

    /// Comme `check`, le rayon circonscrit borné par `slack` h/√3 : sur une carte qui déforme
    /// les longueurs (sphère, tore), l'échelle varie un peu sur un triangle.
    fn check_within(brep: &Brep, h: f64, slack: f64) -> Surface {
        let s = triangulate(brep, h, 0.3).unwrap();
        super::super::tetrahedralize::check(&s).unwrap();
        // conforme : chaque arête bordée par deux triangles, une fois dans chaque sens (pas de
        // sommet posé au milieu d'une arête d'en face)
        let mut sides: HashMap<(usize, usize), i32> = HashMap::new();
        for t in &s.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *sides.entry((a.min(b), a.max(b))).or_default() += if a < b { 1 } else { -1 };
            }
        }
        let mut uses: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &s.triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        assert!(
            uses.values().all(|&n| n == 2),
            "arête bordée une seule fois"
        );
        assert!(sides.values().all(|&n| n == 0), "orientation");
        for p in &s.patches {
            let face = brep.faces.iter().find(|f| f.name == p.name).unwrap();
            let chart = Chart::new(&face.surface, &face.name).unwrap();
            for &t in &p.triangles {
                for v in s.triangles[t] {
                    let gap = chart.gap(s.vertices[v]);
                    assert!(gap < 1e-12, "{} : écart {gap:e}", p.name);
                }
                let [a, b, c] = s.triangles[t].map(|v| s.vertices[v]);
                let l = [norm(sub(b, a)), norm(sub(c, b)), norm(sub(a, c))];
                let radius = l[0] * l[1] * l[2] / (4.0 * triangle_area(&s, s.triangles[t]));
                assert!(radius <= h / 3f64.sqrt() * slack, "rayon {radius}");
            }
        }
        s
    }

    fn patch_area(s: &Surface, name: &str) -> f64 {
        let p = s.patches.iter().find(|p| p.name == name).unwrap();
        p.triangles
            .iter()
            .map(|&t| triangle_area(s, s.triangles[t]))
            .sum()
    }

    #[test]
    fn pave_exact() {
        let mut b = empty();
        block([0.0; 3], [1.0, 2.0, 3.0], [vec![], vec![]], &mut b);
        let s = check(&b, 0.4);
        assert!((s.volume() - 6.0).abs() < 1e-12);
        assert!((patch_area(&s, "x1") - 6.0).abs() < 1e-12);
    }

    #[test]
    fn cylindre_a_couture_converge() {
        // sommets sur le cylindre : surface inscrite, écart en h² / r²
        let exact = PI * 2.0;
        let mut errors = Vec::new();
        for h in [0.4, 0.2] {
            let s = check(&round(1.0, 1.0, 2.0), h);
            let e = (exact - s.volume()) / exact;
            assert!(e > 0.0 && e < 0.5 * h * h, "h = {h} : écart {e:e}");
            errors.push(e);
        }
        assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
    }

    #[test]
    fn couture_ajoutee_a_une_face_fermee() {
        // la face latérale bordée de ses deux cercles seulement : même maillage qu'avec la couture
        let mut b = round(1.0, 1.0, 2.0);
        b.faces[0].loops = vec![vec![(0, false)], vec![(1, true)]];
        b.edges.pop();
        let s = check(&b, 0.3);
        let with_seam = check(&round(1.0, 1.0, 2.0), 0.3);
        assert!((s.volume() - with_seam.volume()).abs() < 1e-3 * TAU);
        // cercle du haut qui commence un quart de tour plus loin : il est coupé à l'angle 0
        b.vertices[1] = [0.0, 1.0, 2.0];
        if let EdgeCurve::Arc { start_angle, .. } = &mut b.edges[1].curve {
            *start_angle = PI / 2.0;
        }
        let s = check(&b, 0.3);
        let e = (TAU - s.volume()) / TAU;
        assert!(e > 0.0 && e < 0.5 * 0.09, "écart {e:e}");
    }

    #[test]
    fn seize_morceaux_par_tour_au_moins() {
        // h bien plus grand que le rayon : les cercles gardent 16 morceaux
        let s = check(&round(1.0, 1.0, 2.0), 10.0);
        let on_bottom = s.vertices.iter().filter(|p| p[2].abs() < 1e-12).count();
        assert!(on_bottom >= 16, "{on_bottom}");
        let e = (PI * 2.0 - s.volume()) / (PI * 2.0);
        // flèche d'un seizième de tour : aire du polygone à 16 côtés
        let polygon = 8.0 * (TAU / 16.0).sin() / PI;
        assert!(e < 1.0 - polygon + 1e-9, "écart {e:e}");
    }

    #[test]
    fn tronc_de_cone() {
        let (r0, r1, h) = (1.0, 2.0, 1.5);
        let exact = PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
        let lateral = PI * (r0 + r1) * (r1 - r0).hypot(h);
        let s = check(&round(r0, r1, h), 0.2);
        let e = (exact - s.volume()) / exact;
        assert!(e > 0.0 && e < 0.02, "volume : écart {e:e}");
        let e = (lateral - patch_area(&s, "cote")) / lateral;
        assert!(e > 0.0 && e < 0.02, "aire : écart {e:e}");
    }

    #[test]
    fn couture_coupee_d_un_seul_cote() {
        // tronc de cône très ouvert : le mailleur 2D coupe la couture d'un seul de ses deux côtés,
        // la face est remaillée pour que l'autre côté ait le point (sinon, sommet au milieu d'une
        // arête d'en face)
        let (r0, r1, h) = (1.0, 3.0, 0.5);
        let s = check(&round(r0, r1, h), 0.05);
        let exact = PI * h / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1);
        assert!((exact - s.volume()) / exact < 1e-3);
    }

    #[test]
    fn pave_perce() {
        let s = check(&holed_block(), 0.3);
        let exact = 32.0 - 2.0 * PI;
        let e = (s.volume() - exact) / exact;
        // le trou inscrit enlève un peu moins que 2π
        assert!(e > 0.0 && e < 1e-2, "écart {e:e}");
        let names: Vec<&str> = s.patches.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["x0", "x1", "y0", "y1", "z0", "z1", "trou"]);
    }

    /// Tronc de cône (rayon 2 en bas, 1 en haut, hauteur 1) coupé par le plan x = 1, tangent au
    /// cercle du haut en V = (1, 0, 1), comme le chanfrein d'un écrou par un de ses pans : le cercle
    /// et l'hyperbole de coupe s'y rejoignent tangents, en pointe, sur la face conique et sur le pan.
    /// Sans les groupes de `Mesher` (une arête découpée part du coin d'un bloc), les coupes vont et
    /// viennent entre les faces vers la pointe jusqu'à des points confondus.
    fn chamfer() -> Brep {
        let mut b = empty();
        let s3 = 3f64.sqrt();
        b.vertices = vec![[1.0, 0.0, 1.0], [1.0, s3, 0.0], [1.0, -s3, 0.0]];
        let top = circle(&mut b, 0, 1.0, 1.0);
        // demi-hyperbole (2 − z)² − y² = 1 de V au bas : conique exacte, poids cosh(T/2) = √1,5
        let mut hyperbola = |end: usize, sign: f64| {
            let curve = BSpline::new(
                2,
                vec![[1.0, 0.0, 1.0], [1.0, sign / s3, 1.0], b.vertices[end]],
                vec![1.0, 1.5f64.sqrt(), 1.0],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            )
            .unwrap();
            b.edges.push(BrepEdge {
                start: 0,
                end,
                curve: EdgeCurve::BSpline {
                    curve,
                    range: [0.0, 1.0],
                },
            });
            b.edges.len() - 1
        };
        let up = hyperbola(1, 1.0);
        let down = hyperbola(2, -1.0);
        b.edges.push(BrepEdge {
            start: 1,
            end: 2,
            curve: EdgeCurve::Arc {
                center: [0.0; 3],
                u: [1.0, 0.0, 0.0],
                v: [0.0, 1.0, 0.0],
                radii: [2.0, 2.0],
                start_angle: PI / 3.0,
                sweep: 4.0 * PI / 3.0,
            },
        });
        let arc = b.edges.len() - 1;
        let chord = line(&mut b, 2, 1);
        face(&mut b, "haut", plane(1.0, true), vec![vec![(top, false)]]);
        face(
            &mut b,
            "bas",
            plane(0.0, false),
            vec![vec![(chord, true), (arc, true)]],
        );
        face(
            &mut b,
            "pan",
            FaceSurface::Plane {
                origin: [1.0, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
            },
            vec![vec![(chord, false), (up, true), (down, false)]],
        );
        face(
            &mut b,
            "chanfrein",
            FaceSurface::Cone {
                origin: [0.0; 3],
                axis: [0.0, 0.0, -1.0],
                radius: 2.0,
                half_angle: PI / 4.0,
            },
            vec![vec![(top, true), (up, false), (arc, false), (down, true)]],
        );
        b
    }

    #[test]
    fn pointe_d_un_chanfrein() {
        // volume : tronc de cône 7π/3 moins, pour chaque rayon ρ de 1 à 2, le segment de disque
        // au-delà de x = 1, ρ² acos(1/ρ) − √(ρ² − 1) (Simpson)
        let n = 20_000;
        let segment = |r: f64| r * r * (1.0 / r).acos() - (r * r - 1.0).max(0.0).sqrt();
        let cut: f64 = (0..=n)
            .map(|i| {
                let w = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * segment(1.0 + i as f64 / n as f64)
            })
            .sum::<f64>()
            / (3.0 * n as f64);
        let exact = 7.0 * PI / 3.0 - cut;
        let mut errors = Vec::new();
        for h in [0.2, 0.1] {
            let s = check(&chamfer(), h);
            let e = (exact - s.volume()) / exact;
            assert!(e > 0.0 && e < 0.5 * h * h, "h = {h} : écart {e:e}");
            errors.push(e);
        }
        assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
        let options = crate::mesher::TetOptions {
            size: Some(0.2),
            ..Default::default()
        };
        crate::mesher::tetrahedralize_brep(&chamfer(), &options).unwrap();
    }

    fn arc(
        b: &mut Brep,
        ends: [usize; 2],
        center: [f64; 3],
        uv: [[f64; 3]; 2],
        r: f64,
        sweep: f64,
    ) -> usize {
        b.edges.push(BrepEdge {
            start: ends[0],
            end: ends[1],
            curve: EdgeCurve::Arc {
                center,
                u: uv[0],
                v: uv[1],
                radii: [r, r],
                start_angle: 0.0,
                sweep,
            },
        });
        b.edges.len() - 1
    }

    /// Sphère de rayon r comme l'écrit OCCT : une face, une couture (demi-méridien d'un pôle à
    /// l'autre) parcourue dans les deux sens, les pôles en contours réduits à un sommet (omis).
    fn sphere(r: f64) -> Brep {
        let mut b = empty();
        b.vertices = vec![[0.0, 0.0, -r], [0.0, 0.0, r]];
        b.edges.push(BrepEdge {
            start: 0,
            end: 1,
            curve: EdgeCurve::Arc {
                center: [0.0; 3],
                u: [1.0, 0.0, 0.0],
                v: [0.0, 0.0, 1.0],
                radii: [r, r],
                start_angle: -PI / 2.0,
                sweep: PI,
            },
        });
        let sphere = FaceSurface::Sphere {
            center: [0.0; 3],
            radius: r,
            convex: true,
        };
        face(&mut b, "boule", sphere, vec![vec![(0, false), (0, true)]]);
        b
    }

    /// Huitième de coque sphérique entre les rayons `ri` et `ro` : faces `interieur`, `exterieur`
    /// et les plans `x0`, `y0`, `z0`.
    pub(in crate::mesher) fn octant(ri: f64, ro: f64) -> Brep {
        let mut b = empty();
        let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        for r in [ri, ro] {
            b.vertices
                .extend([[r, 0.0, 0.0], [0.0, r, 0.0], [0.0, 0.0, r]]);
        }
        let mut arcs = [[0; 3]; 2];
        for (k, r) in [ri, ro].into_iter().enumerate() {
            let [vx, vy, vz] = [3 * k, 3 * k + 1, 3 * k + 2];
            arcs[k] = [
                arc(&mut b, [vx, vy], [0.0; 3], [x, y], r, PI / 2.0),
                arc(&mut b, [vy, vz], [0.0; 3], [y, z], r, PI / 2.0),
                arc(&mut b, [vz, vx], [0.0; 3], [z, x], r, PI / 2.0),
            ];
        }
        let [lx, ly, lz] = [0, 1, 2].map(|k| line(&mut b, k, k + 3));
        let [[xy_i, yz_i, zx_i], [xy_o, yz_o, zx_o]] = arcs;
        let sphere = |radius: f64, convex: bool| FaceSurface::Sphere {
            center: [0.0; 3],
            radius,
            convex,
        };
        let plane = |normal: [f64; 3]| FaceSurface::Plane {
            origin: [0.0; 3],
            normal,
        };
        face(
            &mut b,
            "exterieur",
            sphere(ro, true),
            vec![vec![(xy_o, false), (yz_o, false), (zx_o, false)]],
        );
        face(
            &mut b,
            "interieur",
            sphere(ri, false),
            vec![vec![(zx_i, true), (yz_i, true), (xy_i, true)]],
        );
        face(
            &mut b,
            "z0",
            plane([0.0, 0.0, -1.0]),
            vec![vec![(xy_i, false), (ly, false), (xy_o, true), (lx, true)]],
        );
        face(
            &mut b,
            "x0",
            plane([-1.0, 0.0, 0.0]),
            vec![vec![(yz_i, false), (lz, false), (yz_o, true), (ly, true)]],
        );
        face(
            &mut b,
            "y0",
            plane([0.0, -1.0, 0.0]),
            vec![vec![(zx_i, false), (lx, false), (zx_o, true), (lz, true)]],
        );
        b
    }

    /// Tore complet comme l'écrit OCCT : une face, deux coutures (un parallèle et un méridien, par
    /// le même sommet) parcourues chacune dans les deux sens.
    fn torus(major: f64, minor: f64) -> Brep {
        let mut b = empty();
        b.vertices = vec![[major + minor, 0.0, 0.0]];
        let parallel = arc(
            &mut b,
            [0, 0],
            [0.0; 3],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            major + minor,
            TAU,
        );
        let meridian = arc(
            &mut b,
            [0, 0],
            [major, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            minor,
            TAU,
        );
        let torus = FaceSurface::Torus {
            center: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            major,
            minor,
        };
        face(
            &mut b,
            "tore",
            torus,
            vec![vec![
                (parallel, false),
                (meridian, false),
                (parallel, true),
                (meridian, true),
            ]],
        );
        b
    }

    /// Volume inscrit qui converge à l'ordre 2 vers `exact`, aux tailles `sizes`.
    fn converges(brep: &Brep, exact: f64, sizes: [f64; 2]) -> [f64; 2] {
        let errors = sizes.map(|h| {
            let s = check_within(brep, h, 1.15);
            let e = (exact - s.volume()) / exact;
            assert!(
                e > 0.0,
                "h = {h} : volume {} au-delà de {exact}",
                s.volume()
            );
            e
        });
        assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
        errors
    }

    #[test]
    fn sphere_entiere_coupee_en_deux() {
        let r = 1.5;
        let errors = converges(&sphere(r), 4.0 / 3.0 * PI * r.powi(3), [0.4, 0.2]);
        assert!(errors[1] < 0.5 * (0.2f64 / r).powi(2), "{errors:?}");
        let s = check_within(&sphere(r), 0.4, 1.15);
        assert_eq!(s.patches.len(), 1);
    }

    #[test]
    fn coque_spherique() {
        let (ri, ro) = (1.0, 2.0);
        converges(
            &octant(ri, ro),
            PI / 6.0 * (ro.powi(3) - ri.powi(3)),
            [0.3, 0.15],
        );
    }

    #[test]
    fn bout_arrondi() {
        // cylindre de rayon 1 et de hauteur 2 coiffé d'une demi-sphère : face sphérique bordée par
        // un seul cercle, la couture et le pôle d'OCCT en moins
        let mut b = round(1.0, 1.0, 2.0);
        let top = 1;
        b.faces[2] = BrepFace {
            name: "dome".to_string(),
            surface: FaceSurface::Sphere {
                center: [0.0, 0.0, 2.0],
                radius: 1.0,
                convex: true,
            },
            loops: vec![vec![(top, false)]],
        };
        converges(&b, 2.0 * PI + 2.0 / 3.0 * PI, [0.3, 0.15]);
    }

    #[test]
    fn tore_entier() {
        let (major, minor) = (3.0, 1.0);
        let exact = 2.0 * PI * PI * major * minor * minor;
        converges(&torus(major, minor), exact, [0.4, 0.2]);
        // sans arête, comme OCCT, et l'axe à l'envers : mêmes coutures, ajoutées
        let mut b = torus(major, minor);
        b.edges.clear();
        b.vertices.clear();
        b.faces[0].loops.clear();
        b.faces[0].surface = FaceSurface::Torus {
            center: [0.0; 3],
            axis: [0.0, 0.0, -1.0],
            major,
            minor,
        };
        converges(&b, exact, [0.4, 0.2]);
    }

    #[test]
    fn faces_mal_posees_refusees() {
        // sommet hors de sa surface : cylindre de rayon 1,1 sous des cercles de rayon 1
        let mut b = round(1.0, 1.0, 2.0);
        b.faces[0].surface = FaceSurface::Cylinder {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            radius: 1.1,
        };
        let m = triangulate(&b, 0.5, 0.3).unwrap_err().to_string();
        assert!(m.contains("face « cote » : un sommet est à"), "{m}");
        // cône pointu : la couture va du cercle de base au sommet
        let mut b = empty();
        b.vertices = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 2.0]];
        let base = circle(&mut b, 0, 1.0, 0.0);
        let seam = line(&mut b, 0, 1);
        let cone = FaceSurface::Cone {
            origin: [0.0, 0.0, 2.0],
            axis: [0.0, 0.0, -1.0],
            radius: 0.0,
            half_angle: 0.5f64.atan(),
        };
        face(
            &mut b,
            "cote",
            cone,
            vec![vec![(base, false), (seam, false), (seam, true)]],
        );
        face(&mut b, "fond", plane(0.0, false), vec![vec![(base, true)]]);
        let m = triangulate(&b, 0.5, 0.3).unwrap_err().to_string();
        assert!(m.contains("cône jusqu'à son sommet"), "{m}");
        // demi-angle hors de ]0, π/2[
        let mut b = round(1.0, 2.0, 1.0);
        b.faces[0].surface = FaceSurface::Cone {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            radius: 1.0,
            half_angle: 2.0,
        };
        assert!(matches!(
            triangulate(&b, 0.5, 0.3),
            Err(MailleError::InvalidBrep { .. })
        ));
        // tore croisé
        let mut b = torus(1.0, 0.5);
        b.faces[0].surface = FaceSurface::Torus {
            center: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            major: 1.0,
            minor: 1.0,
        };
        let m = triangulate(&b, 0.5, 0.3).unwrap_err().to_string();
        assert!(m.contains("tore croisé"), "{m}");
    }
}
