//! Remaillage de la surface à la taille h, deuxième étape du mailleur 3D à taille h. Chaque segment
//! (frontière entre deux facettes, voir `facet`) est découpé en morceaux égaux de longueur au plus
//! h, puis chaque facette est maillée dans son plan par le mailleur 2D (bord respecté, Ruppert avec
//! gradation). Quand le mailleur 2D coupe un morceau de bord, le point rejoint le segment et la
//! facette d'en face est remaillée avec lui : la surface reste fermée. Un point de segment est
//! calculé une seule fois, en 3D sur son morceau, pour être le même des deux côtés. Les angles aigus
//! d'une facette (moins de `MIN_CORNER`) sont acceptés : les triangles qui y touchent, ou dont le
//! plus petit côté relie les deux bords du coin, ne sont raffinés que pour leur taille et gardent un
//! angle proche de celui du coin. Mesuré : un prisme dont le coin fait 6° passe, à 5° les coupes ne
//! se stabilisent plus entre les deux faces du dièdre (`InvalidSurface`).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::f64::consts::PI;

use super::facet::{facets, Facets};
use super::geometry::{Boundary, MIN_CORNER};
use super::refine::{Mesher, SizeField};
use super::surface::{cross, dot, norm, sub, Patch, Surface};
use super::{MeshOptions, MAX_NODES};
use crate::error::MailleError;
use crate::material::positive;

/// Maillages de facette au plus, rapportés au nombre de facettes : au-delà, les coupes ne se
/// stabilisent pas d'une facette à l'autre.
const MAX_REMESH: usize = 20;
const STUCK: &str = "la surface ne se laisse pas remailler à la taille h : un coin ou un dièdre de quelques degrés (moins de 6° environ) entre deux faces planes fait couper sans fin ; arrondir ou épaissir ce coin";

/// Surface remaillée à la taille `size` avec la gradation `grading` : mêmes facettes planes, mêmes
/// Patches, triangles de rayon circonscrit au plus h/√3 (celui du triangle équilatéral de côté h,
/// comme en 2D) et d'angle au moins celui de `MeshOptions` (hors coins aigus). `InvalidProperty`
/// si h ou la gradation n'est pas strictement positive, `InvalidSurface` si les coupes ne se
/// stabilisent pas.
pub(crate) fn remesh(surface: &Surface, size: f64, grading: f64) -> Result<Surface, MailleError> {
    let size = positive("taille h", size)?;
    let grading = positive("gradation", grading)?;
    let fs = facets(surface)?;
    let mut points = Points::default();

    // points de chaque segment, du premier coin au second ; taille h de ses morceaux, plus petite
    // que h sur un segment court pour que la gradation parte des petits détails
    let mut chains = Vec::with_capacity(fs.segments.len());
    let mut sizes = Vec::with_capacity(fs.segments.len());
    for seg in &fs.segments {
        let [a, b] =
            [seg.vertices[0], seg.vertices[seg.vertices.len() - 1]].map(|v| surface.vertices[v]);
        let length = norm(sub(b, a));
        let n = (length / size).ceil();
        if n > MAX_NODES as f64 {
            return Err(MailleError::TooManyNodes { max: MAX_NODES });
        }
        let n = (n as usize).max(1);
        let chain = (0..=n)
            .map(|i| match i {
                0 => points.add(a),
                i if i == n => points.add(b),
                i => points.add(lerp(a, b, i as f64 / n as f64)),
            })
            .collect::<Result<Vec<_>, _>>()?;
        chains.push(chain);
        sizes.push(size.min(length));
    }

    let mut meshes: Vec<Vec<[usize; 3]>> = vec![Vec::new(); fs.facets.len()];
    let mut pending: BTreeSet<usize> = (0..fs.facets.len()).collect();
    let mut budget = MAX_REMESH * fs.facets.len();
    let stuck = || MailleError::InvalidSurface {
        problem: STUCK.to_string(),
    };
    while let Some(f) = pending.pop_first() {
        budget = budget.checked_sub(1).ok_or_else(stuck)?;
        let (triangles, origins) = mesh_facet(&fs, f, &chains, &sizes, &points.coords, grading)
            .map_err(|e| match e {
                MailleError::MeshingFailed => stuck(),
                e => e,
            })?;
        let mut splits: BTreeMap<(usize, usize), Vec<(f64, usize)>> = BTreeMap::new();
        let mut global = Vec::with_capacity(origins.len());
        for origin in origins {
            global.push(match origin {
                Origin::Known(i) => i,
                Origin::Interior(p) => points.add(p)?,
                Origin::Split { segment, k, t } => {
                    let c = &chains[segment];
                    let p = lerp(points.coords[c[k]], points.coords[c[k + 1]], t);
                    let i = points.add(p)?;
                    splits.entry((segment, k)).or_default().push((t, i));
                    i
                }
            });
        }
        // du dernier morceau au premier, pour que les rangs k restent justes
        for ((segment, k), mut new) in splits.into_iter().rev() {
            new.sort_by(|a, b| a.0.total_cmp(&b.0));
            chains[segment].splice(k + 1..k + 1, new.into_iter().map(|(_, i)| i));
            pending.extend(fs.segments[segment].facets.into_iter().filter(|&g| g != f));
        }
        meshes[f] = triangles
            .into_iter()
            .map(|t| t.map(|v| global[v]))
            .collect();
    }

    // surface finale, sans les points des maillages remplacés
    let mut number = vec![usize::MAX; points.coords.len()];
    let mut out = Surface {
        vertices: Vec::new(),
        triangles: Vec::new(),
        patches: surface
            .patches
            .iter()
            .map(|p| Patch {
                name: p.name.clone(),
                triangles: Vec::new(),
            })
            .collect(),
    };
    for (f, triangles) in meshes.iter().enumerate() {
        for t in triangles {
            let t = t.map(|v| {
                if number[v] == usize::MAX {
                    number[v] = out.vertices.len();
                    out.vertices.push(points.coords[v]);
                }
                number[v]
            });
            out.patches[fs.facets[f].patch]
                .triangles
                .push(out.triangles.len());
            out.triangles.push(t);
        }
    }
    Ok(out)
}

/// Points de la surface remaillée, chacun une seule fois (repérés par leurs coordonnées exactes).
#[derive(Default)]
struct Points {
    coords: Vec<[f64; 3]>,
    index: HashMap<[u64; 3], usize>,
}

impl Points {
    fn add(&mut self, p: [f64; 3]) -> Result<usize, MailleError> {
        if let Some(&i) = self.index.get(&p.map(|x| (x + 0.0).to_bits())) {
            return Ok(i);
        }
        if self.coords.len() >= MAX_NODES {
            return Err(MailleError::TooManyNodes { max: MAX_NODES });
        }
        self.index
            .insert(p.map(|x| (x + 0.0).to_bits()), self.coords.len());
        self.coords.push(p);
        Ok(self.coords.len() - 1)
    }
}

/// D'où vient un point d'une facette maillée : d'un segment, d'une coupe du morceau `k` du segment
/// `segment` (paramètre t dans le sens du segment), ou de l'intérieur.
enum Origin {
    Known(usize),
    Split { segment: usize, k: usize, t: f64 },
    Interior([f64; 3]),
}

/// Maille la facette `f` dans son plan à partir des points actuels de ses segments ; renvoie les
/// triangles (indices locaux) et l'origine de chaque point local.
fn mesh_facet(
    fs: &Facets,
    f: usize,
    chains: &[Vec<usize>],
    sizes: &[f64],
    points: &[[f64; 3]],
    grading: f64,
) -> Result<(Vec<[usize; 3]>, Vec<Origin>), MailleError> {
    let facet = &fs.facets[f];
    let (e1, e2) = basis(facet.normal);
    let o = facet.origin;
    let flat = |p: [f64; 3]| {
        let d = sub(p, o);
        [dot(d, e1), dot(d, e2)]
    };
    // morceaux de bord dans le sens de la facette : (segment, rang du morceau, à rebours)
    let mut pieces = Vec::new();
    let mut known = Vec::new();
    let mut loops = Vec::new();
    for lp in &facet.loops {
        let mut ring = Vec::new();
        for &(s, reversed) in lp {
            let m = chains[s].len() - 1;
            for j in 0..m {
                let k = if reversed { m - 1 - j } else { j };
                let start = chains[s][if reversed { k + 1 } else { k }];
                ring.push((flat(points[start]), sizes[s]));
                pieces.push((s, k, reversed));
                known.push(start);
            }
        }
        loops.push(ring);
    }
    let boundary = Boundary::polygons(&loops, &format!("facette {f}"));
    let mut m = Mesher::new(&boundary)?;
    m.acute = acute_corners(&loops);
    let options = MeshOptions {
        grading,
        ..MeshOptions::default()
    };
    m.refine(&options, &SizeField::new(&boundary, grading))?;
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
                let (segment, k, reversed) = pieces[piece];
                let t = if reversed { 1.0 - t } else { t };
                Origin::Split { segment, k, t }
            } else {
                let [x, y] = m.d.point(i);
                Origin::Interior(std::array::from_fn(|j| o[j] + x * e1[j] + y * e2[j]))
            }
        })
        .collect();
    Ok((m.d.inside_triangles(), origins))
}

/// Repère orthonormé (e1, e2) du plan de normale unitaire n, avec e1 × e2 = n : un triangle
/// antihoraire dans le plan a la normale n.
fn basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let k = (0..3)
        .min_by(|&a, &b| n[a].abs().total_cmp(&n[b].abs()))
        .unwrap_or(0);
    let mut axis = [0.0; 3];
    axis[k] = 1.0;
    let e1 = cross(n, axis);
    let l = norm(e1);
    let e1 = e1.map(|x| x / l);
    (e1, cross(n, e1))
}

/// Points des contours (numérotés à la suite) où l'angle intérieur est sous `MIN_CORNER`.
fn acute_corners(loops: &[Vec<([f64; 2], f64)>]) -> HashSet<usize> {
    let mut out = HashSet::new();
    let mut first = 0;
    for ring in loops {
        let n = ring.len();
        for i in 0..n {
            let (p, q, r) = (ring[(i + n - 1) % n].0, ring[i].0, ring[(i + 1) % n].0);
            let (u, v) = ([q[0] - p[0], q[1] - p[1]], [r[0] - q[0], r[1] - q[1]]);
            let turn = (u[0] * v[1] - u[1] * v[0]).atan2(u[0] * v[0] + u[1] * v[1]);
            if (PI - turn).to_degrees() < MIN_CORNER {
                out.insert(first + i);
            }
        }
        first += n;
    }
    out
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    std::array::from_fn(|k| a[k] + t * (b[k] - a[k]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesher::testsurface::*;
    use crate::mesher::tetrahedralize;

    /// Plus petit angle (en degrés) et plus grand rayon circonscrit des triangles.
    fn quality(s: &Surface) -> (f64, f64) {
        let (mut angle, mut radius) = (180.0_f64, 0.0_f64);
        for t in &s.triangles {
            let p = t.map(|v| s.vertices[v]);
            let l = [0, 1, 2].map(|k| norm(sub(p[(k + 1) % 3], p[k])));
            let area = 0.5 * norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
            radius = radius.max(l[0] * l[1] * l[2] / (4.0 * area));
            for k in 0..3 {
                let (u, v) = (sub(p[(k + 1) % 3], p[k]), sub(p[(k + 2) % 3], p[k]));
                let a = (dot(u, v) / (norm(u) * norm(v))).acos().to_degrees();
                angle = angle.min(a);
            }
        }
        (angle, radius)
    }

    /// Aire et volume inchangés, surface fermée et orientée, facettes inchangées, rayon
    /// circonscrit au plus h/√3,
    /// angles au moins `min_angle`.
    fn check_remesh(s: &Surface, h: f64, min_angle: f64) -> Surface {
        let r = remesh(s, h, 0.3).unwrap();
        assert!((r.area() - s.area()).abs() < 1e-12 * s.area(), "aire");
        assert!(
            (r.volume() - s.volume()).abs() < 1e-12 * s.volume(),
            "volume"
        );
        assert_eq!(
            facets(&r).unwrap().facets.len(),
            facets(s).unwrap().facets.len()
        );
        let (angle, radius) = quality(&r);
        assert!(
            radius <= h / 3f64.sqrt() * (1.0 + 1e-9),
            "rayon {radius} > h/√3, h = {h}"
        );
        assert!(angle >= min_angle - 1e-9, "angle {angle}");
        r
    }

    #[test]
    fn cube_a_la_taille_h() {
        let s = surface(&cube(1));
        let r = check_remesh(&s, 0.25, 25.0);
        // au moins (1 / h)² triangles par face
        assert!(r.triangles.len() >= 6 * 16 * 2);
        let names: Vec<&str> = r.patches.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["p0", "p1", "p2", "p3", "p4", "p5"]);
        // chaque Patch garde l'aire de sa face
        for p in &r.patches {
            let area: f64 = p
                .triangles
                .iter()
                .map(|&t| {
                    let [a, b, c] = r.triangles[t].map(|v| r.vertices[v]);
                    0.5 * norm(cross(sub(b, a), sub(c, a)))
                })
                .sum();
            assert!((area - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn h_plus_grand_que_la_piece() {
        let s = surface(&cube(2));
        check_remesh(&s, 5.0, 25.0);
    }

    #[test]
    fn plaque_trouee() {
        let s = surface(&holed_plate());
        check_remesh(&s, 0.4, 25.0);
    }

    #[test]
    fn prisme_en_l() {
        let s = surface(&[prism(&L_SHAPE)]);
        check_remesh(&s, 0.3, 25.0);
    }

    #[test]
    fn cylindre_a_bandes_etroites() {
        // bandes de largeur 2 sin(π / 64) ≈ 0,098 < h : la finesse est celle des bandes
        let s = surface(&[cylinder(64)]);
        check_remesh(&s, 0.3, 25.0);
    }

    #[test]
    fn coin_aigu_accepte() {
        // prisme sur un triangle d'angle 20° en (0, 0)
        let a = 20f64.to_radians();
        let s = surface(&[prism(&[[0.0, 0.0], [1.0, 0.0], [a.cos(), a.sin()]])]);
        // les triangles du coin gardent un angle proche du sien (17,9° mesuré)
        check_remesh(&s, 0.2, 15.0);
    }

    #[test]
    fn coin_de_6_degres() {
        // les triangles du coin gardent son petit angle : on contrôle la surface seule
        let a = 6f64.to_radians();
        let s = surface(&[prism(&[[0.0, 0.0], [1.0, 0.0], [a.cos(), a.sin()]])]);
        for h in [0.03, 0.3] {
            let r = remesh(&s, h, 0.3).unwrap();
            assert!((r.volume() - s.volume()).abs() < 1e-12 * s.volume());
            assert_eq!(facets(&r).unwrap().facets.len(), 5);
        }
    }

    #[test]
    fn surface_remaillee_se_laisse_tetraedriser() {
        let s = surface(&cube(1));
        let r = remesh(&s, 0.25, 0.3).unwrap();
        let t = tetrahedralize(&r).unwrap();
        let volume: f64 = t
            .tetrahedra
            .iter()
            .map(|e| {
                let [a, b, c, d] = e.map(|i| {
                    let n = t.nodes[i];
                    [n.x, n.y, n.z]
                });
                dot(sub(b, a), cross(sub(c, a), sub(d, a))) / 6.0
            })
            .sum();
        assert!((volume - 1.0).abs() < 1e-12);
    }

    #[test]
    fn taille_invalide_refusee() {
        let s = surface(&cube(1));
        for h in [0.0, -1.0, f64::NAN] {
            assert!(matches!(
                remesh(&s, h, 0.3),
                Err(MailleError::InvalidProperty { .. })
            ));
        }
    }
}
