//! Facettes planes d'une `Surface`, première étape du remaillage de la surface à la taille h : les
//! triangles voisins, du même Patch et dans le même plan, forment une facette ; les frontières
//! entre facettes sont coupées en segments d'un coin à l'autre. Un coin est un sommet où se
//! croisent plus de deux frontières, ou où la frontière tourne. Les sommets qui ne sont ni sur une
//! frontière ni à un coin disparaîtront au remaillage ; ceux d'un segment, alignés, restent.

use std::collections::{BTreeSet, HashMap};

use super::surface::{cross, dot, norm, sub, Surface};
use super::tetrahedralize::check;
use crate::error::MailleError;

/// Écart au plan toléré, rapporté à la diagonale de la surface : celui de la soudure des sommets
/// d'un STL.
const PLANE_TOLERANCE: f64 = 1e-6;

/// Facettes planes et segments de leurs frontières.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Facets {
    pub facets: Vec<Facet>,
    pub segments: Vec<Segment>,
}

/// Morceau plan et connexe d'un Patch.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Facet {
    /// Indice dans `Surface::patches`.
    pub patch: usize,
    /// Indices dans `Surface::triangles`, croissants.
    pub triangles: Vec<usize>,
    /// Normale unitaire sortante du plan.
    pub normal: [f64; 3],
    /// Un point du plan.
    pub origin: [f64; 3],
    /// Contours fermés, dans le sens des triangles : antihoraire vu de l'extérieur pour le contour
    /// extérieur, horaire pour un trou. Chaque morceau est un segment, parcouru dans son sens
    /// (`false`) ou à rebours (`true`).
    pub loops: Vec<Vec<(usize, bool)>>,
}

/// Frontière entre deux facettes, d'un coin à l'autre.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Segment {
    /// Sommets de la surface, du premier coin au second, alignés entre les deux.
    pub vertices: Vec<usize>,
    /// La facette qui parcourt le segment dans ce sens, puis celle qui le parcourt à rebours.
    pub facets: [usize; 2],
}

/// Facettes planes de la surface. `InvalidSurface` si la surface n'est pas fermée et orientée de
/// façon cohérente, ou si une facette se touche elle-même en un sommet.
pub(crate) fn facets(surface: &Surface) -> Result<Facets, MailleError> {
    let invalid = |problem: &str| MailleError::InvalidSurface {
        problem: problem.to_string(),
    };
    let patch_of = check(surface)?;
    let (vs, ts) = (&surface.vertices, &surface.triangles);

    // triangle qui parcourt chaque arête orientée
    let mut owner: HashMap<(usize, usize), usize> = HashMap::with_capacity(3 * ts.len());
    for (t, tri) in ts.iter().enumerate() {
        for k in 0..3 {
            if owner.insert((tri[k], tri[(k + 1) % 3]), t).is_some() {
                return Err(invalid(
                    "une arête est parcourue deux fois dans le même sens (orientation incohérente ou surface non variété)",
                ));
            }
        }
    }
    if owner.keys().any(|&(a, b)| !owner.contains_key(&(b, a))) {
        return Err(invalid("surface non fermée : une arête n'a qu'un triangle"));
    }
    let twin = |a: usize, b: usize| owner[&(b, a)];

    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in vs {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let tol = PLANE_TOLERANCE * norm(sub(hi, lo));

    // facettes : de proche en proche depuis le plus grand triangle libre, dont le plan sert de
    // référence (le plus précis)
    let normal = |t: usize| {
        let [a, b, c] = ts[t].map(|v| vs[v]);
        cross(sub(b, a), sub(c, a))
    };
    let mut order: Vec<usize> = (0..ts.len()).collect();
    let areas: Vec<f64> = order.iter().map(|&t| norm(normal(t))).collect();
    order.sort_by(|&a, &b| areas[b].total_cmp(&areas[a]));
    let mut facet_of = vec![usize::MAX; ts.len()];
    let mut facets: Vec<Facet> = Vec::new();
    for &seed in &order {
        if facet_of[seed] != usize::MAX {
            continue;
        }
        if areas[seed] == 0.0 {
            return Err(invalid("un triangle est d'aire nulle"));
        }
        let n = normal(seed).map(|x| x / areas[seed]);
        let origin = vs[ts[seed][0]];
        let in_plane = |t: usize| {
            dot(normal(t), n) >= 0.0
                && ts[t]
                    .iter()
                    .all(|&v| dot(sub(vs[v], origin), n).abs() <= tol)
        };
        let f = facets.len();
        facet_of[seed] = f;
        let (mut stack, mut triangles) = (vec![seed], vec![seed]);
        while let Some(t) = stack.pop() {
            for k in 0..3 {
                let o = twin(ts[t][k], ts[t][(k + 1) % 3]);
                if facet_of[o] == usize::MAX && patch_of[o] == patch_of[seed] && in_plane(o) {
                    facet_of[o] = f;
                    stack.push(o);
                    triangles.push(o);
                }
            }
        }
        triangles.sort_unstable();
        facets.push(Facet {
            patch: patch_of[seed],
            triangles,
            normal: n,
            origin,
            loops: Vec::new(),
        });
    }

    // frontières : arêtes entre deux facettes, voisins de chaque sommet le long d'elles
    let facet_on = |a: usize, b: usize| facet_of[owner[&(a, b)]];
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); vs.len()];
    for &(a, b) in owner.keys() {
        if a < b && facet_on(a, b) != facet_on(b, a) {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
    }
    for n in &mut neighbours {
        n.sort_unstable();
    }
    let mut corner: Vec<bool> = (0..vs.len())
        .map(|v| match neighbours[v][..] {
            [] => false,
            [u, w] => {
                let d = sub(vs[w], vs[u]);
                norm(cross(d, sub(vs[v], vs[u]))) > tol * norm(d)
            }
            _ => true,
        })
        .collect();

    // segments : d'un coin au suivant, le long des sommets alignés ; une frontière fermée sans
    // coin (à la tolérance près) prend son premier sommet comme coin
    let mut done: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut segments = Vec::new();
    for pass in 0..2 {
        for v in 0..vs.len() {
            for i in 0..neighbours[v].len() {
                let next = neighbours[v][i];
                if done.contains(&edge(v, next)) || !(corner[v] || pass == 1) {
                    continue;
                }
                corner[v] = true;
                let mut vertices = vec![v, next];
                done.insert(edge(v, next));
                let (mut prev, mut cur) = (v, next);
                while !corner[cur] {
                    // un sommet qui n'est pas un coin a exactement deux voisins
                    let [a, b] = [neighbours[cur][0], neighbours[cur][1]];
                    let following = if a == prev { b } else { a };
                    done.insert(edge(cur, following));
                    vertices.push(following);
                    (prev, cur) = (cur, following);
                }
                segments.push(Segment {
                    facets: [facet_on(v, next), facet_on(next, v)],
                    vertices,
                });
            }
        }
    }

    // contours de chaque facette : ses segments, mis bout à bout dans le sens des triangles
    let mut outgoing: Vec<HashMap<usize, Vec<(usize, bool)>>> = vec![HashMap::new(); facets.len()];
    for (s, seg) in segments.iter().enumerate() {
        let (first, last) = (seg.vertices[0], seg.vertices[seg.vertices.len() - 1]);
        outgoing[seg.facets[0]]
            .entry(first)
            .or_default()
            .push((s, false));
        outgoing[seg.facets[1]]
            .entry(last)
            .or_default()
            .push((s, true));
    }
    let end = |(s, reversed): (usize, bool)| {
        let v = &segments[s].vertices;
        if reversed {
            v[0]
        } else {
            v[v.len() - 1]
        }
    };
    for (f, out) in outgoing.iter_mut().enumerate() {
        if out.values().any(|pieces| pieces.len() > 1) {
            return Err(invalid(
                "une facette plane se touche elle-même en un sommet : remailler la surface autour de ce sommet",
            ));
        }
        let mut starts: Vec<usize> = out.keys().copied().collect();
        starts.sort_unstable();
        for start in starts {
            let mut lp = Vec::new();
            let mut at = start;
            while let Some(piece) = out.get_mut(&at).and_then(Vec::pop) {
                lp.push(piece);
                at = end(piece);
            }
            if !lp.is_empty() {
                facets[f].loops.push(lp);
            }
        }
    }
    Ok(Facets { facets, segments })
}

fn edge(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesher::Patch;

    type Triangle = [[f64; 3]; 3];

    /// Face plane o + s u + t v (s, t dans [0, 1]) en grille mu × mv de carrés coupés en deux,
    /// normale u × v ; `hole` retire une case.
    fn grid(
        o: [f64; 3],
        u: [f64; 3],
        v: [f64; 3],
        [mu, mv]: [usize; 2],
        hole: Option<[usize; 2]>,
    ) -> Vec<Triangle> {
        let p = |i: usize, j: usize| {
            std::array::from_fn(|k| {
                o[k] + u[k] * i as f64 / mu as f64 + v[k] * j as f64 / mv as f64
            })
        };
        let mut out = Vec::new();
        for i in 0..mu {
            for j in 0..mv {
                if hole != Some([i, j]) {
                    out.push([p(i, j), p(i + 1, j), p(i + 1, j + 1)]);
                    out.push([p(i, j), p(i + 1, j + 1), p(i, j + 1)]);
                }
            }
        }
        out
    }

    /// Surface d'une soupe de triangles, sommets soudés à l'identique ; un Patch par groupe.
    fn surface(groups: &[Vec<Triangle>]) -> Surface {
        let mut index: HashMap<[u64; 3], usize> = HashMap::new();
        let mut s = Surface {
            vertices: Vec::new(),
            triangles: Vec::new(),
            patches: Vec::new(),
        };
        for (g, tris) in groups.iter().enumerate() {
            let mut patch = Patch {
                name: format!("p{g}"),
                triangles: Vec::new(),
            };
            for tri in tris {
                let t = tri.map(|p| {
                    *index
                        .entry(p.map(|x| (x + 0.0).to_bits()))
                        .or_insert_with(|| {
                            s.vertices.push(p);
                            s.vertices.len() - 1
                        })
                });
                patch.triangles.push(s.triangles.len());
                s.triangles.push(t);
            }
            s.patches.push(patch);
        }
        s
    }

    const X: [f64; 3] = [1.0, 0.0, 0.0];
    const Y: [f64; 3] = [0.0, 1.0, 0.0];
    const Z: [f64; 3] = [0.0, 0.0, 1.0];
    const O: [f64; 3] = [0.0, 0.0, 0.0];

    /// Cube unité, chaque face en grille m × m, dans l'ordre -x, +x, -y, +y, -z, +z.
    fn cube(m: usize) -> Vec<Vec<Triangle>> {
        vec![
            grid(O, Z, Y, [m, m], None),
            grid(X, Y, Z, [m, m], None),
            grid(O, X, Z, [m, m], None),
            grid(Y, Z, X, [m, m], None),
            grid(O, Y, X, [m, m], None),
            grid(Z, X, Y, [m, m], None),
        ]
    }

    /// Aire signée d'un contour dans le plan de la facette (positive : antihoraire vu de dehors).
    fn loop_area(s: &Surface, fs: &Facets, f: usize, lp: &[(usize, bool)]) -> f64 {
        let mut ring = Vec::new();
        for &(seg, reversed) in lp {
            let mut v = fs.segments[seg].vertices.clone();
            if reversed {
                v.reverse();
            }
            ring.extend_from_slice(&v[..v.len() - 1]);
        }
        let mut c = [0.0; 3];
        for i in 0..ring.len() {
            let x = cross(s.vertices[ring[i]], s.vertices[ring[(i + 1) % ring.len()]]);
            c = std::array::from_fn(|k| c[k] + x[k]);
        }
        0.5 * dot(c, fs.facets[f].normal)
    }

    /// Invariants communs : chaque triangle dans une facette, contours fermés, et la somme des aires
    /// signées des contours égale l'aire de la facette.
    fn check_facets(s: &Surface, fs: &Facets) {
        let mut all: Vec<usize> = fs.facets.iter().flat_map(|f| f.triangles.clone()).collect();
        all.sort_unstable();
        assert_eq!(all, (0..s.triangles.len()).collect::<Vec<_>>());
        for (f, facet) in fs.facets.iter().enumerate() {
            let area: f64 = facet
                .triangles
                .iter()
                .map(|&t| {
                    let [a, b, c] = s.triangles[t].map(|v| s.vertices[v]);
                    0.5 * norm(cross(sub(b, a), sub(c, a)))
                })
                .sum();
            let loops: f64 = facet.loops.iter().map(|lp| loop_area(s, fs, f, lp)).sum();
            assert!(
                (area - loops).abs() < 1e-12,
                "facette {f} : {area} ≠ {loops}"
            );
            for lp in &facet.loops {
                for (k, &(seg, reversed)) in lp.iter().enumerate() {
                    assert_eq!(fs.segments[seg].facets[reversed as usize], f);
                    let v = &fs.segments[seg].vertices;
                    let end = if reversed { v[0] } else { v[v.len() - 1] };
                    let (next, next_rev) = lp[(k + 1) % lp.len()];
                    let w = &fs.segments[next].vertices;
                    assert_eq!(end, if next_rev { w[w.len() - 1] } else { w[0] });
                }
            }
        }
    }

    #[test]
    fn cube_six_facettes_douze_segments() {
        let s = surface(&cube(1));
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), 6);
        assert_eq!(fs.segments.len(), 12);
        assert!(fs.segments.iter().all(|seg| seg.vertices.len() == 2));
        assert!(fs
            .facets
            .iter()
            .all(|f| f.loops.len() == 1 && f.loops[0].len() == 4));
    }

    #[test]
    fn sommets_interieurs_oublies_sommets_alignes_gardes() {
        // grille 2 × 2 : un sommet au centre de chaque face, un au milieu de chaque arête
        let s = surface(&cube(2));
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), 6);
        assert_eq!(fs.segments.len(), 12);
        assert!(fs.segments.iter().all(|seg| seg.vertices.len() == 3));
        let normals: BTreeSet<[i64; 3]> = fs
            .facets
            .iter()
            .map(|f| f.normal.map(|x| x.round() as i64))
            .collect();
        assert_eq!(normals.len(), 6);
    }

    #[test]
    fn deux_patches_dans_un_meme_plan() {
        // face +z coupée en deux Patches par le milieu en x
        let mut groups = cube(2);
        let top = groups.pop().unwrap();
        let (left, right): (Vec<Triangle>, Vec<Triangle>) =
            top.into_iter().partition(|t| t.iter().all(|p| p[0] <= 0.5));
        groups.push(left);
        groups.push(right);
        let s = surface(&groups);
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), 7);
        // les deux arêtes du haut coupées par la frontière, plus la frontière elle-même
        assert_eq!(fs.segments.len(), 15);
        let middle = fs
            .segments
            .iter()
            .filter(|seg| seg.facets.iter().all(|&f| fs.facets[f].normal[2] > 0.5))
            .count();
        assert_eq!(middle, 1);
    }

    #[test]
    fn cylindre_facette() {
        let n = 16;
        let at = |i: usize, z: f64| {
            let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            [a.cos(), a.sin(), z]
        };
        let mut tris = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            tris.push([at(i, 0.0), at(j, 0.0), at(j, 1.0)]);
            tris.push([at(i, 0.0), at(j, 1.0), at(i, 1.0)]);
            tris.push([[0.0, 0.0, 0.0], at(j, 0.0), at(i, 0.0)]);
            tris.push([[0.0, 0.0, 1.0], at(i, 1.0), at(j, 1.0)]);
        }
        let s = surface(&[tris]);
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), n + 2);
        assert_eq!(fs.segments.len(), 3 * n);
        let mut loop_lengths: Vec<usize> = fs.facets.iter().map(|f| f.loops[0].len()).collect();
        loop_lengths.sort_unstable();
        assert_eq!(loop_lengths[..n], vec![4; n][..]);
        assert_eq!(loop_lengths[n..], [n, n]);
    }

    #[test]
    fn prisme_en_l() {
        // L de 2 × 2 moins le carré [1, 2] × [1, 2], épaisseur 1
        let l = |z: f64| {
            [
                [0., 0., z],
                [2., 0., z],
                [2., 1., z],
                [1., 1., z],
                [1., 2., z],
                [0., 2., z],
            ]
        };
        let (b, t) = (l(0.0), l(1.0));
        let mut tris = Vec::new();
        // dessus et dessous en éventail depuis le coin rentrant (indice 3)
        for k in [4, 5, 0, 1] {
            tris.push([t[3], t[k], t[(k + 1) % 6]]);
            tris.push([b[3], b[(k + 1) % 6], b[k]]);
        }
        for k in 0..6 {
            let j = (k + 1) % 6;
            tris.push([b[k], b[j], t[j]]);
            tris.push([b[k], t[j], t[k]]);
        }
        let s = surface(&[tris]);
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), 8);
        assert_eq!(fs.segments.len(), 18);
    }

    #[test]
    fn plaque_trouee_contour_et_trou() {
        // plaque 3 × 3 × 1 percée du carré [1, 2] × [1, 2]
        let s3 = |a: [f64; 3]| a.map(|x| 3.0 * x);
        let groups = vec![
            grid(O, Z, s3(Y), [1, 3], None),
            grid(s3(X), s3(Y), Z, [3, 1], None),
            grid(O, s3(X), Z, [3, 1], None),
            grid(s3(Y), Z, s3(X), [1, 3], None),
            grid(O, s3(Y), s3(X), [3, 3], Some([1, 1])),
            grid(Z, s3(X), s3(Y), [3, 3], Some([1, 1])),
            grid([1., 1., 0.], Y, Z, [1, 1], None),
            grid([2., 1., 0.], Z, Y, [1, 1], None),
            grid([1., 1., 0.], Z, X, [1, 1], None),
            grid([1., 2., 0.], X, Z, [1, 1], None),
        ];
        let s = surface(&groups);
        let fs = facets(&s).unwrap();
        check_facets(&s, &fs);
        assert_eq!(fs.facets.len(), 10);
        // 8 arêtes sur chaque grande face, 8 verticales
        assert_eq!(fs.segments.len(), 24);
        let top = fs.facets.iter().position(|f| f.normal[2] > 0.5).unwrap();
        let mut areas: Vec<f64> = fs.facets[top]
            .loops
            .iter()
            .map(|lp| loop_area(&s, &fs, top, lp))
            .collect();
        areas.sort_by(f64::total_cmp);
        assert_eq!(areas.len(), 2);
        assert!((areas[0] + 1.0).abs() < 1e-12 && (areas[1] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn tolerance_du_plan() {
        // centre de la face +z soulevé : à 1e-9 la face reste plane, à 1e-2 elle ne l'est plus
        for (lift, plane) in [(1e-9, true), (1e-2, false)] {
            let mut s = surface(&cube(2));
            let c = s
                .vertices
                .iter()
                .position(|p| *p == [0.5, 0.5, 1.0])
                .unwrap();
            s.vertices[c][2] += lift;
            let fs = facets(&s).unwrap();
            assert_eq!(fs.facets.len() == 6, plane, "soulèvement {lift}");
        }
    }

    #[test]
    fn surface_ouverte_refusee() {
        let mut s = surface(&cube(1));
        s.triangles.pop();
        if let Some(p) = s.patches.last_mut() {
            p.triangles.pop();
        }
        assert!(matches!(
            facets(&s),
            Err(MailleError::InvalidSurface { problem }) if problem.contains("non fermée")
        ));
    }
}
