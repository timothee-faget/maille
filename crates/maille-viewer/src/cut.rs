//! Plan de coupe : on garde le côté d ≤ 0 du plan (d = n · p − décalage). `section` donne la surface
//! où le plan traverse les tétraèdres, avec le champ interpolé ; `clip_triangles` et
//! `clip_segments` coupent la peau et son contour au même plan. Tout est calculé sur le CPU, sans
//! GPU, et recalculé quand le plan bouge.

use crate::vec3::{dot, lerp, V3};

/// Plan de normale `normal` (unitaire) passant par les points tels que n · p = `offset`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    pub normal: V3,
    pub offset: f64,
}

impl Plane {
    /// Distance signée de `p` au plan : négative ou nulle du côté gardé.
    pub fn distance(&self, p: V3) -> f64 {
        dot(self.normal, p) - self.offset
    }
}

/// Maillage de triangles avec une valeur par sommet (sommets non partagés entre triangles).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Patch {
    pub positions: Vec<V3>,
    pub values: Vec<f64>,
    pub triangles: Vec<[u32; 3]>,
}

impl Patch {
    fn push(&mut self, polygon: &[(V3, f64)]) {
        let first = self.positions.len() as u32;
        for &(p, v) in polygon {
            self.positions.push(p);
            self.values.push(v);
        }
        for k in 1..polygon.len() as u32 - 1 {
            self.triangles.push([first, first + k, first + k + 1]);
        }
    }
}

/// Surface où le plan coupe les tétraèdres de `cells` (leurs 4 premiers Nodes : les sommets, un
/// Tet10 se coupe comme son Tet4), avec le champ `values` interpolé linéairement le long des arêtes.
pub fn section(cells: &[Vec<usize>], positions: &[V3], values: &[f64], plane: Plane) -> Patch {
    let mut patch = Patch::default();
    for cell in cells {
        let v: [usize; 4] = std::array::from_fn(|k| cell[k]);
        let d = v.map(|i| plane.distance(positions[i]));
        let (above, below): (Vec<usize>, Vec<usize>) = (0..4).partition(|&k| d[k] > 0.0);
        // point où l'arête (a au-dessus, b en dessous) traverse le plan
        let cross = |a: usize, b: usize| {
            let t = d[a] / (d[a] - d[b]);
            (
                lerp(positions[v[a]], positions[v[b]], t),
                values[v[a]] + t * (values[v[b]] - values[v[a]]),
            )
        };
        match (above.len(), below.len()) {
            (1, 3) => patch.push(
                &below
                    .iter()
                    .map(|&b| cross(above[0], b))
                    .collect::<Vec<_>>(),
            ),
            (3, 1) => patch.push(
                &above
                    .iter()
                    .map(|&a| cross(a, below[0]))
                    .collect::<Vec<_>>(),
            ),
            // deux sommets de chaque côté : quadrilatère, arêtes prises dans l'ordre du contour
            (2, 2) => patch.push(&[
                cross(above[0], below[0]),
                cross(above[0], below[1]),
                cross(above[1], below[1]),
                cross(above[1], below[0]),
            ]),
            _ => {}
        }
    }
    patch
}

/// Triangles de `triangles` découpés au plan, côté d ≤ 0 gardé, champ interpolé.
pub fn clip_triangles(
    positions: &[V3],
    values: &[f64],
    triangles: &[[u32; 3]],
    plane: Plane,
) -> Patch {
    let mut patch = Patch::default();
    for t in triangles {
        let corners = t.map(|i| (positions[i as usize], values[i as usize]));
        // Sutherland-Hodgman contre un seul plan : au plus un sommet de plus
        let mut polygon: Vec<(V3, f64)> = Vec::with_capacity(4);
        for k in 0..3 {
            let (a, b) = (corners[k], corners[(k + 1) % 3]);
            let (da, db) = (plane.distance(a.0), plane.distance(b.0));
            if da <= 0.0 {
                polygon.push(a);
            }
            if (da <= 0.0) != (db <= 0.0) {
                let s = da / (da - db);
                polygon.push((lerp(a.0, b.0, s), a.1 + s * (b.1 - a.1)));
            }
        }
        if polygon.len() >= 3 {
            patch.push(&polygon);
        }
    }
    patch
}

/// Segments `edges` découpés au plan, côté d ≤ 0 gardé.
pub fn clip_segments(positions: &[V3], edges: &[[u32; 2]], plane: Plane) -> Vec<[V3; 2]> {
    edges
        .iter()
        .filter_map(|e| {
            let (a, b) = (positions[e[0] as usize], positions[e[1] as usize]);
            let (da, db) = (plane.distance(a), plane.distance(b));
            match (da <= 0.0, db <= 0.0) {
                (true, true) => Some([a, b]),
                (false, false) => None,
                (true, false) => Some([a, lerp(a, b, da / (da - db))]),
                (false, true) => Some([lerp(a, b, da / (da - db)), b]),
            }
        })
        .collect()
}

/// Aire d'un ensemble de triangles.
#[cfg(test)]
pub fn area(patch: &Patch) -> f64 {
    patch
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| patch.positions[i as usize]);
            0.5 * crate::vec3::norm(crate::vec3::cross(
                crate::vec3::sub(b, a),
                crate::vec3::sub(c, a),
            ))
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solid::skin;
    use crate::testmesh::{cube, cube_quadratic};

    /// Champ affine : un maillage conforme à P1 le retrouve exactement.
    fn field(p: V3) -> f64 {
        p[0] + 2.0 * p[1] + 3.0 * p[2]
    }

    fn values(positions: &[V3]) -> Vec<f64> {
        positions.iter().copied().map(field).collect()
    }

    fn plane(normal: V3, offset: f64) -> Plane {
        let n = crate::vec3::norm(normal);
        Plane {
            normal: normal.map(|x| x / n),
            offset,
        }
    }

    #[test]
    fn an_axis_aligned_section_of_the_cube_is_a_square_with_the_exact_field() {
        for (cells, positions) in [cube(), cube_quadratic()].map(|(p, c)| (c, p)) {
            let v = values(&positions);
            for c in [-0.7, 0.0, 0.31, 0.9] {
                let s = section(&cells, &positions, &v, plane([1.0, 0.0, 0.0], c));
                // carré de côté 2
                assert!((area(&s) - 4.0).abs() < 1e-12, "x = {c} : {}", area(&s));
                for (p, value) in s.positions.iter().zip(&s.values) {
                    assert!((p[0] - c).abs() < 1e-12);
                    assert!((value - field(*p)).abs() < 1e-12, "champ en {p:?}");
                }
            }
        }
    }

    #[test]
    fn an_oblique_section_of_the_cube_is_a_regular_hexagon() {
        let (positions, cells) = cube();
        let v = values(&positions);
        // x + y + z = 0 passe par le centre : hexagone régulier de côté √2, aire (3√3/2) · 2
        let s = section(&cells, &positions, &v, plane([1.0, 1.0, 1.0], 0.0));
        assert!(
            (area(&s) - 3.0 * 3.0_f64.sqrt()).abs() < 1e-12,
            "{}",
            area(&s)
        );
        for (p, value) in s.positions.iter().zip(&s.values) {
            assert!(((p[0] + p[1] + p[2]) / 3.0_f64.sqrt()).abs() < 1e-12);
            assert!((value - field(*p)).abs() < 1e-12);
        }
    }

    #[test]
    fn a_plane_outside_the_mesh_cuts_nothing() {
        let (positions, cells) = cube();
        let v = values(&positions);
        assert!(section(&cells, &positions, &v, plane([0.0, 0.0, 1.0], 1.5))
            .triangles
            .is_empty());
        assert!(
            section(&cells, &positions, &v, plane([0.0, 0.0, 1.0], -1.5))
                .triangles
                .is_empty()
        );
    }

    #[test]
    fn clipping_the_skin_keeps_the_half_space_with_the_exact_field() {
        let (positions, cells) = cube();
        let s = skin(&cells).unwrap();
        let v = values(&positions);
        // coupe à x = 0 : reste la face x = -1 (4), et la moitié des quatre faces latérales (4 · 2)
        let p = plane([1.0, 0.0, 0.0], 0.0);
        let kept = clip_triangles(&positions, &v, &s.triangles, p);
        assert!((area(&kept) - 12.0).abs() < 1e-12, "{}", area(&kept));
        for (q, value) in kept.positions.iter().zip(&kept.values) {
            assert!(q[0] <= 1e-12);
            assert!((value - field(*q)).abs() < 1e-12);
        }
        // et la coupe complète bien la surface : peau gardée + section = surface fermée du demi-cube
        let cut = section(&cells, &positions, &v, p);
        assert!((area(&kept) + area(&cut) - (12.0 + 4.0)).abs() < 1e-12);
    }

    #[test]
    fn a_plane_that_keeps_everything_or_nothing_changes_nothing_or_removes_everything() {
        let (positions, cells) = cube();
        let s = skin(&cells).unwrap();
        let v = values(&positions);
        let all = clip_triangles(&positions, &v, &s.triangles, plane([1.0, 0.0, 0.0], 5.0));
        assert!((area(&all) - 24.0).abs() < 1e-12);
        let none = clip_triangles(&positions, &v, &s.triangles, plane([1.0, 0.0, 0.0], -5.0));
        assert!(none.triangles.is_empty());
    }

    #[test]
    fn clipped_outline_segments_end_on_the_plane() {
        let (positions, cells) = cube();
        let s = skin(&cells).unwrap();
        let p = plane([0.0, 1.0, 0.0], 0.25);
        let segments = clip_segments(&positions, &s.edges, p);
        assert!(!segments.is_empty());
        for [a, b] in &segments {
            assert!(p.distance(*a) <= 1e-12 && p.distance(*b) <= 1e-12);
        }
        // les arêtes qui traversent le plan s'arrêtent dessus : au moins les 4 arêtes du cube
        // parallèles à y, et aucune ne sort du cube
        let ends: Vec<V3> = segments
            .iter()
            .flatten()
            .copied()
            .filter(|q| p.distance(*q).abs() < 1e-12)
            .collect();
        assert!(ends.len() >= 4, "{} extrémités sur le plan", ends.len());
        assert!(ends
            .iter()
            .all(|q| q[0].abs() <= 1.0 + 1e-12 && q[2].abs() <= 1.0 + 1e-12));
    }
}
