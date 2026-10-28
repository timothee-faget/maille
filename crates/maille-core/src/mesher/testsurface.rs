//! Surfaces fermées construites à la main pour les tests du mailleur 3D.

use std::collections::HashMap;

use super::surface::{Patch, Surface};

pub type Triangle = [[f64; 3]; 3];

pub const X: [f64; 3] = [1.0, 0.0, 0.0];
pub const Y: [f64; 3] = [0.0, 1.0, 0.0];
pub const Z: [f64; 3] = [0.0, 0.0, 1.0];
pub const O: [f64; 3] = [0.0, 0.0, 0.0];

/// Face plane o + s u + t v (s, t dans [0, 1]) en grille mu × mv de carrés coupés en deux,
/// normale u × v ; `hole` retire une case.
pub fn grid(
    o: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
    [mu, mv]: [usize; 2],
    hole: Option<[usize; 2]>,
) -> Vec<Triangle> {
    let p = |i: usize, j: usize| {
        std::array::from_fn(|k| o[k] + u[k] * i as f64 / mu as f64 + v[k] * j as f64 / mv as f64)
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

/// Surface d'une soupe de triangles, sommets soudés à l'identique ; un Patch `p0`, `p1`… par
/// groupe.
pub fn surface(groups: &[Vec<Triangle>]) -> Surface {
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

/// Cube unité, chaque face en grille m × m, dans l'ordre -x, +x, -y, +y, -z, +z.
pub fn cube(m: usize) -> Vec<Vec<Triangle>> {
    vec![
        grid(O, Z, Y, [m, m], None),
        grid(X, Y, Z, [m, m], None),
        grid(O, X, Z, [m, m], None),
        grid(Y, Z, X, [m, m], None),
        grid(O, Y, X, [m, m], None),
        grid(Z, X, Y, [m, m], None),
    ]
}

/// Cylindre facetté de rayon 1 et de hauteur 1, à n côtés, axe z, base en z = 0.
pub fn cylinder(n: usize) -> Vec<Triangle> {
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
    tris
}

/// Prisme droit d'épaisseur 1 sur le polygone `base` (antihoraire, en éventail depuis son premier
/// sommet).
pub fn prism(base: &[[f64; 2]]) -> Vec<Triangle> {
    let n = base.len();
    let (b, t): (Vec<[f64; 3]>, Vec<[f64; 3]>) = base
        .iter()
        .map(|p| ([p[0], p[1], 0.0], [p[0], p[1], 1.0]))
        .unzip();
    let mut tris = Vec::new();
    for k in 1..n - 1 {
        tris.push([t[0], t[k], t[k + 1]]);
        tris.push([b[0], b[k + 1], b[k]]);
    }
    for k in 0..n {
        let j = (k + 1) % n;
        tris.push([b[k], b[j], t[j]]);
        tris.push([b[k], t[j], t[k]]);
    }
    tris
}

/// L de 2 × 2 moins le carré [1, 2] × [1, 2], en partant du coin rentrant (l'éventail reste dedans).
pub const L_SHAPE: [[f64; 2]; 6] = [[1., 1.], [1., 2.], [0., 2.], [0., 0.], [2., 0.], [2., 1.]];

/// Plaque 3 × 3 × 1 percée du carré [1, 2] × [1, 2] : dix faces planes.
pub fn holed_plate() -> Vec<Vec<Triangle>> {
    let s3 = |a: [f64; 3]| a.map(|x| 3.0 * x);
    vec![
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
    ]
}
