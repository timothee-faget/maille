//! Fonctions de forme des triangles Tri3 et Tri6 sur le triangle de référence (ξ, η),
//! et passage aux coordonnées physiques (x, y).
//!
//! Numérotation Tri6 : sommets 0, 1, 2 puis milieux d'arêtes 3 (0-1), 4 (1-2), 5 (2-0).
//!
//! Même chose pour les tétraèdres sur le tétraèdre de référence (ξ, η, ζ), vers (x, y, z).
//! Numérotation Tet10 (comme Abaqus) : sommets 0, 1, 2, 3 puis milieux d'arêtes 4 (0-1), 5 (1-2),
//! 6 (2-0), 7 (0-3), 8 (1-3), 9 (2-3).

/// Valeurs N_i et dérivées [∂N_i/∂ξ, ∂N_i/∂η] du Tri3 (linéaire).
pub fn tri3(xi: f64, eta: f64) -> ([f64; 3], [[f64; 2]; 3]) {
    (
        [1.0 - xi - eta, xi, eta],
        [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]],
    )
}

/// Fonction de forme d'un triangle : (ξ, η) → (valeurs N_i, dérivées [∂N_i/∂ξ, ∂N_i/∂η]).
pub type ShapeFn = fn(f64, f64) -> (Vec<f64>, Vec<[f64; 2]>);

/// `tri3` en `Vec`, pour le code générique sur la forme de l'Element (voir `element::plane_stiffness`).
pub fn tri3_vec(xi: f64, eta: f64) -> (Vec<f64>, Vec<[f64; 2]>) {
    let (n, dn) = tri3(xi, eta);
    (n.to_vec(), dn.to_vec())
}

/// Valeurs N_i et dérivées [∂N_i/∂ξ, ∂N_i/∂η] du Tri6 (quadratique).
pub fn tri6(xi: f64, eta: f64) -> ([f64; 6], [[f64; 2]; 6]) {
    let (l1, l2, l3) = (1.0 - xi - eta, xi, eta);
    (
        [
            l1 * (2.0 * l1 - 1.0),
            l2 * (2.0 * l2 - 1.0),
            l3 * (2.0 * l3 - 1.0),
            4.0 * l1 * l2,
            4.0 * l2 * l3,
            4.0 * l3 * l1,
        ],
        [
            [1.0 - 4.0 * l1, 1.0 - 4.0 * l1],
            [4.0 * l2 - 1.0, 0.0],
            [0.0, 4.0 * l3 - 1.0],
            [4.0 * (l1 - l2), -4.0 * l2],
            [4.0 * l3, 4.0 * l2],
            [-4.0 * l3, 4.0 * (l1 - l3)],
        ],
    )
}

/// `tri6` en `Vec`, voir `tri3_vec`.
pub fn tri6_vec(xi: f64, eta: f64) -> (Vec<f64>, Vec<[f64; 2]>) {
    let (n, dn) = tri6(xi, eta);
    (n.to_vec(), dn.to_vec())
}

/// Transformation du triangle de référence vers l'Element réel, en un point.
#[derive(Debug, Clone, PartialEq)]
pub struct Mapping {
    /// Déterminant du jacobien J = ∂(x, y)/∂(ξ, η). Positif si le triangle est antihoraire.
    pub det: f64,
    /// Dérivées physiques [∂N_i/∂x, ∂N_i/∂y].
    pub dn_dx: Vec<[f64; 2]>,
}

impl Mapping {
    /// `coords` : positions (x, y) des Nodes ; `dn` : dérivées [∂N_i/∂ξ, ∂N_i/∂η] au même point.
    /// Suppose `det ≠ 0` (garanti par `Mesh::new`, qui refuse les triangles dégénérés).
    pub fn new(coords: &[[f64; 2]], dn: &[[f64; 2]]) -> Self {
        // j[i][c] = Σ_k ∂N_k/∂(ξ, η)[i] · coords_k[c]
        let mut j = [[0.0; 2]; 2];
        for (c, d) in coords.iter().zip(dn) {
            for (row, dn_i) in j.iter_mut().zip(d) {
                row[0] += dn_i * c[0];
                row[1] += dn_i * c[1];
            }
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let dn_dx = dn
            .iter()
            .map(|d| {
                [
                    (j[1][1] * d[0] - j[0][1] * d[1]) / det,
                    (-j[1][0] * d[0] + j[0][0] * d[1]) / det,
                ]
            })
            .collect();
        Self { det, dn_dx }
    }
}

/// Fonction de forme d'un tétraèdre : (ξ, η, ζ) → (valeurs N_i, dérivées [∂N_i/∂ξ, ∂N_i/∂η, ∂N_i/∂ζ]).
pub type SolidShapeFn = fn(f64, f64, f64) -> (Vec<f64>, Vec<[f64; 3]>);

/// Coordonnées barycentriques L_i du tétraèdre de référence et leurs dérivées (constantes).
fn tet_barycentric(xi: f64, eta: f64, zeta: f64) -> ([f64; 4], [[f64; 3]; 4]) {
    (
        [1.0 - xi - eta - zeta, xi, eta, zeta],
        [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
}

/// Valeurs N_i = L_i et dérivées du Tet4 (linéaire).
pub fn tet4(xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    let (l, dl) = tet_barycentric(xi, eta, zeta);
    (l.to_vec(), dl.to_vec())
}

/// Arêtes du Tet10 dans l'ordre de ses Nodes milieux 4 à 9.
pub const TET10_EDGES: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];

/// Valeurs et dérivées du Tet10 (quadratique) : sommets L_i (2 L_i − 1), milieux 4 L_a L_b.
pub fn tet10(xi: f64, eta: f64, zeta: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    let (l, dl) = tet_barycentric(xi, eta, zeta);
    let mut n = Vec::with_capacity(10);
    let mut dn = Vec::with_capacity(10);
    for i in 0..4 {
        n.push(l[i] * (2.0 * l[i] - 1.0));
        dn.push(dl[i].map(|d| (4.0 * l[i] - 1.0) * d));
    }
    for (a, b) in TET10_EDGES {
        n.push(4.0 * l[a] * l[b]);
        dn.push(std::array::from_fn(|c| {
            4.0 * (l[a] * dl[b][c] + l[b] * dl[a][c])
        }));
    }
    (n, dn)
}

/// Transformation du tétraèdre de référence vers l'Element réel, en un point.
#[derive(Debug, Clone, PartialEq)]
pub struct Mapping3d {
    /// Déterminant du jacobien J = ∂(x, y, z)/∂(ξ, η, ζ). Positif si le tétraèdre est bien orienté.
    pub det: f64,
    /// Dérivées physiques [∂N_i/∂x, ∂N_i/∂y, ∂N_i/∂z].
    pub dn_dx: Vec<[f64; 3]>,
}

impl Mapping3d {
    /// `coords` : positions (x, y, z) des Nodes ; `dn` : dérivées par rapport à (ξ, η, ζ) au même point.
    /// Suppose `det ≠ 0` (garanti par `Mesh::new`, qui refuse les tétraèdres dégénérés).
    pub fn new(coords: &[[f64; 3]], dn: &[[f64; 3]]) -> Self {
        // j[i][c] = Σ_k ∂N_k/∂(ξ, η, ζ)[i] · coords_k[c]
        let mut j = [[0.0; 3]; 3];
        for (c, d) in coords.iter().zip(dn) {
            for (row, dn_i) in j.iter_mut().zip(d) {
                for (r, x) in row.iter_mut().zip(c) {
                    *r += dn_i * x;
                }
            }
        }
        // J⁻¹ = cof(J)ᵀ / det
        let cof = |r: usize, c: usize| {
            let (r1, r2, c1, c2) = ((r + 1) % 3, (r + 2) % 3, (c + 1) % 3, (c + 2) % 3);
            j[r1][c1] * j[r2][c2] - j[r1][c2] * j[r2][c1]
        };
        let det: f64 = (0..3).map(|c| j[0][c] * cof(0, c)).sum();
        let inv: [[f64; 3]; 3] = std::array::from_fn(|r| std::array::from_fn(|c| cof(c, r) / det));
        let dn_dx = dn
            .iter()
            .map(|d| std::array::from_fn(|r| (0..3).map(|c| inv[r][c] * d[c]).sum()))
            .collect();
        Self { det, dn_dx }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRI3_NODES: [[f64; 2]; 3] = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    const TRI6_NODES: [[f64; 2]; 6] = [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 0.0],
        [0.5, 0.5],
        [0.0, 0.5],
    ];
    const SAMPLES: [[f64; 2]; 4] = [[0.2, 0.3], [0.1, 0.1], [0.6, 0.2], [1.0 / 3.0, 1.0 / 3.0]];

    #[test]
    fn shape_functions_are_one_at_their_node_and_zero_at_the_others() {
        for (j, p) in TRI3_NODES.iter().enumerate() {
            let (n, _) = tri3(p[0], p[1]);
            for (i, v) in n.iter().enumerate() {
                assert!((v - if i == j { 1.0 } else { 0.0 }).abs() < 1e-15);
            }
        }
        for (j, p) in TRI6_NODES.iter().enumerate() {
            let (n, _) = tri6(p[0], p[1]);
            for (i, v) in n.iter().enumerate() {
                assert!((v - if i == j { 1.0 } else { 0.0 }).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn values_sum_to_one_and_derivatives_sum_to_zero() {
        for p in SAMPLES {
            let (n3, d3) = tri3(p[0], p[1]);
            let (n6, d6) = tri6(p[0], p[1]);
            assert!((n3.iter().sum::<f64>() - 1.0).abs() < 1e-15);
            assert!((n6.iter().sum::<f64>() - 1.0).abs() < 1e-15);
            for c in 0..2 {
                assert!(d3.iter().map(|d| d[c]).sum::<f64>().abs() < 1e-15);
                assert!(d6.iter().map(|d| d[c]).sum::<f64>().abs() < 1e-15);
            }
        }
    }

    #[test]
    fn tri3_reproduces_linear_and_tri6_reproduces_quadratic_fields() {
        for p in SAMPLES {
            let (n3, _) = tri3(p[0], p[1]);
            let linear = |x: f64, y: f64| 2.0 + 3.0 * x - 5.0 * y;
            let interp3: f64 = n3
                .iter()
                .zip(TRI3_NODES)
                .map(|(n, c)| n * linear(c[0], c[1]))
                .sum();
            assert!((interp3 - linear(p[0], p[1])).abs() < 1e-14);

            let (n6, _) = tri6(p[0], p[1]);
            let quadratic = |x: f64, y: f64| 1.0 + x - 2.0 * y + 3.0 * x * x - x * y + 4.0 * y * y;
            let interp6: f64 = n6
                .iter()
                .zip(TRI6_NODES)
                .map(|(n, c)| n * quadratic(c[0], c[1]))
                .sum();
            assert!((interp6 - quadratic(p[0], p[1])).abs() < 1e-14);
        }
    }

    #[test]
    fn tri6_derivatives_match_central_differences() {
        let h = 1e-6;
        for p in SAMPLES {
            let (_, d) = tri6(p[0], p[1]);
            let (nxp, _) = tri6(p[0] + h, p[1]);
            let (nxm, _) = tri6(p[0] - h, p[1]);
            let (nyp, _) = tri6(p[0], p[1] + h);
            let (nym, _) = tri6(p[0], p[1] - h);
            for i in 0..6 {
                assert!((d[i][0] - (nxp[i] - nxm[i]) / (2.0 * h)).abs() < 1e-8);
                assert!((d[i][1] - (nyp[i] - nym[i]) / (2.0 * h)).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn tri3_mapping_gives_area_and_constant_gradients() {
        // triangle (0,0), (2,0), (0,3) : aire 3, N0 = 1 − x/2 − y/3, N1 = x/2, N2 = y/3
        let coords = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]];
        let (_, dn) = tri3(0.3, 0.3);
        let m = Mapping::new(&coords, &dn);
        assert!((m.det / 2.0 - 3.0).abs() < 1e-14);
        let expected = [[-0.5, -1.0 / 3.0], [0.5, 0.0], [0.0, 1.0 / 3.0]];
        for (g, e) in m.dn_dx.iter().zip(expected) {
            assert!((g[0] - e[0]).abs() < 1e-14 && (g[1] - e[1]).abs() < 1e-14);
        }
    }

    #[test]
    fn clockwise_triangle_has_negative_determinant() {
        let coords = [[0.0, 0.0], [0.0, 3.0], [2.0, 0.0]];
        let (_, dn) = tri3(0.3, 0.3);
        assert!(Mapping::new(&coords, &dn).det < 0.0);
    }

    #[test]
    fn straight_sided_tri6_recovers_the_gradient_of_x_and_y() {
        // triangle (1,1), (4,2), (2,5), milieux d'arêtes au milieu des côtés
        let c = [[1.0, 1.0], [4.0, 2.0], [2.0, 5.0]];
        let mid = |a: [f64; 2], b: [f64; 2]| [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
        let coords = [
            c[0],
            c[1],
            c[2],
            mid(c[0], c[1]),
            mid(c[1], c[2]),
            mid(c[2], c[0]),
        ];
        for p in SAMPLES {
            let (_, dn) = tri6(p[0], p[1]);
            let m = Mapping::new(&coords, &dn);
            let gx: [f64; 2] =
                std::array::from_fn(|k| m.dn_dx.iter().zip(coords).map(|(g, c)| g[k] * c[0]).sum());
            let gy: [f64; 2] =
                std::array::from_fn(|k| m.dn_dx.iter().zip(coords).map(|(g, c)| g[k] * c[1]).sum());
            assert!((gx[0] - 1.0).abs() < 1e-13 && gx[1].abs() < 1e-13);
            assert!(gy[0].abs() < 1e-13 && (gy[1] - 1.0).abs() < 1e-13);
        }
    }

    const TET10_NODES: [[f64; 3]; 10] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    const TET_SAMPLES: [[f64; 3]; 3] = [[0.1, 0.2, 0.3], [0.25, 0.25, 0.25], [0.6, 0.1, 0.05]];

    #[test]
    fn tet_shape_functions_are_one_at_their_node_and_zero_at_the_others() {
        for (shape, count) in [(tet4 as SolidShapeFn, 4), (tet10, 10)] {
            for (j, p) in TET10_NODES.iter().take(count).enumerate() {
                let (n, _) = shape(p[0], p[1], p[2]);
                for (i, v) in n.iter().enumerate() {
                    assert!((v - if i == j { 1.0 } else { 0.0 }).abs() < 1e-15);
                }
            }
        }
    }

    #[test]
    fn tet10_reproduces_quadratic_fields_and_its_derivatives_match_differences() {
        let q = |x: f64, y: f64, z: f64| {
            1.0 + x - 2.0 * y + z + 3.0 * x * x - x * y + 2.0 * y * z - z * z
        };
        let h = 1e-6;
        for p in TET_SAMPLES {
            let (n, d) = tet10(p[0], p[1], p[2]);
            let interp: f64 = n
                .iter()
                .zip(TET10_NODES)
                .map(|(n, c)| n * q(c[0], c[1], c[2]))
                .sum();
            assert!((interp - q(p[0], p[1], p[2])).abs() < 1e-14);
            assert!((n.iter().sum::<f64>() - 1.0).abs() < 1e-15);
            for c in 0..3 {
                let mut plus = p;
                let mut minus = p;
                plus[c] += h;
                minus[c] -= h;
                let (np, _) = tet10(plus[0], plus[1], plus[2]);
                let (nm, _) = tet10(minus[0], minus[1], minus[2]);
                for i in 0..10 {
                    assert!((d[i][c] - (np[i] - nm[i]) / (2.0 * h)).abs() < 1e-8);
                }
            }
        }
    }

    #[test]
    fn tet_mapping_gives_volume_and_recovers_the_gradient_of_x_y_z() {
        // tétraèdre quelconque, Tet10 à arêtes droites : 6 V = det, ∇x = (1, 0, 0), etc.
        let v = [
            [1.0, 0.5, 0.0],
            [3.0, 1.0, 0.5],
            [1.5, 2.5, 0.0],
            [1.2, 1.0, 2.0],
        ];
        let mut coords = v.to_vec();
        for (a, b) in TET10_EDGES {
            coords.push(std::array::from_fn(|c| (v[a][c] + v[b][c]) / 2.0));
        }
        let edge = |k: usize| -> [f64; 3] { std::array::from_fn(|c| v[k][c] - v[0][c]) };
        let (e1, e2, e3) = (edge(1), edge(2), edge(3));
        let six_v = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
            - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
            + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
        for p in TET_SAMPLES {
            let (_, dn) = tet10(p[0], p[1], p[2]);
            let m = Mapping3d::new(&coords, &dn);
            assert!((m.det - six_v).abs() < 1e-13);
            for x in 0..3 {
                for k in 0..3 {
                    let g: f64 = m.dn_dx.iter().zip(&coords).map(|(g, c)| g[k] * c[x]).sum();
                    assert!((g - if k == x { 1.0 } else { 0.0 }).abs() < 1e-13);
                }
            }
        }
    }

    #[test]
    fn inverted_tetrahedron_has_negative_determinant() {
        let coords = [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let (_, dn) = tet4(0.25, 0.25, 0.25);
        assert!(Mapping3d::new(&coords, &dn).det < 0.0);
    }
}
