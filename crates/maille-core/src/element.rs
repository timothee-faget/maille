use crate::mesh::PlaneState;
use crate::quadrature::{GaussPoint, TetGaussPoint};
use crate::results::Stress;
use crate::shape::{Mapping, Mapping3d, SolidShapeFn};

/// Matrice de rigidité d'une barre 1D : (EA/L) [[1, -1], [-1, 1]].
pub fn bar_stiffness(young: f64, area: f64, length: f64) -> [[f64; 2]; 2] {
    let k = young * area / length;
    [[k, -k], [-k, k]]
}

/// Matrice de rigidité d'une poutre d'Euler-Bernoulli, DOF locaux (v1, θ1, v2, θ2).
pub fn beam_stiffness(young: f64, inertia: f64, l: f64) -> [[f64; 4]; 4] {
    let c = young * inertia / (l * l * l);
    [
        [12.0 * c, 6.0 * l * c, -12.0 * c, 6.0 * l * c],
        [6.0 * l * c, 4.0 * l * l * c, -6.0 * l * c, 2.0 * l * l * c],
        [-12.0 * c, -6.0 * l * c, 12.0 * c, -6.0 * l * c],
        [6.0 * l * c, 2.0 * l * l * c, -6.0 * l * c, 4.0 * l * l * c],
    ]
}

/// Forces nodales équivalentes d'une charge répartie uniforme q sur une barre.
pub fn bar_equivalent_load(q: f64, l: f64) -> [f64; 2] {
    [q * l / 2.0, q * l / 2.0]
}

/// Forces nodales équivalentes d'une charge répartie uniforme q sur une poutre (v1, θ1, v2, θ2).
pub fn beam_equivalent_load(q: f64, l: f64) -> [f64; 4] {
    [
        q * l / 2.0,
        q * l * l / 12.0,
        q * l / 2.0,
        -q * l * l / 12.0,
    ]
}

/// Matrice D (3×3) reliant (εxx, εyy, γxy) à (σxx, σyy, τxy), γxy en déformation d'ingénieur.
pub fn plane_d(young: f64, nu: f64, state: PlaneState) -> [[f64; 3]; 3] {
    match state {
        PlaneState::Stress => {
            let c = young / (1.0 - nu * nu);
            [
                [c, c * nu, 0.0],
                [c * nu, c, 0.0],
                [0.0, 0.0, c * (1.0 - nu) / 2.0],
            ]
        }
        PlaneState::Strain => {
            let c = young / ((1.0 + nu) * (1.0 - 2.0 * nu));
            [
                [c * (1.0 - nu), c * nu, 0.0],
                [c * nu, c * (1.0 - nu), 0.0],
                [0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0],
            ]
        }
    }
}

/// Matrice B (3 × 2n) stockée par colonne : `b[i]` est la colonne i, pour (εxx, εyy, γxy).
/// Les DOF sont ordonnés (ux, uy) par Node.
fn strain_displacement(dn_dx: &[[f64; 2]]) -> Vec<[f64; 3]> {
    let mut b = vec![[0.0; 3]; 2 * dn_dx.len()];
    for (k, grad) in dn_dx.iter().enumerate() {
        b[2 * k] = [grad[0], 0.0, grad[1]];
        b[2 * k + 1] = [0.0, grad[1], grad[0]];
    }
    b
}

/// Rigidité d'un triangle : K = Σ_g w_g t |det J| Bᵀ D B, DOF locaux (ux, uy) par Node.
/// `shape_at(ξ, η)` donne les valeurs N_i et les dérivées [∂N_i/∂ξ, ∂N_i/∂η] de la forme choisie
/// (Tri3, plus tard Tri6).
pub fn plane_stiffness(
    coords: &[[f64; 2]],
    rule: &[GaussPoint],
    shape_at: impl Fn(f64, f64) -> (Vec<f64>, Vec<[f64; 2]>),
    d: &[[f64; 3]; 3],
    thickness: f64,
) -> Vec<Vec<f64>> {
    let m = 2 * coords.len();
    let mut ke = vec![vec![0.0; m]; m];
    for g in rule {
        let map = Mapping::new(coords, &shape_at(g.xi, g.eta).1);
        let b = strain_displacement(&map.dn_dx);
        let factor = g.weight * thickness * map.det.abs();
        for i in 0..m {
            // (D B)_i = D · colonne i de B
            let db: [f64; 3] = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * b[i][c]).sum());
            for j in 0..m {
                ke[j][i] += factor * (0..3).map(|r| b[j][r] * db[r]).sum::<f64>();
            }
        }
    }
    ke
}

/// Contraintes σ = D B u aux points de Gauss d'un triangle ; `u_e` : déplacements locaux (ux, uy) par Node.
pub fn plane_stresses(
    coords: &[[f64; 2]],
    rule: &[GaussPoint],
    shape_at: impl Fn(f64, f64) -> (Vec<f64>, Vec<[f64; 2]>),
    d: &[[f64; 3]; 3],
    u_e: &[f64],
    nu: f64,
    state: PlaneState,
) -> Vec<Stress> {
    rule.iter()
        .map(|g| {
            let (n, dn) = shape_at(g.xi, g.eta);
            let b = strain_displacement(&Mapping::new(coords, &dn).dn_dx);
            let eps: [f64; 3] =
                std::array::from_fn(|r| b.iter().zip(u_e).map(|(col, u)| col[r] * u).sum());
            let sigma: [f64; 3] = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * eps[c]).sum());
            Stress {
                x: n.iter().zip(coords).map(|(n, c)| n * c[0]).sum(),
                y: n.iter().zip(coords).map(|(n, c)| n * c[1]).sum(),
                z: 0.0,
                sigma_xx: sigma[0],
                sigma_yy: sigma[1],
                sigma_zz: match state {
                    PlaneState::Stress => 0.0,
                    PlaneState::Strain => nu * (sigma[0] + sigma[1]),
                },
                tau_xy: sigma[2],
                tau_yz: 0.0,
                tau_xz: 0.0,
            }
        })
        .collect()
}

/// Matrice D (6×6) de l'élasticité isotrope 3D, notation de Voigt (xx, yy, zz, yz, xz, xy),
/// glissements γ en déformation d'ingénieur.
pub fn solid_d(young: f64, nu: f64) -> [[f64; 6]; 6] {
    let c = young / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let g = young / (2.0 * (1.0 + nu));
    let mut d = [[0.0; 6]; 6];
    for (i, row) in d.iter_mut().enumerate().take(3) {
        for (j, v) in row.iter_mut().enumerate().take(3) {
            *v = if i == j { c * (1.0 - nu) } else { c * nu };
        }
    }
    for (i, row) in d.iter_mut().enumerate().skip(3) {
        row[i] = g;
    }
    d
}

/// Matrice B (6 × 3n) stockée par colonne, ordre de Voigt de `solid_d`, DOF (ux, uy, uz) par Node.
fn solid_strain_displacement(dn_dx: &[[f64; 3]]) -> Vec<[f64; 6]> {
    dn_dx
        .iter()
        .flat_map(|&[gx, gy, gz]| {
            [
                [gx, 0.0, 0.0, 0.0, gz, gy],
                [0.0, gy, 0.0, gz, 0.0, gx],
                [0.0, 0.0, gz, gy, gx, 0.0],
            ]
        })
        .collect()
}

/// Rigidité d'un tétraèdre : K = Σ_g w_g |det J| Bᵀ D B, DOF locaux (ux, uy, uz) par Node.
pub fn solid_stiffness(
    coords: &[[f64; 3]],
    rule: &[TetGaussPoint],
    shape_at: SolidShapeFn,
    d: &[[f64; 6]; 6],
) -> Vec<Vec<f64>> {
    let m = 3 * coords.len();
    let mut ke = vec![vec![0.0; m]; m];
    for g in rule {
        let map = Mapping3d::new(coords, &shape_at(g.xi, g.eta, g.zeta).1);
        let b = solid_strain_displacement(&map.dn_dx);
        let factor = g.weight * map.det.abs();
        for i in 0..m {
            let db: [f64; 6] = std::array::from_fn(|r| (0..6).map(|c| d[r][c] * b[i][c]).sum());
            for j in 0..m {
                ke[j][i] += factor * (0..6).map(|r| b[j][r] * db[r]).sum::<f64>();
            }
        }
    }
    ke
}

/// Contraintes σ = D B u aux points de Gauss d'un tétraèdre ; `u_e` : déplacements locaux (ux, uy, uz) par Node.
pub fn solid_stresses(
    coords: &[[f64; 3]],
    rule: &[TetGaussPoint],
    shape_at: SolidShapeFn,
    d: &[[f64; 6]; 6],
    u_e: &[f64],
) -> Vec<Stress> {
    rule.iter()
        .map(|g| {
            let (n, dn) = shape_at(g.xi, g.eta, g.zeta);
            let b = solid_strain_displacement(&Mapping3d::new(coords, &dn).dn_dx);
            let eps: [f64; 6] =
                std::array::from_fn(|r| b.iter().zip(u_e).map(|(col, u)| col[r] * u).sum());
            let sigma: [f64; 6] = std::array::from_fn(|r| (0..6).map(|c| d[r][c] * eps[c]).sum());
            let at = |k: usize| n.iter().zip(coords).map(|(n, c)| n * c[k]).sum();
            Stress {
                x: at(0),
                y: at(1),
                z: at(2),
                sigma_xx: sigma[0],
                sigma_yy: sigma[1],
                sigma_zz: sigma[2],
                tau_yz: sigma[3],
                tau_xz: sigma[4],
                tau_xy: sigma[5],
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bar_stiffness_is_ea_over_l() {
        let k = bar_stiffness(200.0, 2.0, 4.0); // EA/L = 100
        assert_eq!(k, [[100.0, -100.0], [-100.0, 100.0]]);
    }

    #[test]
    fn beam_stiffness_is_symmetric_with_zero_row_sums_for_translation() {
        let k = beam_stiffness(1.0, 1.0, 2.0);
        for (i, row) in k.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!((v - k[j][i]).abs() < 1e-14);
            }
        }
        // translation rigide (v1 = v2 = 1) : aucune force
        for row in &k {
            assert!((row[0] + row[2]).abs() < 1e-14);
        }
    }
}

#[cfg(test)]
mod plane_tests {
    use super::*;
    use crate::quadrature::TRIANGLE_1;
    use crate::shape::tri3_vec;

    const COORDS: [[f64; 2]; 3] = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]]; // aire 3
    const THICKNESS: f64 = 0.5;

    fn tri3_ke(d: &[[f64; 3]; 3]) -> Vec<Vec<f64>> {
        plane_stiffness(&COORDS, &TRIANGLE_1, tri3_vec, d, THICKNESS)
    }

    fn mul(k: &[Vec<f64>], u: &[f64]) -> Vec<f64> {
        k.iter()
            .map(|row| row.iter().zip(u).map(|(a, b)| a * b).sum())
            .collect()
    }

    #[test]
    fn plane_d_matches_the_closed_forms() {
        let (e, nu) = (200.0, 0.25);
        let g = e / (2.0 * (1.0 + nu));

        let s = plane_d(e, nu, PlaneState::Stress);
        assert!((s[0][0] - e / (1.0 - nu * nu)).abs() < 1e-12);
        assert!((s[0][1] - nu * s[0][0]).abs() < 1e-12);
        assert!((s[2][2] - g).abs() < 1e-12);

        let p = plane_d(e, nu, PlaneState::Strain);
        assert!((p[0][0] - e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))).abs() < 1e-12);
        assert!((p[0][1] - e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))).abs() < 1e-12);
        assert!((p[2][2] - g).abs() < 1e-12);
    }

    #[test]
    fn tri3_stiffness_is_symmetric() {
        let ke = tri3_ke(&plane_d(1.0, 0.3, PlaneState::Stress));
        for (i, row) in ke.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!((v - ke[j][i]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn rigid_body_motions_produce_no_force() {
        let ke = tri3_ke(&plane_d(210e3, 0.3, PlaneState::Stress));
        let scale = ke[0][0];
        let translation_x: Vec<f64> = (0..3).flat_map(|_| [1.0, 0.0]).collect();
        let translation_y: Vec<f64> = (0..3).flat_map(|_| [0.0, 1.0]).collect();
        // petite rotation : u = −θ y, v = θ x
        let rotation: Vec<f64> = COORDS.iter().flat_map(|c| [-c[1], c[0]]).collect();
        for u in [translation_x, translation_y, rotation] {
            for f in mul(&ke, &u) {
                assert!(f.abs() < 1e-12 * scale, "f = {f}");
            }
        }
    }

    #[test]
    fn strain_energy_of_a_uniform_strain_is_t_a_eps_d_eps() {
        // Avec 3 modes rigides nuls (test précédent) et ces énergies exactes sur les 3 déformations
        // uniformes, la forme quadratique uᵀKu est entièrement fixée.
        let area = 3.0;
        for state in [PlaneState::Stress, PlaneState::Strain] {
            let d = plane_d(210e3, 0.3, state);
            let ke = tri3_ke(&d);
            for eps in [
                [1e-3, 0.0, 0.0],
                [0.0, 1e-3, 0.0],
                [0.0, 0.0, 1e-3],
                [1e-3, -2e-3, 3e-3],
            ] {
                // u = εxx x + γ/2 y, v = εyy y + γ/2 x
                let u: Vec<f64> = COORDS
                    .iter()
                    .flat_map(|c| {
                        [
                            eps[0] * c[0] + eps[2] / 2.0 * c[1],
                            eps[1] * c[1] + eps[2] / 2.0 * c[0],
                        ]
                    })
                    .collect();
                let energy: f64 = u.iter().zip(mul(&ke, &u)).map(|(a, b)| a * b).sum();
                let expected: f64 = THICKNESS
                    * area
                    * (0..3)
                        .map(|r| eps[r] * (0..3).map(|c| d[r][c] * eps[c]).sum::<f64>())
                        .sum::<f64>();
                assert!(
                    (energy - expected).abs() < 1e-10 * expected,
                    "{state:?} {eps:?}"
                );
            }
        }
    }
}

#[cfg(test)]
mod tri6_tests {
    use super::*;
    use crate::quadrature::TRIANGLE_3;
    use crate::shape::tri6_vec;

    // mêmes sommets que les tests du Tri3 (aire 3), milieux d'arêtes 0-1, 1-2, 2-0
    const COORDS: [[f64; 2]; 6] = [
        [0.0, 0.0],
        [2.0, 0.0],
        [0.0, 3.0],
        [1.0, 0.0],
        [1.0, 1.5],
        [0.0, 1.5],
    ];
    const THICKNESS: f64 = 0.5;

    fn ke(d: &[[f64; 3]; 3]) -> Vec<Vec<f64>> {
        plane_stiffness(&COORDS, &TRIANGLE_3, tri6_vec, d, THICKNESS)
    }

    fn quad_form(k: &[Vec<f64>], u: &[f64]) -> f64 {
        k.iter()
            .zip(u)
            .map(|(row, ui)| ui * row.iter().zip(u).map(|(a, b)| a * b).sum::<f64>())
            .sum()
    }

    #[test]
    fn tri6_stiffness_is_symmetric_and_ignores_rigid_motions() {
        let k = ke(&plane_d(210e3, 0.3, PlaneState::Stress));
        for (i, row) in k.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!((v - k[j][i]).abs() < 1e-9 * k[0][0]);
            }
        }
        let translation_x: Vec<f64> = COORDS.iter().flat_map(|_| [1.0, 0.0]).collect();
        let translation_y: Vec<f64> = COORDS.iter().flat_map(|_| [0.0, 1.0]).collect();
        let rotation: Vec<f64> = COORDS.iter().flat_map(|c| [-c[1], c[0]]).collect();
        for u in [translation_x, translation_y, rotation] {
            for row in &k {
                let f: f64 = row.iter().zip(&u).map(|(a, b)| a * b).sum();
                assert!(f.abs() < 1e-9 * k[0][0], "f = {f}");
            }
        }
    }

    #[test]
    fn tri6_strain_energy_of_a_uniform_strain_is_t_a_eps_d_eps() {
        let area = 3.0;
        for state in [PlaneState::Stress, PlaneState::Strain] {
            let d = plane_d(210e3, 0.3, state);
            let k = ke(&d);
            for eps in [
                [1e-3, 0.0, 0.0],
                [0.0, 1e-3, 0.0],
                [0.0, 0.0, 1e-3],
                [1e-3, -2e-3, 3e-3],
            ] {
                let u: Vec<f64> = COORDS
                    .iter()
                    .flat_map(|c| {
                        [
                            eps[0] * c[0] + eps[2] / 2.0 * c[1],
                            eps[1] * c[1] + eps[2] / 2.0 * c[0],
                        ]
                    })
                    .collect();
                let expected: f64 = THICKNESS
                    * area
                    * (0..3)
                        .map(|r| eps[r] * (0..3).map(|c| d[r][c] * eps[c]).sum::<f64>())
                        .sum::<f64>();
                assert!(
                    (quad_form(&k, &u) - expected).abs() < 1e-9 * expected,
                    "{state:?} {eps:?}"
                );
            }
        }
    }
}

#[cfg(test)]
mod solid_tests {
    use super::*;
    use crate::quadrature::{TETRA_1, TETRA_4};
    use crate::shape::{tet10, tet4, TET10_EDGES};

    const VERTICES: [[f64; 3]; 4] = [
        [1.0, 0.5, 0.0],
        [3.0, 1.0, 0.5],
        [1.5, 2.5, 0.0],
        [1.2, 1.0, 2.0],
    ];

    /// Sommets puis, pour le Tet10, milieux d'arêtes.
    fn coords(quadratic: bool) -> Vec<[f64; 3]> {
        let mut c = VERTICES.to_vec();
        if quadratic {
            for (a, b) in TET10_EDGES {
                c.push(std::array::from_fn(|k| {
                    (VERTICES[a][k] + VERTICES[b][k]) / 2.0
                }));
            }
        }
        c
    }

    fn volume() -> f64 {
        let e = |k: usize| -> [f64; 3] { std::array::from_fn(|c| VERTICES[k][c] - VERTICES[0][c]) };
        let (a, b, c) = (e(1), e(2), e(3));
        (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0
    }

    /// (coordonnées, rigidité, D) du Tet4 puis du Tet10.
    type Kernel = (Vec<[f64; 3]>, Vec<Vec<f64>>, [[f64; 6]; 6]);

    fn kernels() -> [Kernel; 2] {
        let d = solid_d(210e3, 0.3);
        let (c4, c10) = (coords(false), coords(true));
        let k4 = solid_stiffness(&c4, &TETRA_1, tet4, &d);
        let k10 = solid_stiffness(&c10, &TETRA_4, tet10, &d);
        [(c4, k4, d), (c10, k10, d)]
    }

    /// Déplacement u = ε·x d'une déformation uniforme (Voigt, glissements d'ingénieur).
    fn uniform_strain(coords: &[[f64; 3]], eps: [f64; 6]) -> Vec<f64> {
        coords
            .iter()
            .flat_map(|&[x, y, z]| {
                [
                    eps[0] * x + eps[5] / 2.0 * y + eps[4] / 2.0 * z,
                    eps[5] / 2.0 * x + eps[1] * y + eps[3] / 2.0 * z,
                    eps[4] / 2.0 * x + eps[3] / 2.0 * y + eps[2] * z,
                ]
            })
            .collect()
    }

    #[test]
    fn solid_d_matches_lame_coefficients() {
        let (e, nu) = (210e3, 0.3);
        let d = solid_d(e, nu);
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        assert!((d[0][0] - (lambda + 2.0 * mu)).abs() < 1e-9);
        assert!((d[1][2] - lambda).abs() < 1e-9);
        assert!((d[4][4] - mu).abs() < 1e-9);
        assert!(d[0][3].abs() < 1e-15 && d[3][4].abs() < 1e-15);
    }

    #[test]
    fn tet_stiffness_is_symmetric_and_ignores_the_six_rigid_motions() {
        for (c, k, _) in kernels() {
            let scale = k[0][0];
            for (i, row) in k.iter().enumerate() {
                for (j, v) in row.iter().enumerate() {
                    assert!((v - k[j][i]).abs() < 1e-9 * scale);
                }
            }
            let translations = (0..3).map(|axis| {
                c.iter()
                    .flat_map(|_| {
                        std::array::from_fn::<f64, 3, _>(|a| if a == axis { 1.0 } else { 0.0 })
                    })
                    .collect::<Vec<f64>>()
            });
            // petites rotations ω × x autour de x, y, z
            let rotations = (0..3).map(|axis| {
                c.iter()
                    .flat_map(|&[x, y, z]| match axis {
                        0 => [0.0, -z, y],
                        1 => [z, 0.0, -x],
                        _ => [-y, x, 0.0],
                    })
                    .collect::<Vec<f64>>()
            });
            for u in translations.chain(rotations) {
                for row in &k {
                    let f: f64 = row.iter().zip(&u).map(|(a, b)| a * b).sum();
                    assert!(f.abs() < 1e-9 * scale, "f = {f}");
                }
            }
        }
    }

    #[test]
    fn tet_strain_energy_of_a_uniform_strain_is_v_eps_d_eps() {
        for (c, k, d) in kernels() {
            for eps in [
                [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1e-3, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1e-3, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 1e-3],
                [1e-3, -2e-3, 0.5e-3, 3e-3, -1e-3, 2e-3],
            ] {
                let u = uniform_strain(&c, eps);
                let energy: f64 = k
                    .iter()
                    .zip(&u)
                    .map(|(row, ui)| ui * row.iter().zip(&u).map(|(a, b)| a * b).sum::<f64>())
                    .sum();
                let expected = volume()
                    * (0..6)
                        .map(|r| eps[r] * (0..6).map(|s| d[r][s] * eps[s]).sum::<f64>())
                        .sum::<f64>();
                assert!((energy - expected).abs() < 1e-9 * expected, "{eps:?}");
            }
        }
    }

    #[test]
    fn tet_stresses_of_a_uniform_strain_are_d_eps_at_every_gauss_point() {
        let eps = [1e-3, -2e-3, 0.5e-3, 3e-3, -1e-3, 2e-3];
        for (c, _, d) in kernels() {
            let (rule, shape): (&[TetGaussPoint], SolidShapeFn) = if c.len() == 4 {
                (&TETRA_1, tet4)
            } else {
                (&TETRA_4, tet10)
            };
            let sigma: [f64; 6] = std::array::from_fn(|r| (0..6).map(|s| d[r][s] * eps[s]).sum());
            for p in solid_stresses(&c, rule, shape, &d, &uniform_strain(&c, eps)) {
                let got = [
                    p.sigma_xx, p.sigma_yy, p.sigma_zz, p.tau_yz, p.tau_xz, p.tau_xy,
                ];
                for (g, e) in got.iter().zip(sigma) {
                    assert!((g - e).abs() < 1e-9 * d[0][0] * 1e-3, "{got:?}");
                }
            }
        }
    }
}
