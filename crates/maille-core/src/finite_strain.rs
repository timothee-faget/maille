//! Grandes déformations : Lagrangien total, loi de Saint-Venant–Kirchhoff ou néo-Hookéen.
//!
//! Tout est intégré sur la configuration de départ (points de Gauss et fonctions de forme inchangés).
//! Gradient de la transformation F = I + ∇₀u, déformation de Green-Lagrange E = ½ (FᵀF − I),
//! second tenseur de Piola-Kirchhoff S = ∂W/∂E de la loi (`LargeStrainLaw`). Forces internes f = ∫ P ∇₀N dV₀ avec
//! P = F S, tangente K = K_mat + K_géo (la raideur géométrique est le terme en S). Les contraintes
//! sorties sont celles de Cauchy σ = J⁻¹ F S Fᵀ, sur la configuration déformée.
//!
//! Les triangles sont en déformations planes : F₃₃ = 1, et σzz = S₃₃ / J.

use crate::element::{ElementResponse, GaussData};
use crate::plasticity::PlasticState;
use crate::results::Stress;

type Mat3 = [[f64; 3]; 3];

/// Loi de grandes déformations, de coefficients de Lamé λ et μ.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LargeStrainLaw {
    /// Saint-Venant–Kirchhoff : S = λ tr(E) I + 2μ E (Hooke écrite en E et S).
    SaintVenantKirchhoff { lambda: f64, mu: f64 },
    /// Néo-Hookéen compressible : S = μ (I − C⁻¹) + λ ln J C⁻¹.
    NeoHookean { lambda: f64, mu: f64 },
}

/// ∂S_IJ / ∂E_KL.
type Modulus = [[[[f64; 3]; 3]; 3]; 3];

/// (λ, μ) pour E et ν : λ = E ν / ((1 + ν)(1 − 2ν)), μ = E / (2 (1 + ν)).
pub fn lame(young: f64, poisson: f64) -> (f64, f64) {
    (
        young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
        young / (2.0 * (1.0 + poisson)),
    )
}

impl LargeStrainLaw {
    /// Second tenseur de Piola-Kirchhoff S et tenseur d'élasticité C_IJKL = ∂S_IJ / ∂E_KL en F. Non
    /// fini quand le néo-Hookéen reçoit J ≤ 0 (Element retourné) : Newton-Raphson le refuse.
    fn evaluate(&self, f: &Mat3) -> (Mat3, Modulus) {
        let delta = |a: usize, b: usize| f64::from(u8::from(a == b));
        let c = mul(&transpose(f), f);
        match *self {
            Self::SaintVenantKirchhoff { lambda, mu } => {
                let e: Mat3 =
                    std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (c[i][j] - delta(i, j))));
                let trace = e[0][0] + e[1][1] + e[2][2];
                let s = std::array::from_fn(|i| {
                    std::array::from_fn(|j| 2.0 * mu * e[i][j] + lambda * trace * delta(i, j))
                });
                let modulus = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        std::array::from_fn(|k| {
                            std::array::from_fn(|l| {
                                lambda * delta(i, j) * delta(k, l)
                                    + mu * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k))
                            })
                        })
                    })
                });
                (s, modulus)
            }
            Self::NeoHookean { lambda, mu } => {
                let log_j = det(f).ln();
                let inv = inverse(&c);
                let s = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        mu * (delta(i, j) - inv[i][j]) + lambda * log_j * inv[i][j]
                    })
                });
                // C = λ C⁻¹ ⊗ C⁻¹ + (μ − λ ln J) (C⁻¹_IK C⁻¹_JL + C⁻¹_IL C⁻¹_JK)
                let modulus = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        std::array::from_fn(|k| {
                            std::array::from_fn(|l| {
                                lambda * inv[i][j] * inv[k][l]
                                    + (mu - lambda * log_j)
                                        * (inv[i][k] * inv[j][l] + inv[i][l] * inv[j][k])
                            })
                        })
                    })
                });
                (s, modulus)
            }
        }
    }
}

fn mul(a: &Mat3, b: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}

fn transpose(a: &Mat3) -> Mat3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

fn det(a: &Mat3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn inverse(a: &Mat3) -> Mat3 {
    let d = det(a);
    let cofactor = |i: usize, j: usize| {
        let (r0, r1) = ((i + 1) % 3, (i + 2) % 3);
        let (c0, c1) = ((j + 1) % 3, (j + 2) % 3);
        a[r0][c0] * a[r1][c1] - a[r0][c1] * a[r1][c0]
    };
    // A⁻¹ = cofacteurs transposés / det
    std::array::from_fn(|i| std::array::from_fn(|j| cofactor(j, i) / d))
}

/// Réponse d'un Element de dimension `D` (2 : triangle en déformations planes, 3 : tétraèdre) aux
/// déplacements locaux `u_e` ((ux, uy[, uz]) par Node). Les `stresses` sont de Cauchy, positionnées
/// sur la configuration déformée ; `states` est vide de plasticité (un état nul par point).
// indices de tenseurs (i, J, k, L) : la forme indexée est la plus lisible
#[allow(clippy::needless_range_loop)]
pub fn large_strain_response<const D: usize>(
    points: &[GaussData<D>],
    law: &LargeStrainLaw,
    u_e: &[f64],
) -> ElementResponse {
    let m = u_e.len();
    let mut out = ElementResponse {
        internal_force: vec![0.0; m],
        tangent: vec![vec![0.0; m]; m],
        ..ElementResponse::default()
    };
    for p in points {
        // F = I + Σ_a u_a ⊗ ∇₀N_a (F₃₃ = 1 en 2D)
        let mut f: Mat3 =
            std::array::from_fn(|i| std::array::from_fn(|j| f64::from(u8::from(i == j))));
        for (a, grad) in p.dn_dx.iter().enumerate() {
            for i in 0..D {
                for j in 0..D {
                    f[i][j] += u_e[D * a + i] * grad[j];
                }
            }
        }
        let ft = transpose(&f);
        let (s, modulus) = law.evaluate(&f);
        let pk1 = mul(&f, &s);

        // A_iJkL = ∂P_iJ / ∂F_kL = δik S_JL + F_iI C_IJKL F_kK, pour i, k, J, L dans les D directions
        // (en 2D, F_iI = 0 quand I = 3 et i < 3 : les sommes sur I, K restent dans le plan)
        let mut a = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..D {
            for j in 0..D {
                for k in 0..D {
                    for l in 0..D {
                        let mut material = 0.0;
                        for ii in 0..D {
                            for kk in 0..D {
                                material += f[i][ii] * modulus[ii][j][kk][l] * f[k][kk];
                            }
                        }
                        let geometric = if i == k { s[j][l] } else { 0.0 };
                        a[i][j][k][l] = geometric + material;
                    }
                }
            }
        }

        for (na, ga) in p.dn_dx.iter().enumerate() {
            for i in 0..D {
                let row = D * na + i;
                out.internal_force[row] +=
                    p.weight * (0..D).map(|j| pk1[i][j] * ga[j]).sum::<f64>();
                for (nb, gb) in p.dn_dx.iter().enumerate() {
                    for k in 0..D {
                        let mut v = 0.0;
                        for j in 0..D {
                            for l in 0..D {
                                v += ga[j] * a[i][j][k][l] * gb[l];
                            }
                        }
                        out.tangent[row][D * nb + k] += p.weight * v;
                    }
                }
            }
        }

        // Cauchy σ = F S Fᵀ / J, au point déformé x = X + Σ N u
        let j = det(&f);
        let sigma = mul(&pk1, &ft);
        let at = |c: usize| {
            let x = if c < D { p.position[c] } else { 0.0 };
            if c < D {
                x + p
                    .n
                    .iter()
                    .enumerate()
                    .map(|(a, n)| n * u_e[D * a + c])
                    .sum::<f64>()
            } else {
                x
            }
        };
        out.stresses.push(Stress {
            x: at(0),
            y: at(1),
            z: at(2),
            sigma_xx: sigma[0][0] / j,
            sigma_yy: sigma[1][1] / j,
            sigma_zz: sigma[2][2] / j,
            tau_xy: sigma[0][1] / j,
            tau_yz: sigma[1][2] / j,
            tau_xz: sigma[0][2] / j,
        });
        out.states.push(PlasticState::default());
    }
    out
}
