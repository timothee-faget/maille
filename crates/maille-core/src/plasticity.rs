//! Loi de comportement élasto-plastique de von Mises en un point de Gauss, en petites déformations,
//! à écrouissage isotrope linéaire. Pure : ni Mesh ni Element.
//!
//! Notation de Voigt de `solid_d` : (xx, yy, zz, yz, xz, xy), glissements γ = 2 ε en déformation
//! d'ingénieur, contraintes de cisaillement τ telles quelles.

/// État interne d'un point de Gauss, nul au départ.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PlasticState {
    /// Déformation plastique εp (Voigt, glissements d'ingénieur), de trace nulle.
    pub plastic_strain: [f64; 6],
    /// Déformation plastique cumulée p = ∫ √(2/3 ε̇p : ε̇p) dt.
    pub p: f64,
}

/// Réponse du point de Gauss à une déformation totale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialPoint {
    pub stress: [f64; 6],
    /// État à la fin du pas ; à ne garder que si l'incrément converge.
    pub state: PlasticState,
    /// Tangente cohérente dσ/dε, symétrique.
    pub tangent: [[f64; 6]; 6],
    /// Le point a plastifié pendant ce pas.
    pub yielding: bool,
}

/// Matériau de von Mises : σy(p) = σy0 + H p. `yield_stress` infini : élastique.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VonMises {
    pub young: f64,
    pub poisson: f64,
    pub yield_stress: f64,
    pub hardening: f64,
}

impl VonMises {
    /// Retour radial en forme fermée : prédiction élastique depuis l'état `previous`, puis, si la
    /// contrainte équivalente q dépasse σy(p), retour sur la surface le long du déviateur.
    /// `strain` : déformation totale à la fin du pas.
    pub fn update(&self, strain: &[f64; 6], previous: &PlasticState) -> MaterialPoint {
        let g = self.young / (2.0 * (1.0 + self.poisson));
        let bulk = self.young / (3.0 * (1.0 - 2.0 * self.poisson));
        let elastic: [f64; 6] = std::array::from_fn(|k| strain[k] - previous.plastic_strain[k]);
        let volume = elastic[0] + elastic[1] + elastic[2];
        // déviateur d'essai s = 2G e (le tenseur e a γ/2 hors diagonale, d'où G γ en cisaillement)
        let trial: [f64; 6] = std::array::from_fn(|k| {
            if k < 3 {
                2.0 * g * (elastic[k] - volume / 3.0)
            } else {
                g * elastic[k]
            }
        });
        let root_3_2 = 1.5_f64.sqrt();
        let norm = contract(&trial, &trial).sqrt();
        let q = root_3_2 * norm;
        let excess = q - (self.yield_stress + self.hardening * previous.p);

        // dp : incrément de p qui ramène q sur σy(p + dp) ; θ : facteur de réduction du déviateur
        // (s = θ s_essai) ; θ̄ : raidissement perdu dans la direction n quand dp varie avec ε, nul en
        // élastique, égal à θ si H = 0
        let (theta, theta_bar, dp) = if excess > 0.0 {
            let dp = excess / (3.0 * g + self.hardening);
            let theta = 1.0 - 3.0 * g * dp / q;
            (
                theta,
                3.0 * g / (3.0 * g + self.hardening) - (1.0 - theta),
                dp,
            )
        } else {
            (1.0, 0.0, 0.0)
        };
        // n = s/‖s‖, direction d'écoulement ; ε̇p = √(3/2) ṗ n
        let n: [f64; 6] = std::array::from_fn(|k| if dp > 0.0 { trial[k] / norm } else { 0.0 });
        let mut state = *previous;
        for (k, e) in state.plastic_strain.iter_mut().enumerate() {
            let tensor = root_3_2 * dp * n[k];
            *e += if k < 3 { tensor } else { 2.0 * tensor };
        }
        state.p += dp;

        let stress =
            std::array::from_fn(|k| theta * trial[k] + if k < 3 { bulk * volume } else { 0.0 });
        // C = K 1⊗1 + 2Gθ I_dev − 2Gθ̄ n⊗n, appliquée à des glissements d'ingénieur : I_dev vaut
        // δ − 1/3 sur les composantes normales et 1/2 sur les glissements ; n:dε = Σ nₖ dεₖ.
        let tangent = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let delta = if i == j { 1.0 } else { 0.0 };
                let deviatoric = match (i < 3, j < 3) {
                    (true, true) => delta - 1.0 / 3.0,
                    (false, false) if i == j => 0.5,
                    _ => 0.0,
                };
                let volumetric = if i < 3 && j < 3 { bulk } else { 0.0 };
                volumetric + 2.0 * g * theta * deviatoric - 2.0 * g * theta_bar * n[i] * n[j]
            })
        });
        MaterialPoint {
            stress,
            state,
            tangent,
            yielding: dp > 0.0,
        }
    }
}

/// a : b pour deux tenseurs symétriques en Voigt dont les cisaillements sont des composantes du
/// tenseur (contraintes, pas glissements).
fn contract(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    (0..6)
        .map(|k| {
            if k < 3 {
                a[k] * b[k]
            } else {
                2.0 * a[k] * b[k]
            }
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::solid_d;

    const E: f64 = 210e3;
    const NU: f64 = 0.3;
    const SY: f64 = 250.0;
    const H: f64 = 2000.0;

    fn law(hardening: f64) -> VonMises {
        VonMises {
            young: E,
            poisson: NU,
            yield_stress: SY,
            hardening,
        }
    }

    fn von_mises(s: &[f64; 6]) -> f64 {
        let (a, b, c) = (s[0], s[1], s[2]);
        (a * a + b * b + c * c - a * b - b * c - c * a
            + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt()
    }

    /// Déformation totale d'une traction uniaxiale σ avec une déformation plastique cumulée p :
    /// élastique σ/E (1, −ν, −ν) plus plastique p (1, −½, −½).
    fn uniaxial(sigma: f64, p: f64) -> [f64; 6] {
        let lateral = -NU * sigma / E - p / 2.0;
        [sigma / E + p, lateral, lateral, 0.0, 0.0, 0.0]
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn elastic_tangent_is_solid_d() {
        let point = law(H).update(&uniaxial(0.5 * SY, 0.0), &PlasticState::default());
        assert!(!point.yielding);
        let d = solid_d(E, NU);
        for (row, d_row) in point.tangent.iter().zip(&d) {
            for (c, d) in row.iter().zip(d_row) {
                assert!(close(*c, *d, 1e-9 * E), "{row:?}");
            }
        }
    }

    #[test]
    fn uniaxial_tension_follows_the_hardening_slope_then_unloads_elastically() {
        // σ = σy0 + H p : la pente apparente après la limite est E H / (E + H)
        let p = 2e-3;
        let sigma = SY + H * p;
        let loaded = law(H).update(&uniaxial(sigma, p), &PlasticState::default());
        assert!(loaded.yielding);
        assert!(
            close(loaded.stress[0], sigma, 1e-9 * SY),
            "{:?}",
            loaded.stress
        );
        for s in &loaded.stress[1..] {
            assert!(s.abs() < 1e-9 * SY, "{:?}", loaded.stress);
        }
        assert!(close(loaded.state.p, p, 1e-12));
        let ep = loaded.state.plastic_strain;
        assert!(close(ep[0], p, 1e-12) && close(ep[1], -p / 2.0, 1e-12));
        assert!((ep[0] + ep[1] + ep[2]).abs() < 1e-15);

        let total = uniaxial(sigma, p);
        let tangent_slope = E * H / (E + H);
        let further = law(H).update(&uniaxial(sigma + 10.0, p + 10.0 / H), &loaded.state);
        let de = uniaxial(sigma + 10.0, p + 10.0 / H)[0] - total[0];
        assert!(close(
            further.stress[0] - sigma,
            tangent_slope * de,
            1e-9 * SY
        ));

        // décharge complète : σ = 0, la déformation plastique reste
        let unloaded = law(H).update(&uniaxial(0.0, p), &loaded.state);
        assert!(!unloaded.yielding);
        assert!(unloaded.stress.iter().all(|s| s.abs() < 1e-9 * SY));
        assert_eq!(unloaded.state, loaded.state);
    }

    #[test]
    fn perfect_plasticity_stays_on_the_yield_surface() {
        // traction puis cisaillement : q = σy0 à chaque pas plastique
        let mut state = PlasticState::default();
        for k in 1..=10 {
            let s = k as f64;
            let strain = [2e-3 * s, -1e-3 * s, -0.5e-3 * s, 0.0, 0.0, 3e-3 * s];
            let point = law(0.0).update(&strain, &state);
            if point.yielding {
                assert!(close(von_mises(&point.stress), SY, 1e-9 * SY), "pas {k}");
            }
            state = point.state;
        }
        assert!(state.p > 0.0);
    }

    #[test]
    fn pure_shear_yields_at_sigma_y_over_root_3() {
        let g = E / (2.0 * (1.0 + NU));
        let gamma = 10.0 * SY / g;
        let point = law(0.0).update(&[0.0, 0.0, 0.0, 0.0, 0.0, gamma], &PlasticState::default());
        assert!(close(point.stress[5], SY / 3.0_f64.sqrt(), 1e-9 * SY));
    }

    #[test]
    fn pressure_does_not_cause_yielding() {
        let e = 1e-2; // σm = K 3e, bien au-delà de σy
        let point = law(H).update(&[e, e, e, 0.0, 0.0, 0.0], &PlasticState::default());
        assert!(!point.yielding);
        let bulk = E / (3.0 * (1.0 - 2.0 * NU));
        assert!(close(point.stress[0], 3.0 * bulk * e, 1e-9 * bulk * e));
    }

    #[test]
    fn consistent_tangent_matches_finite_differences() {
        let start = law(H)
            .update(
                &[1.5e-3, -0.4e-3, -0.2e-3, 0.3e-3, -0.2e-3, 0.8e-3],
                &PlasticState::default(),
            )
            .state;
        let strain = [3e-3, -1e-3, -0.6e-3, 1e-3, -0.5e-3, 2e-3];
        for hardening in [H, 0.0] {
            let m = law(hardening);
            for (previous, label) in [(PlasticState::default(), "vierge"), (start, "écroui")] {
                let point = m.update(&strain, &previous);
                assert!(point.yielding, "{label}");
                let h = 1e-8;
                for j in 0..6 {
                    let (mut plus, mut minus) = (strain, strain);
                    plus[j] += h;
                    minus[j] -= h;
                    let (sp, sm) = (
                        m.update(&plus, &previous).stress,
                        m.update(&minus, &previous).stress,
                    );
                    for i in 0..6 {
                        let numeric = (sp[i] - sm[i]) / (2.0 * h);
                        assert!(
                            close(point.tangent[i][j], numeric, 1e-7 * E),
                            "H = {hardening}, {label}, C[{i}][{j}] = {} au lieu de {numeric}",
                            point.tangent[i][j]
                        );
                        assert!(close(point.tangent[i][j], point.tangent[j][i], 1e-9 * E));
                    }
                }
            }
        }
    }

    #[test]
    fn infinite_yield_stress_is_elastic() {
        let m = VonMises {
            yield_stress: f64::INFINITY,
            ..law(H)
        };
        let point = m.update(&uniaxial(10.0 * SY, 0.0), &PlasticState::default());
        assert!(!point.yielding);
        assert!(close(point.stress[0], 10.0 * SY, 1e-9 * SY));
    }
}
