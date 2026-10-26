//! Intégration de Gauss sur le triangle de référence (0,0), (1,0), (0,1), d'aire 1/2, sur le
//! tétraèdre de référence (0,0,0), (1,0,0), (0,1,0), (0,0,1), de volume 1/6, et sur le segment [0, 1].

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussPoint {
    pub xi: f64,
    pub eta: f64,
    pub weight: f64,
}

/// 1 point (centre de gravité) : exacte pour les polynômes de degré ≤ 1.
pub const TRIANGLE_1: [GaussPoint; 1] = [GaussPoint {
    xi: 1.0 / 3.0,
    eta: 1.0 / 3.0,
    weight: 0.5,
}];

/// 3 points intérieurs : exacte pour les polynômes de degré ≤ 2.
pub const TRIANGLE_3: [GaussPoint; 3] = [
    GaussPoint {
        xi: 1.0 / 6.0,
        eta: 1.0 / 6.0,
        weight: 1.0 / 6.0,
    },
    GaussPoint {
        xi: 2.0 / 3.0,
        eta: 1.0 / 6.0,
        weight: 1.0 / 6.0,
    },
    GaussPoint {
        xi: 1.0 / 6.0,
        eta: 2.0 / 3.0,
        weight: 1.0 / 6.0,
    },
];

/// 6 points (Dunavant) : exacte pour les polynômes de degré ≤ 4. Sert aux faces des Tet10.
pub const TRIANGLE_6: [GaussPoint; 6] = {
    const A: f64 = 0.445_948_490_915_965;
    const WA: f64 = 0.223_381_589_678_011 / 2.0;
    const B: f64 = 0.091_576_213_509_771;
    const WB: f64 = 0.109_951_743_655_322 / 2.0;
    const fn p(xi: f64, eta: f64, weight: f64) -> GaussPoint {
        GaussPoint { xi, eta, weight }
    }
    [
        p(A, A, WA),
        p(1.0 - 2.0 * A, A, WA),
        p(A, 1.0 - 2.0 * A, WA),
        p(B, B, WB),
        p(1.0 - 2.0 * B, B, WB),
        p(B, 1.0 - 2.0 * B, WB),
    ]
};

/// Point de Gauss sur le segment [0, 1].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinePoint {
    pub s: f64,
    pub weight: f64,
}

/// 2 points : exacte pour les polynômes de degré ≤ 3. 1/2 ∓ 1/(2√3).
pub const LINE_2: [LinePoint; 2] = [
    LinePoint {
        s: 0.211_324_865_405_187_1,
        weight: 0.5,
    },
    LinePoint {
        s: 0.788_675_134_594_812_9,
        weight: 0.5,
    },
];

/// 3 points : exacte pour les polynômes de degré ≤ 5. 1/2 ∓ √(3/5)/2 (poids 5/18), 1/2 (poids 4/9).
pub const LINE_3: [LinePoint; 3] = [
    LinePoint {
        s: 0.112_701_665_379_258_3,
        weight: 5.0 / 18.0,
    },
    LinePoint {
        s: 0.5,
        weight: 4.0 / 9.0,
    },
    LinePoint {
        s: 0.887_298_334_620_741_7,
        weight: 5.0 / 18.0,
    },
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetGaussPoint {
    pub xi: f64,
    pub eta: f64,
    pub zeta: f64,
    pub weight: f64,
}

/// 1 point (centre de gravité) : exacte pour les polynômes de degré ≤ 1.
pub const TETRA_1: [TetGaussPoint; 1] = [TetGaussPoint {
    xi: 0.25,
    eta: 0.25,
    zeta: 0.25,
    weight: 1.0 / 6.0,
}];

/// (5 + 3√5)/20 et (5 − √5)/20 : coordonnées barycentriques des points de `TETRA_4`.
const TET_A: f64 = 0.585_410_196_624_968_5;
const TET_B: f64 = 0.138_196_601_125_010_5;

/// 4 points intérieurs : exacte pour les polynômes de degré ≤ 2. Le point g est près du sommet g
/// (coordonnée barycentrique `TET_A` sur ce sommet, `TET_B` sur les trois autres).
pub const TETRA_4: [TetGaussPoint; 4] = [
    TetGaussPoint {
        xi: TET_B,
        eta: TET_B,
        zeta: TET_B,
        weight: 1.0 / 24.0,
    },
    TetGaussPoint {
        xi: TET_A,
        eta: TET_B,
        zeta: TET_B,
        weight: 1.0 / 24.0,
    },
    TetGaussPoint {
        xi: TET_B,
        eta: TET_A,
        zeta: TET_B,
        weight: 1.0 / 24.0,
    },
    TetGaussPoint {
        xi: TET_B,
        eta: TET_B,
        zeta: TET_A,
        weight: 1.0 / 24.0,
    },
];

/// 14 points (Walkington), tous les poids positifs : exacte pour les polynômes de degré ≤ 5. Sert à
/// la capacité des Tet10 (N Nᵀ de degré 4). Deux groupes de 4 points (coordonnées barycentriques
/// a, a, a, 1 − 3a) et un groupe de 6 (b, b, 1/2 − b, 1/2 − b), constantes résolues à l'arrondi près.
pub const TETRA_14: [TetGaussPoint; 14] = {
    const A1: f64 = 0.092_735_250_310_891_42;
    const W1: f64 = 0.012_248_840_519_393_72;
    const A2: f64 = 0.310_885_919_263_300_8;
    const W2: f64 = 0.018_781_320_953_002_83;
    const B: f64 = 0.045_503_704_125_648_51;
    const W3: f64 = 0.007_091_003_462_846_74;
    const C1: f64 = 1.0 - 3.0 * A1;
    const C2: f64 = 1.0 - 3.0 * A2;
    const D: f64 = 0.5 - B;
    const fn p(xi: f64, eta: f64, zeta: f64, weight: f64) -> TetGaussPoint {
        TetGaussPoint {
            xi,
            eta,
            zeta,
            weight,
        }
    }
    [
        p(A1, A1, A1, W1),
        p(C1, A1, A1, W1),
        p(A1, C1, A1, W1),
        p(A1, A1, C1, W1),
        p(A2, A2, A2, W2),
        p(C2, A2, A2, W2),
        p(A2, C2, A2, W2),
        p(A2, A2, C2, W2),
        p(B, B, D, W3),
        p(B, D, B, W3),
        p(D, B, B, W3),
        p(D, D, B, W3),
        p(D, B, D, W3),
        p(B, D, D, W3),
    ]
};

#[cfg(test)]
mod tests {
    use super::*;

    fn factorial(n: u32) -> f64 {
        (1..=n).map(f64::from).product()
    }

    /// ∫ ξ^a η^b sur le triangle de référence = a! b! / (a+b+2)!
    fn exact(a: u32, b: u32) -> f64 {
        factorial(a) * factorial(b) / factorial(a + b + 2)
    }

    fn integrate(rule: &[GaussPoint], a: u32, b: u32) -> f64 {
        rule.iter()
            .map(|g| g.weight * g.xi.powi(a as i32) * g.eta.powi(b as i32))
            .sum()
    }

    #[test]
    fn weights_sum_to_the_reference_area() {
        assert!((TRIANGLE_1.iter().map(|g| g.weight).sum::<f64>() - 0.5).abs() < 1e-15);
        assert!((TRIANGLE_3.iter().map(|g| g.weight).sum::<f64>() - 0.5).abs() < 1e-15);
    }

    #[test]
    fn one_point_rule_is_exact_up_to_degree_1() {
        for (a, b) in [(0, 0), (1, 0), (0, 1)] {
            assert!((integrate(&TRIANGLE_1, a, b) - exact(a, b)).abs() < 1e-15);
        }
    }

    #[test]
    fn one_point_rule_is_not_exact_at_degree_2() {
        assert!((integrate(&TRIANGLE_1, 2, 0) - exact(2, 0)).abs() > 1e-3);
    }

    /// ∫ ξ^a η^b ζ^c sur le tétraèdre de référence = a! b! c! / (a+b+c+3)!
    fn exact_tet(a: u32, b: u32, c: u32) -> f64 {
        factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3)
    }

    fn integrate_tet(rule: &[TetGaussPoint], a: u32, b: u32, c: u32) -> f64 {
        rule.iter()
            .map(|g| g.weight * g.xi.powi(a as i32) * g.eta.powi(b as i32) * g.zeta.powi(c as i32))
            .sum()
    }

    #[test]
    fn tet_constants_are_the_closed_forms() {
        let r5 = 5.0_f64.sqrt();
        assert!((TET_A - (5.0 + 3.0 * r5) / 20.0).abs() < 1e-16);
        assert!((TET_B - (5.0 - r5) / 20.0).abs() < 1e-16);
    }

    #[test]
    fn tet_rules_are_exact_up_to_their_degree() {
        let linear = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)];
        let quadratic = [
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1),
        ];
        for (a, b, c) in linear {
            assert!((integrate_tet(&TETRA_1, a, b, c) - exact_tet(a, b, c)).abs() < 1e-15);
        }
        assert!((integrate_tet(&TETRA_1, 2, 0, 0) - exact_tet(2, 0, 0)).abs() > 1e-3);
        for (a, b, c) in linear.into_iter().chain(quadratic) {
            assert!(
                (integrate_tet(&TETRA_4, a, b, c) - exact_tet(a, b, c)).abs() < 1e-15,
                "ξ^{a} η^{b} ζ^{c}"
            );
        }
    }

    #[test]
    fn fourteen_point_tet_rule_is_exact_up_to_degree_5() {
        for d in 0..=5 {
            for a in 0..=d {
                for b in 0..=d - a {
                    let c = d - a - b;
                    assert!(
                        (integrate_tet(&TETRA_14, a, b, c) - exact_tet(a, b, c)).abs() < 1e-15,
                        "ξ^{a} η^{b} ζ^{c}"
                    );
                }
            }
        }
        assert!((integrate_tet(&TETRA_14, 6, 0, 0) - exact_tet(6, 0, 0)).abs() > 1e-7);
        assert!(TETRA_14.iter().all(|g| g.weight > 0.0));
    }

    #[test]
    fn six_point_rule_is_exact_up_to_degree_4() {
        for a in 0..=4 {
            for b in 0..=4 - a {
                assert!(
                    (integrate(&TRIANGLE_6, a, b) - exact(a, b)).abs() < 1e-14,
                    "ξ^{a} η^{b}"
                );
            }
        }
        assert!((integrate(&TRIANGLE_6, 6, 0) - exact(6, 0)).abs() > 1e-6);
    }

    #[test]
    fn line_rules_are_exact_up_to_their_degree() {
        let integrate = |rule: &[LinePoint], d: i32| -> f64 {
            rule.iter().map(|p| p.weight * p.s.powi(d)).sum()
        };
        for d in 0..=5 {
            let exact = 1.0 / f64::from(d + 1);
            if d <= 3 {
                assert!((integrate(&LINE_2, d) - exact).abs() < 1e-15, "s^{d}");
            }
            assert!((integrate(&LINE_3, d) - exact).abs() < 1e-15, "s^{d}");
        }
        assert!((integrate(&LINE_2, 4) - 0.2).abs() > 1e-4);
        let r = |x: f64| x.sqrt();
        assert!((LINE_2[0].s - (0.5 - 0.5 / r(3.0))).abs() < 1e-16);
        assert!((LINE_3[0].s - (0.5 - r(0.6) / 2.0)).abs() < 1e-16);
    }

    #[test]
    fn three_point_rule_is_exact_up_to_degree_2() {
        for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)] {
            assert!(
                (integrate(&TRIANGLE_3, a, b) - exact(a, b)).abs() < 1e-15,
                "ξ^{a} η^{b}"
            );
        }
    }
}
