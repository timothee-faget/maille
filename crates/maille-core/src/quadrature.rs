//! Intégration de Gauss sur le triangle de référence (0,0), (1,0), (0,1), d'aire 1/2, et sur le
//! tétraèdre de référence (0,0,0), (1,0,0), (0,1,0), (0,0,1), de volume 1/6.

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
    fn three_point_rule_is_exact_up_to_degree_2() {
        for (a, b) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)] {
            assert!(
                (integrate(&TRIANGLE_3, a, b) - exact(a, b)).abs() < 1e-15,
                "ξ^{a} η^{b}"
            );
        }
    }
}
