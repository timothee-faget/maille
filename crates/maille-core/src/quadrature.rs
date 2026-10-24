//! Intégration de Gauss sur le triangle de référence (0,0), (1,0), (0,1), d'aire 1/2.

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
