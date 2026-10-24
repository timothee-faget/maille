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
