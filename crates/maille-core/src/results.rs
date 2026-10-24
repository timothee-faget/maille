/// Contraintes planes en un point de Gauss d'un triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneStress {
    /// Position (x, y) du point.
    pub x: f64,
    pub y: f64,
    pub sigma_xx: f64,
    pub sigma_yy: f64,
    pub tau_xy: f64,
    /// Nulle en contraintes planes, ν (σxx + σyy) en déformations planes.
    pub sigma_zz: f64,
}

impl PlaneStress {
    /// Contrainte équivalente de von Mises.
    pub fn von_mises(&self) -> f64 {
        let (a, b, c) = (self.sigma_xx, self.sigma_yy, self.sigma_zz);
        (a * a + b * b + c * c - a * b - b * c - c * a + 3.0 * self.tau_xy * self.tau_xy).sqrt()
    }
}

/// Champs de résultats, une valeur par DOF.
#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    pub displacements: Vec<f64>,
    /// Forces de réaction aux DOF contraints, zéro ailleurs.
    pub reactions: Vec<f64>,
    /// Contraintes aux points de Gauss de chaque Element (indexé comme `Mesh::elements`).
    /// Vide pour les Bar et Beam.
    pub stresses: Vec<Vec<PlaneStress>>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stress(sxx: f64, syy: f64, txy: f64, szz: f64) -> PlaneStress {
        PlaneStress {
            x: 0.0,
            y: 0.0,
            sigma_xx: sxx,
            sigma_yy: syy,
            tau_xy: txy,
            sigma_zz: szz,
        }
    }

    #[test]
    fn von_mises_of_uniaxial_stress_is_that_stress() {
        assert!((stress(100.0, 0.0, 0.0, 0.0).von_mises() - 100.0).abs() < 1e-12);
    }

    #[test]
    fn von_mises_of_pure_shear_is_sqrt3_tau() {
        assert!((stress(0.0, 0.0, 50.0, 0.0).von_mises() - 50.0 * 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn von_mises_of_hydrostatic_stress_is_zero() {
        assert!(stress(80.0, 80.0, 0.0, 80.0).von_mises().abs() < 1e-12);
    }

    #[test]
    fn von_mises_accounts_for_sigma_zz_in_plane_strain() {
        // σxx = σ, σzz = νσ (déformations planes) : σ √(1 − ν + ν²)
        let nu: f64 = 0.3;
        let expected = 100.0 * (1.0 - nu + nu * nu).sqrt();
        assert!((stress(100.0, 0.0, 0.0, nu * 100.0).von_mises() - expected).abs() < 1e-12);
    }
}
