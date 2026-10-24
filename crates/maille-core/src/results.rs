use crate::error::MailleError;
use crate::mesh::{Element, Mesh};

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

/// Extrapolation du Tri3 : son unique point de Gauss donne la valeur de ses 3 Nodes.
const TRI3_EXTRAPOLATION: [&[f64]; 3] = [&[1.0]; 3];

/// Extrapolation du Tri6 : la ligne k donne les poids des 3 points de Gauss (`TRIANGLE_3`, le point g
/// est près du sommet g) pour la valeur au Node k. Les points de Gauss forment le triangle de
/// référence réduit de moitié autour de son centre : un point de coordonnées barycentriques λ y a
/// les coordonnées 2λ − 1/3. Exacte pour une contrainte linéaire (flexion pure).
const TRI6_EXTRAPOLATION: [&[f64]; 6] = [
    &[5.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0],
    &[-1.0 / 3.0, 5.0 / 3.0, -1.0 / 3.0],
    &[-1.0 / 3.0, -1.0 / 3.0, 5.0 / 3.0],
    &[2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0], // milieu 0-1
    &[-1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0], // milieu 1-2
    &[2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0], // milieu 2-0
];

/// Contraintes lissées aux Nodes (indexées comme `Mesh::nodes`, `x`, `y` = position du Node) :
/// extrapolées des points de Gauss vers les Nodes de chaque Element, puis moyennées composante par
/// composante entre les Elements qui partagent un Node. Un Node isolé reçoit des contraintes nulles.
pub fn nodal_stresses(mesh: &Mesh, field: &Field) -> Result<Vec<PlaneStress>, MailleError> {
    if field.stresses.len() != mesh.elements.len() {
        return Err(MailleError::FieldMeshMismatch);
    }
    // somme de (σxx, σyy, τxy, σzz) et nombre d'Elements, par Node
    let mut sum = vec![[0.0; 4]; mesh.nodes.len()];
    let mut count = vec![0_usize; mesh.nodes.len()];
    for (e, points) in mesh.elements.iter().zip(&field.stresses) {
        let extrapolation: &[&[f64]] = match e {
            Element::Tri3(_) => &TRI3_EXTRAPOLATION,
            Element::Tri6(_) => &TRI6_EXTRAPOLATION,
            Element::Bar(_) | Element::Beam(_) => return Err(MailleError::NoPlaneStresses),
        };
        if points.len() != extrapolation[0].len() {
            return Err(MailleError::FieldMeshMismatch);
        }
        for (&n, weights) in e.nodes().iter().zip(extrapolation) {
            for (p, w) in points.iter().zip(weights.iter()) {
                sum[n][0] += w * p.sigma_xx;
                sum[n][1] += w * p.sigma_yy;
                sum[n][2] += w * p.tau_xy;
                sum[n][3] += w * p.sigma_zz;
            }
            count[n] += 1;
        }
    }
    Ok(mesh
        .nodes
        .iter()
        .zip(sum.iter().zip(&count))
        .map(|(node, (s, &c))| {
            let c = c.max(1) as f64;
            PlaneStress {
                x: node.x,
                y: node.y,
                sigma_xx: s[0] / c,
                sigma_yy: s[1] / c,
                tau_xy: s[2] / c,
                sigma_zz: s[3] / c,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{Material, Section};
    use crate::mesh::{Bar, Node, PlaneState, Tri3};
    use crate::quadrature::TRIANGLE_3;

    #[test]
    fn tri6_extrapolation_reproduces_a_linear_function_at_the_six_nodes() {
        // f(ξ, η) = 2 + 3ξ − 5η, connue aux points de Gauss, retrouvée aux Nodes de référence
        let f = |xi: f64, eta: f64| 2.0 + 3.0 * xi - 5.0 * eta;
        let at_gauss: Vec<f64> = TRIANGLE_3.iter().map(|g| f(g.xi, g.eta)).collect();
        let nodes = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (0.5, 0.0),
            (0.5, 0.5),
            (0.0, 0.5),
        ];
        for (k, (weights, (xi, eta))) in TRI6_EXTRAPOLATION.iter().zip(nodes).enumerate() {
            let got: f64 = weights.iter().zip(&at_gauss).map(|(w, v)| w * v).sum();
            assert!((got - f(xi, eta)).abs() < 1e-12, "Node {k}");
        }
    }

    #[test]
    fn nodal_stresses_refuse_1d_meshes_and_foreign_fields() {
        let nodes = vec![Node { x: 0.0, y: 0.0 }, Node { x: 1.0, y: 0.0 }];
        let bar = Bar {
            nodes: [0, 1],
            material: Material::new(1.0).unwrap(),
            section: Section::new(1.0).unwrap(),
        };
        let bars = Mesh::new(nodes, vec![Element::Bar(bar)]).unwrap();
        let empty = Field {
            displacements: vec![0.0; 2],
            reactions: vec![0.0; 2],
            stresses: vec![Vec::new()],
        };
        assert_eq!(
            nodal_stresses(&bars, &empty),
            Err(MailleError::NoPlaneStresses)
        );

        let tri = Tri3::new(
            [0, 1, 2],
            Material::with_poisson(1.0, 0.3).unwrap(),
            1.0,
            PlaneState::Stress,
        )
        .unwrap();
        let triangle = Mesh::new(
            vec![
                Node { x: 0.0, y: 0.0 },
                Node { x: 1.0, y: 0.0 },
                Node { x: 0.0, y: 1.0 },
            ],
            vec![Element::Tri3(tri)],
        )
        .unwrap();
        // Field d'un autre Mesh : pas de point de Gauss pour le triangle
        assert_eq!(
            nodal_stresses(&triangle, &empty),
            Err(MailleError::FieldMeshMismatch)
        );
    }

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
