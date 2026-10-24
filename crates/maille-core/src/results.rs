use crate::error::MailleError;
use crate::mesh::{Element, Mesh};

/// Contraintes en un point (point de Gauss ou Node).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stress {
    /// Position (x, y, z) du point.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub sigma_xx: f64,
    pub sigma_yy: f64,
    /// Triangles : nulle en contraintes planes, ν (σxx + σyy) en déformations planes.
    pub sigma_zz: f64,
    pub tau_xy: f64,
    /// Nulle pour les triangles.
    pub tau_yz: f64,
    /// Nulle pour les triangles.
    pub tau_xz: f64,
}

impl Stress {
    /// Contrainte équivalente de von Mises.
    pub fn von_mises(&self) -> f64 {
        let (a, b, c) = (self.sigma_xx, self.sigma_yy, self.sigma_zz);
        let shear =
            self.tau_xy * self.tau_xy + self.tau_yz * self.tau_yz + self.tau_xz * self.tau_xz;
        // ≥ 0 en exact ; max(0) évite un NaN quand l'arrondi rend négatif un état presque hydrostatique
        (a * a + b * b + c * c - a * b - b * c - c * a + 3.0 * shear)
            .max(0.0)
            .sqrt()
    }
}

/// Flux de chaleur q = −k ∇T en un point de Gauss (puissance par unité de surface).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatFlux {
    /// Position (x, y, z) du point.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub qx: f64,
    pub qy: f64,
    /// Nul pour les triangles.
    pub qz: f64,
}

/// Champs de résultats. Selon l'`Analysis`, seuls les champs mécaniques (statique) ou thermiques
/// sont remplis ; les autres restent vides.
#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    /// Statique : un déplacement par DOF.
    pub displacements: Vec<f64>,
    /// Aux DOF contraints, zéro ailleurs. Statique : forces de réaction ; thermique : puissance
    /// fournie au Mesh par la température imposée (négative si la chaleur sort).
    pub reactions: Vec<f64>,
    /// Statique : contraintes aux points de Gauss de chaque Element (indexé comme `Mesh::elements`).
    /// Vide pour les Bar et Beam.
    pub stresses: Vec<Vec<Stress>>,
    /// Thermique : une température par Node.
    pub temperatures: Vec<f64>,
    /// Thermique : flux de chaleur aux points de Gauss de chaque Element.
    pub heat_fluxes: Vec<Vec<HeatFlux>>,
}

/// Extrapolation du Tri3 : son unique point de Gauss donne la valeur de ses 3 Nodes.
const TRI3_EXTRAPOLATION: [&[f64]; 3] = [&[1.0]; 3];

/// Extrapolation du Tet4 : même principe que le Tri3.
const TET4_EXTRAPOLATION: [&[f64]; 4] = [&[1.0]; 4];

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

/// Extrapolation du Tet10 depuis les 4 points de Gauss de `TETRA_4` (le point g est près du sommet g,
/// coordonnée barycentrique a sur ce sommet, b sur les trois autres). Pour un champ linéaire de
/// valeurs c_k aux sommets, f_g = (a − b) c_g + b Σc et Σf = Σc, donc c_k = Σ_g (δ_kg − b) f_g / (a − b),
/// avec 1/(a − b) = √5 : poids (3√5 + 1)/4 sur son point, −(√5 − 1)/4 sur les autres. Un milieu
/// d'arête prend la moyenne de ses deux sommets : φ/2 = (√5 + 1)/4 sur les deux points proches.
const TET10_VERTEX: f64 = 1.927_050_983_124_842_4;
const TET10_MIDDLE: f64 = 0.809_016_994_374_947_5;
const TET10_OTHER: f64 = -0.309_016_994_374_947_45;
const TET10_EXTRAPOLATION: [&[f64]; 10] = [
    &[TET10_VERTEX, TET10_OTHER, TET10_OTHER, TET10_OTHER],
    &[TET10_OTHER, TET10_VERTEX, TET10_OTHER, TET10_OTHER],
    &[TET10_OTHER, TET10_OTHER, TET10_VERTEX, TET10_OTHER],
    &[TET10_OTHER, TET10_OTHER, TET10_OTHER, TET10_VERTEX],
    &[TET10_MIDDLE, TET10_MIDDLE, TET10_OTHER, TET10_OTHER], // milieu 0-1
    &[TET10_OTHER, TET10_MIDDLE, TET10_MIDDLE, TET10_OTHER], // milieu 1-2
    &[TET10_MIDDLE, TET10_OTHER, TET10_MIDDLE, TET10_OTHER], // milieu 2-0
    &[TET10_MIDDLE, TET10_OTHER, TET10_OTHER, TET10_MIDDLE], // milieu 0-3
    &[TET10_OTHER, TET10_MIDDLE, TET10_OTHER, TET10_MIDDLE], // milieu 1-3
    &[TET10_OTHER, TET10_OTHER, TET10_MIDDLE, TET10_MIDDLE], // milieu 2-3
];

/// Contraintes lissées aux Nodes (indexées comme `Mesh::nodes`, `x`, `y`, `z` = position du Node) :
/// extrapolées des points de Gauss vers les Nodes de chaque Element, puis moyennées composante par
/// composante entre les Elements qui partagent un Node. Un Node isolé reçoit des contraintes nulles.
pub fn nodal_stresses(mesh: &Mesh, field: &Field) -> Result<Vec<Stress>, MailleError> {
    if !field.temperatures.is_empty() {
        return Err(MailleError::NoStresses);
    }
    if field.stresses.len() != mesh.elements.len() {
        return Err(MailleError::FieldMeshMismatch);
    }
    // somme de (σxx, σyy, σzz, τxy, τyz, τxz) et nombre d'Elements, par Node
    let mut sum = vec![[0.0; 6]; mesh.nodes.len()];
    let mut count = vec![0_usize; mesh.nodes.len()];
    for (e, points) in mesh.elements.iter().zip(&field.stresses) {
        let extrapolation: &[&[f64]] = match e {
            Element::Tri3(_) => &TRI3_EXTRAPOLATION,
            Element::Tri6(_) => &TRI6_EXTRAPOLATION,
            Element::Tet4(_) => &TET4_EXTRAPOLATION,
            Element::Tet10(_) => &TET10_EXTRAPOLATION,
            Element::Bar(_) | Element::Beam(_) => return Err(MailleError::NoStresses),
        };
        if points.len() != extrapolation[0].len() {
            return Err(MailleError::FieldMeshMismatch);
        }
        for (&n, weights) in e.nodes().iter().zip(extrapolation) {
            for (p, w) in points.iter().zip(weights.iter()) {
                let components = [
                    p.sigma_xx, p.sigma_yy, p.sigma_zz, p.tau_xy, p.tau_yz, p.tau_xz,
                ];
                for (s, c) in sum[n].iter_mut().zip(components) {
                    *s += w * c;
                }
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
            Stress {
                x: node.x,
                y: node.y,
                z: node.z,
                sigma_xx: s[0] / c,
                sigma_yy: s[1] / c,
                sigma_zz: s[2] / c,
                tau_xy: s[3] / c,
                tau_yz: s[4] / c,
                tau_xz: s[5] / c,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::{Material, Section};
    use crate::mesh::{Bar, Node, PlaneState, Tri3};
    use crate::quadrature::{TETRA_4, TRIANGLE_3};

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
    fn tet10_extrapolation_reproduces_a_linear_function_at_the_ten_nodes() {
        let f = |xi: f64, eta: f64, zeta: f64| 2.0 + 3.0 * xi - 5.0 * eta + 7.0 * zeta;
        let at_gauss: Vec<f64> = TETRA_4.iter().map(|g| f(g.xi, g.eta, g.zeta)).collect();
        let nodes = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.5, 0.0, 0.0),
            (0.5, 0.5, 0.0),
            (0.0, 0.5, 0.0),
            (0.0, 0.0, 0.5),
            (0.5, 0.0, 0.5),
            (0.0, 0.5, 0.5),
        ];
        for (k, (weights, (xi, eta, zeta))) in TET10_EXTRAPOLATION.iter().zip(nodes).enumerate() {
            let got: f64 = weights.iter().zip(&at_gauss).map(|(w, v)| w * v).sum();
            assert!((got - f(xi, eta, zeta)).abs() < 1e-12, "Node {k}");
        }
    }

    #[test]
    fn nodal_stresses_refuse_1d_meshes_and_foreign_fields() {
        let nodes = vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)];
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
            temperatures: Vec::new(),
            heat_fluxes: Vec::new(),
        };
        assert_eq!(nodal_stresses(&bars, &empty), Err(MailleError::NoStresses));

        let tri = Tri3::new(
            [0, 1, 2],
            Material::with_poisson(1.0, 0.3).unwrap(),
            1.0,
            PlaneState::Stress,
        )
        .unwrap();
        let triangle = Mesh::new(
            vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.0, 1.0),
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

    fn stress(sxx: f64, syy: f64, txy: f64, szz: f64) -> Stress {
        Stress {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            sigma_xx: sxx,
            sigma_yy: syy,
            sigma_zz: szz,
            tau_xy: txy,
            tau_yz: 0.0,
            tau_xz: 0.0,
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
    fn von_mises_of_out_of_plane_shear_is_sqrt3_tau() {
        let s = Stress {
            tau_yz: 30.0,
            tau_xz: 40.0,
            ..stress(0.0, 0.0, 0.0, 0.0)
        };
        assert!((s.von_mises() - 50.0 * 3.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn von_mises_accounts_for_sigma_zz_in_plane_strain() {
        // σxx = σ, σzz = νσ (déformations planes) : σ √(1 − ν + ν²)
        let nu: f64 = 0.3;
        let expected = 100.0 * (1.0 - nu + nu * nu).sqrt();
        assert!((stress(100.0, 0.0, 0.0, nu * 100.0).von_mises() - expected).abs() < 1e-12);
    }
}
