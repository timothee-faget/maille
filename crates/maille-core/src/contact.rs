//! Contact sans frottement d'un corps déformable contre une surface rigide analytique, par
//! pénalité : un Node esclave qui traverse la surface (interstice g < 0) est repoussé par la force
//! k (−g) n, et la tangente gagne k (n ⊗ n + g ∂n/∂x). La raideur k est automatique : voir
//! `contact_terms`.

use crate::assembly::face_convection;
use crate::error::MailleError;
use crate::mesh::{Element, Mesh};

/// Facteur c de la raideur de pénalité ε = c E / h (E : module d'Young, h : taille des Elements du
/// bord). Grand devant 1 pour que la pénétration soit négligeable, assez petit pour garder une
/// matrice bien conditionnée.
const PENALTY_FACTOR: f64 = 1000.0;

/// Surface rigide analytique. Le corps déformable reste du côté où l'interstice est positif.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RigidSurface {
    /// Demi-espace : `point` est sur le plan, `normal` (non nulle, normalisée par le calcul) pointe
    /// vers le corps déformable. g = (x − point) · normal.
    Plane { point: [f64; 3], normal: [f64; 3] },
    /// Cylindre plein de rayon `radius` autour de la droite par `point` de direction `axis` ; le corps
    /// est à l'extérieur. g = distance à l'axe − radius. En 2D, c'est un cercle de centre
    /// (`point[0]`, `point[1]`) : l'axe est z et `axis` est ignoré.
    Cylinder {
        point: [f64; 3],
        axis: [f64; 3],
        radius: f64,
    },
}

/// Contact entre les Nodes esclaves d'un bord et une surface rigide.
#[derive(Debug, Clone, PartialEq)]
pub struct Contact {
    /// Nodes esclaves : un bord du Mesh, dont les faces ont tous leurs Nodes ici (comme le Set d'un
    /// bord). Chaque Node doit être sur au moins une de ces faces.
    pub nodes: Vec<usize>,
    pub surface: RigidSurface,
    /// Déplacement de la surface rigide quand λ(t) = 1 ; elle se déplace de λ(t) × `motion`.
    pub motion: [f64; 3],
}

/// Interstice g d'un point et normale n (vers le corps), avec ∂n/∂x.
pub(crate) struct Gap {
    pub g: f64,
    pub n: [f64; 3],
    pub dn: [[f64; 3]; 3],
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let norm = dot(&a, &a).sqrt();
    a.map(|v| v / norm)
}

impl RigidSurface {
    /// Refuse une surface non finie, un plan sans normale, un cylindre sans axe ni rayon. `dim` : 2
    /// (triangles) ou 3 (tétraèdres).
    pub(crate) fn check(&self, dim: usize) -> Result<(), MailleError> {
        let finite = |v: &[f64; 3]| v.iter().all(|x| x.is_finite());
        let bad = |problem: &'static str| Err(MailleError::InvalidRigidSurface { problem });
        match self {
            Self::Plane { point, normal } => {
                if !finite(point) || !finite(normal) {
                    return bad("le plan a un point ou une normale non fini");
                }
                if normal[..dim].iter().all(|v| *v == 0.0) {
                    return bad("la normale du plan est nulle dans le plan du maillage");
                }
            }
            Self::Cylinder {
                point,
                axis,
                radius,
            } => {
                if !finite(point) || !finite(axis) || !radius.is_finite() {
                    return bad("le cylindre a une valeur non finie");
                }
                if *radius <= 0.0 {
                    return bad("le rayon du cylindre doit être strictement positif");
                }
                if dim == 3 && axis.iter().all(|v| *v == 0.0) {
                    return bad("l'axe du cylindre est nul");
                }
            }
        }
        Ok(())
    }

    /// Interstice du point `x` (z = 0 en 2D).
    pub(crate) fn gap(&self, x: [f64; 3], dim: usize) -> Gap {
        match *self {
            Self::Plane { point, normal } => {
                let mut normal = normal;
                if dim == 2 {
                    normal[2] = 0.0;
                }
                let n = unit(normal);
                let d = [x[0] - point[0], x[1] - point[1], x[2] - point[2]];
                Gap {
                    g: dot(&d, &n),
                    n,
                    dn: [[0.0; 3]; 3],
                }
            }
            Self::Cylinder {
                point,
                axis,
                radius,
            } => {
                let a = if dim == 2 {
                    [0.0, 0.0, 1.0]
                } else {
                    unit(axis)
                };
                let d = [x[0] - point[0], x[1] - point[1], x[2] - point[2]];
                let along = dot(&d, &a);
                let r: [f64; 3] = std::array::from_fn(|k| d[k] - along * a[k]);
                let distance = dot(&r, &r).sqrt();
                if distance < 1e-300 {
                    // sur l'axe : pas de direction de poussée
                    return Gap {
                        g: -radius,
                        n: [0.0; 3],
                        dn: [[0.0; 3]; 3],
                    };
                }
                let n = r.map(|v| v / distance);
                // ∂n/∂x = (I − a ⊗ a − n ⊗ n) / distance
                let dn = std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        let delta = if i == j { 1.0 } else { 0.0 };
                        (delta - a[i] * a[j] - n[i] * n[j]) / distance
                    })
                });
                Gap {
                    g: distance - radius,
                    n,
                    dn,
                }
            }
        }
    }
}

/// Par Node esclave d'un `Contact` : raideur de pénalité k = ε × aire tributaire, et cette aire
/// (longueur en 2D, épaisseur comprise).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ContactTerms {
    pub stiffness: Vec<f64>,
    pub area: Vec<f64>,
}

/// Taille d'un Element de surface : racine de l'aire d'un triangle, racine cubique du volume d'un
/// tétraèdre (sommets seulement).
fn element_size(mesh: &Mesh, i: usize) -> f64 {
    match &mesh.elements[i] {
        Element::Tet4(_) | Element::Tet10(_) => {
            let c = mesh.coords_3d(&mesh.elements[i].nodes()[..4]);
            let edge = |k: usize| -> [f64; 3] { std::array::from_fn(|c3| c[k][c3] - c[0][c3]) };
            let (a, b, d) = (edge(1), edge(2), edge(3));
            let det = a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
                + a[2] * (b[0] * d[1] - b[1] * d[0]);
            (det.abs() / 6.0).cbrt()
        }
        _ => {
            let c = mesh.coords(&mesh.elements[i].nodes()[..3]);
            let area = 0.5
                * ((c[1][0] - c[0][0]) * (c[2][1] - c[0][1])
                    - (c[2][0] - c[0][0]) * (c[1][1] - c[0][1]))
                    .abs();
            area.sqrt()
        }
    }
}

/// Aire tributaire et raideur de pénalité de chaque Node esclave. L'aire vient des faces de Mesh dont
/// tous les Nodes sont esclaves (∫ N dS, épaisseur comprise) ; ε = c E / h, E et h (taille) étant
/// ceux de l'Element de chaque face, pondérés par l'aire.
pub(crate) fn contact_terms(mesh: &Mesh, contact: &Contact) -> Result<ContactTerms, MailleError> {
    let slave: std::collections::HashMap<usize, usize> = contact
        .nodes
        .iter()
        .enumerate()
        .map(|(k, &n)| (n, k))
        .collect();
    let mut area = vec![0.0; contact.nodes.len()];
    let mut weighted = vec![0.0; contact.nodes.len()];
    for (i, e) in mesh.elements.iter().enumerate() {
        let (young, size) = (e.material().young.unwrap_or(0.0), element_size(mesh, i));
        for face in 0.. {
            let Some(positions) = e.face_positions(face) else {
                break;
            };
            let global: Vec<usize> = positions.iter().map(|&p| e.nodes()[p]).collect();
            if !global.iter().all(|n| slave.contains_key(n)) {
                continue;
            }
            let Some(terms) = face_convection(mesh, i, face, 1.0, 1.0) else {
                continue;
            };
            for (&p, share) in terms.positions.iter().zip(&terms.load) {
                let k = slave[&e.nodes()[p]];
                area[k] += share;
                weighted[k] += share * PENALTY_FACTOR * young / size;
            }
        }
    }
    if let Some(k) = area.iter().position(|a| *a <= 0.0) {
        return Err(MailleError::ContactNodeWithoutFace {
            node: contact.nodes[k],
        });
    }
    Ok(ContactTerms {
        stiffness: weighted,
        area,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_plane_gap_is_the_signed_distance_along_the_normal() {
        let plane = RigidSurface::Plane {
            point: [1.0, 2.0, 3.0],
            normal: [0.0, 0.0, 4.0],
        };
        let gap = plane.gap([5.0, 5.0, 3.5], 3);
        assert!((gap.g - 0.5).abs() < 1e-15);
        assert_eq!(gap.n, [0.0, 0.0, 1.0]);
        // en 2D, la composante z de la normale est ignorée
        let slanted = RigidSurface::Plane {
            point: [0.0; 3],
            normal: [3.0, 4.0, 9.0],
        };
        let gap = slanted.gap([1.0, 1.0, 0.0], 2);
        assert!((gap.g - 1.4).abs() < 1e-15);
    }

    #[test]
    fn the_cylinder_gap_is_the_distance_to_the_axis_minus_the_radius() {
        let cylinder = RigidSurface::Cylinder {
            point: [1.0, 0.0, 0.0],
            axis: [0.0, 0.0, 2.0],
            radius: 2.0,
        };
        // 3D : le point (4, 4, 7) est à 5 de l'axe vertical par (1, 0)
        let gap = cylinder.gap([4.0, 4.0, 7.0], 3);
        assert!((gap.g - 3.0).abs() < 1e-15);
        assert!((gap.n[0] - 0.6).abs() < 1e-15 && (gap.n[1] - 0.8).abs() < 1e-15);
        assert!(gap.n[2].abs() < 1e-15);
        // en 2D, l'axe est z quel que soit `axis`
        let tilted = RigidSurface::Cylinder {
            point: [1.0, 0.0, 0.0],
            axis: [1.0, 0.0, 0.0],
            radius: 2.0,
        };
        assert!((tilted.gap([4.0, 4.0, 0.0], 2).g - 3.0).abs() < 1e-15);
    }

    #[test]
    fn the_normal_derivative_matches_finite_differences() {
        let cylinder = RigidSurface::Cylinder {
            point: [0.5, -0.5, 0.0],
            axis: [1.0, 2.0, 2.0],
            radius: 1.0,
        };
        let x = [2.0, 1.0, 0.3];
        let gap = cylinder.gap(x, 3);
        let h = 1e-6;
        for j in 0..3 {
            let (mut plus, mut minus) = (x, x);
            plus[j] += h;
            minus[j] -= h;
            let (np, nm) = (cylinder.gap(plus, 3).n, cylinder.gap(minus, 3).n);
            for i in 0..3 {
                let numeric = (np[i] - nm[i]) / (2.0 * h);
                assert!((gap.dn[i][j] - numeric).abs() < 1e-8, "dn[{i}][{j}]");
            }
        }
    }
}
