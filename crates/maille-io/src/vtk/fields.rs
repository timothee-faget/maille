//! Grandeurs aux Nodes d'un Field du cœur.

use maille_core::mesh::Mesh;
use maille_core::results::{nodal_equivalent_plastic_strains, nodal_stresses, Field};

use super::{vtk_error, PointArray};
use crate::IoError;

fn array(name: &str, components: usize, values: Vec<f64>) -> PointArray {
    PointArray {
        name: name.to_string(),
        components,
        values,
    }
}

/// Une valeur par DOF, `dofs` DOF par Node, complétée par des zéros jusqu'à 3 composantes.
fn vector(name: &str, values: &[f64], dofs: usize) -> PointArray {
    let padded = values
        .chunks(dofs)
        .flat_map(|n| (0..3).map(|k| n.get(k).copied().unwrap_or(0.0)))
        .collect();
    array(name, 3, padded)
}

/// Grandeurs écrites pour un Field.
///
/// - thermique : `temperature` et `reaction` (puissance fournie aux Nodes imposés) ;
/// - statique : `displacement` et `reaction` (3 composantes), `stress` (6 composantes dans
///   l'ordre de VTK : xx, yy, zz, xy, yz, xz) et `von_mises` lissés aux Nodes (`nodal_stresses`),
///   et `plastic_strain` (p lissée) si le calcul est plastique.
///
/// 2D et 3D seulement : les Bar et Beam n'ont pas de contraintes aux Nodes.
pub fn result_arrays(mesh: &Mesh, field: &Field) -> Result<Vec<PointArray>, IoError> {
    let n = mesh.nodes.len();
    if !field.temperatures.is_empty() {
        return Ok(vec![
            array("temperature", 1, field.temperatures.clone()),
            array("reaction", 1, field.reactions.clone()),
        ]);
    }
    if n == 0 || !field.displacements.len().is_multiple_of(n) {
        return Err(vtk_error(
            "le Field ne correspond pas au Mesh : même maillage requis".to_string(),
        ));
    }
    let dofs = field.displacements.len() / n;
    if !(2..=3).contains(&dofs) {
        return Err(vtk_error(
            "seuls les maillages 2D et 3D (triangles, tétraèdres) s'écrivent en VTK".to_string(),
        ));
    }
    let stresses = nodal_stresses(mesh, field)?;
    let mut arrays = vec![
        vector("displacement", &field.displacements, dofs),
        vector("reaction", &field.reactions, dofs),
        array(
            "stress",
            6,
            stresses
                .iter()
                .flat_map(|s| {
                    [
                        s.sigma_xx, s.sigma_yy, s.sigma_zz, s.tau_xy, s.tau_yz, s.tau_xz,
                    ]
                })
                .collect(),
        ),
        array(
            "von_mises",
            1,
            stresses.iter().map(|s| s.von_mises()).collect(),
        ),
    ];
    if !field.equivalent_plastic_strains.is_empty() {
        arrays.push(array(
            "plastic_strain",
            1,
            nodal_equivalent_plastic_strains(mesh, field)?,
        ));
    }
    Ok(arrays)
}
