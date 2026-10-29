//! Lecture des résultats pour l'affichage : un `.pvd` ou un `.vtu` devient ce que le viewer montre.
//! Partagé par `maille view` et l'espace Résultats de l'interface.

use std::error::Error;
use std::fs;
use std::path::Path;

use maille_core::results::Stress;
use maille_io::vtk::{self, VtuResults};
use maille_viewer::NodalResults;

fn read(path: &Path) -> Result<String, Box<dyn Error>> {
    fs::read_to_string(path)
        .map_err(|e| format!("impossible de lire {} : {e}", path.display()).into())
}

/// Le `.vtu` demandé : directement, ou l'Increment `index` d'un `.pvd` (le dernier par défaut).
pub fn open(path: &Path, index: Option<usize>) -> Result<VtuResults, Box<dyn Error>> {
    let vtu = match path.extension().and_then(|e| e.to_str()) {
        Some("pvd") => {
            let entries = vtk::read_pvd(&read(path)?)?;
            let k = index.unwrap_or(entries.len() - 1);
            let entry = entries.get(k).ok_or_else(|| {
                format!(
                    "incrément {k} : {} en contient {}, numérotés depuis 0",
                    path.display(),
                    entries.len()
                )
            })?;
            path.parent().unwrap_or(Path::new(".")).join(&entry.file)
        }
        Some("vtu") if index.is_none() => path.to_path_buf(),
        Some("vtu") => return Err("un .vtu n'a qu'un incrément : retirer l'indice".into()),
        _ => return Err(format!("{} : attendu un fichier .pvd ou .vtu", path.display()).into()),
    };
    Ok(vtk::read_results(&read(&vtu)?)?)
}

/// Résultats relus → ce que le viewer affiche (triangles et tétraèdres).
pub fn to_nodal(r: &VtuResults) -> Result<NodalResults, Box<dyn Error>> {
    let array = |name: &str, components: usize| {
        r.arrays
            .iter()
            .find(|a| a.name == name && a.components == components)
            .map(|a| a.values.as_slice())
    };
    let temperatures = array("temperature", 1).unwrap_or_default().to_vec();
    let stresses: Vec<Stress> = array("stress", 6)
        .map(|s| {
            s.chunks(6)
                .zip(&r.nodes)
                .map(|(c, p)| Stress {
                    x: p.x,
                    y: p.y,
                    z: p.z,
                    sigma_xx: c[0],
                    sigma_yy: c[1],
                    sigma_zz: c[2],
                    tau_xy: c[3],
                    tau_yz: c[4],
                    tau_xz: c[5],
                })
                .collect()
        })
        .unwrap_or_default();
    let displacements: Vec<[f64; 3]> = array("displacement", 3)
        .map(|d| d.chunks(3).map(|u| [u[0], u[1], u[2]]).collect())
        .unwrap_or_default();
    if temperatures.is_empty() && (stresses.is_empty() || displacements.is_empty()) {
        return Err(
            "ce fichier n'a ni « temperature », ni « displacement » et « stress » : rien à afficher"
                .into(),
        );
    }
    Ok(NodalResults {
        nodes: r.nodes.iter().map(|p| [p.x, p.y, p.z]).collect(),
        cells: r.cells.clone(),
        displacements,
        stresses,
        temperatures,
        plastic_strains: array("plastic_strain", 1).unwrap_or_default().to_vec(),
        large_strain: r.large_strain,
    })
}

#[cfg(test)]
mod tests {
    use maille_core::mesh::Node;
    use maille_io::msh::CellKind;
    use maille_io::vtk::PointArray;

    use super::*;

    fn results(kind: CellKind, arrays: Vec<PointArray>) -> VtuResults {
        let nodes = vec![
            Node::new_2d(0.0, 0.0),
            Node::new_2d(1.0, 0.0),
            Node::new_2d(0.0, 1.0),
        ];
        VtuResults {
            nodes,
            kind,
            cells: vec![vec![0, 1, 2]],
            arrays,
            units: None,
            large_strain: false,
        }
    }

    fn array(name: &str, components: usize) -> PointArray {
        PointArray {
            name: name.to_string(),
            components,
            values: (0..3 * components).map(|i| i as f64).collect(),
        }
    }

    #[test]
    fn mechanical_results_become_nodal_results() {
        let r = results(
            CellKind::Tri3,
            vec![array("displacement", 3), array("stress", 6)],
        );
        let n = to_nodal(&r).unwrap();
        assert_eq!(
            n.displacements,
            [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]]
        );
        assert_eq!(n.stresses[1].sigma_xx, 6.0);
        assert_eq!((n.stresses[1].x, n.stresses[2].y), (1.0, 1.0));
        assert!(n.temperatures.is_empty() && n.plastic_strains.is_empty());
    }

    #[test]
    fn thermal_results_need_only_the_temperature() {
        let n = to_nodal(&results(CellKind::Tri3, vec![array("temperature", 1)])).unwrap();
        assert_eq!(n.temperatures, [0.0, 1.0, 2.0]);
    }

    #[test]
    fn tetrahedra_are_shown_and_empty_files_are_refused() {
        let solid = to_nodal(&results(CellKind::Tet4, vec![array("temperature", 1)])).unwrap();
        assert_eq!(solid.temperatures.len(), 3);
        let empty = to_nodal(&results(CellKind::Tri3, vec![array("stress", 6)]));
        assert!(empty.unwrap_err().to_string().contains("rien à afficher"));
    }
}
