//! L'espace Maillage sans fenêtre : importer un fichier CAO (DXF, STL, STEP) ou un `.msh`, le mailler
//! avec la même recette que `maille mesh`, regarder la qualité. `app` ne fait que l'afficher.

use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use maille_io::meshing::{self, Recipe, Source};
use maille_io::msh::{self, MshFile, DOMAIN_SET};
use maille_io::IoError;
use maille_viewer::NodalResults;

use crate::quality::{self, Quality};

/// Un fichier CAO lu, prêt à être maillé.
#[derive(Debug, Clone, PartialEq)]
pub struct Cad {
    pub path: PathBuf,
    pub source: Source,
    bytes: Vec<u8>,
    /// Taille h proposée au départ, d'après l'étendue de la pièce.
    pub suggested: f64,
}

/// Ce que l'import a trouvé.
#[derive(Debug, Clone, PartialEq)]
pub enum Imported {
    /// Une géométrie à mailler.
    Cad(Cad),
    /// Un `.msh` : déjà maillé, avec ses Sets, mais sans géométrie pour le remailler.
    Mesh(PathBuf, Box<Built>),
}

/// Maillage obtenu, avec sa qualité.
#[derive(Debug, Clone, PartialEq)]
pub struct Built {
    pub file: MshFile,
    /// Unité de longueur annoncée par un STEP.
    pub unit: Option<String>,
    pub quality: Quality,
}

impl Built {
    pub fn new(file: MshFile, unit: Option<String>) -> Self {
        let quality = quality::of(&file);
        Self {
            file,
            unit,
            quality,
        }
    }

    /// Ce que le viewer montre : le maillage, colorié par son défaut de forme (1 − qualité : rouge
    /// pour un Element écrasé).
    pub fn nodal(&self) -> NodalResults {
        NodalResults {
            nodes: self.file.nodes.iter().map(|p| [p.x, p.y, p.z]).collect(),
            cells: self.file.cells.clone(),
            distortion: self.quality.per_node.iter().map(|q| 1.0 - q).collect(),
            ..NodalResults::default()
        }
    }

    /// Noms des Sets de bord (tous sauf le domaine).
    pub fn sets(&self) -> Vec<&str> {
        self.file
            .sets
            .keys()
            .filter(|name| *name != DOMAIN_SET)
            .map(String::as_str)
            .collect()
    }
}

/// Lit un fichier : un `.msh`, ou un fichier CAO dont le format vient de l'extension.
pub fn import(path: &Path) -> Result<Imported, Box<dyn Error>> {
    let read = |what: &str| {
        fs::read(path).map_err(|e| format!("impossible de lire {what} {} : {e}", path.display()))
    };
    let is_msh = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("msh"));
    if is_msh {
        let text = String::from_utf8_lossy(&read("le maillage")?).into_owned();
        let file = msh::read(&text)?;
        return Ok(Imported::Mesh(
            path.to_path_buf(),
            Box::new(Built::new(file, None)),
        ));
    }
    let source = Source::of(path).ok_or_else(|| {
        format!(
            "{} : format non pris en charge ; importer un dessin .dxf, une surface .stl, une pièce .step ou un maillage .msh",
            path.display()
        )
    })?;
    let bytes = read("le fichier")?;
    let suggested = meshing::suggested_size(source, meshing::extent(source, &bytes)?);
    Ok(Imported::Cad(Cad {
        path: path.to_path_buf(),
        source,
        bytes,
        suggested,
    }))
}

/// Maille la géométrie à la taille `size` : la même recette que `maille mesh`, donc le même maillage.
pub fn mesh(cad: &Cad, size: f64, quadratic: bool) -> Result<Built, IoError> {
    let recipe = Recipe {
        size: Some(size),
        quadratic,
        ..Recipe::default()
    };
    let meshed = meshing::mesh(cad.source, &cad.bytes, &recipe)?;
    Ok(Built::new(meshed.file, meshed.unit))
}

/// Écrit le maillage à côté de la source, comme `maille mesh` (`<source>.msh`) ; rend le chemin.
pub fn save(cad: &Cad, built: &Built) -> Result<PathBuf, Box<dyn Error>> {
    let path = cad.path.with_extension("msh");
    fs::write(&path, msh::write(&built.file)?)
        .map_err(|e| format!("impossible d'écrire {} : {e}", path.display()))?;
    Ok(path)
}

/// Une phrase sur le maillage : Nodes, Elements, ordre.
pub fn summary(built: &Built) -> String {
    use maille_io::msh::CellKind::{Tet10, Tet4, Tri3, Tri6};
    let (what, order) = match built.file.kind {
        Tri3 => ("triangles", "ordre 1"),
        Tri6 => ("triangles", "ordre 2"),
        Tet4 => ("tétraèdres", "ordre 1"),
        Tet10 => ("tétraèdres", "ordre 2"),
    };
    format!(
        "{} Nodes, {} {what} ({order})",
        built.file.nodes.len(),
        built.file.cells.len()
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_io::msh::CellKind;

    fn example(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../examples")
            .join(name)
    }

    fn cad(name: &str) -> Cad {
        match import(&example(name)).unwrap() {
            Imported::Cad(c) => c,
            other => panic!("{name} : {other:?}"),
        }
    }

    #[test]
    fn a_dxf_is_imported_then_meshed_into_a_colored_mesh() {
        let c = cad("plaque.dxf");
        assert_eq!(c.source, Source::Dxf);
        assert!(c.suggested > 0.0);
        let built = mesh(&c, c.suggested, false).unwrap();
        assert_eq!(built.file.kind, CellKind::Tri3);
        assert_eq!(
            built.sets().len(),
            4,
            "un Set par calque : {:?}",
            built.sets()
        );
        let nodal = built.nodal();
        assert_eq!(nodal.distortion.len(), nodal.nodes.len());
        assert_eq!(nodal.cells.len(), built.file.cells.len());
        assert!(nodal.distortion.iter().all(|d| (0.0..=1.0).contains(d)));
        assert!(built.quality.min > 0.3, "{:?}", built.quality);
        assert!(summary(&built).contains("triangles (ordre 1)"));
    }

    #[test]
    fn the_quadratic_recipe_gives_second_order_elements() {
        let c = cad("plaque.dxf");
        let built = mesh(&c, c.suggested, true).unwrap();
        assert_eq!(built.file.kind, CellKind::Tri6);
        assert!(summary(&built).contains("ordre 2"));
    }

    #[test]
    fn a_msh_is_opened_as_it_is_with_its_sets() {
        match import(&example("tube.msh")).unwrap() {
            Imported::Mesh(path, built) => {
                assert!(path.ends_with("tube.msh"));
                assert!(built.sets().contains(&"interieur"), "{:?}", built.sets());
                assert!(built.quality.elements > 0);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn an_stl_and_a_step_are_imported_with_a_suggested_size() {
        for name in ["cube.stl", "sphere.step"] {
            let c = cad(name);
            assert!(c.suggested > 0.0, "{name}");
        }
        let cube = cad("cube.stl");
        let built = mesh(&cube, cube.suggested, false).unwrap();
        assert_eq!(built.file.kind, CellKind::Tet4);
        assert!(built.quality.elements > 0 && built.quality.min > 0.0);
    }

    #[test]
    fn unknown_formats_and_missing_files_say_what_to_do() {
        let e = import(&example("README.md")).unwrap_err().to_string();
        assert!(e.contains("format non pris en charge"), "{e}");
        let e = import(&example("absent.dxf")).unwrap_err().to_string();
        assert!(e.contains("impossible de lire"), "{e}");
    }

    #[test]
    fn saving_writes_the_msh_next_to_the_source() {
        let dir = std::env::temp_dir().join(format!("maille-gui-save-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        fs::copy(example("plaque.dxf"), dir.join("plaque.dxf")).unwrap();
        let c = match import(&dir.join("plaque.dxf")).unwrap() {
            Imported::Cad(c) => c,
            other => panic!("{other:?}"),
        };
        let built = mesh(&c, c.suggested, false).unwrap();
        let path = save(&c, &built).unwrap();
        assert_eq!(path, dir.join("plaque.msh"));
        let again = msh::read(&fs::read_to_string(&path).unwrap()).unwrap();
        assert_eq!(again, built.file);
        let _ = fs::remove_dir_all(&dir);
    }
}
