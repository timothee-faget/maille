//! L'espace Maillage sans fenêtre : importer un fichier CAO (DXF, STL, STEP) ou un `.msh`, le mailler
//! avec la même recette que `maille mesh`, regarder la qualité. `app` ne fait que l'afficher.

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use maille_io::meshing::{self, Recipe, Source};
use maille_io::msh::{self, MshFile, Set, DOMAIN_SET};
use maille_io::IoError;
use maille_viewer::NodalResults;

use crate::files::File;
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

/// Maillage obtenu, avec sa qualité et ses Sets.
#[derive(Debug, Clone, PartialEq)]
pub struct Built {
    pub file: MshFile,
    /// Unité de longueur annoncée par un STEP.
    pub unit: Option<String>,
    pub quality: Quality,
    /// Entités de la géométrie que le mailleur a nommées (calque, `solid` ou `surface_N`, face) : ce
    /// qu'on peut désigner à l'écran, triées.
    pub entities: Vec<String>,
    /// Sets utilisateur : nom → entités regroupées, ce qui s'écrit dans `[meshing] sets`.
    pub user: BTreeMap<String, Vec<String>>,
    /// Sets tirés sur le Mesh (cadre, angle) : ils n'ont pas d'entité de géométrie, vivent dans le
    /// `.msh` et ne survivent pas à un remaillage.
    pub mesh_sets: BTreeSet<String>,
}

impl Built {
    /// Maillage tel que le mailleur ou un `.msh` le donne : tous ses Sets sont des entités.
    pub fn new(file: MshFile, unit: Option<String>) -> Self {
        let quality = quality::of(&file);
        let entities = file
            .sets
            .keys()
            .filter(|name| *name != DOMAIN_SET)
            .cloned()
            .collect();
        Self {
            file,
            unit,
            quality,
            entities,
            user: BTreeMap::new(),
            mesh_sets: BTreeSet::new(),
        }
    }

    /// Ajoute un Set tiré sur le Mesh (Nodes, Elements ou faces). Erreur si le nom est vide ou déjà
    /// pris, ou si le Set est vide.
    pub fn add_mesh_set(&mut self, name: &str, set: Set) -> Result<(), IoError> {
        let name = name.trim().to_string();
        if name.is_empty() {
            return Err(IoError::InvalidModel {
                problem: "le Set n'a pas de nom : en écrire un, par exemple « encastrement »"
                    .to_string(),
            });
        }
        if self.file.sets.contains_key(&name) || name == DOMAIN_SET {
            return Err(IoError::DuplicateSetName { name });
        }
        if set.nodes.is_empty()
            && set.faces.is_empty()
            && set.edges.is_empty()
            && set.elements.is_empty()
        {
            return Err(IoError::InvalidModel {
                problem: "rien n'est sélectionné : trace un cadre ou clique une face d'abord"
                    .to_string(),
            });
        }
        self.file.sets.insert(name.clone(), set);
        self.mesh_sets.insert(name);
        Ok(())
    }

    /// Crée un Set qui regroupe `members` (des entités de la géométrie). Erreur si le nom est vide ou
    /// déjà pris, ou si une entité n'existe pas.
    pub fn add_set(&mut self, name: &str, members: &[String]) -> Result<(), IoError> {
        if name.trim().is_empty() {
            return Err(IoError::InvalidModel {
                problem: "le Set n'a pas de nom : en écrire un, par exemple « encastrement »"
                    .to_string(),
            });
        }
        // seules les entités se regroupent : un Set utilisateur n'est pas une entité
        let entities: BTreeMap<String, Vec<String>> =
            BTreeMap::from([(name.trim().to_string(), members.to_vec())]);
        let mut grouped = self.file.clone();
        grouped
            .sets
            .retain(|k, _| self.entities.contains(k) || k == DOMAIN_SET);
        meshing::combine_sets(&mut grouped, &entities)?;
        let name = name.trim().to_string();
        if self.file.sets.contains_key(&name) {
            return Err(IoError::DuplicateSetName { name });
        }
        if let Some(set) = grouped.sets.remove(&name) {
            self.file.sets.insert(name.clone(), set);
        }
        self.user.insert(name, members.to_vec());
        Ok(())
    }

    /// Retire un Set utilisateur ; faux si ce n'en est pas un (une entité ne se supprime pas).
    pub fn remove_set(&mut self, name: &str) -> bool {
        if self.user.remove(name).is_none() && !self.mesh_sets.remove(name) {
            return false;
        }
        self.file.sets.remove(name);
        true
    }

    /// Renomme un Set utilisateur. Erreur si `name` n'en est pas un ou si le nouveau nom est vide ou
    /// déjà pris.
    pub fn rename_set(&mut self, name: &str, new_name: &str) -> Result<(), IoError> {
        let new_name = new_name.trim();
        if new_name.is_empty() {
            return Err(IoError::InvalidModel {
                problem: "le Set n'a pas de nom : en écrire un".to_string(),
            });
        }
        if new_name == name {
            return Ok(());
        }
        if self.file.sets.contains_key(new_name) || new_name == DOMAIN_SET {
            return Err(IoError::DuplicateSetName {
                name: new_name.to_string(),
            });
        }
        if !self.user.contains_key(name) && !self.mesh_sets.contains(name) {
            return Err(IoError::InvalidModel {
                problem: format!(
                    "« {name} » n'est pas un Set créé ici : seuls ceux-là se renomment"
                ),
            });
        }
        let set = self.file.sets.remove(name).unwrap_or_default();
        if self.mesh_sets.remove(name) {
            self.mesh_sets.insert(new_name.to_string());
        } else {
            let members = self.user.remove(name).unwrap_or_default();
            self.user.insert(new_name.to_string(), members);
        }
        self.file.sets.insert(new_name.to_string(), set);
        Ok(())
    }

    /// Entités de la géométrie que désigne un Set : lui-même pour une entité, ses membres pour un Set
    /// utilisateur.
    pub fn entities_of(&self, set: &str) -> Vec<String> {
        match self.user.get(set) {
            Some(members) => members.clone(),
            None if self.entities.iter().any(|e| e == set) => vec![set.to_string()],
            None => Vec::new(),
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

    /// Noms des Sets d'Elements, le domaine compris.
    pub fn element_sets(&self) -> Vec<&str> {
        self.file
            .sets
            .iter()
            .filter(|(_, set)| !set.elements.is_empty())
            .map(|(name, _)| name.as_str())
            .collect()
    }

    /// Noms des Sets qu'une force linéique peut viser : des bords (faces) en 2D, des arêtes en 3D.
    pub fn edge_sets(&self) -> Vec<&str> {
        let dim = self.file.kind.dim();
        self.file
            .sets
            .iter()
            .filter(|(_, set)| {
                if dim == 2 {
                    !set.faces.is_empty()
                } else {
                    !set.edges.is_empty()
                }
            })
            .map(|(name, _)| name.as_str())
            .collect()
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

/// Maille la géométrie à la taille `size` : la même recette que `maille mesh`, donc le même maillage,
/// puis regroupe les entités des Sets utilisateur `user` (ils survivent au remaillage).
pub fn mesh(
    cad: &Cad,
    size: f64,
    quadratic: bool,
    user: &BTreeMap<String, Vec<String>>,
) -> Result<Built, IoError> {
    let recipe = Recipe {
        size: Some(size),
        quadratic,
        ..Recipe::default()
    };
    let meshed = meshing::mesh(cad.source, &cad.bytes, &recipe)?;
    let mut built = Built::new(meshed.file, meshed.unit);
    meshing::combine_sets(&mut built.file, user)?;
    built.user = user.clone();
    Ok(built)
}

/// La recette du maillage telle que le TOML du modèle l'écrit : `[meshing]` avec les Sets
/// utilisateur. La source est donnée relativement au dossier du modèle.
pub fn recipe(
    source: &str,
    size: f64,
    quadratic: bool,
    user: &BTreeMap<String, Vec<String>>,
) -> maille_io::model::MeshingDef {
    maille_io::model::MeshingDef {
        source: source.to_string(),
        size: Some(size),
        quadratic: quadratic.then_some(true),
        feature_angle: None,
        sets: user.clone(),
    }
}

/// Le maillage à écrire à côté de la source, comme `maille mesh` (`<source>.msh`).
pub fn mesh_file(cad: &Cad, built: &Built) -> Result<File, Box<dyn Error>> {
    Ok((
        cad.path.with_extension("msh"),
        msh::write(&built.file)?.into_bytes(),
    ))
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
        let built = mesh(&c, c.suggested, false, &BTreeMap::new()).unwrap();
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
        let built = mesh(&c, c.suggested, true, &BTreeMap::new()).unwrap();
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
        let built = mesh(&cube, cube.suggested, false, &BTreeMap::new()).unwrap();
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

    fn names(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn a_user_set_is_added_renamed_and_removed() {
        let c = cad("plaque.dxf");
        let mut built = mesh(&c, c.suggested, false, &BTreeMap::new()).unwrap();
        let before = built.file.sets.len();
        built
            .add_set("verticaux", &names(&["gauche", "droite"]))
            .unwrap();
        assert_eq!(built.file.sets.len(), before + 1);
        assert_eq!(built.user["verticaux"], names(&["gauche", "droite"]));
        // une entité de plus dans les Sets de bord : l'ensemble des entités ne change pas
        assert_eq!(built.entities.len(), before - 1);
        assert_eq!(built.entities_of("verticaux"), names(&["gauche", "droite"]));
        assert_eq!(built.entities_of("bas"), names(&["bas"]));
        assert!(built.entities_of("inconnu").is_empty());

        built.rename_set("verticaux", "cotes").unwrap();
        assert!(
            built.file.sets.contains_key("cotes") && !built.file.sets.contains_key("verticaux")
        );
        assert!(built.user.contains_key("cotes"));
        assert!(
            built.rename_set("bas", "autre").is_err(),
            "une entité ne se renomme pas"
        );

        assert!(built.remove_set("cotes"));
        assert!(!built.file.sets.contains_key("cotes") && built.user.is_empty());
        assert!(!built.remove_set("bas"), "une entité ne se supprime pas");
        assert!(built.file.sets.contains_key("bas"));
    }

    #[test]
    fn a_bad_set_is_refused_and_leaves_the_mesh_alone() {
        let c = cad("plaque.dxf");
        let mut built = mesh(&c, c.suggested, false, &BTreeMap::new()).unwrap();
        let before = built.clone();
        assert!(built.add_set("  ", &names(&["bas"])).is_err());
        assert!(built.add_set("x", &[]).is_err());
        assert!(built.add_set("x", &names(&["fantome"])).is_err());
        assert!(
            built.add_set("bas", &names(&["haut"])).is_err(),
            "nom d'entité déjà pris"
        );
        built.add_set("a", &names(&["bas"])).unwrap();
        assert!(
            built.add_set("a", &names(&["haut"])).is_err(),
            "nom déjà pris"
        );
        assert!(
            built.add_set("b", &names(&["a"])).is_err(),
            "un Set ne regroupe que des entités"
        );
        assert!(built.rename_set("a", "bas").is_err());
        assert!(built.rename_set("a", " ").is_err());
        built.remove_set("a");
        assert_eq!(built, before);
    }

    #[test]
    fn user_sets_survive_remeshing_and_write_into_the_recipe() {
        let c = cad("plaque.dxf");
        let user: BTreeMap<String, Vec<String>> =
            BTreeMap::from([("verticaux".to_string(), names(&["gauche", "droite"]))]);
        let coarse = mesh(&c, c.suggested, false, &user).unwrap();
        let fine = mesh(&c, c.suggested / 2.0, true, &user).unwrap();
        for built in [&coarse, &fine] {
            let (g, d) = (&built.file.sets["gauche"], &built.file.sets["droite"]);
            assert_eq!(
                built.file.sets["verticaux"].faces.len(),
                g.faces.len() + d.faces.len()
            );
            assert_eq!(built.user, user);
        }
        assert!(fine.file.nodes.len() > coarse.file.nodes.len());
        // la recette qui s'écrira dans le TOML du modèle
        let def = recipe("plaque.dxf", c.suggested, true, &user);
        assert_eq!(
            (def.source.as_str(), def.quadratic),
            ("plaque.dxf", Some(true))
        );
        assert_eq!(def.sets, user);
        assert_eq!(recipe("plaque.dxf", 1.0, false, &user).quadratic, None);
    }

    #[test]
    fn the_msh_goes_next_to_the_source() {
        let c = match import(&example("plaque.dxf")).unwrap() {
            Imported::Cad(c) => c,
            other => panic!("{other:?}"),
        };
        let built = mesh(&c, c.suggested, false, &BTreeMap::new()).unwrap();
        let (path, bytes) = mesh_file(&c, &built).unwrap();
        assert_eq!(path, example("plaque.msh"));
        let again = msh::read(std::str::from_utf8(&bytes).unwrap()).unwrap();
        assert_eq!(again, built.file);
    }

    #[test]
    fn a_set_picked_on_the_mesh_is_added_renamed_removed_and_written_in_the_msh() {
        let Imported::Mesh(_, mut built) = import(&example("plaque.msh")).unwrap() else {
            panic!("plaque.msh");
        };
        let nodes = crate::mesh_pick::node_set(&[0, 1, 2]);
        built.add_mesh_set("  appui ", nodes.clone()).unwrap();
        assert!(built.mesh_sets.contains("appui") && built.sets().contains(&"appui"));
        assert!(matches!(
            built.add_mesh_set("appui", nodes.clone()),
            Err(IoError::DuplicateSetName { .. })
        ));
        assert!(built.add_mesh_set("gauche", nodes.clone()).is_err());
        assert!(built.add_mesh_set("   ", nodes).is_err());
        assert!(built
            .add_mesh_set("vide", crate::mesh_pick::node_set(&[]))
            .is_err());

        // le Set suit le fichier `.msh`
        let again = msh::read(&msh::write(&built.file).unwrap()).unwrap();
        assert_eq!(again.sets["appui"].nodes, vec![0, 1, 2]);

        built.rename_set("appui", "base").unwrap();
        assert!(built.mesh_sets.contains("base") && !built.file.sets.contains_key("appui"));
        assert!(
            built.rename_set("gauche", "x").is_err(),
            "une entité ne se renomme pas"
        );
        assert!(built.remove_set("base") && !built.file.sets.contains_key("base"));
        assert!(!built.remove_set("gauche"));
    }
}
