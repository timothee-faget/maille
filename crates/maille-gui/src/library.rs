//! La bibliothèque de Materials dans l'interface : la livrée plus celle de l'utilisateur
//! (`~/.config/maille/materials.toml`), et le passage d'un Material de la bibliothèque (en SI) au
//! formulaire du Modèle, dans le système d'unités choisi. Sans fenêtre.

use std::fs;
use std::path::{Path, PathBuf};

use maille_io::materials::{Library, LibraryMaterial, UnitSystem as LibraryUnits};

use crate::model_space::{MaterialForm, UnitSystem};

/// Le système d'unités du brouillon, tel que la bibliothèque sait convertir. Un système libre en
/// texte (« mm, N, s ») marche ; un autre texte ne donne rien.
pub fn library_units(units: &UnitSystem) -> Option<LibraryUnits> {
    match units {
        UnitSystem::Si => Some(LibraryUnits::SI),
        UnitSystem::MmNMpa => Some(LibraryUnits::MM_N_S),
        UnitSystem::Custom(text) => LibraryUnits::parse(text).ok(),
    }
}

/// Le formulaire d'un Material de la bibliothèque, et sa limite d'élasticité dans le système `units`.
/// Le Material est élastique et thermique : la plasticité reste à la demande de l'utilisateur, car
/// elle n'existe pas en contraintes planes.
pub fn material_form(
    name: &str,
    material: &LibraryMaterial,
    units: LibraryUnits,
) -> (MaterialForm, Option<f64>) {
    let def = material.to_def(units);
    let form = MaterialForm {
        name: name.to_string(),
        young: def.young.unwrap_or(0.0),
        poisson: def.poisson.unwrap_or(0.0),
        conductivity: def.conductivity,
        density: def.density,
        specific_heat: def.specific_heat,
        expansion: def.expansion,
        yield_stress: None,
        hardening: None,
        neo_hookean: false,
    };
    (form, material.yield_stress_in(units))
}

/// Le fichier de Materials de l'utilisateur : `$XDG_CONFIG_HOME/maille/materials.toml`, sinon
/// `$HOME/.config/maille/materials.toml`.
pub fn user_path(xdg_config_home: Option<&str>, home: Option<&str>) -> Option<PathBuf> {
    let base = match (xdg_config_home, home) {
        (Some(xdg), _) if !xdg.is_empty() => PathBuf::from(xdg),
        (_, Some(home)) if !home.is_empty() => Path::new(home).join(".config"),
        _ => return None,
    };
    Some(base.join("maille").join("materials.toml"))
}

/// Le fichier de l'utilisateur d'après l'environnement.
pub fn default_user_path() -> Option<PathBuf> {
    user_path(
        std::env::var("XDG_CONFIG_HOME").ok().as_deref(),
        std::env::var("HOME").ok().as_deref(),
    )
}

/// Le fichier de l'utilisateur seul : vide s'il n'existe pas encore.
fn read_user(path: &Path) -> Result<Library, String> {
    match fs::read_to_string(path) {
        Ok(text) => Library::parse(&text).map_err(|e| format!("{} : {e}", path.display())),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(Library::default()),
        Err(e) => Err(format!("{} : {e}", path.display())),
    }
}

/// La bibliothèque livrée, complétée par le fichier de l'utilisateur. Si ce fichier est illisible, on
/// garde la livrée et le message dit pourquoi.
pub fn load(path: Option<&Path>) -> (Library, Option<String>) {
    let builtin = Library::builtin();
    let Some(path) = path else {
        return (builtin, None);
    };
    match read_user(path) {
        Ok(user) => (builtin.merged(user), None),
        Err(e) => (
            builtin,
            Some(format!("Ma bibliothèque n'a pas pu être lue ({e}).")),
        ),
    }
}

/// Ajoute `material` sous `name` au fichier de l'utilisateur (créé au besoin) et rend la
/// bibliothèque rechargée. Un fichier existant mais illisible n'est jamais écrasé.
pub fn save(
    path: &Path,
    name: &str,
    material: LibraryMaterial,
) -> Result<(Library, Option<String>), String> {
    let mut user = read_user(path)?;
    user.insert(name, material).map_err(|e| e.to_string())?;
    let text = user.to_toml().map_err(|e| e.to_string())?;
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir).map_err(|e| format!("{} : {e}", dir.display()))?;
    }
    fs::write(path, text).map_err(|e| format!("{} : {e}", path.display()))?;
    Ok(load(Some(path)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_file(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("maille-lib-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        dir.join("maille").join("materials.toml")
    }

    #[test]
    fn the_units_of_the_draft_map_to_the_library() {
        assert_eq!(library_units(&UnitSystem::Si), Some(LibraryUnits::SI));
        assert_eq!(
            library_units(&UnitSystem::MmNMpa),
            Some(LibraryUnits::MM_N_S)
        );
        let custom = UnitSystem::Custom("cm, kN, s".into());
        assert!(library_units(&custom).is_some());
        assert_eq!(
            library_units(&UnitSystem::Custom("n'importe quoi".into())),
            None
        );
    }

    #[test]
    fn a_library_steel_fills_the_form_in_the_chosen_units_without_plasticity() {
        let library = Library::builtin();
        let steel = library.get("acier").unwrap();
        let (form, yield_stress) = material_form("acier", steel, LibraryUnits::MM_N_S);
        assert!((form.young - 210e3).abs() < 1e-6);
        assert_eq!(form.poisson, 0.3);
        assert!(form.density.is_some_and(|d| (d - 7.85e-9).abs() < 1e-18));
        assert_eq!(form.yield_stress, None);
        assert!(yield_stress.is_some_and(|s| (s - 235.0).abs() < 1e-9));
        let (si, _) = material_form("acier", steel, LibraryUnits::SI);
        assert!((si.young - 210e9).abs() < 1.0);
    }

    #[test]
    fn a_library_steel_is_the_hand_made_steel_and_keeps_the_model_linear() {
        let mut draft = crate::model_space::Draft::default();
        let library = Library::builtin();
        let units = library_units(&draft.units).unwrap();
        let (form, _) = material_form("acier", library.get("acier").unwrap(), units);
        let by_hand = draft.steel();
        assert!((form.young - by_hand.young).abs() < 1e-6 && form.poisson == by_hand.poisson);
        draft.material = Some(form);
        // pas de plasticité d'office : un calcul élastique en contraintes planes reste possible
        assert!(!draft.is_nonlinear());
    }

    #[test]
    fn the_user_file_lives_under_the_config_directory() {
        assert_eq!(
            user_path(Some("/cfg"), Some("/home/tim")),
            Some(PathBuf::from("/cfg/maille/materials.toml"))
        );
        assert_eq!(
            user_path(None, Some("/home/tim")),
            Some(PathBuf::from("/home/tim/.config/maille/materials.toml"))
        );
        assert_eq!(user_path(Some(""), None), None);
    }

    #[test]
    fn saving_creates_the_file_and_the_material_comes_back_on_loading() {
        let path = temp_file("save");
        let (builtin, note) = load(Some(&path));
        assert!(note.is_none() && builtin.get("maison").is_none());

        let mut mine = Library::builtin().get("acier").unwrap().clone();
        mine.young = Some(200e9);
        let (library, note) = save(&path, "maison", mine).unwrap();
        assert!(note.is_none());
        assert_eq!(library.get("maison").unwrap().young, Some(200e9));
        assert!(library.get("cuivre").is_some());
        // un second enregistrement garde le premier
        save(
            &path,
            "autre",
            Library::builtin().get("verre").unwrap().clone(),
        )
        .unwrap();
        let (library, _) = load(Some(&path));
        assert!(library.get("maison").is_some() && library.get("autre").is_some());
        let _ = fs::remove_dir_all(path.parent().unwrap().parent().unwrap());
    }

    #[test]
    fn a_broken_user_file_is_reported_and_never_overwritten() {
        let path = temp_file("broken");
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, "[materials.x]\nyoung = -1\n").unwrap();
        let (library, note) = load(Some(&path));
        assert!(library.get("acier").is_some());
        assert!(note.is_some_and(|n| n.contains("young")));

        let steel = Library::builtin().get("acier").unwrap().clone();
        assert!(save(&path, "y", steel).is_err());
        assert!(fs::read_to_string(&path).unwrap().contains("young = -1"));
        let _ = fs::remove_dir_all(path.parent().unwrap().parent().unwrap());
    }
}
