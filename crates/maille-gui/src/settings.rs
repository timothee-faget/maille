//! Réglages de l'interface qui survivent à la fermeture : le thème, dans
//! `~/.config/maille/settings.toml` (à côté de `materials.toml`).

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use crate::library;

/// Le thème voulu : un des deux, ou celui du système.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ThemeChoice {
    #[default]
    System,
    Light,
    Dark,
}

impl ThemeChoice {
    /// Le choix suivant, dans l'ordre Système, Clair, Sombre.
    pub fn next(self) -> Self {
        match self {
            Self::System => Self::Light,
            Self::Light => Self::Dark,
            Self::Dark => Self::System,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::System => "Système",
            Self::Light => "Clair",
            Self::Dark => "Sombre",
        }
    }

    fn key(self) -> &'static str {
        match self {
            Self::System => "system",
            Self::Light => "light",
            Self::Dark => "dark",
        }
    }

    /// Lit `theme = "…"` ; tout ce qui n'est pas reconnu donne `System`.
    pub fn parse(text: &str) -> Self {
        text.lines()
            .filter_map(|line| line.split_once('='))
            .find(|(name, _)| name.trim() == "theme")
            .map(|(_, value)| value.trim().trim_matches('"'))
            .map_or(Self::System, |value| match value {
                "light" => Self::Light,
                "dark" => Self::Dark,
                _ => Self::System,
            })
    }

    fn to_text(self) -> String {
        format!("theme = \"{}\"\n", self.key())
    }
}

/// Le fichier de réglages d'après l'environnement.
pub fn default_path() -> Option<PathBuf> {
    library::default_user_path().map(|p| p.with_file_name("settings.toml"))
}

/// Le thème enregistré ; `System` si le fichier manque ou ne se lit pas.
pub fn load(path: Option<&Path>) -> ThemeChoice {
    path.and_then(|p| fs::read_to_string(p).ok())
        .map_or(ThemeChoice::System, |text| ThemeChoice::parse(&text))
}

/// Enregistre le thème, en créant le dossier au besoin.
pub fn save(path: &Path, choice: ThemeChoice) -> io::Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir)?;
    }
    fs::write(path, choice.to_text())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_saved_theme_is_read_back() {
        let dir = std::env::temp_dir().join(format!("maille-settings-{}", std::process::id()));
        let path = dir.join("maille").join("settings.toml");
        for choice in [ThemeChoice::System, ThemeChoice::Light, ThemeChoice::Dark] {
            save(&path, choice).unwrap();
            assert_eq!(load(Some(&path)), choice);
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn a_missing_or_unreadable_file_means_the_system_theme() {
        assert_eq!(load(None), ThemeChoice::System);
        assert_eq!(
            load(Some(Path::new("/nulle/part/settings.toml"))),
            ThemeChoice::System
        );
        assert_eq!(
            ThemeChoice::parse("theme = \"violet\""),
            ThemeChoice::System
        );
        assert_eq!(ThemeChoice::parse("n'importe quoi"), ThemeChoice::System);
        assert_eq!(ThemeChoice::parse("theme = \"dark\"\n"), ThemeChoice::Dark);
    }

    #[test]
    fn the_choices_cycle() {
        let mut choice = ThemeChoice::System;
        for _ in 0..3 {
            choice = choice.next();
        }
        assert_eq!(choice, ThemeChoice::System);
    }
}
