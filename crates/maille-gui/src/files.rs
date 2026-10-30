//! Fichiers écrits par l'interface (maillage, modèle) sans écraser le travail de l'utilisateur : un
//! fichier qui existe déjà n'est remplacé sans question que si son contenu ne change pas, ou si
//! l'interface l'a écrit elle-même pendant cette session et que personne ne l'a modifié depuis.

use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::fs;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;

/// Un fichier à écrire et son contenu.
pub type File = (PathBuf, Vec<u8>);

/// Empreinte de ce que l'interface a écrit à chaque chemin pendant la session.
#[derive(Debug, Default)]
pub struct Written(HashMap<PathBuf, u64>);

impl Written {
    /// Les fichiers de `files` qu'écrire ferait perdre : ils existent, leur contenu change, et ce
    /// n'est pas l'interface qui les a écrits tels quels.
    pub fn conflicts(&self, files: &[File]) -> Vec<PathBuf> {
        files
            .iter()
            .filter(|(path, bytes)| match fs::read(path) {
                // absent ou illisible : rien à perdre, l'écriture dira l'erreur s'il y en a une
                Err(_) => false,
                Ok(existing) => existing != *bytes && self.0.get(path) != Some(&digest(&existing)),
            })
            .map(|(path, _)| path.clone())
            .collect()
    }

    /// Écrit les fichiers, dans l'ordre, et retient leur empreinte.
    pub fn write(&mut self, files: &[File]) -> Result<(), String> {
        for (path, bytes) in files {
            fs::write(path, bytes)
                .map_err(|e| format!("impossible d'écrire {} : {e}", path.display()))?;
            self.0.insert(path.clone(), digest(bytes));
        }
        Ok(())
    }
}

fn digest(bytes: &[u8]) -> u64 {
    let mut h = DefaultHasher::new();
    bytes.hash(&mut h);
    h.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("maille-gui-files-{name}"));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn a_missing_file_or_the_same_content_is_written_without_asking() {
        let d = dir("absent");
        let path = d.join("plaque_modele.toml");
        let written = Written::default();
        let file = (path.clone(), b"a = 1\n".to_vec());
        assert!(written.conflicts(std::slice::from_ref(&file)).is_empty());
        fs::write(&path, b"a = 1\n").unwrap();
        assert!(written.conflicts(&[file]).is_empty());
    }

    #[test]
    fn a_file_the_user_wrote_is_a_conflict() {
        let d = dir("utilisateur");
        let path = d.join("plaque.msh");
        fs::write(&path, b"maillage de Gmsh").unwrap();
        let written = Written::default();
        assert_eq!(
            written.conflicts(&[(path.clone(), b"maillage de Maille".to_vec())]),
            [path]
        );
    }

    #[test]
    fn a_file_the_interface_wrote_is_rewritten_until_someone_edits_it() {
        let d = dir("interface");
        let (msh, toml) = (d.join("plaque.msh"), d.join("plaque_modele.toml"));
        let mut written = Written::default();
        written
            .write(&[
                (msh.clone(), b"v1".to_vec()),
                (toml.clone(), b"v1".to_vec()),
            ])
            .unwrap();
        let next = [
            (msh.clone(), b"v2".to_vec()),
            (toml.clone(), b"v2".to_vec()),
        ];
        assert!(written.conflicts(&next).is_empty());
        // retouché à la main : seul ce fichier est en conflit
        fs::write(&toml, b"v1 retouche").unwrap();
        assert_eq!(written.conflicts(&next), std::slice::from_ref(&toml));
        // une fois écrasé avec l'accord de l'utilisateur, il redevient celui de l'interface
        written.write(&next).unwrap();
        assert_eq!(fs::read(&toml).unwrap(), b"v2");
        assert!(written
            .conflicts(&[(toml, b"v3".to_vec()), (msh, b"v3".to_vec())])
            .is_empty());
    }

    #[test]
    fn a_write_error_names_the_file() {
        let d = dir("erreur");
        let path = d.join("absent").join("plaque.msh");
        let e = Written::default()
            .write(&[(path.clone(), b"x".to_vec())])
            .unwrap_err();
        assert!(e.contains(&path.display().to_string()), "{e}");
    }
}
