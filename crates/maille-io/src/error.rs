use std::fmt;

use crate::msh::CellKind;

/// Erreur de lecture ou d'écriture d'un fichier. Les messages disent quoi corriger.
#[derive(Debug, Clone, PartialEq)]
pub enum IoError {
    /// Fichier mal formé à la ligne `line` (numérotée depuis 1).
    Syntax {
        line: usize,
        message: String,
    },
    /// Version autre que 4.1 (`found` est la version lue).
    UnsupportedVersion {
        line: usize,
        found: String,
    },
    /// Fichier binaire.
    BinaryFormat {
        line: usize,
    },
    /// Section obligatoire absente.
    MissingSection {
        name: &'static str,
    },
    /// Type d'élément Gmsh que Maille ne sait pas lire.
    UnsupportedElement {
        line: usize,
        gmsh_type: u32,
    },
    /// Aucun triangle ni tétraèdre.
    NoCells,
    /// Deux types d'Elements différents (Tri3 avec Tri6, Tet4 avec Tet10…).
    MixedCellKinds {
        first: CellKind,
        found: CellKind,
    },
    /// Un maillage 2D doit être dans le plan z = 0.
    NonPlanar {
        z: f64,
    },
    /// Deux groupes physiques portent le même nom.
    DuplicateSetName {
        name: String,
    },
    DuplicateNodeTag {
        line: usize,
        tag: i64,
    },
    UnknownNode {
        line: usize,
        tag: i64,
    },
    /// Un élément de bord du groupe `set` ne correspond à aucune face d'Element.
    BoundaryFaceNotFound {
        line: usize,
        set: String,
    },
    /// `MshFile` incohérent à l'écriture (indice hors limites, Set mélangé…).
    Invalid {
        problem: String,
    },
}

fn gmsh_name(t: u32) -> String {
    match t {
        3 => "quadrangle".to_string(),
        5 => "hexaèdre".to_string(),
        6 => "prisme".to_string(),
        7 => "pyramide".to_string(),
        10 => "quadrangle à 9 Nodes".to_string(),
        t => format!("de type {t}"),
    }
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Syntax { line, message } => write!(f, "ligne {line} : {message}"),
            Self::UnsupportedVersion { line, found } => write!(
                f,
                "ligne {line} : format {found} non pris en charge, Maille lit le .msh 4.1. Ré-exporter depuis Gmsh avec `Mesh.MshFileVersion = 4.1`"
            ),
            Self::BinaryFormat { line } => write!(
                f,
                "ligne {line} : fichier binaire, Maille lit le .msh en ASCII. Ré-exporter depuis Gmsh avec `Mesh.Binary = 0`"
            ),
            Self::MissingSection { name } => write!(f, "section ${name} absente du fichier"),
            Self::UnsupportedElement { line, gmsh_type } => write!(
                f,
                "ligne {line} : élément {} non pris en charge. Maille lit les triangles (Tri3, Tri6) et les tétraèdres (Tet4, Tet10) : mailler en triangles ou en tétraèdres",
                gmsh_name(*gmsh_type)
            ),
            Self::NoCells => write!(
                f,
                "le fichier ne contient ni triangle ni tétraèdre dans un groupe physique ou une entité : pensez à `Mesh.SaveAll = 1` ou à un groupe physique de surface (2D) ou de volume (3D)"
            ),
            Self::MixedCellKinds { first, found } => write!(
                f,
                "le fichier mélange des Elements {first:?} et {found:?} : utiliser un seul type (ordre 1 ou ordre 2)"
            ),
            Self::NonPlanar { z } => write!(
                f,
                "maillage de triangles hors du plan z = 0 (z = {z}) : Maille résout le 2D dans le plan xy"
            ),
            Self::DuplicateSetName { name } => write!(
                f,
                "deux groupes physiques portent le nom « {name} » : leur donner des noms différents"
            ),
            Self::DuplicateNodeTag { line, tag } => {
                write!(f, "ligne {line} : le Node {tag} est défini deux fois")
            }
            Self::UnknownNode { line, tag } => {
                write!(f, "ligne {line} : le Node {tag} n'est pas défini dans $Nodes")
            }
            Self::BoundaryFaceNotFound { line, set } => write!(
                f,
                "ligne {line} : un élément du groupe « {set} » ne correspond à aucune face d'Element : le bord doit être maillé avec le domaine"
            ),
            Self::Invalid { problem } => write!(f, "MshFile invalide : {problem}"),
        }
    }
}

impl std::error::Error for IoError {}
