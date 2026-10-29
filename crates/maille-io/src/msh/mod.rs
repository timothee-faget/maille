//! Maillage Gmsh `.msh` 4.1 ASCII.
//!
//! `read` produit un `MshFile` : Nodes, connectivité (numérotation de Maille, orientation
//! corrigée) et Sets nommés (les groupes physiques de Gmsh). Ce n'est pas un `Mesh` : un `.msh` n'a
//! que la géométrie, et les Elements de Maille portent leur Material. Le chargeur de modèle crée
//! les Elements à partir des Sets, puis appelle `Mesh::new`. `write` fait le chemin inverse.

mod numbering;
mod reader;
mod writer;

use std::collections::BTreeMap;

use maille_core::mesh::Node;
use maille_core::mesher::{BoundarySet, Tetrahedralization, Triangulation};

use crate::IoError;

pub use numbering::face_positions;
pub use reader::read;
pub use writer::write;

/// Nom du Set d'Elements que reçoit le domaine d'une `Triangulation`.
pub const DOMAIN_SET: &str = "domaine";

/// Type des Elements d'un maillage (un seul par fichier).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellKind {
    Tri3,
    Tri6,
    Tet4,
    Tet10,
}

impl CellKind {
    /// Dimension de l'espace : 2 pour les triangles, 3 pour les tétraèdres.
    pub fn dim(self) -> usize {
        match self {
            Self::Tri3 | Self::Tri6 => 2,
            Self::Tet4 | Self::Tet10 => 3,
        }
    }

    pub fn node_count(self) -> usize {
        match self {
            Self::Tri3 => 3,
            Self::Tri6 => 6,
            Self::Tet4 => 4,
            Self::Tet10 => 10,
        }
    }

    fn gmsh_type(self) -> u32 {
        match self {
            Self::Tri3 => 2,
            Self::Tri6 => 9,
            Self::Tet4 => 4,
            Self::Tet10 => 11,
        }
    }

    /// Type Gmsh de l'élément de bord (ligne en 2D, triangle en 3D).
    fn face_gmsh_type(self) -> u32 {
        match self {
            Self::Tri3 => 1,
            Self::Tri6 => 8,
            Self::Tet4 => 2,
            Self::Tet10 => 9,
        }
    }

    fn from_gmsh_type(t: u32) -> Option<Self> {
        match t {
            2 => Some(Self::Tri3),
            9 => Some(Self::Tri6),
            4 => Some(Self::Tet4),
            11 => Some(Self::Tet10),
            _ => None,
        }
    }
}

/// Ensemble nommé : un groupe physique de Gmsh. Les indices sont ceux de `MshFile` (depuis 0).
///
/// - groupe de la dimension du maillage : `elements` ;
/// - groupe d'une dimension de moins (bord) : `faces` et les `nodes` de ces faces ;
/// - groupe de dimension inférieure (points, et lignes en 3D) : `nodes` seulement.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Set {
    /// Nodes triés, sans doublon.
    pub nodes: Vec<usize>,
    /// Faces (élément, face) triées, numérotées comme `Element::face_nodes`.
    pub faces: Vec<(usize, usize)>,
    /// Elements triés, sans doublon.
    pub elements: Vec<usize>,
}

/// Contenu d'un `.msh` : géométrie et Sets, sans Material.
#[derive(Debug, Clone, PartialEq)]
pub struct MshFile {
    pub nodes: Vec<Node>,
    pub kind: CellKind,
    /// Connectivité, `kind.node_count()` Nodes par Element, dans l'ordre de `Tri3`, `Tri6`, `Tet4`
    /// ou `Tet10`. Triangles antihoraires et tétraèdres positifs.
    pub cells: Vec<Vec<usize>>,
    pub sets: BTreeMap<String, Set>,
}

impl MshFile {
    /// Fichier d'une `Triangulation` : un Set par bord nommé, plus `DOMAIN_SET` pour tous les
    /// triangles. `DuplicateSetName` si un bord s'appelle déjà `domaine`.
    pub fn from_triangulation<const N: usize>(t: &Triangulation<N>) -> Result<Self, IoError> {
        Ok(Self {
            nodes: t.nodes.clone(),
            kind: if N == 3 {
                CellKind::Tri3
            } else {
                CellKind::Tri6
            },
            cells: t.triangles.iter().map(|c| c.to_vec()).collect(),
            sets: mesher_sets(
                t.boundary_names(),
                |n| t.boundary(n).ok(),
                t.triangles.len(),
            )?,
        })
    }

    /// Fichier d'une `Tetrahedralization` : un Set par bord nommé (un patch de la surface), plus
    /// `DOMAIN_SET` pour tous les tétraèdres. `DuplicateSetName` si un bord s'appelle déjà
    /// `domaine`.
    pub fn from_tetrahedralization<const N: usize>(
        t: &Tetrahedralization<N>,
    ) -> Result<Self, IoError> {
        Ok(Self {
            nodes: t.nodes.clone(),
            kind: if N == 4 {
                CellKind::Tet4
            } else {
                CellKind::Tet10
            },
            cells: t.tetrahedra.iter().map(|c| c.to_vec()).collect(),
            sets: mesher_sets(
                t.boundary_names(),
                |n| t.boundary(n).ok(),
                t.tetrahedra.len(),
            )?,
        })
    }
}

/// Sets d'un mailleur : un par bord nommé, plus `DOMAIN_SET` pour les `count` Elements.
fn mesher_sets<'a>(
    names: impl Iterator<Item = &'a str>,
    boundary: impl Fn(&str) -> Option<&'a BoundarySet>,
    count: usize,
) -> Result<BTreeMap<String, Set>, IoError> {
    let mut sets = BTreeMap::new();
    for name in names {
        if name == DOMAIN_SET {
            return Err(IoError::DuplicateSetName {
                name: name.to_string(),
            });
        }
        // le nom vient de boundary_names : `boundary` ne peut pas échouer
        if let Some(b) = boundary(name) {
            sets.insert(
                name.to_string(),
                Set {
                    nodes: b.nodes.clone(),
                    faces: b.faces.clone(),
                    elements: Vec::new(),
                },
            );
        }
    }
    sets.insert(
        DOMAIN_SET.to_string(),
        Set {
            elements: (0..count).collect(),
            ..Set::default()
        },
    );
    Ok(sets)
}
