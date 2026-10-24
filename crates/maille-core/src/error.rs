use std::fmt;

use crate::mesh::Dof;

/// Erreur de modèle ou de calcul. Les messages disent quoi corriger.
#[derive(Debug, Clone, PartialEq)]
pub enum MailleError {
    EmptyMesh,
    InvalidProperty {
        name: &'static str,
        value: f64,
    },
    UnknownNode {
        node: usize,
    },
    ZeroLengthElement {
        element: usize,
    },
    DuplicateBc {
        node: usize,
        dof: Dof,
    },
    /// Bar et Beam dans le même Mesh : les DOF ne sont pas compatibles (pas de portique en M0).
    MixedElementTypes,
    /// Beam dont le Node 2 est à gauche du Node 1 : inverser l'ordre des Nodes.
    ReversedBeam {
        element: usize,
    },
    /// DOF qui n'existe pas pour ce type d'Element (Bar : U ; Beam : V, Theta).
    InvalidDof {
        node: usize,
        dof: Dof,
    },
    UnknownElement {
        element: usize,
    },
    /// Triangle d'aire nulle : Nodes alignés ou confondus.
    DegenerateTriangle {
        element: usize,
    },
    /// Triangle dont les Nodes sont dans le sens horaire : inverser l'ordre de deux Nodes.
    ClockwiseTriangle {
        element: usize,
    },
    /// Un Element 2D a besoin de ν : construire le Material avec `Material::with_poisson`.
    MissingPoisson,
    UnsupportedDistributedLoad {
        element: usize,
    },
    /// Matrice non définie positive : structure non contrainte ou mécanisme.
    Singular,
}

impl fmt::Display for MailleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMesh => write!(f, "le Mesh ne contient aucun Node ou aucun Element"),
            Self::InvalidProperty { name, value } => {
                write!(f, "{name} = {value} : doit être fini et strictement positif")
            }
            Self::UnknownNode { node } => {
                write!(f, "le Node {node} n'existe pas dans le Mesh")
            }
            Self::ZeroLengthElement { element } => write!(
                f,
                "l'Element {element} a une longueur nulle : séparer ses deux Nodes"
            ),
            Self::DuplicateBc { node, dof } => write!(
                f,
                "deux BC sur le DOF {dof:?} du Node {node} : n'en garder qu'une"
            ),
            Self::MixedElementTypes => write!(
                f,
                "le Mesh mélange des Elements de familles différentes (Bar, Beam, triangles) : n'en garder qu'une"
            ),
            Self::ReversedBeam { element } => write!(
                f,
                "la Beam {element} va vers les x décroissants : inverser l'ordre de ses deux Nodes"
            ),
            Self::InvalidDof { node, dof } => write!(
                f,
                "le DOF {dof:?} du Node {node} n'existe pas pour ce Mesh (Bar : U ; Beam : V, Theta ; triangles : Ux, Uy)"
            ),
            Self::UnknownElement { element } => {
                write!(f, "l'Element {element} n'existe pas dans le Mesh")
            }
            Self::DegenerateTriangle { element } => write!(
                f,
                "le triangle {element} a une aire nulle : ses Nodes sont alignés ou confondus"
            ),
            Self::ClockwiseTriangle { element } => write!(
                f,
                "le triangle {element} est dans le sens horaire : échanger deux de ses Nodes pour l'orienter dans le sens antihoraire"
            ),
            Self::MissingPoisson => write!(
                f,
                "le Material d'un triangle n'a pas de coefficient de Poisson : le créer avec Material::with_poisson"
            ),
            Self::UnsupportedDistributedLoad { element } => write!(
                f,
                "l'Element {element} ne supporte pas de charge répartie (Bar et Beam seulement) : utiliser des Load nodales"
            ),
            Self::Singular => write!(
                f,
                "système singulier : la structure n'est pas assez contrainte (mouvement de corps rigide) \
                 ou contient un mécanisme. Ajouter une BC ou vérifier la connectivité des Elements"
            ),
        }
    }
}

impl std::error::Error for MailleError {}
