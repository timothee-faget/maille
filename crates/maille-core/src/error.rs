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
    /// Tri6 dont un Node milieu n'est pas au milieu de son arête (côtés courbes non supportés).
    CurvedTriangle {
        element: usize,
    },
    /// Un Element 2D a besoin de ν : construire le Material avec `Material::with_poisson`.
    MissingPoisson,
    UnsupportedDistributedLoad {
        element: usize,
    },
    /// Tétraèdre de volume nul : Nodes coplanaires ou confondus.
    DegenerateTetrahedron {
        element: usize,
    },
    /// Tétraèdre retourné (volume négatif) : échanger deux de ses Nodes.
    InvertedTetrahedron {
        element: usize,
    },
    /// Tet10 dont un Node milieu n'est pas au milieu de son arête (arêtes courbes non supportées).
    CurvedTetrahedron {
        element: usize,
    },
    /// Matrice non définie positive : structure non contrainte ou mécanisme.
    Singular,
    /// Contraintes aux Nodes demandées sur un Mesh de Bar ou de Beam.
    NoStresses,
    /// Traitement réservé aux triangles (le viewer est 2D).
    PlaneMeshRequired,
    /// Le Field ne correspond pas au Mesh (nombre d'Elements ou de points de Gauss).
    FieldMeshMismatch,
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
                "le Mesh mélange des Elements de familles différentes (Bar, Beam, triangles, tétraèdres) : n'en garder qu'une"
            ),
            Self::ReversedBeam { element } => write!(
                f,
                "la Beam {element} va vers les x décroissants : inverser l'ordre de ses deux Nodes"
            ),
            Self::InvalidDof { node, dof } => write!(
                f,
                "le DOF {dof:?} du Node {node} n'existe pas pour ce Mesh (Bar : U ; Beam : V, Theta ; triangles : Ux, Uy ; tétraèdres : Ux, Uy, Uz)"
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
            Self::CurvedTriangle { element } => write!(
                f,
                "le Tri6 {element} a un Node milieu qui n'est pas au milieu de son arête : seuls les côtés droits sont supportés"
            ),
            Self::MissingPoisson => write!(
                f,
                "le Material d'un triangle ou d'un tétraèdre n'a pas de coefficient de Poisson : le créer avec Material::with_poisson"
            ),
            Self::UnsupportedDistributedLoad { element } => write!(
                f,
                "l'Element {element} ne supporte pas de charge répartie (Bar et Beam seulement) : utiliser des Load nodales"
            ),
            Self::DegenerateTetrahedron { element } => write!(
                f,
                "le tétraèdre {element} a un volume nul : ses Nodes sont coplanaires ou confondus"
            ),
            Self::InvertedTetrahedron { element } => write!(
                f,
                "le tétraèdre {element} est retourné : vu depuis son Node 3, la face 0-1-2 doit tourner dans le sens antihoraire. Échanger les sommets 1 et 2 (pour un Tet10, échanger aussi les milieux 4 et 6, puis 8 et 9)"
            ),
            Self::CurvedTetrahedron { element } => write!(
                f,
                "le Tet10 {element} a un Node milieu qui n'est pas au milieu de son arête : seules les arêtes droites sont supportées"
            ),
            Self::Singular => write!(
                f,
                "système singulier : la structure n'est pas assez contrainte (mouvement de corps rigide) \
                 ou contient un mécanisme. Ajouter une BC ou vérifier la connectivité des Elements"
            ),
            Self::NoStresses => write!(
                f,
                "les contraintes aux Nodes n'existent que pour les triangles et les tétraèdres, pas pour les Bar et Beam"
            ),
            Self::PlaneMeshRequired => write!(
                f,
                "seuls les Mesh de triangles (Tri3, Tri6) sont pris en charge ici : le viewer est 2D pour l'instant"
            ),
            Self::FieldMeshMismatch => write!(
                f,
                "le Field ne vient pas de ce Mesh : passer le Field renvoyé par la résolution du Model construit sur ce Mesh"
            ),
        }
    }
}

impl std::error::Error for MailleError {}
