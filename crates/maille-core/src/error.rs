use std::fmt;

use crate::mesh::Dof;
use crate::mesher::CurveId;

/// Erreur de modèle ou de calcul. Les messages disent quoi corriger.
#[derive(Debug, Clone, PartialEq)]
pub enum MailleError {
    EmptyMesh,
    InvalidProperty {
        name: &'static str,
        value: f64,
    },
    /// Valeur qui doit seulement être finie (pas NaN ni infinie), négative ou nulle permise.
    NonFinite {
        name: &'static str,
        value: f64,
    },
    /// Valeur qui doit être finie et positive ou nulle.
    NegativeProperty {
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
    /// Elements de familles différentes (Bar, Beam, triangles, tétraèdres) dans le même Mesh :
    /// les DOF ne sont pas compatibles (pas de portique, pas de mélange 2D/3D).
    MixedElementTypes,
    /// Beam dont le Node 2 est à gauche du Node 1 : inverser l'ordre des Nodes.
    ReversedBeam {
        element: usize,
    },
    /// DOF qui n'existe pas pour ce type d'Element ou cette analyse (Bar : U ; Beam : V, Theta ;
    /// triangles : Ux, Uy ; tétraèdres : Ux, Uy, Uz ; thermique : T).
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
    /// Un triangle ou un tétraèdre a besoin de ν en statique : construire le Material avec
    /// `Material::with_poisson`.
    MissingPoisson {
        element: usize,
    },
    /// Une analyse statique a besoin de E.
    MissingYoung {
        element: usize,
    },
    /// Une analyse thermique a besoin de k.
    MissingConductivity {
        element: usize,
    },
    /// Une analyse thermique transitoire a besoin de ρ.
    MissingDensity {
        element: usize,
    },
    /// Une analyse thermique transitoire a besoin de c.
    MissingSpecificHeat {
        element: usize,
    },
    /// Step sans aucun incrément.
    NoIncrements,
    /// Point `point` d'une Amplitude invalide (ou Amplitude sans point).
    InvalidAmplitude {
        point: usize,
    },
    /// θ hors de [½, 1] : en dessous de ½, le schéma n'est stable que pour de petits Δt.
    ThetaOutOfRange {
        theta: f64,
    },
    /// Champ de températures (condition initiale d'une Step, températures d'un calcul statique) qui
    /// n'a pas une valeur par Node.
    TemperatureCount {
        expected: usize,
        got: usize,
    },
    /// Température non finie (NaN ou infinie) au Node `node`.
    NonFiniteTemperature {
        node: usize,
    },
    /// Champ de températures (dilatation) ou calcul quasi-statique sur un Model qui n'est pas statique.
    StaticOnly,
    /// Plasticité, ou calcul quasi-statique, sur un Element qui n'est ni un triangle en déformations
    /// planes ni un tétraèdre.
    UnsupportedPlasticElement {
        element: usize,
    },
    /// `Model::solve` (élastique linéaire) sur un Material plastique.
    PlasticityNeedsQuasiStatic {
        element: usize,
    },
    /// Champ de températures dans un calcul quasi-statique (thermo-plasticité non gérée).
    QuasiStaticWithTemperatures,
    /// Grandes déformations sur un Element qui n'est ni un triangle en déformations planes ni un
    /// tétraèdre.
    UnsupportedLargeStrainElement {
        element: usize,
    },
    /// Plasticité en grandes déformations.
    LargeStrainPlasticity {
        element: usize,
    },
    /// Champ de températures en grandes déformations.
    LargeStrainWithTemperatures,
    /// Loi hyperélastique (néo-Hookéen) avec `Kinematics::Small`.
    HyperelasticityNeedsLargeStrain {
        element: usize,
    },
    /// `Model::solve` (linéaire) avec `Kinematics::Large`.
    LargeStrainNeedsQuasiStatic,
    /// Newton-Raphson n'a pas trouvé l'équilibre de l'incrément `increment` (temps `time`).
    NotConverged {
        increment: usize,
        time: f64,
    },
    /// Tangente singulière alors que des points de Gauss plastifient : la charge dépasse la charge
    /// limite de la structure (mécanisme plastique).
    PlasticCollapse {
        increment: usize,
        time: f64,
    },
    /// La dilatation d'un champ de températures demande α.
    MissingExpansion {
        element: usize,
    },
    /// La Beam n'a pas de dilatation thermique.
    UnsupportedThermalStrain {
        element: usize,
    },
    /// Face inexistante (triangle : 0 à 2, tétraèdre : 0 à 3 ; aucune pour Bar et Beam).
    UnknownFace {
        element: usize,
        face: usize,
    },
    /// Donnée thermique (source, convection) ou Step transitoire sur un Model qui n'est pas thermique.
    ThermalOnly,
    /// Pression demandée sur un Model thermique.
    PressureStaticOnly,
    /// Pression sur une Bar ou une Beam, qui n'ont pas de face.
    UnsupportedPressureElement {
        element: usize,
    },
    /// Bar et Beam n'ont pas de version thermique.
    UnsupportedThermalElement {
        element: usize,
    },
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
    /// Les gradients conjugués n'ont pas atteint la tolérance (résidu relatif) en `iterations`.
    IterativeNotConverged {
        iterations: usize,
        residual: f64,
    },
    /// Contraintes aux Nodes demandées sur un Mesh de Bar ou de Beam, ou sur un Field thermique.
    NoStresses,
    /// Le Field ne correspond pas au Mesh (nombre d'Elements ou de points de Gauss).
    FieldMeshMismatch,
    /// Déformation plastique cumulée p aux Nodes demandée sur un Field qui n'en a pas.
    NoEquivalentPlasticStrains,
    /// Le mailleur bute sur un cas dégénéré à l'arrondi près (points de la géométrie confondus ou
    /// presque cocirculaires).
    MeshingFailed,
    /// Le maillage dépasserait `max` Nodes (tailles h trop petites pour la géométrie).
    TooManyNodes {
        max: usize,
    },
    /// Surface 3D inutilisable par le mailleur (indice hors limites, coordonnée non finie…).
    InvalidSurface {
        problem: String,
    },
    /// Pièce B-rep ou courbe d'arête inutilisable (nœuds d'une B-spline mal rangés, poids négatif…).
    InvalidBrep {
        problem: String,
    },
    /// Taille h d'une Curve nulle, négative ou non finie.
    InvalidCurveSize {
        curve: CurveId,
        value: f64,
    },
    /// Coordonnée non finie (NaN ou infinie) du point `what` (« départ », « centre ») d'une Curve.
    NonFiniteCurve {
        curve: CurveId,
        what: &'static str,
    },
    /// Contour sans aucune Curve (`hole` vaut `None` pour le contour extérieur).
    EmptyBoundary {
        hole: Option<usize>,
    },
    /// Segment dont le point d'arrivée (le début de la Curve suivante) est confondu avec son départ.
    ZeroLengthCurve {
        curve: CurveId,
    },
    /// Arc dont le départ et l'arrivée ne sont pas à la même distance du centre, ou de rayon nul.
    ArcRadiusMismatch {
        curve: CurveId,
    },
    /// Deux Curves qui se coupent ou se touchent ailleurs qu'à leur raccord.
    SelfIntersectingBoundary {
        first: CurveId,
        second: CurveId,
    },
    /// Trou hors du contour extérieur ou à l'intérieur d'un autre trou.
    HoleOutsideBoundary {
        hole: usize,
    },
    /// Angle intérieur trop aigu (en degrés) au raccord de deux Curves : Ruppert ne termine pas.
    SharpCorner {
        first: CurveId,
        second: CurveId,
        angle: f64,
    },
    /// Angle minimal du mailleur hors de ]0, 33°].
    MinAngleOutOfRange {
        angle: f64,
    },
    /// Bord nommé inexistant ; `known` liste les noms portés par les Curves.
    UnknownBoundary {
        name: String,
        known: Vec<String>,
    },
}

impl fmt::Display for MailleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMesh => write!(f, "le Mesh ne contient aucun Node ou aucun Element"),
            Self::InvalidProperty { name, value } => {
                write!(f, "{name} = {value} : doit être fini et strictement positif")
            }
            Self::NonFinite { name, value } => write!(f, "{name} = {value} : doit être fini"),
            Self::NegativeProperty { name, value } => {
                write!(f, "{name} = {value} : doit être fini et positif ou nul")
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
                "le DOF {dof:?} du Node {node} n'existe pas pour ce Mesh (Bar : U ; Beam : V, Theta ; triangles : Ux, Uy ; tétraèdres : Ux, Uy, Uz ; thermique : T)"
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
            Self::MissingPoisson { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de coefficient de Poisson, requis en statique pour un triangle ou un tétraèdre : le créer avec Material::with_poisson"
            ),
            Self::MissingConductivity { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de conductivité thermique, requise en thermique : le créer avec Material::conductive ou ajouter with_conductivity"
            ),
            Self::MissingDensity { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de masse volumique, requise en thermique transitoire : ajouter with_density"
            ),
            Self::MissingSpecificHeat { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de chaleur massique, requise en thermique transitoire : ajouter with_specific_heat"
            ),
            Self::NoIncrements => write!(
                f,
                "la Step n'a aucun incrément : en demander au moins un"
            ),
            Self::InvalidAmplitude { point } => write!(
                f,
                "le point {point} de l'Amplitude est invalide : le premier point doit être (0, 0), les temps strictement croissants et toutes les valeurs finies"
            ),
            Self::ThetaOutOfRange { theta } => write!(
                f,
                "θ = {theta} : doit être entre 0,5 (Crank-Nicolson) et 1 (Euler implicite) ; en dessous de 0,5, le calcul diverge si Δt est trop grand"
            ),
            Self::TemperatureCount { expected, got } => write!(
                f,
                "le champ de températures a {got} valeurs pour {expected} Nodes : en donner une par Node, par exemple vec![T0; nombre de Nodes] ou les temperatures d'un Field thermique"
            ),
            Self::NonFiniteTemperature { node } => write!(
                f,
                "la température du Node {node} n'est pas finie (NaN ou infinie)"
            ),
            Self::StaticOnly => write!(
                f,
                "le champ de températures (dilatation thermique) et le calcul quasi-statique ne servent qu'en statique : créer le Model avec Model::new (pour la dilatation, le champ vient d'un calcul thermique sur le même Mesh)"
            ),
            Self::UnsupportedPlasticElement { element } => write!(
                f,
                "l'Element {element} n'a pas de version élasto-plastique : la plasticité et le calcul quasi-statique demandent des triangles en déformations planes (PlaneState::Strain) ou des tétraèdres"
            ),
            Self::PlasticityNeedsQuasiStatic { element } => write!(
                f,
                "le Material de l'Element {element} est plastique : Model::solve est élastique linéaire, utiliser Model::solve_quasi_static avec une Step::quasi_static et une Amplitude"
            ),
            Self::QuasiStaticWithTemperatures => write!(
                f,
                "le calcul quasi-statique ne prend pas de champ de températures (la thermo-plasticité n'est pas gérée) : retirer with_temperatures"
            ),
            Self::UnsupportedLargeStrainElement { element } => write!(
                f,
                "l'Element {element} ne supporte pas les grandes déformations, qui demandent des triangles en déformations planes (les contraintes planes ne sont pas gérées) ou des tétraèdres : utiliser Kinematics::Small, ou changer d'Element"
            ),
            Self::LargeStrainPlasticity { element } => write!(
                f,
                "le Material de l'Element {element} est plastique : la plasticité en grandes déformations n'est pas gérée, utiliser Kinematics::Small ou retirer with_plasticity"
            ),
            Self::LargeStrainWithTemperatures => write!(
                f,
                "les grandes déformations ne prennent pas de champ de températures : retirer with_temperatures ou utiliser Kinematics::Small"
            ),
            Self::HyperelasticityNeedsLargeStrain { element } => write!(
                f,
                "le Material de l'Element {element} est hyperélastique (néo-Hookéen) : la loi n'existe qu'en grandes déformations, utiliser Model::with_kinematics(Kinematics::Large)"
            ),
            Self::LargeStrainNeedsQuasiStatic => write!(
                f,
                "Kinematics::Large demande un calcul non linéaire : utiliser Model::solve_quasi_static avec une Step::quasi_static et une Amplitude, au lieu de Model::solve"
            ),
            Self::NotConverged { increment, time } => write!(
                f,
                "Newton-Raphson ne converge pas à l'incrément {increment} (t = {time}) : augmenter le nombre d'incréments de la Step, ou réduire la charge. Si la plasticité est parfaite (H = 0), la charge dépasse peut-être la charge limite de la structure ; en grandes déformations, la structure flambe peut-être (ou ses Elements se retournent : mailler plus fin)"
            ),
            Self::PlasticCollapse { increment, time } => write!(
                f,
                "ruine plastique à l'incrément {increment} (t = {time}) : la charge dépasse la charge limite de la structure, qui n'a plus de rigidité. Réduire l'Amplitude, imposer un déplacement plutôt qu'une force, ou ajouter de l'écrouissage (H > 0)"
            ),
            Self::MissingExpansion { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de coefficient de dilatation, requis pour la dilatation d'un champ de températures : ajouter with_expansion"
            ),
            Self::UnsupportedThermalStrain { element } => write!(
                f,
                "l'Element {element} est une Beam, qui n'a pas de dilatation thermique (il faudrait un gradient dans l'épaisseur) : utiliser des Bar, des triangles ou des tétraèdres"
            ),
            Self::UnknownFace { element, face } => write!(
                f,
                "l'Element {element} n'a pas de face {face} : un triangle a les faces 0 à 2 (arêtes 0-1, 1-2, 2-0), un tétraèdre les faces 0 à 3 (0-1-2, 0-3-1, 1-3-2, 2-3-0)"
            ),
            Self::PressureStaticOnly => write!(
                f,
                "la pression n'existe qu'en statique : en thermique, utiliser la convection ou une source de chaleur"
            ),
            Self::UnsupportedPressureElement { element } => write!(
                f,
                "l'Element {element} est une Bar ou une Beam, qui n'a pas de face : la pression s'applique aux triangles et aux tétraèdres (pour une Bar ou une Beam, utiliser une charge répartie)"
            ),
            Self::ThermalOnly => write!(
                f,
                "les sources de chaleur, la convection et le calcul transitoire n'existent qu'en thermique : créer le Model avec Model::with_analysis(mesh, Analysis::Thermal, ...)"
            ),
            Self::UnsupportedThermalElement { element } => write!(
                f,
                "l'Element {element} est une Bar ou une Beam, qui n'existent pas en thermique : mailler la pièce en triangles ou en tétraèdres"
            ),
            Self::MissingYoung { element } => write!(
                f,
                "le Material de l'Element {element} n'a pas de module d'Young, requis en statique : le créer avec Material::new ou Material::with_poisson"
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
                "le tétraèdre {element} est retourné : vu depuis son Node 3, la face 0-1-2 doit tourner dans le sens antihoraire. Échanger ses Nodes 1 et 2 (pour un Tet10, échanger aussi les Nodes 4 et 6, puis 8 et 9)"
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
            Self::IterativeNotConverged {
                iterations,
                residual,
            } => write!(
                f,
                "le solveur itératif n'a pas convergé en {iterations} itérations (résidu relatif {residual:.1e}) : \
                 modèle mal conditionné ou presque singulier (Elements très aplatis, structure peu contrainte). \
                 Vérifier les BC, ou utiliser SolverKind::Direct"
            ),
            Self::NoStresses => write!(
                f,
                "les contraintes aux Nodes n'existent qu'en statique, pour les triangles et les tétraèdres : pas pour les Bar et Beam (lire les efforts dans les réactions) ni pour un Field thermique (lire temperatures et heat_fluxes)"
            ),
            Self::FieldMeshMismatch => write!(
                f,
                "le Field ne vient pas de ce Mesh : passer le Field renvoyé par la résolution du Model construit sur ce Mesh"
            ),
            Self::NoEquivalentPlasticStrains => write!(
                f,
                "le Field n'a pas de déformation plastique cumulée p : elle n'existe que dans les Increments de Model::solve_quasi_static"
            ),
            Self::MeshingFailed => write!(
                f,
                "le mailleur n'a pas pu trianguler la géométrie (points de la géométrie confondus ou presque cocirculaires à l'arrondi près) : déplacer légèrement un point de la géométrie ou changer les tailles h"
            ),
            Self::InvalidSurface { problem } => write!(f, "surface 3D inutilisable : {problem}"),
            Self::InvalidBrep { problem } => write!(f, "pièce B-rep inutilisable : {problem}"),
            Self::TooManyNodes { max } => write!(
                f,
                "le maillage dépasserait {max} Nodes : augmenter les tailles h des Curves"
            ),
            Self::InvalidCurveSize { curve, value } => write!(
                f,
                "{curve} a une taille h = {value} : elle doit être finie et strictement positive"
            ),
            Self::NonFiniteCurve { curve, what } => write!(
                f,
                "{curve} a un point de {what} dont une coordonnée n'est pas finie (NaN ou infinie)"
            ),
            Self::EmptyBoundary { hole: None } => {
                write!(f, "le contour extérieur ne contient aucune Curve")
            }
            Self::EmptyBoundary { hole: Some(h) } => {
                write!(f, "le trou {h} ne contient aucune Curve : le retirer ou le décrire")
            }
            Self::ZeroLengthCurve { curve } => write!(
                f,
                "{curve} a une longueur nulle : elle finit au début de la Curve suivante, qui est confondu avec son propre début. La retirer"
            ),
            Self::ArcRadiusMismatch { curve } => write!(
                f,
                "{curve} est un arc dont le début et la fin (le début de la Curve suivante) ne sont pas à la même distance du centre, ou de rayon nul : corriger le centre ou les points"
            ),
            Self::SelfIntersectingBoundary { first, second } => write!(
                f,
                "{first} et {second} se coupent ou se touchent : un contour ne doit pas se recouper, et les trous ne doivent toucher ni le contour extérieur ni les autres trous"
            ),
            Self::HoleOutsideBoundary { hole } => write!(
                f,
                "le trou {hole} est hors du contour extérieur ou dans un autre trou : le déplacer dans le domaine"
            ),
            Self::SharpCorner {
                first,
                second,
                angle,
            } => write!(
                f,
                "angle de {angle:.1}° entre {first} et {second} : le mailleur demande au moins {}° entre deux Curves qui se suivent. Arrondir ou couper l'angle",
                crate::mesher::MIN_CORNER
            ),
            Self::MinAngleOutOfRange { angle } => write!(
                f,
                "angle minimal du mailleur = {angle}° : doit être strictement positif et au plus {}° (au-delà, le raffinement peut ne pas terminer)",
                crate::mesher::MAX_MIN_ANGLE
            ),
            Self::UnknownBoundary { name, known } if known.is_empty() => write!(
                f,
                "aucune Curve ne porte le nom « {name} » : la géométrie n'a aucun bord nommé"
            ),
            Self::UnknownBoundary { name, known } => write!(
                f,
                "aucune Curve ne porte le nom « {name} » : les bords nommés sont {}",
                known
                    .iter()
                    .map(|n| format!("« {n} »"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        }
    }
}

impl std::error::Error for MailleError {}
