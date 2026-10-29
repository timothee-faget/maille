//! Ce qu'on dessine, calculé sans GPU : triangles et arêtes à tracer, valeur de chaque champ aux
//! Nodes, positions déformées. Un Mesh 2D est dans le plan z = 0 ; un Mesh de tétraèdres est
//! dessiné par sa peau (voir `solid`), et ses cellules servent au plan de coupe.

use std::collections::HashSet;
use std::fmt;

use maille_core::mesh::{Element, Mesh};
use maille_core::results::{nodal_equivalent_plastic_strains, nodal_stresses, Field, Stress};
use maille_core::MailleError;

use crate::solid::skin;

/// Part de la taille du Mesh occupée par le plus grand déplacement, dans la déformée amplifiée.
const DEFORMATION_RATIO: f64 = 0.1;

/// Champ affiché, une valeur par Node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantity {
    VonMises,
    SigmaXx,
    SigmaYy,
    SigmaZz,
    TauXy,
    TauYz,
    TauXz,
    /// Norme du déplacement |u|.
    Displacement,
    Temperature,
    /// Déformation plastique cumulée p.
    PlasticStrain,
    /// Défaut de forme (0 : parfait), pour regarder un Mesh sans résultats de calcul.
    Distortion,
}

impl Quantity {
    /// Champs d'un Field 2D élasto-plastique, dans l'ordre des touches 1 à 6 : les champs
    /// mécaniques, puis p.
    pub const PLASTIC: [Self; 6] = [
        Self::VonMises,
        Self::SigmaXx,
        Self::SigmaYy,
        Self::TauXy,
        Self::Displacement,
        Self::PlasticStrain,
    ];

    /// Champs d'un Field 2D statique (touches 1 à 5) : ceux d'un Field élasto-plastique sans p.
    pub const MECHANICAL: &'static [Self] = Self::PLASTIC.split_at(5).0;

    /// Champs d'un Field statique de tétraèdres (touches 1 à 8), puis p (touche 9) si plastique.
    pub const SOLID_PLASTIC: [Self; 9] = [
        Self::VonMises,
        Self::SigmaXx,
        Self::SigmaYy,
        Self::SigmaZz,
        Self::TauXy,
        Self::TauYz,
        Self::TauXz,
        Self::Displacement,
        Self::PlasticStrain,
    ];

    pub const SOLID: &'static [Self] = Self::SOLID_PLASTIC.split_at(8).0;

    pub fn label(self) -> &'static str {
        match self {
            Self::VonMises => "von Mises",
            Self::SigmaXx => "σxx",
            Self::SigmaYy => "σyy",
            Self::SigmaZz => "σzz",
            Self::TauXy => "τxy",
            Self::TauYz => "τyz",
            Self::TauXz => "τxz",
            Self::Displacement => "|u|",
            Self::Temperature => "T",
            Self::PlasticStrain => "p (déformation plastique cumulée)",
            Self::Distortion => "Défaut de forme (0 : parfait)",
        }
    }
}

pub struct Scene {
    nodes: Vec<[f64; 3]>,
    /// (ux, uy, uz) par Node.
    displacements: Vec<[f64; 3]>,
    /// Vide pour un Field thermique.
    stresses: Vec<Stress>,
    /// Vide pour un Field statique.
    temperatures: Vec<f64>,
    /// Déformation plastique cumulée p aux Nodes, vide hors calcul élasto-plastique.
    plastic_strains: Vec<f64>,
    /// Défaut de forme aux Nodes d'un Mesh montré sans résultats, vide sinon.
    distortion: Vec<f64>,
    /// Grandes déformations (`NodalResults::large_strain`).
    large_strain: bool,
    /// Triangles à remplir (indices de Nodes) : 1 par Tri3, 4 par Tri6 (antihoraires) ; pour des
    /// tétraèdres, la peau, normales sortantes.
    pub triangles: Vec<[u32; 3]>,
    /// Segments du contour des Elements (ou de la peau), sans doublon : les arêtes internes des
    /// Tri6 et des faces de Tet10 découpés n'y sont pas, on voit le vrai Mesh.
    pub edges: Vec<[u32; 2]>,
    /// Cellules des tétraèdres (4 ou 10 Nodes), pour le plan de coupe ; vide en 2D.
    pub solids: Vec<Vec<usize>>,
}

/// Résultat lissé aux Nodes, tel que le relit un `.vtu` : une valeur par Node, pas de points de
/// Gauss. Un champ statique remplit `displacements` et `stresses`, un champ thermique
/// `temperatures` ; `plastic_strains` reste vide hors calcul élasto-plastique.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct NodalResults {
    /// Positions (x, y, z) ; z = 0 pour un Mesh 2D.
    pub nodes: Vec<[f64; 3]>,
    /// Cellules : triangles de 3 Nodes (Tri3) ou 6 (Tri6, sommets puis milieux 0-1, 1-2, 2-0),
    /// tétraèdres de 4 (Tet4) ou 10 (Tet10).
    pub cells: Vec<Vec<usize>>,
    /// (ux, uy, uz) ; uz = 0 pour un Mesh 2D.
    pub displacements: Vec<[f64; 3]>,
    pub stresses: Vec<Stress>,
    pub temperatures: Vec<f64>,
    pub plastic_strains: Vec<f64>,
    /// Un Mesh sans calcul : son défaut de forme (0 : parfait) à chaque Node ; les autres champs restent vides.
    pub distortion: Vec<f64>,
    /// Grandes déformations : la déformée se montre à l'échelle 1 (aucune amplification).
    pub large_strain: bool,
}

/// Ajoute les triangles à remplir et le contour d'un triangle de 3 (Tri3) ou 6 Nodes (Tri6).
/// Faux, sans rien ajouter, pour un autre nombre de Nodes.
fn add_cell(triangles: &mut Vec<[u32; 3]>, edges: &mut HashSet<[u32; 2]>, n: &[u32]) -> bool {
    let (tris, outline): (&[[usize; 3]], &[[usize; 2]]) = match n.len() {
        3 => (&[[0, 1, 2]], &[[0, 1], [1, 2], [2, 0]]),
        // sommets 0, 1, 2 puis milieux 3 (0-1), 4 (1-2), 5 (2-0)
        6 => (
            &[[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]],
            &[[0, 3], [3, 1], [1, 4], [4, 2], [2, 5], [5, 0]],
        ),
        _ => return false,
    };
    triangles.extend(tris.iter().map(|t| t.map(|k| n[k])));
    edges.extend(
        outline
            .iter()
            .map(|&[a, b]| [n[a].min(n[b]), n[a].max(n[b])]),
    );
    true
}

/// Triangles, contour et cellules de tétraèdres de `cells` (3 ou 6 Nodes : triangles ; 4 ou 10 :
/// tétraèdres). Erreur si une cellule est d'une autre taille ou si les deux sortes se mélangent.
fn geometry(cells: &[Vec<usize>]) -> Result<Geometry, SceneError> {
    let is_triangle = |c: &Vec<usize>| matches!(c.len(), 3 | 6);
    let is_solid = |c: &Vec<usize>| matches!(c.len(), 4 | 10);
    if let Some(c) = cells.iter().find(|c| !is_triangle(c) && !is_solid(c)) {
        return Err(SceneError::UnsupportedCell(c.len()));
    }
    if cells.iter().any(is_triangle) && cells.iter().any(is_solid) {
        return Err(SceneError::Mixed);
    }
    if cells.iter().any(is_solid) {
        let skin = skin(cells).ok_or(SceneError::UnsupportedCell(0))?;
        return Ok(Geometry {
            triangles: skin.triangles,
            edges: skin.edges,
            solids: cells.to_vec(),
        });
    }
    let mut triangles = Vec::new();
    let mut edges = HashSet::new();
    for cell in cells {
        let ids: Vec<u32> = cell.iter().map(|&i| i as u32).collect();
        add_cell(&mut triangles, &mut edges, &ids);
    }
    let mut edges: Vec<[u32; 2]> = edges.into_iter().collect();
    edges.sort_unstable();
    Ok(Geometry {
        triangles,
        edges,
        solids: Vec::new(),
    })
}

struct Geometry {
    triangles: Vec<[u32; 3]>,
    edges: Vec<[u32; 2]>,
    solids: Vec<Vec<usize>>,
}

/// Raison pour laquelle un Mesh ne peut pas être affiché.
#[derive(Debug, Clone, PartialEq)]
pub enum SceneError {
    /// Refus du cœur : Mesh de Bar ou de Beam, ou Field d'un autre Mesh.
    Core(MailleError),
    /// Cellule de ce nombre de Nodes : le viewer affiche des triangles (3 ou 6 Nodes) et des
    /// tétraèdres (4 ou 10).
    UnsupportedCell(usize),
    /// Triangles et tétraèdres dans le même Mesh.
    Mixed,
}

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Core(e) => e.fmt(f),
            Self::UnsupportedCell(n) => write!(
                f,
                "le viewer affiche des triangles (3 ou 6 Nodes) et des tétraèdres (4 ou 10), pas des cellules de {n} Nodes"
            ),
            Self::Mixed => write!(
                f,
                "le viewer n'affiche pas un Mesh qui mélange triangles et tétraèdres"
            ),
        }
    }
}

impl std::error::Error for SceneError {}

impl From<MailleError> for SceneError {
    fn from(e: MailleError) -> Self {
        Self::Core(e)
    }
}

impl Scene {
    /// Erreur si le Mesh n'est pas fait de triangles ou de tétraèdres (Bar, Beam, mélange) ou si le
    /// Field vient d'un autre Mesh.
    pub fn new(mesh: &Mesh, field: &Field) -> Result<Self, SceneError> {
        let thermal = !field.temperatures.is_empty();
        if thermal && field.temperatures.len() != mesh.nodes.len() {
            return Err(MailleError::FieldMeshMismatch.into());
        }
        let stresses = if thermal {
            Vec::new()
        } else {
            nodal_stresses(mesh, field)?
        };
        let cells: Vec<Vec<usize>> = mesh
            .elements
            .iter()
            .map(|e| match e {
                // déjà refusé par `nodal_stresses`
                Element::Bar(_) | Element::Beam(_) => Err(MailleError::NoStresses.into()),
                _ => Ok(e.nodes().to_vec()),
            })
            .collect::<Result<_, SceneError>>()?;
        let Geometry {
            triangles,
            edges,
            solids,
        } = geometry(&cells)?;
        let n = mesh.nodes.len();
        let displacements = if thermal {
            vec![[0.0; 3]; n]
        } else if solids.is_empty() {
            field
                .displacements
                .as_chunks::<2>()
                .0
                .iter()
                .map(|u| [u[0], u[1], 0.0])
                .collect()
        } else {
            field.displacements.as_chunks::<3>().0.to_vec()
        };
        Ok(Self {
            nodes: mesh.nodes.iter().map(|p| [p.x, p.y, p.z]).collect(),
            displacements,
            stresses,
            temperatures: field.temperatures.clone(),
            plastic_strains: if field.equivalent_plastic_strains.is_empty() {
                Vec::new()
            } else {
                nodal_equivalent_plastic_strains(mesh, field)?
            },
            distortion: Vec::new(),
            large_strain: false,
            triangles,
            edges,
            solids,
        })
    }

    /// Scène d'un résultat déjà lissé aux Nodes (relu d'un `.vtu`). Erreur si une cellule n'est ni
    /// un triangle ni un tétraèdre ou si un champ n'a pas une valeur par Node.
    pub fn from_nodal(r: &NodalResults) -> Result<Self, SceneError> {
        let n = r.nodes.len();
        let mesh_only = !r.distortion.is_empty();
        let thermal = !r.temperatures.is_empty();
        let sized = |len: usize, wanted: usize| len == 0 || len == wanted;
        if (mesh_only && r.distortion.len() != n)
            || (thermal && r.temperatures.len() != n)
            || (!mesh_only && !thermal && (r.stresses.len() != n || r.displacements.len() != n))
            || !sized(r.plastic_strains.len(), n)
            || r.cells.iter().flatten().any(|&i| i >= n)
        {
            return Err(MailleError::FieldMeshMismatch.into());
        }
        let Geometry {
            triangles,
            edges,
            solids,
        } = geometry(&r.cells)?;
        Ok(Self {
            nodes: r.nodes.clone(),
            displacements: if thermal || mesh_only {
                vec![[0.0; 3]; n]
            } else {
                r.displacements.clone()
            },
            stresses: r.stresses.clone(),
            temperatures: r.temperatures.clone(),
            plastic_strains: r.plastic_strains.clone(),
            distortion: r.distortion.clone(),
            large_strain: r.large_strain,
            triangles,
            edges,
            solids,
        })
    }

    /// Vrai pour un Mesh de tétraèdres.
    pub fn is_solid(&self) -> bool {
        !self.solids.is_empty()
    }

    /// Champs affichables : la température pour un Field thermique, sinon les champs mécaniques,
    /// plus p pour un Field élasto-plastique.
    pub fn quantities(&self) -> &'static [Quantity] {
        if !self.distortion.is_empty() {
            return &[Quantity::Distortion];
        }
        match (
            !self.temperatures.is_empty(),
            self.is_solid(),
            !self.plastic_strains.is_empty(),
        ) {
            (true, ..) => &[Quantity::Temperature],
            (_, false, true) => &Quantity::PLASTIC,
            (_, false, false) => Quantity::MECHANICAL,
            (_, true, true) => &Quantity::SOLID_PLASTIC,
            (_, true, false) => Quantity::SOLID,
        }
    }

    /// Vrai si le Field a des déplacements à amplifier (statique).
    pub fn has_displacements(&self) -> bool {
        self.temperatures.is_empty() && self.distortion.is_empty()
    }

    /// Valeur du champ à chaque Node (vide si le champ n'est pas dans `quantities`).
    pub fn values(&self, quantity: Quantity) -> Vec<f64> {
        let stress = |f: fn(&Stress) -> f64| self.stresses.iter().map(f).collect();
        match quantity {
            Quantity::VonMises => stress(Stress::von_mises),
            Quantity::SigmaXx => stress(|s| s.sigma_xx),
            Quantity::SigmaYy => stress(|s| s.sigma_yy),
            Quantity::SigmaZz => stress(|s| s.sigma_zz),
            Quantity::TauXy => stress(|s| s.tau_xy),
            Quantity::TauYz => stress(|s| s.tau_yz),
            Quantity::TauXz => stress(|s| s.tau_xz),
            Quantity::Displacement => self
                .displacements
                .iter()
                .map(|u| crate::vec3::norm(*u))
                .collect(),
            Quantity::Temperature => self.temperatures.clone(),
            Quantity::PlasticStrain => self.plastic_strains.clone(),
            Quantity::Distortion => self.distortion.clone(),
        }
    }

    /// Coins (min, max) de la boîte englobante des positions x + scale·u (scale = 0 : non déformé).
    pub fn bounds(&self, scale: f64) -> ([f64; 3], [f64; 3]) {
        self.positions(scale).iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(lo, hi), p| {
                (
                    std::array::from_fn(|k| lo[k].min(p[k])),
                    std::array::from_fn(|k| hi[k].max(p[k])),
                )
            },
        )
    }

    /// Facteur d'amplification tel que le plus grand déplacement vaille 10 % de la taille du Mesh.
    /// Nul si rien ne bouge. 1 en grandes déformations : la déformée est déjà celle du calcul.
    pub fn deformation_scale(&self) -> f64 {
        if self.large_strain {
            return 1.0;
        }
        let (lo, hi) = self.bounds(0.0);
        let size = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        let max_u = self
            .values(Quantity::Displacement)
            .into_iter()
            .fold(0.0, f64::max);
        if max_u > 0.0 {
            DEFORMATION_RATIO * size / max_u
        } else {
            0.0
        }
    }

    /// Positions x + scale·u des Nodes (scale = 0 : Mesh non déformé).
    pub fn positions(&self, scale: f64) -> Vec<[f64; 3]> {
        self.nodes
            .iter()
            .zip(&self.displacements)
            .map(|(p, u)| std::array::from_fn(|k| p[k] + scale * u[k]))
            .collect()
    }
}

/// (min, max) des valeurs ; un champ constant donne un intervalle nul, à gérer par l'appelant.
pub fn range(values: &[f64]) -> (f64, f64) {
    values
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_core::material::{Material, Section};
    use maille_core::mesh::{Bar, Node, PlaneState, Tet4, Tri3, Tri6};

    fn material() -> Material {
        Material::with_poisson(1.0, 0.3).unwrap()
    }

    /// Field fabriqué à la main : contrainte uniforme (σxx, σyy, τxy) et déplacements donnés.
    fn field(mesh: &Mesh, sigma: [f64; 3], displacements: Vec<f64>) -> Field {
        let stress = Stress {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            sigma_xx: sigma[0],
            sigma_yy: sigma[1],
            sigma_zz: 0.0,
            tau_xy: sigma[2],
            tau_yz: 0.0,
            tau_xz: 0.0,
        };
        let stresses = mesh
            .elements
            .iter()
            .map(|e| match e {
                Element::Tri6(_) => vec![stress; 3],
                _ => vec![stress],
            })
            .collect();
        Field {
            reactions: vec![0.0; displacements.len()],
            displacements,
            stresses,
            equivalent_plastic_strains: Vec::new(),
            temperatures: Vec::new(),
            heat_fluxes: Vec::new(),
            contact_forces: Vec::new(),
            rigid_forces: Vec::new(),
            contact_pressures: Vec::new(),
        }
    }

    /// Carré unité coupé en deux Tri3 qui partagent l'arête 0-2.
    fn two_tri3() -> Mesh {
        let nodes = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
            .map(|(x, y)| Node::new_2d(x, y))
            .to_vec();
        let tri = |n| Element::Tri3(Tri3::new(n, material(), 1.0, PlaneState::Stress).unwrap());
        Mesh::new(nodes, vec![tri([0, 1, 2]), tri([0, 2, 3])]).unwrap()
    }

    fn one_tri6() -> Mesh {
        let nodes = [
            (0.0, 0.0),
            (2.0, 0.0),
            (0.0, 2.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        .map(|(x, y)| Node::new_2d(x, y))
        .to_vec();
        let tri6 = Tri6::new([0, 1, 2, 3, 4, 5], material(), 1.0, PlaneState::Stress).unwrap();
        Mesh::new(nodes, vec![Element::Tri6(tri6)]).unwrap()
    }

    fn signed_area(p: &[[f64; 3]], t: [u32; 3]) -> f64 {
        let [a, b, c] = t.map(|i| p[i as usize]);
        ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2.0
    }

    #[test]
    fn shared_edges_are_drawn_once() {
        let mesh = two_tri3();
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], vec![0.0; 8])).unwrap();
        assert_eq!(scene.triangles, vec![[0, 1, 2], [0, 2, 3]]);
        assert_eq!(scene.edges, vec![[0, 1], [0, 2], [0, 3], [1, 2], [2, 3]]);
    }

    #[test]
    fn tri6_is_split_into_four_counterclockwise_triangles_covering_it() {
        let mesh = one_tri6();
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], vec![0.0; 12])).unwrap();
        let p = scene.positions(0.0);
        assert_eq!(scene.triangles.len(), 4);
        for &t in &scene.triangles {
            assert!(
                (signed_area(&p, t) - 0.5).abs() < 1e-12,
                "{t:?} : aire 2/4, antihoraire"
            );
        }
        // contour : 3 arêtes coupées en 2 au Node milieu, pas d'arête interne (3-4, 4-5, 5-3)
        assert_eq!(scene.edges.len(), 6);
        assert!(!scene.edges.contains(&[3, 4]) && !scene.edges.contains(&[4, 5]));
    }

    #[test]
    fn values_are_the_nodal_stresses_and_displacement_norms() {
        let mesh = two_tri3();
        let u = vec![3.0, 4.0, 0.0, 0.0, -6.0, 8.0, 0.0, 1.0];
        let scene = Scene::new(&mesh, &field(&mesh, [100.0, 0.0, 50.0], u)).unwrap();
        let von_mises = (100.0_f64.powi(2) + 3.0 * 50.0_f64.powi(2)).sqrt();
        for v in scene.values(Quantity::VonMises) {
            assert!((v - von_mises).abs() < 1e-9);
        }
        assert!(scene
            .values(Quantity::SigmaXx)
            .iter()
            .all(|v| (v - 100.0).abs() < 1e-9));
        assert!(scene
            .values(Quantity::SigmaYy)
            .iter()
            .all(|v| v.abs() < 1e-9));
        assert!(scene
            .values(Quantity::TauXy)
            .iter()
            .all(|v| (v - 50.0).abs() < 1e-9));
        let norms = scene.values(Quantity::Displacement);
        for (got, expected) in norms.iter().zip([5.0, 0.0, 10.0, 1.0]) {
            assert!((got - expected).abs() < 1e-12, "{norms:?}");
        }
        let (lo, hi) = range(&norms);
        assert!(lo.abs() < 1e-12 && (hi - 10.0).abs() < 1e-12);
    }

    #[test]
    fn a_plastic_field_adds_the_plastic_strain() {
        let mesh = two_tri3();
        let elastic = Scene::new(&mesh, &field(&mesh, [100.0, 0.0, 0.0], vec![0.0; 8])).unwrap();
        assert_eq!(elastic.quantities(), Quantity::MECHANICAL);
        // p = 0,01 dans le premier Tri3, 0 dans le second : moyenne 0,005 sur l'arête partagée 0-2
        let plastic = Field {
            equivalent_plastic_strains: vec![vec![0.01], vec![0.0]],
            ..field(&mesh, [100.0, 0.0, 0.0], vec![0.0; 8])
        };
        let scene = Scene::new(&mesh, &plastic).unwrap();
        assert_eq!(scene.quantities(), &Quantity::PLASTIC);
        let p = scene.values(Quantity::PlasticStrain);
        for (got, expected) in p.iter().zip([0.005, 0.01, 0.005, 0.0]) {
            assert!((got - expected).abs() < 1e-15, "{p:?}");
        }
    }

    #[test]
    fn large_strain_shows_the_deformed_shape_at_scale_one() {
        let nodes = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let zero = |p: &[f64; 3]| Stress {
            x: p[0],
            y: p[1],
            z: p[2],
            sigma_xx: 0.0,
            sigma_yy: 0.0,
            sigma_zz: 0.0,
            tau_xy: 0.0,
            tau_yz: 0.0,
            tau_xz: 0.0,
        };
        let displacements = vec![[0.0; 3], [0.5, 0.0, 0.0], [0.0, -0.25, 0.0]];
        let results = |large_strain| NodalResults {
            stresses: nodes.iter().map(zero).collect(),
            nodes: nodes.clone(),
            cells: vec![vec![0, 1, 2]],
            displacements: displacements.clone(),
            large_strain,
            ..NodalResults::default()
        };
        // petites déformations : le plus grand déplacement vaut un dixième du Mesh
        let small = Scene::from_nodal(&results(false)).unwrap();
        assert!((small.deformation_scale() - 0.4).abs() < 1e-12);
        // grandes déformations : x + u, sans amplification
        let large = Scene::from_nodal(&results(true)).unwrap();
        assert!((large.deformation_scale() - 1.0).abs() < 1e-15);
        let moved = large.positions(large.deformation_scale());
        assert_eq!(moved[1], [2.5, 0.0, 0.0]);
        assert_eq!(moved[2], [0.0, 1.75, 0.0]);
    }

    #[test]
    fn deformed_shape_moves_the_largest_displacement_by_a_tenth_of_the_mesh() {
        let mesh = two_tri3();
        let u = vec![0.0, 0.0, 1e-3, 0.0, 2e-3, -1e-3, 0.0, 0.0];
        let scene = Scene::new(&mesh, &field(&mesh, [1.0, 0.0, 0.0], u)).unwrap();
        let scale = scene.deformation_scale();
        let (before, after) = (scene.positions(0.0), scene.positions(scale));
        let moved = before
            .iter()
            .zip(&after)
            .map(|(a, b)| (b[0] - a[0]).hypot(b[1] - a[1]))
            .fold(0.0, f64::max);
        assert!((moved - 0.1).abs() < 1e-12, "{moved}");
        let (lo, hi) = scene.bounds(0.0);
        assert!(lo.iter().all(|c| c.abs() < 1e-12));
        assert!(hi[..2].iter().all(|c| (c - 1.0).abs() < 1e-12));
        // un Mesh 2D est dans le plan z = 0
        assert_eq!((lo[2], hi[2]), (0.0, 0.0));
        // la boîte de la déformée englobe le Node 2 déplacé de 0,1 (u ∝ (2, −1))
        let (lo, hi) = scene.bounds(scale);
        let moved = [1.0 + 0.2 / 5_f64.sqrt(), 1.0 - 0.1 / 5_f64.sqrt()];
        assert!((hi[0] - moved[0]).abs() < 1e-12 && lo[1].abs() < 1e-12);
    }

    #[test]
    fn a_mesh_that_does_not_move_has_no_deformation() {
        let mesh = two_tri3();
        let scene = Scene::new(&mesh, &field(&mesh, [0.0; 3], vec![0.0; 8])).unwrap();
        assert!(scene.deformation_scale().abs() < 1e-12);
    }

    #[test]
    fn bars_cannot_be_shown() {
        let nodes = vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)];
        let bar = Bar {
            nodes: [0, 1],
            material: Material::new(1.0).unwrap(),
            section: Section::new(1.0).unwrap(),
        };
        let mesh = Mesh::new(nodes, vec![Element::Bar(bar)]).unwrap();
        let field = Field {
            displacements: vec![0.0; 2],
            reactions: vec![0.0; 2],
            stresses: vec![Vec::new()],
            equivalent_plastic_strains: Vec::new(),
            temperatures: Vec::new(),
            heat_fluxes: Vec::new(),
            contact_forces: Vec::new(),
            rigid_forces: Vec::new(),
            contact_pressures: Vec::new(),
        };
        assert_eq!(
            Scene::new(&mesh, &field).err(),
            Some(SceneError::Core(MailleError::NoStresses))
        );
    }

    #[test]
    fn a_thermal_field_shows_only_the_temperature_without_deformation() {
        let mesh = two_tri3();
        let temperatures = vec![10.0, 20.0, 30.0, 40.0];
        let field = Field {
            displacements: Vec::new(),
            reactions: vec![0.0; 4],
            stresses: Vec::new(),
            equivalent_plastic_strains: Vec::new(),
            temperatures: temperatures.clone(),
            heat_fluxes: Vec::new(),
            contact_forces: Vec::new(),
            rigid_forces: Vec::new(),
            contact_pressures: Vec::new(),
        };
        let scene = Scene::new(&mesh, &field).unwrap();
        assert_eq!(scene.quantities(), &[Quantity::Temperature]);
        assert!(!scene.has_displacements());
        for (v, t) in scene
            .values(Quantity::Temperature)
            .iter()
            .zip(&temperatures)
        {
            assert!((v - t).abs() < 1e-12);
        }
        assert!(scene.deformation_scale().abs() < 1e-12);

        // Field thermique d'un autre Mesh (3 températures pour 4 Nodes)
        let foreign = Field {
            temperatures: vec![1.0; 3],
            ..field
        };
        assert_eq!(
            Scene::new(&mesh, &foreign).err(),
            Some(SceneError::Core(MailleError::FieldMeshMismatch))
        );
    }

    #[test]
    fn a_tetrahedron_shows_its_four_faces_and_six_edges() {
        let nodes = vec![
            Node::new(0.0, 0.0, 0.0),
            Node::new(1.0, 0.0, 0.0),
            Node::new(0.0, 1.0, 0.0),
            Node::new(0.0, 0.0, 1.0),
        ];
        let tet = Tet4::new([0, 1, 2, 3], material());
        let mesh = Mesh::new(nodes, vec![Element::Tet4(tet)]).unwrap();
        // contrainte uniforme σzz = 7 et déplacement (0, 0, 2) partout
        let mut f = field(&mesh, [0.0; 3], vec![0.0; 12]);
        f.stresses[0][0].sigma_zz = 7.0;
        for u in f.displacements.chunks_mut(3) {
            u[2] = 2.0;
        }
        let scene = Scene::new(&mesh, &f).unwrap();
        assert!(scene.is_solid());
        assert_eq!(scene.triangles.len(), 4);
        assert_eq!(scene.edges.len(), 6);
        assert_eq!(scene.quantities(), Quantity::SOLID);
        assert!(scene
            .values(Quantity::SigmaZz)
            .iter()
            .all(|v| (v - 7.0).abs() < 1e-9));
        assert!(scene
            .values(Quantity::Displacement)
            .iter()
            .all(|v| (v - 2.0).abs() < 1e-12));
        let (lo, hi) = scene.bounds(0.0);
        assert_eq!((lo, hi), ([0.0; 3], [1.0; 3]));
    }

    #[test]
    fn triangles_and_tetrahedra_do_not_mix() {
        let mut mixed = nodal_triangle();
        mixed.cells = vec![vec![0, 1, 2], vec![0, 1, 2, 3]];
        assert_eq!(Scene::from_nodal(&mixed).err(), Some(SceneError::Mixed));
    }

    /// Un Tri6 à champ de contrainte uniforme, comme `Scene::new` le voit et comme un `.vtu` le donne.
    fn nodal_triangle() -> NodalResults {
        let stress = Stress {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            sigma_xx: 3.0,
            sigma_yy: 0.0,
            sigma_zz: 0.0,
            tau_xy: 0.0,
            tau_yz: 0.0,
            tau_xz: 0.0,
        };
        NodalResults {
            nodes: vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            cells: vec![vec![0, 1, 2, 3, 4, 5]],
            displacements: vec![[0.5, 0.0, 0.0]; 6],
            stresses: vec![stress; 6],
            ..NodalResults::default()
        }
    }

    #[test]
    fn nodal_results_build_the_same_scene_as_a_field() {
        let nodes: Vec<Node> = nodal_triangle()
            .nodes
            .iter()
            .map(|p| Node::new_2d(p[0], p[1]))
            .collect();
        let tri = Tri6::new([0, 1, 2, 3, 4, 5], material(), 1.0, PlaneState::Stress).unwrap();
        let mesh = Mesh::new(nodes, vec![Element::Tri6(tri)]).unwrap();
        let from_field =
            Scene::new(&mesh, &field(&mesh, [3.0, 0.0, 0.0], [0.5, 0.0].repeat(6))).unwrap();
        let from_nodal = Scene::from_nodal(&nodal_triangle()).unwrap();
        assert_eq!(from_nodal.triangles, from_field.triangles);
        assert_eq!(from_nodal.edges, from_field.edges);
        for q in [
            Quantity::VonMises,
            Quantity::SigmaXx,
            Quantity::Displacement,
        ] {
            let (a, b) = (from_nodal.values(q), from_field.values(q));
            assert!(
                a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-12),
                "{q:?}"
            );
        }
    }

    #[test]
    fn a_mesh_without_results_shows_only_its_distortion() {
        let (nodes, cells) = crate::testmesh::cube();
        let distortion = vec![0.5; nodes.len()];
        let scene = Scene::from_nodal(&NodalResults {
            nodes,
            cells,
            distortion: distortion.clone(),
            ..NodalResults::default()
        })
        .unwrap();
        assert_eq!(scene.quantities(), &[Quantity::Distortion]);
        assert_eq!(scene.values(Quantity::Distortion), distortion);
        assert!(!scene.has_displacements() && scene.is_solid());
        // un défaut qui n'a pas une valeur par Node est refusée
        let wrong = NodalResults {
            nodes: vec![[0.0; 3]; 4],
            cells: vec![vec![0, 1, 2, 3]],
            distortion: vec![1.0],
            ..NodalResults::default()
        };
        assert!(Scene::from_nodal(&wrong).is_err());
    }

    #[test]
    fn nodal_results_refuse_unknown_cells_and_wrong_sizes() {
        let mut odd = nodal_triangle();
        odd.cells = vec![vec![0, 1, 2, 3, 4]];
        assert_eq!(
            Scene::from_nodal(&odd).err(),
            Some(SceneError::UnsupportedCell(5))
        );

        let mut short = nodal_triangle();
        short.stresses.pop();
        assert_eq!(
            Scene::from_nodal(&short).err(),
            Some(SceneError::Core(MailleError::FieldMeshMismatch))
        );

        let mut outside = nodal_triangle();
        outside.cells = vec![vec![0, 1, 9]];
        assert_eq!(
            Scene::from_nodal(&outside).err(),
            Some(SceneError::Core(MailleError::FieldMeshMismatch))
        );
    }
}
