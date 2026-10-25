//! Mailleur 2D : triangulation de Delaunay (Bowyer-Watson), bords respectés, raffinement de Ruppert.
//!
//! Le mailleur ne produit que la connectivité : l'utilisateur crée ses Elements (Tri3, Tri6) avec
//! son Material à partir de `Triangulation::triangles`, et pose ses BC et ses Convection sur les
//! bords nommés (`Triangulation::boundary`).

mod delaunay;
mod geometry;
mod predicates;
mod refine;

use std::collections::{BTreeMap, HashMap};

use crate::error::MailleError;
use crate::mesh::Node;
pub use geometry::{Curve, CurveId, CurveKind, Geometry, MIN_CORNER};

/// Plafond de Nodes produits par le mailleur : au-delà, `MeshingFailed` plutôt qu'un calcul sans fin.
pub const MAX_NODES: usize = 1_000_000;

/// Maillage triangulaire produit par le mailleur : Nodes, connectivité et bords nommés.
#[derive(Debug, Clone, PartialEq)]
pub struct Triangulation<const N: usize> {
    pub nodes: Vec<Node>,
    /// Connectivité des triangles, sommets dans le sens antihoraire (comme `Tri3::nodes`).
    pub triangles: Vec<[usize; N]>,
    boundaries: BTreeMap<String, BoundarySet>,
}

/// Ce que le mailleur sait d'un bord nommé (toutes les Curves qui portent ce nom).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BoundarySet {
    /// Nodes posés sur le bord, triés.
    pub nodes: Vec<usize>,
    /// Faces de triangles sur le bord, (triangle, face) triés, faces numérotées comme
    /// `Element::face_nodes` : 0 pour l'arête 0-1, 1 pour 1-2, 2 pour 2-0.
    pub faces: Vec<(usize, usize)>,
}

impl<const N: usize> Triangulation<N> {
    /// Bord nommé `name`. `UnknownBoundary` (avec la liste des noms) si aucune Curve ne le porte.
    pub fn boundary(&self, name: &str) -> Result<&BoundarySet, MailleError> {
        self.boundaries
            .get(name)
            .ok_or_else(|| MailleError::UnknownBoundary {
                name: name.to_string(),
                known: self.boundaries.keys().cloned().collect(),
            })
    }
}

/// Réglages du mailleur.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshOptions {
    /// Angle minimal garanti des triangles, en degrés, de 0 (exclu) à `MAX_MIN_ANGLE`.
    pub min_angle: f64,
    /// Gradation g : de combien la taille des triangles peut croître par unité de longueur en
    /// s'éloignant d'un bord fin, strictement positive.
    pub grading: f64,
}

impl Default for MeshOptions {
    fn default() -> Self {
        Self {
            min_angle: 25.0,
            grading: 0.3,
        }
    }
}

/// Plus grand angle minimal accepté, en degrés : Ruppert est démontré jusqu'à 20,7° et termine en
/// pratique jusque vers 33°.
pub const MAX_MIN_ANGLE: f64 = 33.0;

/// Maille la géométrie en triangles. Chaque Curve est découpée en morceaux de longueur au plus sa
/// taille h ; tout morceau absent de la triangulation de Delaunay est coupé en son milieu (sur la
/// Curve) jusqu'à ce que le bord soit respecté. Puis Ruppert raffine jusqu'à l'angle minimal et
/// à la taille locale : celle du bord le plus proche, croissant de `grading` par unité de longueur.
pub fn triangulate(
    geometry: &Geometry,
    options: &MeshOptions,
) -> Result<Triangulation<3>, MailleError> {
    if !(options.min_angle > 0.0 && options.min_angle <= MAX_MIN_ANGLE) {
        return Err(MailleError::MinAngleOutOfRange {
            angle: options.min_angle,
        });
    }
    let grading = crate::material::positive("gradation", options.grading)?;
    let boundary = geometry.discretize()?;
    let mut m = refine::Mesher::new(&boundary)?;
    m.refine(options, &refine::SizeField::new(&boundary, grading))?;
    m.carve()?;

    // Nodes renumérotés dans l'ordre de création, sans ceux qu'aucun triangle n'utilise
    let mut triangles = m.d.inside_triangles();
    let mut number = vec![usize::MAX; m.d.point_count()];
    for t in &triangles {
        for &v in t {
            number[v] = 0;
        }
    }
    let mut nodes = Vec::new();
    for (i, n) in number.iter_mut().enumerate() {
        if *n == 0 {
            *n = nodes.len();
            let p = m.d.point(i);
            nodes.push(Node::new_2d(p[0], p[1]));
        }
    }
    for t in &mut triangles {
        *t = t.map(|v| number[v]);
    }

    // faces de bord : le triangle du domaine qui parcourt a → b
    let mut faces = HashMap::new();
    for (e, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            faces.insert((t[k], t[(k + 1) % 3]), (e, k));
        }
    }
    let mut boundaries: BTreeMap<String, BoundarySet> = BTreeMap::new();
    for s in &m.segments {
        let (a, b) = (number[s.a], number[s.b]);
        let face = *faces.get(&(a, b)).ok_or(MailleError::MeshingFailed)?;
        let set = boundaries
            .entry(boundary.pieces[s.piece].id.name.clone())
            .or_default();
        set.nodes.extend([a, b]);
        set.faces.push(face);
    }
    for set in boundaries.values_mut() {
        set.nodes.sort_unstable();
        set.nodes.dedup();
        set.faces.sort_unstable();
    }
    Ok(Triangulation {
        nodes,
        triangles,
        boundaries,
    })
}
