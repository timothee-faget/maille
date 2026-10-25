//! Mailleur 2D : triangulation de Delaunay (Bowyer-Watson), bords respectés, raffinement de Ruppert.
//!
//! Le mailleur ne produit que la connectivité : l'utilisateur crée ses Elements (Tri3, Tri6) avec
//! son Material à partir de `Triangulation::triangles`, et pose ses BC et ses Convection sur les
//! bords nommés (`Triangulation::boundary`).

mod delaunay;
mod geometry;
mod predicates;

use std::collections::BTreeMap;

use crate::error::MailleError;
use crate::mesh::Node;
use delaunay::Delaunay;
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

/// Maille la géométrie en triangles. Chaque Curve est découpée en morceaux de longueur au plus sa
/// taille h, puis tout morceau absent de la triangulation de Delaunay est coupé en son milieu (sur
/// la Curve) jusqu'à ce que le bord soit respecté.
pub fn triangulate(geometry: &Geometry) -> Result<Triangulation<3>, MailleError> {
    let boundary = geometry.discretize()?;
    let (mut min, mut max) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for p in &boundary.points {
        for j in 0..2 {
            min[j] = min[j].min(p[j]);
            max[j] = max[j].max(p[j]);
        }
    }
    let mut d = Delaunay::new(min, max);
    for &p in &boundary.points {
        d.insert(p)?; // indice dans la triangulation = indice dans boundary.points
    }

    // bord respecté : un morceau présent devient contraint, un morceau absent est coupé en deux
    let mut segments = boundary.segments.clone();
    let mut pending: Vec<usize> = (0..segments.len()).collect();
    while let Some(s) = pending.pop() {
        let seg = segments[s];
        if d.has_edge(seg.a, seg.b) {
            d.set_constrained(seg.a, seg.b, true);
            continue;
        }
        if d.point_count() >= MAX_NODES {
            return Err(MailleError::MeshingFailed);
        }
        let t = (seg.t[0] + seg.t[1]) / 2.0;
        let m = d.insert(boundary.pieces[seg.piece].point(t))?;
        segments[s].b = m;
        segments[s].t[1] = t;
        segments.push(geometry::Segment {
            a: m,
            t: [t, seg.t[1]],
            ..seg
        });
        pending.extend([s, segments.len() - 1]);
    }

    let triangles = d.domain(segments.iter().map(|s| (s.a, s.b)))?;
    let nodes = (0..d.point_count())
        .map(|i| {
            let p = d.point(i);
            Node::new_2d(p[0], p[1])
        })
        .collect();

    // faces de bord : le triangle du domaine qui parcourt a → b
    let mut faces = std::collections::HashMap::new();
    for (e, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            faces.insert((t[k], t[(k + 1) % 3]), (e, k));
        }
    }
    let mut boundaries: BTreeMap<String, BoundarySet> = BTreeMap::new();
    for s in &segments {
        let face = *faces.get(&(s.a, s.b)).ok_or(MailleError::MeshingFailed)?;
        let set = boundaries
            .entry(boundary.pieces[s.piece].id.name.clone())
            .or_default();
        set.nodes.extend([s.a, s.b]);
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
