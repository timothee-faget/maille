//! Cas d'exemple : plaque trouée de Kirsch, traction uniaxiale σ selon x, trou de rayon a,
//! contraintes planes, Tri6. On modélise un quart de plaque (symétries), maillé par le mailleur
//! (fin au trou, grossier au bord), et on impose sur le bord extérieur r = R les déplacements
//! exacts de Kirsch : la plaque finie n'introduit aucune erreur de modèle.
//!
//! Unités cohérentes (N, mm, MPa) : le cœur ne convertit rien.

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Bc, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

use crate::report;

const YOUNG: f64 = 210e3; // MPa
const POISSON: f64 = 0.3;
const SIGMA: f64 = 100.0; // MPa
const HOLE: f64 = 10.0; // mm, rayon a
const OUTER: f64 = 40.0; // mm, rayon R
const THICKNESS: f64 = 1.0; // mm
const SIZE_HOLE: f64 = 1.0; // mm, taille des triangles au trou
const SIZE_OUTER: f64 = 5.0; // mm, ailleurs

/// Déplacements de Kirsch (ux, uy) en contraintes planes, κ = (3 − ν)/(1 + ν).
fn kirsch_displacement(x: f64, y: f64) -> (f64, f64) {
    let (r, t) = (x.hypot(y), y.atan2(x));
    let a = HOLE;
    let kappa = (3.0 - POISSON) / (1.0 + POISSON);
    let mu = YOUNG / (2.0 * (1.0 + POISSON));
    let f = SIGMA / (8.0 * mu);
    let ur = f
        * ((kappa - 1.0) * r
            + 2.0 * a * a / r
            + 2.0 * (r + (kappa + 1.0) * a * a / r - a.powi(4) / r.powi(3)) * (2.0 * t).cos());
    let ut = -f * 2.0 * (r + (kappa - 1.0) * a * a / r + a.powi(4) / r.powi(3)) * (2.0 * t).sin();
    (ur * t.cos() - ut * t.sin(), ur * t.sin() + ut * t.cos())
}

/// Quart de plaque maillé : bords nommés pour les BC.
pub fn solve() -> Result<(Mesh, Field), MailleError> {
    let geometry = Geometry {
        outer: vec![
            Curve::line([HOLE, 0.0], SIZE_OUTER, "sym_y"),
            Curve::arc([OUTER, 0.0], [0.0, 0.0], false, SIZE_OUTER, "exterieur"),
            Curve::line([0.0, OUTER], SIZE_OUTER, "sym_x"),
            Curve::arc([0.0, HOLE], [0.0, 0.0], true, SIZE_HOLE, "trou"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())?.to_quadratic();
    let material = Material::with_poisson(YOUNG, POISSON)?;
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, THICKNESS, PlaneState::Stress).map(Element::Tri6))
        .collect::<Result<_, _>>()?;

    let outer = &t.boundary("exterieur")?.nodes;
    let mut bcs = Vec::new();
    for &n in outer {
        let (u, v) = kirsch_displacement(t.nodes[n].x, t.nodes[n].y);
        bcs.push(Bc {
            node: n,
            dof: Dof::Ux,
            value: u,
        });
        bcs.push(Bc {
            node: n,
            dof: Dof::Uy,
            value: v,
        });
    }
    // symétries : Ux = 0 sur l'axe y, Uy = 0 sur l'axe x, sauf aux coins déjà imposés
    for (name, dof) in [("sym_x", Dof::Ux), ("sym_y", Dof::Uy)] {
        for &n in t
            .boundary(name)?
            .nodes
            .iter()
            .filter(|n| !outer.contains(n))
        {
            bcs.push(Bc {
                node: n,
                dof,
                value: 0.0,
            });
        }
    }
    let mesh = Mesh::new(t.nodes, elements)?;
    let field = Model::new(mesh.clone(), bcs, vec![])?.solve()?;
    Ok((mesh, field))
}

/// Indice du Node posé en (x, y) : les débuts de Curves sont toujours des Nodes du maillage.
fn node_at(mesh: &Mesh, x: f64, y: f64) -> Result<usize, String> {
    mesh.nodes
        .iter()
        .position(|p| (p.x - x).hypot(p.y - y) < 1e-9 * OUTER)
        .ok_or_else(|| format!("aucun Node en ({x}, {y})"))
}

/// Affiche le rapport face à Kirsch, puis ouvre le viewer.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, field) = solve()?;
    let nodal = nodal_stresses(&mesh, &field)?;

    println!(
        "Plaque trouée de Kirsch : a = {HOLE} mm, R = {OUTER} mm, σ = {SIGMA} MPa, E = {YOUNG} MPa, ν = {POISSON}"
    );
    println!(
        "Quart de plaque maillé (h = {SIZE_HOLE} mm au trou, {SIZE_OUTER} mm ailleurs), {} Tri6, {} Nodes, {} DOF, contraintes lissées aux Nodes",
        mesh.elements.len(),
        mesh.nodes.len(),
        field.displacements.len()
    );
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    // au bord du trou : σθ = σ (1 − 2 cos 2θ), soit 3σ en θ = π/2 (σxx) et −σ en θ = 0 (σyy)
    report(
        "σxx en (0, a) (MPa)",
        nodal[node_at(&mesh, 0.0, HOLE)?].sigma_xx,
        3.0 * SIGMA,
    );
    report(
        "σyy en (a, 0) (MPa)",
        nodal[node_at(&mesh, HOLE, 0.0)?].sigma_yy,
        -SIGMA,
    );
    maille_viewer::show(&mesh, &field)
}
