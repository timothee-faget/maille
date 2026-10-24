//! Cas d'exemple : plaque trouée de Kirsch, traction uniaxiale σ selon x, trou de rayon a,
//! contraintes planes, Tri6. On modélise un quart de plaque (symétries) et on impose sur le bord
//! extérieur r = R les déplacements exacts de Kirsch : la plaque finie n'introduit aucune erreur de modèle.
//!
//! Unités cohérentes (N, mm, MPa) : le cœur ne convertit rien.

use std::f64::consts::FRAC_PI_2;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::model::{Bc, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

use crate::grid::{rect_grid, to_quadratic};
use crate::report;

const YOUNG: f64 = 210e3; // MPa
const POISSON: f64 = 0.3;
const SIGMA: f64 = 100.0; // MPa
const HOLE: f64 = 10.0; // mm, rayon a
const OUTER: f64 = 40.0; // mm, rayon R
const THICKNESS: f64 = 1.0; // mm
const DIVISIONS: usize = 16; // en rayon et en angle

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

/// Indice du sommet (i, j) de la grille polaire : i en rayon (0 au trou), j en angle (0 sur l'axe x).
fn idx(i: usize, j: usize) -> usize {
    i + j * (DIVISIONS + 1)
}

/// Quart de plaque : grille polaire, rayons en progression géométrique de a à R (plus fins au trou).
pub fn solve() -> Result<(Mesh, Field), MailleError> {
    let (grid, tris) = rect_grid(DIVISIONS, DIVISIONS, 1.0, 1.0);
    let corners: Vec<Node> = grid
        .iter()
        .map(|g| {
            let (r, t) = (HOLE * (OUTER / HOLE).powf(g.x), FRAC_PI_2 * g.y);
            Node {
                x: r * t.cos(),
                y: r * t.sin(),
            }
        })
        .collect();
    let q = to_quadratic(&corners, &tris);
    let material = Material::with_poisson(YOUNG, POISSON)?;
    let elements = q
        .elements
        .iter()
        .map(|&n| Tri6::new(n, material, THICKNESS, PlaneState::Stress).map(Element::Tri6))
        .collect::<Result<_, _>>()?;

    // bord extérieur : sommets i = DIVISIONS et milieux de leurs arêtes
    let mut outer: Vec<usize> = (0..=DIVISIONS).map(|j| idx(DIVISIONS, j)).collect();
    outer.extend((0..DIVISIONS).map(|j| q.mid(idx(DIVISIONS, j), idx(DIVISIONS, j + 1))));
    outer.sort_unstable();

    let mut bcs = Vec::new();
    for (n, p) in q.nodes.iter().enumerate() {
        if outer.binary_search(&n).is_ok() {
            let (u, v) = kirsch_displacement(p.x, p.y);
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
        } else {
            // symétries : Ux = 0 sur l'axe y, Uy = 0 sur l'axe x
            if p.x.abs() < 1e-9 * OUTER {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Ux,
                    value: 0.0,
                });
            }
            if p.y.abs() < 1e-9 * OUTER {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Uy,
                    value: 0.0,
                });
            }
        }
    }
    let mesh = Mesh::new(q.nodes, elements)?;
    let field = Model::new(mesh.clone(), bcs, vec![])?.solve()?;
    Ok((mesh, field))
}

/// Affiche le rapport face à Kirsch, puis ouvre le viewer.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, field) = solve()?;
    let nodal = nodal_stresses(&mesh, &field)?;

    println!(
        "Plaque trouée de Kirsch : a = {HOLE} mm, R = {OUTER} mm, σ = {SIGMA} MPa, E = {YOUNG} MPa, ν = {POISSON}"
    );
    println!(
        "Quart de plaque, {} Tri6, {} Nodes, {} DOF, contraintes lissées aux Nodes",
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
        nodal[idx(0, DIVISIONS)].sigma_xx,
        3.0 * SIGMA,
    );
    report("σyy en (a, 0) (MPa)", nodal[idx(0, 0)].sigma_yy, -SIGMA);
    maille_viewer::show(&mesh, &field)
}
