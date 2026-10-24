//! Cas d'exemple : conduction dans un cylindre creux (tube épais), régime stationnaire, Tri6.
//! Températures imposées Tᵢ sur le rayon intérieur rᵢ et Tₑ sur le rayon extérieur rₑ :
//! T(r) = Tᵢ + (Tₑ − Tᵢ) ln(r/rᵢ) / ln(rₑ/rᵢ). On modélise un quart de section (les symétries sont
//! adiabatiques, condition naturelle) et on impose sur les bords la solution exacte aux Nodes, qui
//! ne sont pas tous sur le cercle (côtés droits) : la géométrie approchée n'introduit pas d'erreur
//! de modèle.
//!
//! Unités cohérentes (W, mm, K) : le cœur ne convertit rien.

use std::f64::consts::FRAC_PI_2;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::model::{Analysis, Bc, Model};
use maille_core::results::Field;
use maille_core::MailleError;

use crate::grid::{rect_grid, to_quadratic};
use crate::report;

const CONDUCTIVITY: f64 = 0.05; // W/(mm·K), acier
const INNER: f64 = 10.0; // mm, rᵢ
const OUTER: f64 = 40.0; // mm, rₑ
const T_INNER: f64 = 200.0; // °C
const T_OUTER: f64 = 20.0; // °C
const THICKNESS: f64 = 1.0; // mm (longueur de tube modélisée)
const DIVISIONS: usize = 16; // en rayon et en angle

fn exact(p: &Node) -> f64 {
    let r = p.x.hypot(p.y);
    T_INNER + (T_OUTER - T_INNER) * (r / INNER).ln() / (OUTER / INNER).ln()
}

/// Indice du sommet (i, j) de la grille polaire : i en rayon (0 au bord intérieur), j en angle.
fn idx(i: usize, j: usize) -> usize {
    i + j * (DIVISIONS + 1)
}

/// Quart de section, grille polaire à rayons réguliers. Renvoie aussi les Nodes du bord intérieur.
pub fn solve() -> Result<(Mesh, Field, Vec<usize>), MailleError> {
    let (grid, tris) = rect_grid(DIVISIONS, DIVISIONS, 1.0, 1.0);
    let corners: Vec<Node> = grid
        .iter()
        .map(|g| {
            let (r, t) = (INNER + (OUTER - INNER) * g.x, FRAC_PI_2 * g.y);
            Node::new_2d(r * t.cos(), r * t.sin())
        })
        .collect();
    let q = to_quadratic(&corners, &tris);
    let material = Material::conductive(CONDUCTIVITY)?;
    let elements = q
        .elements
        .iter()
        .map(|&n| Tri6::new(n, material, THICKNESS, PlaneState::Stress).map(Element::Tri6))
        .collect::<Result<_, _>>()?;

    // bords : sommets i = 0 et i = DIVISIONS, et milieux de leurs arêtes
    let edge = |i: usize| {
        let mut nodes: Vec<usize> = (0..=DIVISIONS).map(|j| idx(i, j)).collect();
        nodes.extend((0..DIVISIONS).map(|j| q.mid(idx(i, j), idx(i, j + 1))));
        nodes
    };
    let inner = edge(0);
    let bcs = inner
        .iter()
        .chain(&edge(DIVISIONS))
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: exact(&q.nodes[node]),
        })
        .collect();
    let mesh = Mesh::new(q.nodes, elements)?;
    let field = Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, vec![])?.solve()?;
    Ok((mesh, field, inner))
}

/// Affiche le rapport face à la solution analytique, puis ouvre le viewer.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, field, inner) = solve()?;

    println!(
        "Cylindre creux : rᵢ = {INNER} mm, rₑ = {OUTER} mm, Tᵢ = {T_INNER} °C, Tₑ = {T_OUTER} °C, k = {CONDUCTIVITY} W/(mm·K)"
    );
    println!(
        "Quart de section, {} Tri6, {} Nodes, {} DOF",
        mesh.elements.len(),
        mesh.nodes.len(),
        field.temperatures.len()
    );
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    let middle = idx(DIVISIONS / 2, DIVISIONS / 2);
    report(
        "T à r = 25 mm (°C)",
        field.temperatures[middle],
        exact(&mesh.nodes[middle]),
    );
    // puissance par mm de tube sur le quart : (π/2) k t (Tᵢ − Tₑ) / ln(rₑ/rᵢ)
    let power: f64 = inner.iter().map(|&n| field.reactions[n]).sum();
    report(
        "puissance, quart (W)",
        power,
        FRAC_PI_2 * CONDUCTIVITY * THICKNESS * (T_INNER - T_OUTER) / (OUTER / INNER).ln(),
    );
    maille_viewer::show(&mesh, &field)
}
