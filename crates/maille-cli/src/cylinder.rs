//! Cas d'exemple : conduction dans un cylindre creux (tube épais), régime stationnaire, Tri6.
//! Températures imposées Tᵢ sur le rayon intérieur rᵢ et Tₑ sur le rayon extérieur rₑ :
//! T(r) = Tᵢ + (Tₑ − Tᵢ) ln(r/rᵢ) / ln(rₑ/rᵢ). On modélise un quart de section, maillé par le
//! mailleur (les symétries sont adiabatiques, condition naturelle), et on impose sur les bords la
//! solution exacte aux Nodes, qui ne sont pas tous sur le cercle (côtés droits) : la géométrie
//! approchée n'introduit pas d'erreur de modèle.
//!
//! Unités cohérentes (W, mm, K) : le cœur ne convertit rien.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Analysis, Bc, Model};
use maille_core::results::Field;
use maille_core::MailleError;

use crate::report;

const CONDUCTIVITY: f64 = 0.05; // W/(mm·K), acier
const INNER: f64 = 10.0; // mm, rᵢ
const OUTER: f64 = 40.0; // mm, rₑ
const T_INNER: f64 = 200.0; // °C
const T_OUTER: f64 = 20.0; // °C
const THICKNESS: f64 = 1.0; // mm (longueur de tube modélisée)
const SIZE: f64 = 2.0; // mm, taille des triangles

fn exact(p: &Node) -> f64 {
    let r = p.x.hypot(p.y);
    T_INNER + (T_OUTER - T_INNER) * (r / INNER).ln() / (OUTER / INNER).ln()
}

/// Quart de section maillé. Renvoie aussi les Nodes du bord intérieur.
pub fn solve() -> Result<(Mesh, Field, Vec<usize>), MailleError> {
    let geometry = Geometry {
        outer: vec![
            Curve::line([INNER, 0.0], SIZE, "sym"),
            Curve::arc([OUTER, 0.0], [0.0, 0.0], false, SIZE, "exterieur"),
            Curve::line([0.0, OUTER], SIZE, "sym"),
            Curve::arc([0.0, INNER], [0.0, 0.0], true, SIZE, "interieur"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())?.to_quadratic();
    let material = Material::conductive(CONDUCTIVITY)?;
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, THICKNESS, PlaneState::Stress).map(Element::Tri6))
        .collect::<Result<_, _>>()?;

    let inner = t.boundary("interieur")?.nodes.clone();
    let bcs = inner
        .iter()
        .chain(&t.boundary("exterieur")?.nodes)
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: exact(&t.nodes[node]),
        })
        .collect();
    let mesh = Mesh::new(t.nodes, elements)?;
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
        "Quart de section maillé (h = {SIZE} mm), {} Tri6, {} Nodes, {} DOF",
        mesh.elements.len(),
        mesh.nodes.len(),
        field.temperatures.len()
    );
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    // Node le plus proche de r = 25 mm, θ = 45° : la solution exacte y est prise à sa position
    let target = [25.0 * FRAC_PI_4.cos(), 25.0 * FRAC_PI_4.sin()];
    let middle = (0..mesh.nodes.len())
        .min_by(|&a, &b| {
            let d = |n: usize| (mesh.nodes[n].x - target[0]).hypot(mesh.nodes[n].y - target[1]);
            d(a).total_cmp(&d(b))
        })
        .unwrap_or(0);
    report(
        "T vers r = 25 mm (°C)",
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
