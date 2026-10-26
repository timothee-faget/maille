//! Cas d'exemple : trempe d'une plaque d'acier, thermique transitoire, Tri6. La plaque, d'épaisseur
//! L, est à T₀ ; ses deux faces x = 0 et x = L passent à Tₛ à t = 0 (trempe idéale : température de
//! surface imposée, sans résistance d'échange). Ses bords y = 0 et y = H sont adiabatiques, donc
//! T ne dépend que de x : T(x, t) = Tₛ + (T₀ − Tₛ) Σ_{n impair} 4/(nπ) sin(nπx/L) exp(−α (nπ/L)² t),
//! avec α = k/(ρc).
//!
//! Unités cohérentes (W, mm, K, s, kg) : le cœur ne convertit rien.

use std::f64::consts::PI;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Analysis, Bc, Model, Step};
use maille_core::results::Increment;
use maille_core::MailleError;

use crate::report;

const CONDUCTIVITY: f64 = 0.05; // W/(mm·K), acier
const DENSITY: f64 = 7.85e-6; // kg/mm³
const SPECIFIC_HEAT: f64 = 460.0; // J/(kg·K)
const WIDTH: f64 = 20.0; // mm, épaisseur L de la plaque (en x)
const HEIGHT: f64 = 10.0; // mm, hauteur modélisée (en y)
const THICKNESS: f64 = 1.0; // mm (profondeur modélisée)
const T_START: f64 = 850.0; // °C, T₀
const T_SURFACE: f64 = 20.0; // °C, Tₛ
const DURATION: f64 = 10.0; // s
const INCREMENTS: usize = 200;
const SIZE: f64 = 1.0; // mm, taille des triangles

fn exact(p: &Node, t: f64) -> f64 {
    let alpha = CONDUCTIVITY / (DENSITY * SPECIFIC_HEAT);
    let sum: f64 = (0..500)
        .map(|i| {
            let n = (2 * i + 1) as f64;
            let k = n * PI / WIDTH;
            4.0 / (n * PI) * (k * p.x).sin() * (-alpha * k * k * t).exp()
        })
        .sum();
    T_SURFACE + (T_START - T_SURFACE) * sum
}

/// Plaque maillée et intégrée en Euler implicite (le départ est discontinu).
pub fn solve() -> Result<(Mesh, Vec<Increment>), MailleError> {
    let geometry = Geometry {
        outer: vec![
            Curve::line([0.0, 0.0], SIZE, "isole"),
            Curve::line([WIDTH, 0.0], SIZE, "trempe"),
            Curve::line([WIDTH, HEIGHT], SIZE, "isole"),
            Curve::line([0.0, HEIGHT], SIZE, "trempe"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())?.to_quadratic();
    let material = Material::conductive(CONDUCTIVITY)?
        .with_density(DENSITY)?
        .with_specific_heat(SPECIFIC_HEAT)?;
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, THICKNESS, PlaneState::Stress).map(Element::Tri6))
        .collect::<Result<_, _>>()?;
    let bcs = t
        .boundary("trempe")?
        .nodes
        .iter()
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: T_SURFACE,
        })
        .collect();
    let mesh = Mesh::new(t.nodes, elements)?;
    let step = Step::transient(DURATION, INCREMENTS)?;
    let history = Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, vec![])?
        .solve_step(&step, &vec![T_START; mesh.nodes.len()])?;
    Ok((mesh, history))
}

/// Affiche la température au cœur face à la solution analytique, puis ouvre le viewer sur le
/// dernier incrément.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, history) = solve()?;

    println!(
        "Trempe d'une plaque : L = {WIDTH} mm, T₀ = {T_START} °C, Tₛ = {T_SURFACE} °C, k = {CONDUCTIVITY} W/(mm·K), ρ = {DENSITY} kg/mm³, c = {SPECIFIC_HEAT} J/(kg·K)"
    );
    println!(
        "Maillée (h = {SIZE} mm), {} Tri6, {} Nodes ; {DURATION} s en {INCREMENTS} incréments d'Euler implicite",
        mesh.elements.len(),
        mesh.nodes.len(),
    );
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    // Node le plus proche du cœur (x = L/2, y = H/2) : la solution exacte y est prise à sa position
    let core = (0..mesh.nodes.len())
        .min_by(|&a, &b| {
            let d =
                |n: usize| (mesh.nodes[n].x - WIDTH / 2.0).hypot(mesh.nodes[n].y - HEIGHT / 2.0);
            d(a).total_cmp(&d(b))
        })
        .unwrap_or(0);
    let increment = DURATION / INCREMENTS as f64;
    for time in [1.0, 2.0, 5.0, 10.0] {
        let inc = &history[(time / increment).round() as usize];
        report(
            &format!("T au cœur, t = {time} s (°C)"),
            inc.field.temperatures[core],
            exact(&mesh.nodes[core], inc.time),
        );
    }
    let last = history.last().ok_or("aucun incrément")?;
    maille_viewer::show(&mesh, &last.field)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temperatures_follow_the_fourier_series_at_the_reported_times() {
        // au tout début, le front est plus raide que h : la comparaison commence à t = 1 s
        // (longueur de diffusion √(αt) ≈ 3,7 mm)
        let (mesh, history) = solve().unwrap();
        for s in [20, 40, 100, 200] {
            let inc = &history[s];
            for (p, t) in mesh.nodes.iter().zip(&inc.field.temperatures) {
                let error = (t - exact(p, inc.time)).abs();
                assert!(
                    error < 1e-2 * (T_START - T_SURFACE),
                    "t = {} : écart {error:.2} °C en ({}, {})",
                    inc.time,
                    p.x,
                    p.y
                );
            }
        }
    }
}
