//! Cas d'exemple : contraintes thermiques dans un tube épais, couplage thermo-mécanique séquentiel,
//! Tri6. BC de température Tᵢ sur le rayon intérieur rᵢ et Tₑ sur le rayon extérieur rₑ : on
//! calcule d'abord le champ T (profil en ln r), puis la dilatation en déformations planes (tube long,
//! εz = 0), surfaces libres. Solution de Timoshenko et Goodier, avec ΔT = T − T_ref et
//! I(r) = ∫ᵣᵢʳ ΔT r dr :
//! σr = αE / ((1 − ν) r²) [(r² − rᵢ²)/(rₑ² − rᵢ²) I(rₑ) − I(r)],
//! σθ = αE / ((1 − ν) r²) [(r² + rᵢ²)/(rₑ² − rᵢ²) I(rₑ) + I(r) − ΔT r²],
//! σz = ν (σr + σθ) − E α ΔT.
//! On modélise un quart de section, maillé par le mailleur, avec des symétries en appui glissant.
//! Les Tri6 ont des côtés droits : près des rayons, la surface libre est celle des cordes, ce qui
//! dévie localement la contrainte d'environ 1 % (l'erreur moyenne est de l'ordre de 0,1 %).
//!
//! Unités cohérentes, une par calcul : W, mm, K en thermique ; N, mm, MPa en mécanique, α en 1/K.
//! Le cœur ne convertit rien.

use std::f64::consts::FRAC_PI_4;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Analysis, Bc, Model};
use maille_core::results::{Field, Stress};
use maille_core::MailleError;

use crate::report;

const YOUNG: f64 = 210e3; // MPa
const POISSON: f64 = 0.3;
const EXPANSION: f64 = 1.2e-5; // 1/K
const CONDUCTIVITY: f64 = 0.05; // W/(mm·K)
const INNER: f64 = 10.0; // mm, rᵢ
const OUTER: f64 = 20.0; // mm, rₑ
const T_INNER: f64 = 120.0; // °C
const T_OUTER: f64 = 20.0; // °C
const T_REF: f64 = 20.0; // °C, température sans contrainte
const SIZE: f64 = 1.0; // mm, taille des triangles

fn beta() -> f64 {
    (T_OUTER - T_INNER) / (OUTER / INNER).ln()
}

/// ΔT(r) = T(r) − T_ref.
fn delta(r: f64) -> f64 {
    T_INNER - T_REF + beta() * (r / INNER).ln()
}

/// I(r) = ∫ᵣᵢʳ ΔT r dr, avec ∫ ln(r/rᵢ) r dr = r²/2 ln(r/rᵢ) − r²/4.
fn integral(r: f64) -> f64 {
    (T_INNER - T_REF) * (r * r - INNER * INNER) / 2.0
        + beta() * (r * r / 2.0 * (r / INNER).ln() - r * r / 4.0 + INNER * INNER / 4.0)
}

/// (σr, σθ, σz) exactes en r.
fn exact(r: f64) -> [f64; 3] {
    let c = EXPANSION * YOUNG / ((1.0 - POISSON) * r * r);
    let span = OUTER * OUTER - INNER * INNER;
    let sr = c * ((r * r - INNER * INNER) / span * integral(OUTER) - integral(r));
    let st =
        c * ((r * r + INNER * INNER) / span * integral(OUTER) + integral(r) - delta(r) * r * r);
    [sr, st, POISSON * (sr + st) - YOUNG * EXPANSION * delta(r)]
}

/// (σr, σθ, σz) d'un point de Gauss.
fn polar(s: &Stress) -> [f64; 3] {
    let r = s.x.hypot(s.y);
    let (c, sn) = (s.x / r, s.y / r);
    [
        s.sigma_xx * c * c + s.sigma_yy * sn * sn + 2.0 * s.tau_xy * c * sn,
        s.sigma_xx * sn * sn + s.sigma_yy * c * c - 2.0 * s.tau_xy * c * sn,
        s.sigma_zz,
    ]
}

/// Quart de tube maillé : calcul thermique, puis statique sur le même Mesh. Renvoie le Mesh et le
/// Field mécanique.
pub fn solve() -> Result<(Mesh, Field), MailleError> {
    let geometry = Geometry {
        outer: vec![
            Curve::line([INNER, 0.0], SIZE, "sym_x"),
            Curve::arc([OUTER, 0.0], [0.0, 0.0], false, SIZE, "exterieur"),
            Curve::line([0.0, OUTER], SIZE, "sym_y"),
            Curve::arc([0.0, INNER], [0.0, 0.0], true, SIZE, "interieur"),
        ],
        holes: vec![],
    };
    let t = triangulate(&geometry, &MeshOptions::default())?.to_quadratic();
    let material = Material::with_poisson(YOUNG, POISSON)?
        .with_expansion(EXPANSION)?
        .with_conductivity(CONDUCTIVITY)?;
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, 1.0, PlaneState::Strain).map(Element::Tri6))
        .collect::<Result<_, _>>()?;
    let mesh = Mesh::new(t.nodes.clone(), elements)?;

    // thermique : solution exacte imposée sur les deux rayons (Nodes pas tous sur le cercle)
    let bcs = t
        .boundary("interieur")?
        .nodes
        .iter()
        .chain(&t.boundary("exterieur")?.nodes)
        .map(|&node| Bc {
            node,
            dof: Dof::T,
            value: T_REF + delta(t.nodes[node].x.hypot(t.nodes[node].y)),
        })
        .collect();
    let thermal = Model::with_analysis(mesh.clone(), Analysis::Thermal, bcs, vec![])?.solve()?;

    // mécanique : symétries en appui glissant, surfaces libres, températures du calcul thermique
    let mut bcs = Vec::new();
    for (name, dof) in [("sym_x", Dof::Uy), ("sym_y", Dof::Ux)] {
        bcs.extend(t.boundary(name)?.nodes.iter().map(|&node| Bc {
            node,
            dof,
            value: 0.0,
        }));
    }
    let mechanical = Model::new(mesh.clone(), bcs, vec![])?
        .with_temperatures(thermal.temperatures, T_REF)?
        .solve()?;
    Ok((mesh, mechanical))
}

/// Affiche les contraintes face à la solution analytique, puis ouvre le viewer sur le Field
/// mécanique.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, mechanical) = solve()?;

    println!(
        "Tube épais chauffé de l'intérieur : rᵢ = {INNER} mm, rₑ = {OUTER} mm, Tᵢ = {T_INNER} °C, Tₑ = {T_OUTER} °C, T_ref = {T_REF} °C"
    );
    println!(
        "E = {YOUNG} MPa, ν = {POISSON}, α = {EXPANSION} /K, déformations planes ; quart de section maillé (h = {SIZE} mm), {} Tri6, {} Nodes",
        mesh.elements.len(),
        mesh.nodes.len()
    );
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    // points de Gauss les plus proches de θ = 45° près de rᵢ, au milieu et près de rₑ
    let points: Vec<&Stress> = mechanical.stresses.iter().flatten().collect();
    let nearest = |r: f64| {
        let target = [r * FRAC_PI_4.cos(), r * FRAC_PI_4.sin()];
        points
            .iter()
            .min_by(|a, b| {
                let d = |s: &Stress| (s.x - target[0]).hypot(s.y - target[1]);
                d(a).total_cmp(&d(b))
            })
            .copied()
    };
    for (label, r, k) in [
        ("σθ près de rᵢ (MPa)", INNER, 1),
        ("σz près de rᵢ (MPa)", INNER, 2),
        ("σr au milieu (MPa)", (INNER + OUTER) / 2.0, 0),
        ("σθ près de rₑ (MPa)", OUTER, 1),
    ] {
        let s = nearest(r).ok_or("aucun point de Gauss")?;
        report(label, polar(s)[k], exact(s.x.hypot(s.y))[k]);
    }
    maille_viewer::show(&mesh, &mechanical)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stresses_follow_timoshenko() {
        // erreur RMS sur (σr, σθ, σz) aux points de Gauss sous 0,2 % du plus grand |σθ| ; point par
        // point sous 2 % : près des rayons, le bord en cordes (côtés droits) dévie localement
        let (_, mechanical) = solve().unwrap();
        let scale = exact(INNER)[1].abs();
        let (mut sum, mut count) = (0.0, 0);
        for s in mechanical.stresses.iter().flatten() {
            let (got, want) = (polar(s), exact(s.x.hypot(s.y)));
            for k in 0..3 {
                let error = (got[k] - want[k]).abs();
                assert!(error < 2e-2 * scale, "({}, {}) composante {k}", s.x, s.y);
                sum += error * error;
                count += 1;
            }
        }
        let rms = (sum / count as f64).sqrt() / scale;
        assert!(rms < 2e-3, "erreur RMS {rms:.2e}");
    }
}
