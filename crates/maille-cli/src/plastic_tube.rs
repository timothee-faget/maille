//! Cas d'exemple : tube épais en acier mis en pression jusqu'à plastifier une partie de son
//! épaisseur, élasto-plasticité de von Mises parfaitement plastique, Tri6 en déformations planes.
//! Le chargement est un déplacement radial δ imposé à l'alésage, augmenté à chaque incrément ; la
//! pression se lit dans les réactions (effort de cerclage du quart de tube, p rᵢ = −Σ R_y sur y = 0).
//!
//! Solution de Hill pour un matériau incompressible, avec k = σy/√3 et c le rayon de la zone
//! plastique :
//! - zone plastique (r < c) : σr = −k (1 − c²/rₑ²) + 2k ln(r/c), σθ = σr + 2k ;
//! - zone élastique (r > c) : σr = −k c²/rₑ² (rₑ²/r² − 1), σθ = k c²/rₑ² (rₑ²/r² + 1) ;
//! - pression p = k (1 − c²/rₑ² + 2 ln(c/rᵢ)), limite 2k ln(rₑ/rᵢ) quand c = rₑ ;
//! - déplacement u = 1,5 k c²/(E r), d'où c pour un déplacement δ de l'alésage.
//!
//! Avant la première plastification, Lamé (élastique). L'écoulement plastique ne change pas le
//! volume ; pour rester proche de l'hypothèse de Hill, ν = 0,49 (et non 0,3 comme l'acier) : il
//! reste un écart de l'ordre de 1 %.
//!
//! Unités cohérentes : N, mm, MPa. Le cœur ne convertit rien.

use std::f64::consts::FRAC_PI_4;

use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Bc, Model, Step};
use maille_core::results::{Increment, Stress};
use maille_core::MailleError;

use crate::report;

const YOUNG: f64 = 210e3; // MPa
const POISSON: f64 = 0.49;
const YIELD: f64 = 250.0; // MPa, σy
const INNER: f64 = 10.0; // mm, rᵢ
const OUTER: f64 = 20.0; // mm, rₑ
const PLASTIC_RADIUS: f64 = 15.0; // mm, rayon c visé en fin de chargement
const SIZE: f64 = 1.0; // mm, taille des triangles
const INCREMENTS: usize = 10;

fn k() -> f64 {
    YIELD / 3.0_f64.sqrt()
}

/// Pression de Hill pour un rayon plastique c.
fn hill_pressure(c: f64) -> f64 {
    k() * (1.0 - c * c / (OUTER * OUTER) + 2.0 * (c / INNER).ln())
}

/// Rayon plastique de Hill pour un déplacement δ de l'alésage (inférieur à rᵢ : encore élastique).
fn hill_radius(bore: f64) -> f64 {
    (bore * YOUNG * INNER / (1.5 * k())).sqrt()
}

/// (σr, σθ) de Hill en r pour un rayon plastique c.
fn hill(r: f64, c: f64) -> [f64; 2] {
    if r <= c {
        let sr = -k() * (1.0 - c * c / (OUTER * OUTER)) + 2.0 * k() * (r / c).ln();
        [sr, sr + 2.0 * k()]
    } else {
        let f = k() * c * c / (OUTER * OUTER);
        [
            -f * (OUTER * OUTER / (r * r) - 1.0),
            f * (OUTER * OUTER / (r * r) + 1.0),
        ]
    }
}

/// Pression de Lamé (élastique) pour un déplacement δ de l'alésage.
fn lame_pressure(bore: f64) -> f64 {
    let span = OUTER * OUTER - INNER * INNER;
    let (a, b) = (INNER * INNER / span, INNER * INNER * OUTER * OUTER / span);
    let unit = (1.0 + POISSON) / YOUNG * ((1.0 - 2.0 * POISSON) * a * INNER + b / INNER);
    bore / unit
}

/// Pression de référence : Lamé tant que le tube est élastique, Hill ensuite.
fn exact_pressure(bore: f64) -> f64 {
    let c = hill_radius(bore);
    if c <= INNER {
        lame_pressure(bore)
    } else {
        hill_pressure(c)
    }
}

/// Déplacement de l'alésage qui plastifie le tube jusqu'au rayon `PLASTIC_RADIUS` (Hill).
fn final_bore() -> f64 {
    1.5 * k() * PLASTIC_RADIUS * PLASTIC_RADIUS / (YOUNG * INNER)
}

/// (σr, σθ) d'un point de Gauss.
fn polar(s: &Stress) -> [f64; 2] {
    let r = s.x.hypot(s.y);
    let (c, sn) = (s.x / r, s.y / r);
    [
        s.sigma_xx * c * c + s.sigma_yy * sn * sn + 2.0 * s.tau_xy * c * sn,
        s.sigma_xx * sn * sn + s.sigma_yy * c * c - 2.0 * s.tau_xy * c * sn,
    ]
}

/// Quart de tube maillé, chargé par un déplacement radial de l'alésage en `INCREMENTS` pas. Renvoie
/// le Mesh, les Increments et les Nodes du plan de symétrie y = 0 (pour lire la pression).
pub fn solve() -> Result<(Mesh, Vec<Increment>, Vec<usize>), MailleError> {
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
    let material = Material::with_poisson(YOUNG, POISSON)?.with_plasticity(YIELD, 0.0)?;
    let elements = t
        .triangles
        .iter()
        .map(|&n| Tri6::new(n, material, 1.0, PlaneState::Strain).map(Element::Tri6))
        .collect::<Result<_, _>>()?;
    let mesh = Mesh::new(t.nodes.clone(), elements)?;

    // symétries en appui glissant, puis déplacement radial δ (x, y)/r à l'alésage sur les DOF que
    // les symétries laissent libres
    let sym_x = t.boundary("sym_x")?.nodes.clone();
    let sym_y = t.boundary("sym_y")?.nodes.clone();
    let mut bcs = Vec::new();
    for (nodes, dof) in [(&sym_x, Dof::Uy), (&sym_y, Dof::Ux)] {
        bcs.extend(nodes.iter().map(|&node| Bc {
            node,
            dof,
            value: 0.0,
        }));
    }
    let bore = final_bore();
    for &node in &t.boundary("interieur")?.nodes {
        let n = mesh.nodes[node];
        let r = n.x.hypot(n.y);
        if !sym_y.contains(&node) {
            bcs.push(Bc {
                node,
                dof: Dof::Ux,
                value: bore * n.x / r,
            });
        }
        if !sym_x.contains(&node) {
            bcs.push(Bc {
                node,
                dof: Dof::Uy,
                value: bore * n.y / r,
            });
        }
    }
    let step = Step::quasi_static(1.0, INCREMENTS)?;
    let history =
        Model::new(mesh.clone(), bcs, vec![])?.solve_quasi_static(&step, &Amplitude::ramp(1.0)?)?;
    Ok((mesh, history, sym_x))
}

/// Pression interne d'un incrément : −Σ R_y sur le plan y = 0, divisé par rᵢ.
fn pressure(increment: &Increment, sym_x: &[usize]) -> f64 {
    -sym_x
        .iter()
        .map(|&n| increment.field.reactions[2 * n + 1])
        .sum::<f64>()
        / INNER
}

/// Rayon du point de Gauss plastifié le plus éloigné de l'axe (0 si aucun).
fn plastic_front(increment: &Increment) -> f64 {
    let field = &increment.field;
    field
        .stresses
        .iter()
        .flatten()
        .zip(field.equivalent_plastic_strains.iter().flatten())
        .filter(|(_, &p)| p > 0.0)
        .map(|(s, _)| s.x.hypot(s.y))
        .fold(0.0, f64::max)
}

/// Affiche la courbe pression-déplacement et les contraintes face à Lamé et Hill, puis ouvre le
/// viewer sur le dernier incrément.
pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let (mesh, history, sym_x) = solve()?;

    println!(
        "Tube épais mis en pression : rᵢ = {INNER} mm, rₑ = {OUTER} mm, E = {YOUNG} MPa, ν = {POISSON}, σy = {YIELD} MPa (parfaitement plastique)"
    );
    println!(
        "Déformations planes ; quart de section maillé (h = {SIZE} mm), {} Tri6, {} Nodes ; déplacement de l'alésage en {INCREMENTS} incréments",
        mesh.elements.len(),
        mesh.nodes.len()
    );
    println!(
        "Pression limite (toute l'épaisseur plastique) : {:.1} MPa",
        hill_pressure(OUTER)
    );
    println!(
        "{:>8}{:>12}{:>14}{:>14}{:>12}{:>12}",
        "λ", "δ (mm)", "p (MPa)", "Lamé/Hill", "c (mm)", "itérations"
    );
    for inc in &history[1..] {
        let bore = inc.time * final_bore();
        let front = plastic_front(inc);
        println!(
            "{:>8.2}{:>12.3e}{:>14.2}{:>14.2}{:>12}{:>12}",
            inc.time,
            bore,
            pressure(inc, &sym_x),
            exact_pressure(bore),
            if front > 0.0 {
                format!("{front:.1}")
            } else {
                "-".to_string()
            },
            inc.residuals.len() - 1
        );
    }

    let last = history.last().ok_or("aucun incrément")?;
    println!();
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "Hill", "écart rel."
    );
    report(
        "pression (MPa)",
        pressure(last, &sym_x),
        hill_pressure(PLASTIC_RADIUS),
    );
    report("rayon plastique (mm)", plastic_front(last), PLASTIC_RADIUS);
    // points de Gauss les plus proches de θ = 45° près de rᵢ, dans la zone plastique et près de rₑ
    let points: Vec<&Stress> = last.field.stresses.iter().flatten().collect();
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
        ("σr près de rᵢ (MPa)", INNER, 0),
        ("σθ près de rᵢ (MPa)", INNER, 1),
        ("σθ en zone plastique (MPa)", 12.5, 1),
        ("σθ près de rₑ (MPa)", OUTER, 1),
    ] {
        let s = nearest(r).ok_or("aucun point de Gauss")?;
        report(label, polar(s)[k], hill(s.x.hypot(s.y), PLASTIC_RADIUS)[k]);
    }
    maille_viewer::show(&mesh, &last.field)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pressure_and_plastic_zone_follow_lame_then_hill() {
        let (_, history, sym_x) = solve().unwrap();
        for inc in &history[1..] {
            let bore = inc.time * final_bore();
            let (p, exact) = (pressure(inc, &sym_x), exact_pressure(bore));
            assert!(
                (p / exact - 1.0).abs() < 1e-2,
                "λ = {} : {p} au lieu de {exact}",
                inc.time
            );
        }
        let front = plastic_front(history.last().unwrap());
        assert!((front - PLASTIC_RADIUS).abs() < SIZE, "c = {front}");
    }
}
