//! Plaque trouée de Kirsch : traction uniaxiale σ (selon x) à l'infini, trou circulaire de rayon a,
//! contraintes planes. On modélise un quart de plaque (symétries) et on impose sur le bord extérieur
//! les déplacements exacts de Kirsch : la plaque finie n'introduit alors aucune erreur de modèle.

mod common;

use common::{rect_grid, to_quadratic};
use maille_core::element::plane_d;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri3, Tri6};
use maille_core::model::{Bc, Model};
use maille_core::results::{nodal_stresses, Field};

const E: f64 = 210e3;
const NU: f64 = 0.3;
const SIGMA: f64 = 100.0;
const A: f64 = 1.0;
const R: f64 = 4.0;

/// Contraintes de Kirsch en coordonnées polaires (σr, σθ, τrθ).
fn kirsch_polar(r: f64, theta: f64) -> [f64; 3] {
    let (q, c, s) = ((A / r).powi(2), (2.0 * theta).cos(), (2.0 * theta).sin());
    [
        SIGMA / 2.0 * (1.0 - q) + SIGMA / 2.0 * (1.0 - 4.0 * q + 3.0 * q * q) * c,
        SIGMA / 2.0 * (1.0 + q) - SIGMA / 2.0 * (1.0 + 3.0 * q * q) * c,
        -SIGMA / 2.0 * (1.0 + 2.0 * q - 3.0 * q * q) * s,
    ]
}

/// Contraintes de Kirsch en coordonnées cartésiennes (σx, σy, τxy) au point (x, y).
fn kirsch_stress(x: f64, y: f64) -> [f64; 3] {
    let (r, t) = (x.hypot(y), y.atan2(x));
    let [sr, st, trt] = kirsch_polar(r, t);
    let (c, s) = (t.cos(), t.sin());
    [
        sr * c * c + st * s * s - 2.0 * trt * s * c,
        sr * s * s + st * c * c + 2.0 * trt * s * c,
        (sr - st) * s * c + trt * (c * c - s * s),
    ]
}

/// Déplacements de Kirsch (ux, uy) en contraintes planes, κ = (3 − ν)/(1 + ν).
fn kirsch_displacement(x: f64, y: f64) -> (f64, f64) {
    let (r, t) = (x.hypot(y), y.atan2(x));
    let kappa = (3.0 - NU) / (1.0 + NU);
    let mu = E / (2.0 * (1.0 + NU));
    let f = SIGMA / (8.0 * mu);
    let ur = f
        * ((kappa - 1.0) * r
            + 2.0 * A * A / r
            + 2.0 * (r + (kappa + 1.0) * A * A / r - A.powi(4) / r.powi(3)) * (2.0 * t).cos());
    let ut = -f * 2.0 * (r + (kappa - 1.0) * A * A / r + A.powi(4) / r.powi(3)) * (2.0 * t).sin();
    (ur * t.cos() - ut * t.sin(), ur * t.sin() + ut * t.cos())
}

#[test]
fn kirsch_stresses_leave_the_hole_free_and_peak_at_three_sigma() {
    for k in 0..=8 {
        let theta = std::f64::consts::FRAC_PI_2 * k as f64 / 8.0;
        let [sr, _, trt] = kirsch_polar(A, theta);
        assert!(
            sr.abs() < 1e-12 * SIGMA && trt.abs() < 1e-12 * SIGMA,
            "θ = {theta}"
        );
    }
    assert!((kirsch_polar(A, std::f64::consts::FRAC_PI_2)[1] - 3.0 * SIGMA).abs() < 1e-12);
    // loin du trou : traction uniaxiale σ selon x
    let far = kirsch_stress(1e4, 0.0);
    assert!((far[0] - SIGMA).abs() < 1e-6 * SIGMA && far[1].abs() < 1e-6 * SIGMA);
}

#[test]
fn kirsch_displacements_derive_to_kirsch_stresses() {
    // Vérification des formules de déplacement avant usage : σ = D ε avec ε par différences
    // finies centrées des déplacements, comparée aux contraintes de Kirsch.
    let d = plane_d(E, NU, PlaneState::Stress);
    let h = 1e-6;
    for (x, y) in [(1.3, 0.2), (0.4, 1.7), (2.5, 2.5), (3.0, 0.1), (0.05, 1.1)] {
        let dx = |c: usize| {
            let (p, m) = (kirsch_displacement(x + h, y), kirsch_displacement(x - h, y));
            [(p.0 - m.0) / (2.0 * h), (p.1 - m.1) / (2.0 * h)][c]
        };
        let dy = |c: usize| {
            let (p, m) = (kirsch_displacement(x, y + h), kirsch_displacement(x, y - h));
            [(p.0 - m.0) / (2.0 * h), (p.1 - m.1) / (2.0 * h)][c]
        };
        let eps = [dx(0), dy(1), dy(0) + dx(1)];
        let sigma: [f64; 3] = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * eps[c]).sum());
        let exact = kirsch_stress(x, y);
        for k in 0..3 {
            assert!(
                (sigma[k] - exact[k]).abs() < 1e-5 * SIGMA,
                "({x}, {y}) composante {k} : {} contre {}",
                sigma[k],
                exact[k]
            );
        }
    }
}

struct Solved {
    mesh: Mesh,
    field: Field,
    dofs: usize,
}

/// Quart de plaque : grille polaire nr × nθ, rayons en progression géométrique de A à R.
fn solve_kirsch(nr: usize, ntheta: usize, quadratic: bool) -> Solved {
    let (grid, tris) = rect_grid(nr, ntheta, 1.0, 1.0);
    let corners: Vec<Node> = grid
        .iter()
        .map(|g| {
            let (r, t) = (A * (R / A).powf(g.x), std::f64::consts::FRAC_PI_2 * g.y);
            Node {
                x: r * t.cos(),
                y: r * t.sin(),
            }
        })
        .collect();
    let material = Material::with_poisson(E, NU).unwrap();
    let idx = |i: usize, j: usize| i + j * (nr + 1);

    let (nodes, elements, mut outer): (Vec<Node>, Vec<Element>, Vec<usize>) = if quadratic {
        let q = to_quadratic(&corners, &tris);
        let mut outer: Vec<usize> = (0..=ntheta).map(|j| idx(nr, j)).collect();
        outer.extend((0..ntheta).map(|j| q.mid(idx(nr, j), idx(nr, j + 1))));
        let elements = q
            .elements
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material, 1.0, PlaneState::Stress).unwrap()))
            .collect();
        (q.nodes, elements, outer)
    } else {
        let elements = tris
            .iter()
            .map(|&n| Element::Tri3(Tri3::new(n, material, 1.0, PlaneState::Stress).unwrap()))
            .collect();
        (
            corners,
            elements,
            (0..=ntheta).map(|j| idx(nr, j)).collect(),
        )
    };
    outer.sort_unstable();

    let mut bcs = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
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
            // symétries : Ux = 0 sur l'axe y (θ = π/2), Uy = 0 sur l'axe x (θ = 0)
            if p.x.abs() < 1e-9 {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Ux,
                    value: 0.0,
                });
            }
            if p.y.abs() < 1e-9 {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Uy,
                    value: 0.0,
                });
            }
        }
    }
    let mesh = Mesh::new(nodes, elements).unwrap();
    let field = Model::new(mesh.clone(), bcs, vec![])
        .unwrap()
        .solve()
        .unwrap();
    let dofs = field.displacements.len();
    Solved { mesh, field, dofs }
}

/// (erreur RMS relative, erreur max / (3σ)) des contraintes aux points de Gauss contre Kirsch.
fn stress_errors(s: &Solved) -> (f64, f64) {
    let (mut diff2, mut exact2, mut max) = (0.0, 0.0, 0.0_f64);
    for p in s.field.stresses.iter().flatten() {
        let exact = kirsch_stress(p.x, p.y);
        let got = [p.sigma_xx, p.sigma_yy, p.tau_xy];
        for k in 0..3 {
            diff2 += (got[k] - exact[k]).powi(2);
            exact2 += exact[k].powi(2);
            max = max.max((got[k] - exact[k]).abs());
        }
    }
    ((diff2 / exact2).sqrt(), max / (3.0 * SIGMA))
}

const TRI3_GRIDS: [usize; 4] = [4, 8, 16, 32];
const TRI6_GRIDS: [usize; 4] = [2, 4, 8, 16];

fn rms_errors(quadratic: bool, grids: [usize; 4]) -> Vec<f64> {
    grids
        .iter()
        .map(|&n| stress_errors(&solve_kirsch(n, n, quadratic)).0)
        .collect()
}

#[test]
fn stress_error_decreases_with_every_refinement() {
    // Le maillage est divisé par 2 à chaque étape. Ordre attendu : stress en O(h) pour le Tri3
    // (rapport ≈ 2), en O(h²) pour le Tri6 (rapport → 4). Seuils à mi-chemin du minimum utile :
    // 1,4 (ordre > 0,5) et 2 (ordre > 1), pour ne pas dépendre de l'erreur de géométrie du trou polygonal.
    for (name, quadratic, grids, min_ratio) in [
        ("Tri3", false, TRI3_GRIDS, 1.4),
        ("Tri6", true, TRI6_GRIDS, 2.0),
    ] {
        let e = rms_errors(quadratic, grids);
        for w in e.windows(2) {
            assert!(w[0] / w[1] > min_ratio, "{name} : erreurs RMS {e:?}");
        }
    }
}

#[test]
fn tri6_reaches_one_percent_and_beats_tri3_with_three_times_fewer_dofs() {
    let tri6 = solve_kirsch(16, 16, true);
    let (rms, max) = stress_errors(&tri6);
    assert!(rms < 0.01, "RMS Tri6 16x16 : {rms}");
    assert!(max < 0.02, "erreur max / 3σ Tri6 16x16 : {max}");

    let tri6_coarse = solve_kirsch(8, 8, true);
    let tri3_fine = solve_kirsch(32, 32, false);
    assert!(tri6_coarse.dofs * 3 < tri3_fine.dofs);
    assert!(stress_errors(&tri6_coarse).0 < stress_errors(&tri3_fine).0);
}

#[test]
fn tri6_finds_the_stress_concentration_at_the_hole() {
    // Pic de σxx sur l'axe y, au bord du trou : Kt = 3 exactement au bord. Les points de Gauss sont à
    // l'intérieur des éléments (r > a) : le maximum calculé est inférieur à 3σ et doit égaler Kirsch
    // au même point.
    let s = solve_kirsch(16, 16, true);
    let peak = s
        .field
        .stresses
        .iter()
        .flatten()
        .max_by(|a, b| a.sigma_xx.total_cmp(&b.sigma_xx))
        .unwrap();
    assert!(
        peak.x < 0.05 * A && (peak.y - A) < 0.05 * A,
        "pic en ({}, {})",
        peak.x,
        peak.y
    );
    let exact = kirsch_stress(peak.x, peak.y)[0];
    assert!(
        (peak.sigma_xx - exact).abs() < 0.01 * 3.0 * SIGMA,
        "{} contre {exact}",
        peak.sigma_xx
    );
    assert!(peak.sigma_xx > 0.95 * 3.0 * SIGMA && peak.sigma_xx < 3.0 * SIGMA);
}

#[test]
fn nodal_smoothing_brings_the_hole_peak_closer_to_three_sigma() {
    // Le Node (0, a) est au bord du trou : σxx lissé y est comparé directement à Kt σ = 3σ,
    // alors que les points de Gauss restent à l'intérieur des Elements.
    let s = solve_kirsch(16, 16, true);
    let nodal = nodal_stresses(&s.mesh, &s.field).unwrap();
    let at_hole = nodal
        .iter()
        .find(|p| p.x.abs() < 1e-9 && (p.y - A).abs() < 1e-9)
        .unwrap();
    let gauss_peak = s
        .field
        .stresses
        .iter()
        .flatten()
        .map(|p| p.sigma_xx)
        .fold(f64::MIN, f64::max);
    let (nodal_gap, gauss_gap) = (
        (at_hole.sigma_xx - 3.0 * SIGMA).abs(),
        (gauss_peak - 3.0 * SIGMA).abs(),
    );
    assert!(
        nodal_gap < gauss_gap,
        "lissé {nodal_gap} contre Gauss {gauss_gap}"
    );
    assert!(
        nodal_gap < 0.01 * 3.0 * SIGMA,
        "écart lissé / 3σ : {}",
        nodal_gap / (3.0 * SIGMA)
    );
}
