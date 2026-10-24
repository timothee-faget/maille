//! Console 2D en contraintes planes, chargée en bout, comparée à la théorie des poutres.

mod common;

use common::{rect_grid, to_quadratic};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, PlaneState, Tri3, Tri6};
use maille_core::model::{Bc, Load, Model};
use maille_core::results::Field;

const E: f64 = 210e3;
const NU: f64 = 0.3;
const L: f64 = 100.0;
const H: f64 = 10.0;
const T: f64 = 5.0;
const P: f64 = 100.0;

/// Flèche de la théorie des poutres de Timoshenko : PL³/3EI + PL/(κGA), κ = 5/6 (section rectangulaire).
fn beam_theory_deflection() -> f64 {
    let (i, a, g) = (T * H.powi(3) / 12.0, T * H, E / (2.0 * (1.0 + NU)));
    P * L.powi(3) / (3.0 * E * i) + P * L / (5.0 / 6.0 * g * a)
}

struct Result {
    mesh: Mesh,
    field: Field,
    /// Flèche du Node central de la section d'extrémité (x = L, y = H/2).
    tip: f64,
}

/// Console nx × ny cases (ny pair), encastrée sur x = 0 (Ux = Uy = 0), effort P réparti sur la section d'extrémité.
fn solve_cantilever(nx: usize, ny: usize, quadratic: bool) -> Result {
    assert!(
        ny.is_multiple_of(2),
        "ny pair pour avoir un Node en y = H/2"
    );
    let (grid_nodes, tris) = rect_grid(nx, ny, L, H);
    let material = Material::with_poisson(E, NU).unwrap();
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let dy = H / ny as f64;

    // forces nodales consistantes d'un effort P uniformément réparti sur la section d'extrémité
    let mut shares: Vec<(usize, f64)> = Vec::new();
    let (nodes, elements): (Vec<Node>, Vec<Element>) = if quadratic {
        let q = to_quadratic(&grid_nodes, &tris);
        for j in 0..ny {
            let (a, b) = (idx(nx, j), idx(nx, j + 1));
            let segment = P * dy / H;
            shares.push((a, segment / 6.0));
            shares.push((q.mid(a, b), segment * 4.0 / 6.0));
            shares.push((b, segment / 6.0));
        }
        let elements = q
            .elements
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material, T, PlaneState::Stress).unwrap()))
            .collect();
        (q.nodes, elements)
    } else {
        for j in 0..=ny {
            let tributary = if j == 0 || j == ny { dy / 2.0 } else { dy };
            shares.push((idx(nx, j), P * tributary / H));
        }
        let elements = tris
            .iter()
            .map(|&n| Element::Tri3(Tri3::new(n, material, T, PlaneState::Stress).unwrap()))
            .collect();
        (grid_nodes, elements)
    };

    let mesh = Mesh::new(nodes.clone(), elements).unwrap();
    let mut bcs = Vec::new();
    for (n, _) in nodes.iter().enumerate().filter(|(_, p)| p.x.abs() < 1e-12) {
        bcs.push(Bc {
            node: n,
            dof: Dof::Ux,
            value: 0.0,
        });
        bcs.push(Bc {
            node: n,
            dof: Dof::Uy,
            value: 0.0,
        });
    }
    let loads = shares
        .into_iter()
        .map(|(node, value)| Load {
            node,
            dof: Dof::Uy,
            value,
        })
        .collect();
    let field = Model::new(mesh.clone(), bcs, loads)
        .unwrap()
        .solve()
        .unwrap();
    let tip_node = nodes
        .iter()
        .position(|p| (p.x - L).abs() < 1e-12 && (p.y - H / 2.0).abs() < 1e-12)
        .unwrap();
    let tip = field.displacements[2 * tip_node + 1];
    Result { mesh, field, tip }
}

/// Écart relatif à la théorie des poutres.
fn error(r: &Result) -> f64 {
    (r.tip / beam_theory_deflection() - 1.0).abs()
}

#[test]
fn tri6_matches_beam_theory_within_two_percent() {
    // L/H = 10 : les variantes de la théorie (κ = 5/6 ou 2/3, gauchissement bloqué à l'encastrement)
    // diffèrent de moins de ~1 % de la flèche, d'où la tolérance de 2 %.
    for (nx, ny) in [(10, 2), (20, 4)] {
        let r = solve_cantilever(nx, ny, true);
        assert!(error(&r) < 0.02, "Tri6 {nx}x{ny} : écart {:.4}", error(&r));
    }
}

#[test]
fn tri3_locks_in_bending_and_converges_from_below() {
    // Le Tri3 (déformation constante) ne sait pas représenter la flexion : il est trop rigide,
    // d'autant plus que le maillage est grossier (verrouillage), et converge par en dessous.
    let ratios: Vec<f64> = [(20, 2), (40, 4), (80, 8)]
        .iter()
        .map(|&(nx, ny)| solve_cantilever(nx, ny, false).tip / beam_theory_deflection())
        .collect();
    assert!(ratios.iter().all(|&r| r < 1.0), "ratios = {ratios:?}");
    assert!(
        ratios.windows(2).all(|w| w[0] < w[1]),
        "ratios = {ratios:?}"
    );
    assert!(
        ratios[0] < 0.8,
        "le verrouillage doit être net sur le maillage grossier : {ratios:?}"
    );
}

#[test]
fn tri6_with_few_dofs_beats_tri3_with_many() {
    let tri6 = solve_cantilever(10, 2, true);
    let tri3 = solve_cantilever(80, 8, false);
    assert!(tri6.field.displacements.len() < tri3.field.displacements.len() / 5);
    assert!(error(&tri6) < error(&tri3));
}

#[test]
fn reactions_balance_force_and_moment_for_both_elements() {
    for quadratic in [false, true] {
        let r = solve_cantilever(20, 2, quadratic);
        let n = r.mesh.nodes.len();
        let fy: f64 = (0..n).map(|i| r.field.reactions[2 * i + 1]).sum();
        // moment des réactions (toutes en x = 0) : −Σ y Rx, qui doit équilibrer P L
        let moment: f64 = (0..n)
            .map(|i| r.mesh.nodes[i].y * r.field.reactions[2 * i])
            .sum();
        assert!(
            (fy + P).abs() < 1e-9 * P,
            "quadratic = {quadratic} : ΣRy = {fy}"
        );
        assert!(
            (moment - P * L).abs() < 1e-9 * P * L,
            "quadratic = {quadratic} : moment = {moment}"
        );
    }
}

#[test]
fn tri6_bending_stress_matches_beam_theory_at_midspan() {
    // Charge vers le haut : courbure concave vers le haut, fibres du bas tendues,
    // σxx = −P (L − x) (y − H/2) / I. Erreur en h² pour des Tri6 de côté H/4 : tolérance 5 % de la
    // contrainte maximale de la section.
    let r = solve_cantilever(20, 4, true);
    let i = T * H.powi(3) / 12.0;
    let mut checked = 0;
    for points in &r.field.stresses {
        for p in points.iter().filter(|p| p.x > 0.3 * L && p.x < 0.7 * L) {
            let sigma_max = P * (L - p.x) * (H / 2.0) / i;
            let expected = -P * (L - p.x) * (p.y - H / 2.0) / i;
            assert!(
                (p.sigma_xx - expected).abs() < 0.05 * sigma_max,
                "x = {}, y = {} : σxx = {}, attendu {}",
                p.x,
                p.y,
                p.sigma_xx,
                expected
            );
            checked += 1;
        }
    }
    assert!(checked > 20, "{checked} points vérifiés");
}
