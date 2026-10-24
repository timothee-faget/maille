//! Élasticité 3D, Tet4 et Tet10 : patch test, traction simple, console comparée à la théorie des poutres.

mod common;

use common::{box_grid, to_quadratic_tet};
use maille_core::element::solid_d;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, Node, Tet10, Tet4};
use maille_core::model::{Bc, Load, Model};
use maille_core::results::{nodal_stresses, Field, Stress};
use maille_core::MailleError;

const E: f64 = 210e3;
const NU: f64 = 0.3;
const DOFS: [Dof; 3] = [Dof::Ux, Dof::Uy, Dof::Uz];

fn material() -> Material {
    Material::with_poisson(E, NU).unwrap()
}

/// Mesh du pavé en Tet4, ou en Tet10 si `quadratic`.
fn box_mesh(n: [usize; 3], size: [f64; 3], quadratic: bool) -> Mesh {
    let (nodes, tets) = box_grid(n, size);
    if quadratic {
        let q = to_quadratic_tet(&nodes, &tets);
        let elements = q
            .elements
            .iter()
            .map(|&t| Element::Tet10(Tet10::new(t, material())))
            .collect();
        Mesh::new(q.nodes, elements).unwrap()
    } else {
        let elements = tets
            .iter()
            .map(|&t| Element::Tet4(Tet4::new(t, material())))
            .collect();
        Mesh::new(nodes, elements).unwrap()
    }
}

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-9 * scale
}

/// Composantes de Voigt (xx, yy, zz, yz, xz, xy) d'une contrainte.
fn voigt(s: &Stress) -> [f64; 6] {
    [
        s.sigma_xx, s.sigma_yy, s.sigma_zz, s.tau_yz, s.tau_xz, s.tau_xy,
    ]
}

#[test]
fn uniform_strain_patch_is_exact_for_tet4_and_tet10() {
    // déplacement linéaire u = ε·x imposé au bord : les Nodes intérieurs et toutes les contraintes
    // doivent être exacts (σ = D ε partout)
    let eps = [1e-3, -2e-3, 0.5e-3, 3e-3, -1e-3, 2e-3];
    let exact = |p: &Node| {
        [
            eps[0] * p.x + eps[5] / 2.0 * p.y + eps[4] / 2.0 * p.z,
            eps[5] / 2.0 * p.x + eps[1] * p.y + eps[3] / 2.0 * p.z,
            eps[4] / 2.0 * p.x + eps[3] / 2.0 * p.y + eps[2] * p.z,
        ]
    };
    let size = [2.0, 1.0, 1.5];
    let d = solid_d(E, NU);
    let sigma: [f64; 6] = std::array::from_fn(|r| (0..6).map(|c| d[r][c] * eps[c]).sum());
    let scale = sigma.iter().fold(0.0_f64, |m, s| m.max(s.abs()));

    for quadratic in [false, true] {
        let mesh = box_mesh([3, 2, 2], size, quadratic);
        let on_boundary = |p: &Node| {
            [p.x, p.y, p.z]
                .iter()
                .zip(size)
                .any(|(&c, l)| c.abs() < 1e-12 || (c - l).abs() < 1e-12)
        };
        let mut bcs = Vec::new();
        for (n, p) in mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| on_boundary(p))
        {
            for (dof, value) in DOFS.into_iter().zip(exact(p)) {
                bcs.push(Bc {
                    node: n,
                    dof,
                    value,
                });
            }
        }
        let interior = mesh.nodes.len() - bcs.len() / 3;
        assert!(interior > 0, "le patch doit avoir des Nodes libres");
        let field = Model::new(mesh.clone(), bcs, Vec::new())
            .unwrap()
            .solve()
            .unwrap();

        for (n, p) in mesh.nodes.iter().enumerate() {
            for (c, e) in exact(p).into_iter().enumerate() {
                let u = field.displacements[3 * n + c];
                assert!((u - e).abs() < 1e-12, "Node {n} : {u} au lieu de {e}");
            }
        }
        let gauss = field.stresses.iter().flatten();
        for s in gauss.chain(&nodal_stresses(&mesh, &field).unwrap()) {
            for (g, e) in voigt(s).iter().zip(sigma) {
                assert!(
                    close(*g, e, scale),
                    "quadratique {quadratic} : {:?}",
                    voigt(s)
                );
            }
        }
    }
}

#[test]
fn uniaxial_traction_gives_e_delta_over_l_and_poisson_contraction() {
    // barre L × 1 × 1, plans de symétrie x = 0, y = 0, z = 0, allongement δ imposé en x = L
    let (l, delta) = (4.0, 0.01);
    for quadratic in [false, true] {
        let mesh = box_mesh([4, 1, 1], [l, 1.0, 1.0], quadratic);
        let mut bcs = Vec::new();
        for (n, p) in mesh.nodes.iter().enumerate() {
            for (dof, c) in DOFS.into_iter().zip([p.x, p.y, p.z]) {
                if c.abs() < 1e-12 {
                    bcs.push(Bc {
                        node: n,
                        dof,
                        value: 0.0,
                    });
                }
            }
            if (p.x - l).abs() < 1e-12 {
                bcs.push(Bc {
                    node: n,
                    dof: Dof::Ux,
                    value: delta,
                });
            }
        }
        let field = Model::new(mesh.clone(), bcs, Vec::new())
            .unwrap()
            .solve()
            .unwrap();

        let strain = delta / l;
        for s in field.stresses.iter().flatten() {
            assert!(
                close(s.sigma_xx, E * strain, E * strain),
                "σxx = {}",
                s.sigma_xx
            );
            assert!(s.sigma_yy.abs() < 1e-9 * E * strain);
        }
        for (n, p) in mesh.nodes.iter().enumerate() {
            let (uy, uz) = (
                field.displacements[3 * n + 1],
                field.displacements[3 * n + 2],
            );
            assert!((uy + NU * strain * p.y).abs() < 1e-12, "Uy du Node {n}");
            assert!((uz + NU * strain * p.z).abs() < 1e-12, "Uz du Node {n}");
        }
        // somme des réactions sur la face tirée = σ A = E δ / L (section unité)
        let pull: f64 = mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| (p.x - l).abs() < 1e-12)
            .map(|(n, _)| field.reactions[3 * n])
            .sum();
        assert!(close(pull, E * strain, E * strain), "réaction {pull}");
    }
}

const L: f64 = 10.0;
const B: f64 = 1.0;
const H: f64 = 1.0;
const P: f64 = 1.0;

/// Flèche de Timoshenko : PL³/3EI + PL/(κGA), κ = 5/6 (section rectangulaire).
fn beam_theory_deflection() -> f64 {
    let (i, a, g) = (B * H.powi(3) / 12.0, B * H, E / (2.0 * (1.0 + NU)));
    P * L.powi(3) / (3.0 * E * i) + P * L / (5.0 / 6.0 * g * a)
}

/// Faces d'un tétraèdre (indices de sommets) et, pour chacune, ses trois arêtes (indices des
/// Nodes milieux du Tet10).
const FACES: [([usize; 3], [usize; 3]); 4] = [
    ([0, 1, 2], [4, 5, 6]),
    ([0, 1, 3], [4, 8, 7]),
    ([1, 2, 3], [5, 9, 8]),
    ([0, 2, 3], [6, 9, 7]),
];

/// Console L × B × H encastrée en x = 0, effort P (selon z) réparti uniformément sur la section
/// x = L. Forces nodales consistantes : sur une face triangulaire d'aire A, A/3 par sommet pour le
/// Tet4 ; A/3 par milieu d'arête et 0 aux sommets pour le Tet10. Renvoie la flèche Uz du Node
/// au centre de la section d'extrémité (ny et nz pairs).
fn cantilever_tip(n: [usize; 3], quadratic: bool) -> f64 {
    let mesh = box_mesh(n, [L, B, H], quadratic);
    let at_tip = |k: usize| (mesh.nodes[k].x - L).abs() < 1e-12;
    let mut loads = Vec::new();
    for e in &mesh.elements {
        let nodes = e.nodes();
        for (vertices, mids) in FACES {
            if !vertices.iter().all(|&v| at_tip(nodes[v])) {
                continue;
            }
            let p: Vec<Node> = vertices.iter().map(|&v| mesh.nodes[nodes[v]]).collect();
            // la face est dans le plan x = L : aire = |Δy₁ Δz₂ − Δz₁ Δy₂| / 2
            let area = ((p[1].y - p[0].y) * (p[2].z - p[0].z)
                - (p[1].z - p[0].z) * (p[2].y - p[0].y))
                .abs()
                / 2.0;
            let share = P * area / (B * H) / 3.0;
            let carriers = if quadratic { mids } else { vertices };
            for k in carriers {
                loads.push(Load {
                    node: nodes[k],
                    dof: Dof::Uz,
                    value: share,
                });
            }
        }
    }
    let total: f64 = loads.iter().map(|l| l.value).sum();
    assert!((total - P).abs() < 1e-12, "effort total {total}");

    let mut bcs = Vec::new();
    for (k, _) in mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, p)| p.x.abs() < 1e-12)
    {
        for dof in DOFS {
            bcs.push(Bc {
                node: k,
                dof,
                value: 0.0,
            });
        }
    }
    let tip = mesh
        .nodes
        .iter()
        .position(|p| {
            (p.x - L).abs() < 1e-12
                && (p.y - B / 2.0).abs() < 1e-12
                && (p.z - H / 2.0).abs() < 1e-12
        })
        .unwrap();
    let field: Field = Model::new(mesh, bcs, loads).unwrap().solve().unwrap();
    field.displacements[3 * tip + 2]
}

#[test]
fn tet10_cantilever_matches_beam_theory() {
    // L/H = 10 : comme en 2D, les variantes de la théorie (κ, gauchissement bloqué à
    // l'encastrement) diffèrent d'environ 1 %, d'où la tolérance de 2 %.
    for n in [[10, 2, 2], [20, 2, 2]] {
        let ratio = cantilever_tip(n, true) / beam_theory_deflection();
        assert!(
            (ratio - 1.0).abs() < 0.02,
            "Tet10 {n:?} : rapport {ratio:.4}"
        );
    }
}

#[test]
fn tet4_cantilever_locks_and_converges_slowly() {
    // déformation constante par Element : la flexion verrouille, la flèche reste très en dessous
    let coarse = cantilever_tip([10, 2, 2], false) / beam_theory_deflection();
    let fine = cantilever_tip([20, 4, 4], false) / beam_theory_deflection();
    assert!(coarse < 0.6, "Tet4 grossier : rapport {coarse:.4}");
    assert!(fine > coarse && fine < 0.95, "Tet4 fin : rapport {fine:.4}");
}

#[test]
fn bad_tetrahedra_are_refused_with_their_index() {
    let nodes = vec![
        Node::new(0.0, 0.0, 0.0),
        Node::new(1.0, 0.0, 0.0),
        Node::new(0.0, 1.0, 0.0),
        Node::new(0.0, 0.0, 1.0),
        Node::new(1.0, 1.0, 0.0),
    ];
    let tet = |n| Element::Tet4(Tet4::new(n, material()));
    assert_eq!(
        Mesh::new(nodes.clone(), vec![tet([0, 1, 2, 3]), tet([0, 2, 1, 3])]),
        Err(MailleError::InvertedTetrahedron { element: 1 })
    );
    assert_eq!(
        Mesh::new(nodes.clone(), vec![tet([0, 1, 2, 4])]),
        Err(MailleError::DegenerateTetrahedron { element: 0 })
    );

    // Tet10 dont le milieu de l'arête 0-1 est déplacé
    let mut q = to_quadratic_tet(&nodes[..4], &[[0, 1, 2, 3]]);
    let m01 = q.mid(0, 1);
    q.nodes[m01].y += 0.1;
    let tet10 = Element::Tet10(Tet10::new(q.elements[0], material()));
    assert_eq!(
        Mesh::new(q.nodes, vec![tet10]),
        Err(MailleError::CurvedTetrahedron { element: 0 })
    );
}

#[test]
fn tetrahedra_have_three_dofs_and_no_rotation() {
    let mesh = box_mesh([1, 1, 1], [1.0, 1.0, 1.0], false);
    assert_eq!(mesh.dofs_per_node(), 3);
    assert_eq!(mesh.dof_index(2, Dof::Uz), Ok(8));
    assert_eq!(
        mesh.dof_index(0, Dof::Theta),
        Err(MailleError::InvalidDof {
            node: 0,
            dof: Dof::Theta
        })
    );
}

#[test]
fn tetrahedron_without_poisson_is_refused_by_the_static_model() {
    let nodes = vec![
        Node::new(0.0, 0.0, 0.0),
        Node::new(1.0, 0.0, 0.0),
        Node::new(0.0, 1.0, 0.0),
        Node::new(0.0, 0.0, 1.0),
    ];
    let tet = Tet4::new([0, 1, 2, 3], Material::new(E).unwrap());
    let mesh = Mesh::new(nodes, vec![Element::Tet4(tet)]).unwrap();
    assert_eq!(
        Model::new(mesh, Vec::new(), Vec::new()),
        Err(MailleError::MissingPoisson { element: 0 })
    );
}
