//! Validation de l'élément poutre d'Euler-Bernoulli contre des solutions analytiques.

use maille_core::material::Material;
use maille_core::mesh::{Bar, Beam, Dof, Element, Mesh, Node};
use maille_core::model::{Bc, DistributedLoad, Load, Model};
use maille_core::results::Field;
use maille_core::MailleError;

const RTOL: f64 = 1e-10;
const E: f64 = 210e3;
const I: f64 = 8.0e4;
const L: f64 = 3.0;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= RTOL * b.abs().max(1e-30)
}

fn beam_mesh(n: usize) -> Mesh {
    let nodes = (0..=n)
        .map(|i| Node::new_2d(L * i as f64 / n as f64, 0.0))
        .collect();
    let m = Material::new(E).unwrap();
    let elements = (0..n)
        .map(|i| Element::Beam(Beam::new([i, i + 1], m, I).unwrap()))
        .collect();
    Mesh::new(nodes, elements).unwrap()
}

fn clamp_left() -> Vec<Bc> {
    vec![
        Bc {
            node: 0,
            dof: Dof::V,
            value: 0.0,
        },
        Bc {
            node: 0,
            dof: Dof::Theta,
            value: 0.0,
        },
    ]
}

fn simple_supports(n: usize) -> Vec<Bc> {
    vec![
        Bc {
            node: 0,
            dof: Dof::V,
            value: 0.0,
        },
        Bc {
            node: n,
            dof: Dof::V,
            value: 0.0,
        },
    ]
}

// DOF globaux : V = 2·node, Theta = 2·node + 1
fn v(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node]
}
fn theta(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node + 1]
}

#[test]
fn cantilever_tip_load() {
    let p = 1500.0;
    for n in [1, 4, 10] {
        let f = Model::new(
            beam_mesh(n),
            clamp_left(),
            vec![Load {
                node: n,
                dof: Dof::V,
                value: p,
            }],
        )
        .unwrap()
        .solve()
        .unwrap();
        assert!(close(v(&f, n), p * L.powi(3) / (3.0 * E * I)), "n = {n}");
        assert!(
            close(theta(&f, n), p * L.powi(2) / (2.0 * E * I)),
            "n = {n}"
        );
        // équilibre : effort tranchant et moment d'encastrement
        assert!(close(f.reactions[0], -p), "n = {n}");
        assert!(close(f.reactions[1], -p * L), "n = {n}");
    }
}

#[test]
fn cantilever_tip_moment() {
    let m = 900.0;
    let n = 3;
    let f = Model::new(
        beam_mesh(n),
        clamp_left(),
        vec![Load {
            node: n,
            dof: Dof::Theta,
            value: m,
        }],
    )
    .unwrap()
    .solve()
    .unwrap();
    assert!(close(v(&f, n), m * L.powi(2) / (2.0 * E * I)));
    assert!(close(theta(&f, n), m * L / (E * I)));
    assert!(close(f.reactions[1], -m));
}

#[test]
fn simply_supported_center_load() {
    let p = 2000.0;
    for n in [2, 6] {
        let f = Model::new(
            beam_mesh(n),
            simple_supports(n),
            vec![Load {
                node: n / 2,
                dof: Dof::V,
                value: p,
            }],
        )
        .unwrap()
        .solve()
        .unwrap();
        assert!(
            close(v(&f, n / 2), p * L.powi(3) / (48.0 * E * I)),
            "n = {n}"
        );
        assert!(
            close(theta(&f, 0), p * L.powi(2) / (16.0 * E * I)),
            "n = {n}"
        );
        assert!(
            close(theta(&f, n), -p * L.powi(2) / (16.0 * E * I)),
            "n = {n}"
        );
        assert!(close(f.reactions[0], -p / 2.0), "n = {n}");
        assert!(close(f.reactions[2 * n], -p / 2.0), "n = {n}");
    }
}

#[test]
fn cantilever_uniform_load() {
    let q = 700.0;
    for n in [1, 5] {
        let f = Model::new(beam_mesh(n), clamp_left(), vec![])
            .unwrap()
            .with_distributed_loads(
                (0..n)
                    .map(|e| DistributedLoad {
                        element: e,
                        value: q,
                    })
                    .collect(),
            )
            .unwrap()
            .solve()
            .unwrap();
        assert!(close(v(&f, n), q * L.powi(4) / (8.0 * E * I)), "n = {n}");
        assert!(
            close(theta(&f, n), q * L.powi(3) / (6.0 * E * I)),
            "n = {n}"
        );
        assert!(close(f.reactions[0], -q * L), "n = {n}");
        assert!(close(f.reactions[1], -q * L * L / 2.0), "n = {n}");
    }
}

#[test]
fn simply_supported_uniform_load() {
    let q = 700.0;
    let n = 4;
    let f = Model::new(beam_mesh(n), simple_supports(n), vec![])
        .unwrap()
        .with_distributed_loads(
            (0..n)
                .map(|e| DistributedLoad {
                    element: e,
                    value: q,
                })
                .collect(),
        )
        .unwrap()
        .solve()
        .unwrap();
    assert!(close(v(&f, n / 2), 5.0 * q * L.powi(4) / (384.0 * E * I)));
}

#[test]
fn beam_without_enough_supports_is_singular() {
    // une seule appui simple laisse la rotation rigide libre
    let result = Model::new(
        beam_mesh(2),
        vec![Bc {
            node: 0,
            dof: Dof::V,
            value: 0.0,
        }],
        vec![Load {
            node: 2,
            dof: Dof::V,
            value: 1.0,
        }],
    )
    .unwrap()
    .solve();
    assert_eq!(result, Err(MailleError::Singular));
}

#[test]
fn invalid_beam_inputs_are_typed_errors() {
    let m = Material::new(E).unwrap();
    assert!(matches!(
        Beam::new([0, 1], m, 0.0),
        Err(MailleError::InvalidProperty { .. })
    ));
    let beam = Element::Beam(Beam::new([1, 0], m, I).unwrap());
    assert_eq!(
        Mesh::new(
            vec![Node::new_2d(0.0, 0.0), Node::new_2d(1.0, 0.0)],
            vec![beam]
        ),
        Err(MailleError::ReversedBeam { element: 0 })
    );
    // Bar et Beam mélangés
    let bar = Element::Bar(Bar {
        nodes: [1, 2],
        material: m,
        section: maille_core::material::Section::new(1.0).unwrap(),
    });
    let beam = Element::Beam(Beam::new([0, 1], m, I).unwrap());
    let nodes = vec![
        Node::new_2d(0.0, 0.0),
        Node::new_2d(1.0, 0.0),
        Node::new_2d(2.0, 0.0),
    ];
    assert_eq!(
        Mesh::new(nodes, vec![beam, bar]),
        Err(MailleError::MixedElementTypes)
    );
    // DOF axial sur une poutre
    assert_eq!(
        Model::new(
            beam_mesh(1),
            vec![Bc {
                node: 0,
                dof: Dof::U,
                value: 0.0
            }],
            vec![]
        ),
        Err(MailleError::InvalidDof {
            node: 0,
            dof: Dof::U
        })
    );
    assert_eq!(
        Model::new(beam_mesh(1), clamp_left(), vec![])
            .unwrap()
            .with_distributed_loads(vec![DistributedLoad {
                element: 3,
                value: 1.0
            }]),
        Err(MailleError::UnknownElement { element: 3 })
    );
}
