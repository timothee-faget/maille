//! Validation de l'élément barre contre des solutions analytiques.

use maille_core::material::{Material, Section};
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node};
use maille_core::model::{Bc, Load, Model};
use maille_core::MailleError;

const RTOL: f64 = 1e-10;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= RTOL * b.abs().max(1.0)
}

fn bar(n1: usize, n2: usize, e: f64, a: f64) -> Element {
    Element::Bar(Bar {
        nodes: [n1, n2],
        material: Material::new(e).unwrap(),
        section: Section::new(a).unwrap(),
    })
}

/// Barre homogène de longueur `l` découpée en `n` éléments.
fn uniform_bar(n: usize, l: f64, e: f64, a: f64) -> Mesh {
    let nodes = (0..=n)
        .map(|i| Node {
            x: l * i as f64 / n as f64,
        })
        .collect();
    let elements = (0..n).map(|i| bar(i, i + 1, e, a)).collect();
    Mesh::new(nodes, elements).unwrap()
}

#[test]
fn fixed_free_bar_end_load_is_pl_over_ea() {
    let (l, e, a, p) = (2.0, 210e3, 50.0, 1000.0);
    // le résultat nodal est exact quel que soit le nombre d'Elements
    for n in [1, 2, 10, 50] {
        let model = Model::new(
            uniform_bar(n, l, e, a),
            vec![Bc {
                dof: Dof::U,
                node: 0,
                value: 0.0,
            }],
            vec![Load {
                dof: Dof::U,
                node: n,
                value: p,
            }],
        )
        .unwrap();
        let field = model.solve().unwrap();
        assert!(close(field.displacements[n], p * l / (e * a)), "n = {n}");
        // équilibre : la réaction s'oppose à la charge
        assert!(close(field.reactions[0], -p), "n = {n}");
    }
}

#[test]
fn two_bars_in_series_add_their_elongations() {
    let (e1, a1, l1) = (200e3, 40.0, 1.5);
    let (e2, a2, l2) = (70e3, 25.0, 0.5);
    let p = 800.0;
    let mesh = Mesh::new(
        vec![Node { x: 0.0 }, Node { x: l1 }, Node { x: l1 + l2 }],
        vec![bar(0, 1, e1, a1), bar(1, 2, e2, a2)],
    )
    .unwrap();
    let field = Model::new(
        mesh,
        vec![Bc {
            dof: Dof::U,
            node: 0,
            value: 0.0,
        }],
        vec![Load {
            dof: Dof::U,
            node: 2,
            value: p,
        }],
    )
    .unwrap()
    .solve()
    .unwrap();
    assert!(close(field.displacements[1], p * l1 / (e1 * a1)));
    assert!(close(
        field.displacements[2],
        p * l1 / (e1 * a1) + p * l2 / (e2 * a2)
    ));
}

#[test]
fn bar_fixed_at_both_ends_splits_the_reactions() {
    // charge P au nœud intermédiaire : R_gauche = −P·L2/L, R_droite = −P·L1/L
    let (e, a, l1, l2, p) = (100e3, 10.0, 1.0, 3.0, 400.0);
    let mesh = Mesh::new(
        vec![Node { x: 0.0 }, Node { x: l1 }, Node { x: l1 + l2 }],
        vec![bar(0, 1, e, a), bar(1, 2, e, a)],
    )
    .unwrap();
    let field = Model::new(
        mesh,
        vec![
            Bc {
                dof: Dof::U,
                node: 0,
                value: 0.0,
            },
            Bc {
                dof: Dof::U,
                node: 2,
                value: 0.0,
            },
        ],
        vec![Load {
            dof: Dof::U,
            node: 1,
            value: p,
        }],
    )
    .unwrap()
    .solve()
    .unwrap();
    let l = l1 + l2;
    assert!(close(field.reactions[0], -p * l2 / l));
    assert!(close(field.reactions[2], -p * l1 / l));
    assert!(close(field.reactions[0] + field.reactions[2], -p));
}

#[test]
fn imposed_displacement_gives_uniform_strain() {
    let (l, e, a, d) = (4.0, 210e3, 10.0, 0.01);
    let field = Model::new(
        uniform_bar(4, l, e, a),
        vec![
            Bc {
                dof: Dof::U,
                node: 0,
                value: 0.0,
            },
            Bc {
                dof: Dof::U,
                node: 4,
                value: d,
            },
        ],
        vec![],
    )
    .unwrap()
    .solve()
    .unwrap();
    assert!(close(field.displacements[2], d / 2.0));
    assert!(close(field.reactions[4], e * a * d / l));
}

#[test]
fn unconstrained_bar_is_reported_as_singular() {
    let result = Model::new(
        uniform_bar(2, 1.0, 1.0, 1.0),
        vec![],
        vec![Load {
            dof: Dof::U,
            node: 2,
            value: 1.0,
        }],
    )
    .unwrap()
    .solve();
    assert_eq!(result, Err(MailleError::Singular));
}

#[test]
fn invalid_inputs_are_typed_errors() {
    assert!(matches!(
        Material::new(0.0),
        Err(MailleError::InvalidProperty { .. })
    ));
    assert!(matches!(
        Section::new(f64::NAN),
        Err(MailleError::InvalidProperty { .. })
    ));
    assert_eq!(Mesh::new(vec![], vec![]), Err(MailleError::EmptyMesh));
    assert_eq!(
        Mesh::new(
            vec![Node { x: 0.0 }, Node { x: 0.0 }],
            vec![bar(0, 1, 1.0, 1.0)]
        ),
        Err(MailleError::ZeroLengthElement { element: 0 })
    );
    assert_eq!(
        Mesh::new(vec![Node { x: 0.0 }], vec![bar(0, 1, 1.0, 1.0)]),
        Err(MailleError::UnknownNode { node: 1 })
    );
    let mesh = uniform_bar(1, 1.0, 1.0, 1.0);
    assert_eq!(
        Model::new(
            mesh.clone(),
            vec![Bc {
                dof: Dof::U,
                node: 5,
                value: 0.0
            }],
            vec![]
        ),
        Err(MailleError::UnknownNode { node: 5 })
    );
    assert_eq!(
        Model::new(
            mesh,
            vec![
                Bc {
                    dof: Dof::U,
                    node: 0,
                    value: 0.0
                },
                Bc {
                    dof: Dof::U,
                    node: 0,
                    value: 1.0
                }
            ],
            vec![]
        ),
        Err(MailleError::DuplicateBc {
            dof: Dof::U,
            node: 0
        })
    );
}
