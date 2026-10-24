//! Validation de l'élément Tri3 : patch test et traction uniaxiale.

use maille_core::element::plane_d;
use maille_core::material::Material;
use maille_core::mesh::{Bar, Dof, Element, Mesh, Node, PlaneState, Tri3};
use maille_core::model::{Bc, DistributedLoad, Load, Model};
use maille_core::results::{nodal_stresses, Field};
use maille_core::MailleError;

const RTOL: f64 = 1e-10;

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= RTOL * scale
}

fn tri(nodes: [usize; 3], nu: f64, state: PlaneState) -> Element {
    Element::Tri3(
        Tri3::new(
            nodes,
            Material::with_poisson(210e3, nu).unwrap(),
            2.0,
            state,
        )
        .unwrap(),
    )
}

fn ux(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node]
}
fn uy(f: &Field, node: usize) -> f64 {
    f.displacements[2 * node + 1]
}

/// Rectangle w × h de nx × ny cases, chaque case coupée en 2 triangles antihoraires.
/// Le Node (i, j) a pour indice i + j (nx + 1).
fn rect_mesh(nx: usize, ny: usize, w: f64, h: f64, nu: f64, state: PlaneState) -> Mesh {
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push(Node {
                x: w * i as f64 / nx as f64,
                y: h * j as f64 / ny as f64,
            });
        }
    }
    let mut elements = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            elements.push(tri(
                [idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)],
                nu,
                state,
            ));
            elements.push(tri(
                [idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)],
                nu,
                state,
            ));
        }
    }
    Mesh::new(nodes, elements).unwrap()
}

/// Traction uniaxiale σ sur le bord droit : symétries à gauche (Ux = 0) et en bas (Uy = 0).
fn uniaxial_tension(
    nx: usize,
    ny: usize,
    w: f64,
    h: f64,
    nu: f64,
    state: PlaneState,
) -> (Field, f64) {
    let sigma = 100.0;
    let thickness = 2.0;
    let mesh = rect_mesh(nx, ny, w, h, nu, state);
    let idx = |i: usize, j: usize| i + j * (nx + 1);
    let mut bcs = Vec::new();
    for j in 0..=ny {
        bcs.push(Bc {
            node: idx(0, j),
            dof: Dof::Ux,
            value: 0.0,
        });
    }
    for i in 0..=nx {
        bcs.push(Bc {
            node: idx(i, 0),
            dof: Dof::Uy,
            value: 0.0,
        });
    }
    // forces nodales consistantes d'une traction uniforme : longueur tributaire de chaque Node du bord
    let dy = h / ny as f64;
    let loads = (0..=ny)
        .map(|j| {
            let tributary = if j == 0 || j == ny { dy / 2.0 } else { dy };
            Load {
                node: idx(nx, j),
                dof: Dof::Ux,
                value: sigma * thickness * tributary,
            }
        })
        .collect();
    (
        Model::new(mesh, bcs, loads).unwrap().solve().unwrap(),
        sigma,
    )
}

#[test]
fn uniaxial_tension_plane_stress_is_exact() {
    let (nx, ny, w, h, nu) = (3, 2, 6.0, 4.0, 0.3);
    let (f, sigma) = uniaxial_tension(nx, ny, w, h, nu, PlaneState::Stress);
    let e = 210e3;
    for j in 0..=ny {
        for i in 0..=nx {
            let n = i + j * (nx + 1);
            let (x, y) = (w * i as f64 / nx as f64, h * j as f64 / ny as f64);
            assert!(close(ux(&f, n), sigma * x / e, sigma * w / e), "ux {i},{j}");
            assert!(
                close(uy(&f, n), -nu * sigma * y / e, sigma * h / e),
                "uy {i},{j}"
            );
        }
    }
}

#[test]
fn uniaxial_tension_plane_strain_uses_the_reduced_moduli() {
    // déformations planes = contraintes planes avec E' = E/(1−ν²) et ν' = ν/(1−ν)
    let (nx, ny, w, h, nu) = (2, 3, 5.0, 3.0, 0.3);
    let (f, sigma) = uniaxial_tension(nx, ny, w, h, nu, PlaneState::Strain);
    let e = 210e3;
    let (e_eq, nu_eq) = (e / (1.0 - nu * nu), nu / (1.0 - nu));
    for j in 0..=ny {
        for i in 0..=nx {
            let n = i + j * (nx + 1);
            let (x, y) = (w * i as f64 / nx as f64, h * j as f64 / ny as f64);
            assert!(
                close(ux(&f, n), sigma * x / e_eq, sigma * w / e_eq),
                "ux {i},{j}"
            );
            assert!(
                close(uy(&f, n), -nu_eq * sigma * y / e_eq, sigma * h / e_eq),
                "uy {i},{j}"
            );
        }
    }
}

#[test]
fn reactions_balance_the_applied_traction() {
    let (nx, ny, w, h) = (3, 3, 6.0, 4.0);
    let (f, sigma) = uniaxial_tension(nx, ny, w, h, 0.3, PlaneState::Stress);
    let total: f64 = (0..=ny).map(|j| f.reactions[2 * (j * (nx + 1))]).sum();
    assert!(close(total, -sigma * 2.0 * h, sigma * 2.0 * h));
}

/// Maillage irrégulier de [0,2]² avec 2 Nodes intérieurs (4 et 5), 6 triangles antihoraires.
fn patch_mesh(state: PlaneState) -> Mesh {
    let nodes = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (0.7, 0.6),
        (1.3, 1.2),
    ]
    .map(|(x, y)| Node { x, y })
    .to_vec();
    let elements = [
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 5],
        [3, 0, 4],
        [4, 1, 5],
        [4, 5, 3],
    ]
    .map(|n| tri(n, 0.3, state))
    .to_vec();
    Mesh::new(nodes, elements).unwrap()
}

#[test]
fn patch_test_linear_displacement_field_is_reproduced_exactly() {
    // champ linéaire quelconque (les trois déformations sont non nulles), imposé aux 4 coins
    let (a, b) = ([1e-3, 2e-3, -1.5e-3], [-5e-4, 3e-3, 2.5e-3]);
    let field = |x: f64, y: f64| (a[0] + a[1] * x + a[2] * y, b[0] + b[1] * x + b[2] * y);
    for state in [PlaneState::Stress, PlaneState::Strain] {
        let mesh = patch_mesh(state);
        let bcs: Vec<Bc> = (0..4)
            .flat_map(|n| {
                let (u, v) = field(mesh.nodes[n].x, mesh.nodes[n].y);
                [
                    Bc {
                        node: n,
                        dof: Dof::Ux,
                        value: u,
                    },
                    Bc {
                        node: n,
                        dof: Dof::Uy,
                        value: v,
                    },
                ]
            })
            .collect();
        let nodes = mesh.nodes.clone();
        let f = Model::new(mesh.clone(), bcs, vec![])
            .unwrap()
            .solve()
            .unwrap();
        for n in [4, 5] {
            let (u, v) = field(nodes[n].x, nodes[n].y);
            assert!(close(ux(&f, n), u, 1e-3), "{state:?} ux {n}");
            assert!(close(uy(&f, n), v, 1e-3), "{state:?} uy {n}");
        }
        // patch test fort : contrainte constante σ = D ε dans les 6 triangles irréguliers
        let eps = [a[1], b[2], a[2] + b[1]];
        let d = plane_d(210e3, 0.3, state);
        let sigma: [f64; 3] = std::array::from_fn(|r| (0..3).map(|c| d[r][c] * eps[c]).sum());
        assert_eq!(f.stresses.len(), 6);
        for (e, points) in f.stresses.iter().enumerate() {
            assert_eq!(points.len(), 1);
            let p = points[0];
            let scale = sigma.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            assert!(
                close(p.sigma_xx, sigma[0], scale),
                "{state:?} σxx élément {e}"
            );
            assert!(
                close(p.sigma_yy, sigma[1], scale),
                "{state:?} σyy élément {e}"
            );
            assert!(
                close(p.tau_xy, sigma[2], scale),
                "{state:?} τxy élément {e}"
            );
        }
        // lissage aux Nodes : une contrainte constante reste exacte à chaque Node
        let nodal = nodal_stresses(&mesh, &f).unwrap();
        assert_eq!(nodal.len(), nodes.len());
        for (n, p) in nodal.iter().enumerate() {
            let scale = sigma.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            assert!(close(p.x, nodes[n].x, 1.0) && close(p.y, nodes[n].y, 1.0));
            assert!(close(p.sigma_xx, sigma[0], scale), "{state:?} σxx Node {n}");
            assert!(close(p.sigma_yy, sigma[1], scale), "{state:?} σyy Node {n}");
            assert!(close(p.tau_xy, sigma[2], scale), "{state:?} τxy Node {n}");
        }
    }
}

#[test]
fn unconstrained_plate_is_reported_as_singular() {
    let mesh = rect_mesh(1, 1, 1.0, 1.0, 0.3, PlaneState::Stress);
    let result = Model::new(
        mesh,
        vec![],
        vec![Load {
            node: 3,
            dof: Dof::Uy,
            value: 1.0,
        }],
    )
    .unwrap()
    .solve();
    assert_eq!(result, Err(MailleError::Singular));
}

#[test]
fn invalid_plane_inputs_are_typed_errors() {
    let m = Material::with_poisson(210e3, 0.3).unwrap();
    let unit = vec![
        Node { x: 0.0, y: 0.0 },
        Node { x: 1.0, y: 0.0 },
        Node { x: 0.0, y: 1.0 },
    ];
    let make =
        |nodes: [usize; 3]| Element::Tri3(Tri3::new(nodes, m, 1.0, PlaneState::Stress).unwrap());

    // Poisson manquant ou hors bornes, épaisseur invalide
    assert_eq!(
        Tri3::new(
            [0, 1, 2],
            Material::new(210e3).unwrap(),
            1.0,
            PlaneState::Stress
        ),
        Err(MailleError::MissingPoisson)
    );
    for nu in [0.5, -1.0, f64::NAN] {
        assert!(matches!(
            Material::with_poisson(210e3, nu),
            Err(MailleError::InvalidProperty { .. })
        ));
    }
    assert!(matches!(
        Tri3::new([0, 1, 2], m, 0.0, PlaneState::Stress),
        Err(MailleError::InvalidProperty { .. })
    ));

    // géométrie
    assert_eq!(
        Mesh::new(unit.clone(), vec![make([0, 2, 1])]),
        Err(MailleError::ClockwiseTriangle { element: 0 })
    );
    let aligned = vec![
        Node { x: 0.0, y: 0.0 },
        Node { x: 1.0, y: 1.0 },
        Node { x: 2.0, y: 2.0 },
    ];
    assert_eq!(
        Mesh::new(aligned, vec![make([0, 1, 2])]),
        Err(MailleError::DegenerateTriangle { element: 0 })
    );
    assert_eq!(
        Mesh::new(unit.clone(), vec![make([0, 1, 7])]),
        Err(MailleError::UnknownNode { node: 7 })
    );

    // familles et DOF
    let bar = Element::Bar(Bar {
        nodes: [0, 1],
        material: Material::new(1.0).unwrap(),
        section: maille_core::material::Section::new(1.0).unwrap(),
    });
    assert_eq!(
        Mesh::new(unit.clone(), vec![make([0, 1, 2]), bar]),
        Err(MailleError::MixedElementTypes)
    );
    let mesh = Mesh::new(unit.clone(), vec![make([0, 1, 2])]).unwrap();
    for dof in [Dof::U, Dof::V, Dof::Theta] {
        assert_eq!(
            Model::new(
                mesh.clone(),
                vec![Bc {
                    node: 0,
                    dof,
                    value: 0.0
                }],
                vec![]
            ),
            Err(MailleError::InvalidDof { node: 0, dof })
        );
    }
    assert_eq!(
        Model::new(
            mesh,
            vec![
                Bc {
                    node: 0,
                    dof: Dof::Ux,
                    value: 0.0
                },
                Bc {
                    node: 0,
                    dof: Dof::Ux,
                    value: 1.0
                },
            ],
            vec![]
        ),
        Err(MailleError::DuplicateBc {
            node: 0,
            dof: Dof::Ux
        })
    );
}

#[test]
fn distributed_load_on_a_triangle_is_refused() {
    let mesh = rect_mesh(1, 1, 1.0, 1.0, 0.3, PlaneState::Stress);
    let result = Model::new(
        mesh,
        vec![
            Bc {
                node: 0,
                dof: Dof::Ux,
                value: 0.0,
            },
            Bc {
                node: 0,
                dof: Dof::Uy,
                value: 0.0,
            },
            Bc {
                node: 1,
                dof: Dof::Uy,
                value: 0.0,
            },
        ],
        vec![],
    )
    .unwrap()
    .with_distributed_loads(vec![DistributedLoad {
        element: 0,
        value: 1.0,
    }])
    .unwrap()
    .solve();
    assert_eq!(
        result,
        Err(MailleError::UnsupportedDistributedLoad { element: 0 })
    );
}

#[test]
fn uniaxial_tension_gives_uniform_stress_in_every_triangle() {
    let (nx, ny, nu) = (3, 2, 0.3);
    for state in [PlaneState::Stress, PlaneState::Strain] {
        let (f, sigma) = uniaxial_tension(nx, ny, 6.0, 4.0, nu, state);
        assert_eq!(f.stresses.len(), 2 * nx * ny);
        for (e, points) in f.stresses.iter().enumerate() {
            let p = points[0];
            assert!(close(p.sigma_xx, sigma, sigma), "{state:?} σxx élément {e}");
            assert!(close(p.sigma_yy, 0.0, sigma), "{state:?} σyy élément {e}");
            assert!(close(p.tau_xy, 0.0, sigma), "{state:?} τxy élément {e}");
            let (szz, vm) = match state {
                PlaneState::Stress => (0.0, sigma),
                PlaneState::Strain => (nu * sigma, sigma * (1.0 - nu + nu * nu).sqrt()),
            };
            assert!(close(p.sigma_zz, szz, sigma), "{state:?} σzz élément {e}");
            assert!(
                close(p.von_mises(), vm, sigma),
                "{state:?} von Mises élément {e}"
            );
        }
    }
}

#[test]
fn tri3_gauss_point_is_the_centroid() {
    // premier triangle de la grille : (0,0), (2,0), (2,2) pour 3 × 2 cases sur 6 × 4
    let (f, _) = uniaxial_tension(3, 2, 6.0, 4.0, 0.3, PlaneState::Stress);
    let p = f.stresses[0][0];
    assert!(close(p.x, 4.0 / 3.0, 1.0) && close(p.y, 2.0 / 3.0, 1.0));
}
