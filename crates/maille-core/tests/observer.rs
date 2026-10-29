//! Rappel de progression et d'annulation des calculs longs : mêmes résultats avec et sans rappel,
//! événements dans l'ordre, arrêt propre sur `Cancelled`.

mod common;

use common::rect_grid;
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tri3};
use maille_core::model::{Amplitude, Analysis, Bc, Event, Flow, Model, Step};
use maille_core::MailleError;

const L: f64 = 2.0;
const W: f64 = 1.0;

fn bar(material: Material) -> Mesh {
    let (nodes, tris) = rect_grid(4, 2, L, W);
    let elements = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, material, 1.0, PlaneState::Strain).unwrap()))
        .collect();
    Mesh::new(nodes, elements).unwrap()
}

/// Traction plastique selon x : symétries en x = 0 et y = 0, déplacement imposé en x = L.
fn tension(delta: f64) -> Model {
    let material = Material::with_poisson(210e3, 0.3)
        .and_then(|m| m.with_plasticity(250.0, 2000.0))
        .unwrap();
    let mesh = bar(material);
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if n.x.abs() < 1e-12 {
            bcs.push(Bc {
                node: i,
                dof: Dof::Ux,
                value: 0.0,
            });
        }
        if n.y.abs() < 1e-12 {
            bcs.push(Bc {
                node: i,
                dof: Dof::Uy,
                value: 0.0,
            });
        }
        if (n.x - L).abs() < 1e-12 {
            bcs.push(Bc {
                node: i,
                dof: Dof::Ux,
                value: delta,
            });
        }
    }
    Model::new(mesh, bcs, Vec::new()).unwrap()
}

/// Barre chauffée à gauche, froide à droite, en transitoire.
fn heated() -> Model {
    let material = Material::conductive(1.0)
        .and_then(|m| m.with_density(2.0))
        .and_then(|m| m.with_specific_heat(0.5))
        .unwrap();
    let mesh = bar(material);
    let mut bcs = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if n.x.abs() < 1e-12 {
            bcs.push(Bc {
                node: i,
                dof: Dof::T,
                value: 100.0,
            });
        }
        if (n.x - L).abs() < 1e-12 {
            bcs.push(Bc {
                node: i,
                dof: Dof::T,
                value: 0.0,
            });
        }
    }
    Model::with_analysis(mesh, Analysis::Thermal, bcs, Vec::new()).unwrap()
}

#[test]
fn the_quasi_static_results_are_the_same_with_and_without_an_observer() {
    let model = tension(0.02);
    let step = Step::quasi_static(1.0, 5).unwrap();
    let ramp = Amplitude::ramp(1.0).unwrap();
    let plain = model.solve_quasi_static(&step, &ramp).unwrap();
    let mut events = Vec::new();
    let watched = model
        .solve_quasi_static_with(&step, &ramp, &mut |e| {
            events.push(e.clone());
            Flow::Continue
        })
        .unwrap();
    assert_eq!(plain, watched);

    // un Event::Increment par incrément, dans l'ordre, et les itérations qui les précèdent
    let done: Vec<usize> = events
        .iter()
        .filter_map(|e| match e {
            Event::Increment { done, total, .. } => {
                assert_eq!(*total, 5);
                Some(*done)
            }
            _ => None,
        })
        .collect();
    assert_eq!(done, vec![1, 2, 3, 4, 5]);
    // les résidus annoncés sont ceux de l'Increment (sans coupe d'incrément ici)
    for (k, inc) in watched.iter().enumerate().skip(1) {
        let residuals: Vec<f64> = events
            .iter()
            .filter_map(|e| match e {
                Event::Iteration {
                    increment,
                    residual,
                    ..
                } if *increment == k => Some(*residual),
                _ => None,
            })
            .collect();
        assert_eq!(residuals, inc.residuals, "incrément {k}");
    }
    assert!(!events.iter().any(|e| matches!(e, Event::Cut { .. })));
}

#[test]
fn stopping_a_quasi_static_calculation_returns_cancelled() {
    let model = tension(0.02);
    let step = Step::quasi_static(1.0, 5).unwrap();
    let ramp = Amplitude::ramp(1.0).unwrap();
    // arrêt au troisième incrément fini
    let result = model.solve_quasi_static_with(&step, &ramp, &mut |e| match e {
        Event::Increment { done: 3, .. } => Flow::Stop,
        _ => Flow::Continue,
    });
    assert_eq!(result, Err(MailleError::Cancelled));
    // arrêt dès la première itération
    let result = model.solve_quasi_static_with(&step, &ramp, &mut |_| Flow::Stop);
    assert_eq!(result, Err(MailleError::Cancelled));
}

#[test]
fn the_transient_results_are_the_same_with_and_without_an_observer() {
    let model = heated();
    let step = Step::transient(4.0, 8).unwrap();
    let initial = vec![0.0; 15];
    let plain = model.solve_step(&step, &initial).unwrap();
    let mut times = Vec::new();
    let watched = model
        .solve_step_with(&step, &initial, &mut |e| {
            if let Event::Increment { done, total, time } = e {
                assert_eq!(*total, 8);
                times.push((*done, *time));
            }
            Flow::Continue
        })
        .unwrap();
    assert_eq!(plain, watched);
    assert_eq!(times.len(), 8);
    assert_eq!(times[7], (8, 4.0));

    let result = model.solve_step_with(&step, &initial, &mut |e| match e {
        Event::Increment { done: 2, .. } => Flow::Stop,
        _ => Flow::Continue,
    });
    assert_eq!(result, Err(MailleError::Cancelled));
}

#[test]
fn cancelled_is_said_in_plain_words() {
    assert_eq!(MailleError::Cancelled.to_string(), "calcul annulé");
}
