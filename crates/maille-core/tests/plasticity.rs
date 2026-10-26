//! Élasto-plasticité de von Mises : Step quasi-statique, Amplitude, Newton-Raphson.

use maille_core::model::Amplitude;
use maille_core::MailleError;

#[test]
fn amplitude_is_piecewise_linear_then_constant() {
    let a = Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (3.0, -1.0)]).unwrap();
    for (t, expected) in [
        (0.0, 0.0),
        (0.25, 0.25),
        (1.0, 1.0),
        (2.0, 0.0),
        (2.5, -0.5),
        (3.0, -1.0),
        (10.0, -1.0),
    ] {
        assert!((a.at(t) - expected).abs() < 1e-15, "λ({t}) = {}", a.at(t));
    }
    let ramp = Amplitude::ramp(4.0).unwrap();
    assert!((ramp.at(1.0) - 0.25).abs() < 1e-15 && (ramp.at(5.0) - 1.0).abs() < 1e-15);
}

#[test]
fn amplitude_inputs_are_typed_errors() {
    for (points, bad) in [
        (vec![], 0),
        (vec![(0.0, 0.5), (1.0, 1.0)], 0),
        (vec![(0.1, 0.0), (1.0, 1.0)], 0),
        (vec![(0.0, 0.0), (1.0, 1.0), (1.0, 2.0)], 2),
        (vec![(0.0, 0.0), (1.0, f64::NAN)], 1),
        (vec![(0.0, 0.0), (f64::INFINITY, 1.0)], 1),
    ] {
        assert_eq!(
            Amplitude::new(points),
            Err(MailleError::InvalidAmplitude { point: bad })
        );
    }
    assert!(matches!(
        Amplitude::ramp(0.0),
        Err(MailleError::InvalidProperty { .. })
    ));
}
