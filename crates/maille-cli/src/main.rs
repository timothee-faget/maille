//! Cas d'exemple comparés à la solution analytique :
//! - `console` (par défaut) : console en acier chargée en bout ;
//! - `kirsch` : plaque trouée de Kirsch, voir `kirsch.rs`.
//!
//! Unités cohérentes (N, mm, MPa) : le cœur ne convertit rien.

mod grid;
mod kirsch;

use maille_core::material::Material;
use maille_core::mesh::{Beam, Dof, Element, Mesh, Node};
use maille_core::model::{Bc, Load, Model};

const YOUNG: f64 = 210e3; // MPa
const INERTIA: f64 = 8.0e4; // mm⁴
const LENGTH: f64 = 1000.0; // mm
const TIP_LOAD: f64 = 100.0; // N
const ELEMENTS: usize = 4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    match std::env::args().nth(1).as_deref() {
        None | Some("console") => console(),
        Some("kirsch") => kirsch::run(),
        Some(other) => Err(format!("cas inconnu « {other} » : choisir console ou kirsch").into()),
    }
}

fn console() -> Result<(), Box<dyn std::error::Error>> {
    let material = Material::new(YOUNG)?;
    let nodes = (0..=ELEMENTS)
        .map(|i| Node {
            x: LENGTH * i as f64 / ELEMENTS as f64,
            y: 0.0,
        })
        .collect();
    let elements = (0..ELEMENTS)
        .map(|i| Beam::new([i, i + 1], material, INERTIA).map(Element::Beam))
        .collect::<Result<_, _>>()?;
    let mesh = Mesh::new(nodes, elements)?;

    let bcs = vec![
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
    ];
    let loads = vec![Load {
        node: ELEMENTS,
        dof: Dof::V,
        value: TIP_LOAD,
    }];
    let field = Model::new(mesh, bcs, loads)?.solve()?;

    let deflection = field.displacements[2 * ELEMENTS];
    let rotation = field.displacements[2 * ELEMENTS + 1];
    let exact_deflection = TIP_LOAD * LENGTH.powi(3) / (3.0 * YOUNG * INERTIA);
    let exact_rotation = TIP_LOAD * LENGTH.powi(2) / (2.0 * YOUNG * INERTIA);

    println!("Console : L = {LENGTH} mm, E = {YOUNG} MPa, I = {INERTIA} mm4, P = {TIP_LOAD} N, {ELEMENTS} Elements");
    println!(
        "{:<26}{:>14}{:>14}{:>12}",
        "", "calcul", "analytique", "écart rel."
    );
    report("flèche en bout (mm)", deflection, exact_deflection);
    report("rotation en bout (rad)", rotation, exact_rotation);
    report("réaction verticale (N)", field.reactions[0], -TIP_LOAD);
    report(
        "moment d'encastr. (N.mm)",
        field.reactions[1],
        -TIP_LOAD * LENGTH,
    );
    Ok(())
}

fn report(label: &str, computed: f64, exact: f64) {
    let error = (computed - exact).abs() / exact.abs();
    println!("{label:<26}{computed:>14.6e}{exact:>14.6e}{error:>12.1e}");
}
