//! Le binaire `maille` :
//! - `maille solve modele.toml` résout un modèle (fichiers `.msh` et `.toml`) et écrit les
//!   résultats (`.vtu` et `.pvd`), voir `files.rs` ;
//! - `maille view resultats.pvd` affiche des résultats, voir `files.rs` ;
//! - `maille mesh dessin.dxf --size h` maille un dessin DXF 2D, `maille mesh piece.stl [--size h]`
//!   une surface STL 3D, `maille mesh piece.step --size h` une pièce STEP, voir `mesh.rs` ;
//! - les cas d'exemple ci-dessous, comparés à la solution analytique.
//!
//! Cas d'exemple :
//! - `console` : console en acier chargée en bout ;
//! - `kirsch` : plaque trouée de Kirsch, voir `kirsch.rs` ;
//! - `cylindre` : conduction dans un cylindre creux, voir `cylinder.rs` ;
//! - `trempe` : trempe d'une plaque, thermique transitoire, voir `quench.rs` ;
//! - `cylindre-thermo` : contraintes thermiques dans un tube épais, voir `thermal_stress.rs` ;
//! - `cylindre-plastique` : tube épais qui plastifie sous pression, voir `plastic_tube.rs`.
//!
//! Unités cohérentes (N, mm, MPa en mécanique ; W, mm, K en thermique, plus s et kg en transitoire) :
//! le cœur ne convertit rien.

mod cylinder;
mod files;
mod kirsch;
mod mesh;
mod plastic_tube;
mod quench;
mod thermal_stress;

use std::f64::consts::FRAC_PI_4;

use maille_core::material::Material;
use maille_core::mesh::{Beam, Dof, Element, Mesh, Node};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions, Triangulation};
use maille_core::model::{Bc, Load, Model};
use maille_core::results::Stress;
use maille_core::MailleError;

const YOUNG: f64 = 210e3; // MPa
const INERTIA: f64 = 8.0e4; // mm⁴
const LENGTH: f64 = 1000.0; // mm
const TIP_LOAD: f64 = 100.0; // N
const ELEMENTS: usize = 4;

const USAGE: &str = "maille : calcul par éléments finis, mécanique et thermique

Usage :
  maille [resultats.pvd|.vtu]                 ouvre l'interface graphique (sur ces résultats, s'il y en a)
  maille solve <modele.toml> [-o <préfixe>]   résout un modèle, écrit <préfixe>.pvd et les .vtu
  maille view <resultats.pvd|.vtu> [indice]   affiche un incrément (le dernier par défaut)
  maille view <résultats> --png <image.png>   écrit une image au lieu d'ouvrir une fenêtre
                                              (--size LxH, --angles lacet,tangage, --field k, --deformed, --cut x=0.5)
  maille mesh <dessin.dxf> --size <h> [--quadratic] [-o <maillage.msh>]
                                              maille un dessin DXF 2D (calque = nom de bord)
  maille mesh <piece.stl> [--size <h>] [--feature-angle <degrés>] [--quadratic] [-o <maillage.msh>]
                                              maille une surface STL 3D (solid = nom de bord) ;
                                              sans --size, à la finesse des triangles du STL
  maille mesh <piece.step> --size <h> [--quadratic] [-o <maillage.msh>]
                                              maille une pièce STEP (nom de face = nom de bord)

Cas d'exemple, comparés à la solution analytique :
  maille console | kirsch | cylindre | trempe | cylindre-thermo | cylindre-plastique

Sans -o, les résultats sont écrits à côté du modèle, sous son nom. Les chemins du modèle
(`mesh`, `from`) sont relatifs à son dossier.";

fn main() {
    if let Err(e) = run() {
        eprintln!("erreur : {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        None => maille_gui::run(None),
        Some(path) if path.ends_with(".pvd") || path.ends_with(".vtu") => {
            maille_gui::run(Some(path.into()))
        }
        Some("help" | "--help" | "-h") => {
            println!("{USAGE}");
            Ok(())
        }
        Some("solve") => files::solve(&args[1..]),
        Some("view") => files::view(&args[1..]),
        Some("mesh") => mesh::mesh(&args[1..]),
        Some("console") => console(),
        Some("kirsch") => kirsch::run(),
        Some("cylindre") => cylinder::run(),
        Some("trempe") => quench::run(),
        Some("cylindre-thermo") => thermal_stress::run(),
        Some("cylindre-plastique") => plastic_tube::run(),
        Some(other) => Err(format!(
            "commande inconnue « {other} » : choisir solve, view, mesh, console, kirsch, cylindre, trempe, cylindre-thermo ou cylindre-plastique (maille help pour le détail)"
        )
        .into()),
    }
}

fn console() -> Result<(), Box<dyn std::error::Error>> {
    let material = Material::new(YOUNG)?;
    let nodes = (0..=ELEMENTS)
        .map(|i| Node::new_2d(LENGTH * i as f64 / ELEMENTS as f64, 0.0))
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

/// Quart de tube rᵢ < r < rₑ maillé en Tri6 (taille h), bords nommés sym_x (y = 0), exterieur,
/// sym_y (x = 0) et interieur.
fn quarter_tube(inner: f64, outer: f64, size: f64) -> Result<Triangulation<6>, MailleError> {
    let geometry = Geometry {
        outer: vec![
            Curve::line([inner, 0.0], size, "sym_x"),
            Curve::arc([outer, 0.0], [0.0, 0.0], false, size, "exterieur"),
            Curve::line([0.0, outer], size, "sym_y"),
            Curve::arc([0.0, inner], [0.0, 0.0], true, size, "interieur"),
        ],
        holes: vec![],
    };
    Ok(triangulate(&geometry, &MeshOptions::default())?.to_quadratic())
}

/// (σr, σθ, σz) d'un point de Gauss.
fn polar(s: &Stress) -> [f64; 3] {
    let r = s.x.hypot(s.y);
    let (c, sn) = (s.x / r, s.y / r);
    [
        s.sigma_xx * c * c + s.sigma_yy * sn * sn + 2.0 * s.tau_xy * c * sn,
        s.sigma_xx * sn * sn + s.sigma_yy * c * c - 2.0 * s.tau_xy * c * sn,
        s.sigma_zz,
    ]
}

/// Point de Gauss le plus proche du point de rayon r à θ = 45°.
fn nearest(stresses: &[Vec<Stress>], r: f64) -> Option<&Stress> {
    let target = [r * FRAC_PI_4.cos(), r * FRAC_PI_4.sin()];
    let d = |s: &Stress| (s.x - target[0]).hypot(s.y - target[1]);
    stresses
        .iter()
        .flatten()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
}

fn report(label: &str, computed: f64, exact: f64) {
    let error = (computed - exact).abs() / exact.abs();
    println!("{label:<26}{computed:>14.6e}{exact:>14.6e}{error:>12.1e}");
}
