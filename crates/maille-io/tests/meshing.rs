//! Du fichier CAO au maillage : les exemples de `examples/` donnent les maillages attendus, et un
//! message par erreur.

use std::path::{Path, PathBuf};

use maille_io::meshing::{extent, mesh, suggested_size, Recipe, Source};
use maille_io::msh::{CellKind, DOMAIN_SET};
use maille_io::IoError;

fn example(name: &str) -> (PathBuf, Vec<u8>) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(name);
    let bytes = std::fs::read(&path).unwrap();
    (path, bytes)
}

#[test]
fn the_format_comes_from_the_extension() {
    assert_eq!(Source::of(Path::new("a/plaque.DXF")), Some(Source::Dxf));
    assert_eq!(Source::of(Path::new("cube.stl")), Some(Source::Stl));
    assert_eq!(Source::of(Path::new("pièce.stp")), Some(Source::Step));
    assert_eq!(Source::of(Path::new("pièce.step")), Some(Source::Step));
    assert_eq!(Source::of(Path::new("maillage.msh")), None);
    assert_eq!(Source::of(Path::new("sans_extension")), None);
    assert!(Source::Dxf.needs_size() && Source::Step.needs_size() && !Source::Stl.needs_size());
}

#[test]
fn a_dxf_drawing_is_meshed_into_named_triangles() {
    let (_, bytes) = example("plaque.dxf");
    let recipe = Recipe {
        size: Some(5.0),
        ..Recipe::default()
    };
    let linear = mesh(Source::Dxf, &bytes, &recipe).unwrap().file;
    assert_eq!(linear.kind, CellKind::Tri3);
    assert!(linear.sets.len() > 1 && linear.sets.contains_key(DOMAIN_SET));
    let quadratic = mesh(
        Source::Dxf,
        &bytes,
        &Recipe {
            quadratic: true,
            ..recipe
        },
    )
    .unwrap()
    .file;
    assert_eq!(quadratic.kind, CellKind::Tri6);
    assert_eq!(quadratic.cells.len(), linear.cells.len());
    // même recette, même maillage
    assert_eq!(mesh(Source::Dxf, &bytes, &recipe).unwrap().file, linear);
}

#[test]
fn an_stl_surface_is_meshed_without_a_size() {
    let (_, bytes) = example("cube.stl");
    let meshed = mesh(Source::Stl, &bytes, &Recipe::default()).unwrap();
    assert_eq!(meshed.file.kind, CellKind::Tet4);
    assert!(meshed.unit.is_none() && !meshed.file.cells.is_empty());
}

#[test]
fn a_step_piece_reports_its_unit() {
    let (_, bytes) = example("sphere.step");
    let size = suggested_size(Source::Step, extent(Source::Step, &bytes).unwrap());
    let meshed = mesh(
        Source::Step,
        &bytes,
        &Recipe {
            size: Some(size),
            ..Recipe::default()
        },
    )
    .unwrap();
    assert_eq!(meshed.file.kind, CellKind::Tet4);
    assert!(meshed.unit.is_some());
}

#[test]
fn a_missing_size_is_refused_with_a_message_that_says_what_to_do() {
    for (source, name) in [(Source::Dxf, "plaque.dxf"), (Source::Step, "sphere.step")] {
        let (_, bytes) = example(name);
        match mesh(source, &bytes, &Recipe::default()) {
            Err(IoError::InvalidModel { problem }) => {
                assert!(problem.contains("taille h"), "{problem}")
            }
            other => panic!("{name} : {other:?}"),
        }
    }
}

#[test]
fn the_extent_is_the_diagonal_of_the_box_of_the_vertices() {
    let (_, cube) = example("cube.stl");
    // cube unité : diagonale √3
    assert!((extent(Source::Stl, &cube).unwrap() - 3f64.sqrt()).abs() < 1e-9);
    let (_, plaque) = example("plaque.dxf");
    assert!(extent(Source::Dxf, &plaque).unwrap() > 0.0);
}

#[test]
fn the_suggested_size_is_a_round_fraction_of_the_extent() {
    let close = |got: f64, want: f64| assert!((got - want).abs() < 1e-12, "{got} contre {want}");
    // 2D : un vingtième ; 3D : un dixième ; deux chiffres significatifs
    close(suggested_size(Source::Dxf, 100.0), 5.0);
    close(suggested_size(Source::Stl, 100.0), 10.0);
    close(suggested_size(Source::Dxf, 3.0_f64.sqrt()), 0.087);
    close(suggested_size(Source::Dxf, 0.0), 1.0);
}
