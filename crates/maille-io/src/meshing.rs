//! Du fichier CAO au maillage : dessin DXF 2D, surface STL 3D ou pièce STEP → `MshFile`. Partagé par
//! `maille mesh` et l'interface, qui obtiennent ainsi le même maillage. Le calque d'une entité DXF, le
//! nom d'un `solid` STL (ou `surface_N`), le nom d'une face STEP (ou `face_N`) nomme le bord.

use std::path::Path;

use maille_core::mesher::{
    tetrahedralize, tetrahedralize_brep, triangulate, MeshOptions, TetOptions,
};

use crate::msh::MshFile;
use crate::{dxf, step, stl, IoError};

/// Angle vif par défaut du STL, en degrés : au-delà, deux triangles voisins sont deux surfaces.
pub const FEATURE_ANGLE: f64 = 30.0;

/// Format du fichier CAO.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Dxf,
    Stl,
    Step,
}

impl Source {
    /// Le format d'après l'extension du fichier (sans tenir compte de la casse) ; `None` s'il n'est
    /// pas pris en charge.
    pub fn of(path: &Path) -> Option<Self> {
        match path.extension()?.to_str()?.to_ascii_lowercase().as_str() {
            "dxf" => Some(Self::Dxf),
            "stl" => Some(Self::Stl),
            "step" | "stp" => Some(Self::Step),
            _ => None,
        }
    }

    /// Le STL se maille sans taille (finesse de ses triangles) ; le DXF et le STEP l'exigent.
    pub fn needs_size(self) -> bool {
        self != Self::Stl
    }
}

/// Réglages du maillage.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Recipe {
    /// Taille h des Elements, dans les unités du fichier.
    pub size: Option<f64>,
    /// Angle vif du STL, en degrés.
    pub feature_angle: f64,
    /// Elements du second ordre (Tri6, Tet10).
    pub quadratic: bool,
}

impl Default for Recipe {
    fn default() -> Self {
        Self {
            size: None,
            feature_angle: FEATURE_ANGLE,
            quadratic: false,
        }
    }
}

/// Maillage obtenu.
#[derive(Debug, Clone, PartialEq)]
pub struct Meshed {
    pub file: MshFile,
    /// Unité de longueur annoncée par un fichier STEP, sans rien convertir.
    pub unit: Option<String>,
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn missing_size(source: Source) -> IoError {
    IoError::InvalidModel {
        problem: match source {
            Source::Dxf => "la taille h manque : la donner, dans les unités du dessin",
            _ => "la taille h manque : la donner, dans l'unité de longueur du fichier STEP",
        }
        .to_string(),
    }
}

/// Maille le contenu d'un fichier CAO. Erreur si le format exige une taille et que `recipe.size`
/// manque, ou si le mailleur refuse la géométrie (le message dit quoi corriger).
pub fn mesh(source: Source, bytes: &[u8], recipe: &Recipe) -> Result<Meshed, IoError> {
    let tet_options = TetOptions {
        size: recipe.size,
        ..TetOptions::default()
    };
    match source {
        Source::Dxf => {
            let size = recipe.size.ok_or_else(|| missing_size(source))?;
            let geometry = dxf::read_geometry(&text(bytes), size)?;
            let t = triangulate(&geometry, &MeshOptions::default())?;
            let file = if recipe.quadratic {
                MshFile::from_triangulation(&t.to_quadratic())?
            } else {
                MshFile::from_triangulation(&t)?
            };
            Ok(Meshed { file, unit: None })
        }
        Source::Stl => {
            let surface = stl::read_surface(bytes, recipe.feature_angle)?;
            let t = tetrahedralize(&surface, &tet_options)?;
            let file = if recipe.quadratic {
                MshFile::from_tetrahedralization(&t.to_quadratic())?
            } else {
                MshFile::from_tetrahedralization(&t)?
            };
            Ok(Meshed { file, unit: None })
        }
        Source::Step => {
            recipe.size.ok_or_else(|| missing_size(source))?;
            let (brep, unit) = step::read_brep(&text(bytes))?;
            let t = tetrahedralize_brep(&brep, &tet_options)?;
            let file = if recipe.quadratic {
                MshFile::from_tetrahedralization(&t.to_quadratic())?
            } else {
                MshFile::from_tetrahedralization(&t)?
            };
            Ok(Meshed {
                file,
                unit: Some(unit.to_string()),
            })
        }
    }
}

/// Longueur caractéristique de la pièce : la diagonale de la boîte de ses sommets, dans les unités du
/// fichier. Sert à proposer une taille h de départ.
pub fn extent(source: Source, bytes: &[u8]) -> Result<f64, IoError> {
    let points: Vec<[f64; 3]> = match source {
        Source::Dxf => {
            let g = dxf::read_geometry(&text(bytes), 1.0)?;
            g.outer
                .iter()
                .chain(g.holes.iter().flatten())
                .map(|c| [c.start[0], c.start[1], 0.0])
                .collect()
        }
        Source::Stl => stl::read_surface(bytes, FEATURE_ANGLE)?.vertices,
        Source::Step => step::read_brep(&text(bytes))?.0.vertices,
    };
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    Ok(if points.is_empty() {
        0.0
    } else {
        (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
    })
}

/// Taille h de départ pour une pièce d'étendue `extent` : un vingtième en 2D, un dixième en 3D (un
/// maillage 3D grossit vite), arrondie à deux chiffres significatifs.
pub fn suggested_size(source: Source, extent: f64) -> f64 {
    let raw = extent / if source == Source::Dxf { 20.0 } else { 10.0 };
    if raw <= 0.0 || !raw.is_finite() {
        return 1.0;
    }
    let scale = 10f64.powf(raw.log10().floor() - 1.0);
    (raw / scale).round() * scale
}
