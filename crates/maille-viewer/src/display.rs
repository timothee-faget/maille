//! Ce qu'on montre et comment : champ, colormap, déformée, caméra, plan de coupe. Partagé par la
//! fenêtre (`app`) et la capture hors écran : chacun lui donne un `Renderer` à remplir.

use crate::camera::Camera;
use crate::colormap::Colormap;
use crate::cut::{clip_segments, clip_triangles, section, Patch, Plane};
use crate::legend;
use crate::render::Renderer;
use crate::scene::{range, Quantity, Scene};
use crate::vec3::V3;

/// Plan de coupe perpendiculaire à un axe, à une fraction de la boîte du Mesh (0 : au minimum,
/// 1 : au maximum). On garde le côté des petites coordonnées, ou des grandes si `flipped`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cut {
    pub axis: usize,
    pub fraction: f64,
    pub flipped: bool,
}

pub struct Display {
    pub scene: Scene,
    pub quantity: Quantity,
    pub colormap: Colormap,
    pub deformed: bool,
    pub camera: Camera,
    /// Vrai tant que l'utilisateur n'a ni zoomé ni déplacé : la vue se recadre avec la fenêtre.
    pub follow_window: bool,
    pub cut: Option<Cut>,
}

impl Display {
    pub fn new(scene: Scene) -> Self {
        Self {
            camera: Camera::fit(scene.bounds(0.0), [1.0, 1.0], 0.0),
            quantity: scene.quantities()[0],
            scene,
            colormap: Colormap::Rainbow,
            deformed: false,
            follow_window: true,
            cut: None,
        }
    }

    pub fn deformation_scale(&self) -> f64 {
        if self.deformed {
            self.scene.deformation_scale()
        } else {
            0.0
        }
    }

    /// Le plan de coupe courant, avec sa position le long de l'axe, dans le Mesh affiché.
    fn plane(&self) -> Option<(Plane, f64)> {
        let cut = self.cut.filter(|_| self.scene.is_solid())?;
        let (lo, hi) = self.scene.bounds(self.deformation_scale());
        let position = lo[cut.axis] + cut.fraction * (hi[cut.axis] - lo[cut.axis]);
        let sign = if cut.flipped { -1.0 } else { 1.0 };
        let mut normal = [0.0; 3];
        normal[cut.axis] = sign;
        Some((
            Plane {
                normal,
                offset: sign * position,
            },
            position,
        ))
    }

    /// Recadre la caméra sur le Mesh affiché, sans changer son orientation.
    pub fn refit(&mut self, renderer: &Renderer) {
        let reserved = legend::reserved_width(renderer.scale_factor());
        let bounds = self.scene.bounds(self.deformation_scale());
        self.camera.refit(bounds, renderer.size(), reserved);
        self.follow_window = true;
    }

    /// Envoie au GPU la géométrie, le champ, la caméra et la légende courants.
    pub fn refresh(&self, renderer: &mut Renderer) {
        let values = self.scene.values(self.quantity);
        let (lo, hi) = range(&values);
        let scale = self.deformation_scale();
        let positions = self.scene.positions(scale);

        let mut surface = Patch {
            positions: positions.clone(),
            values,
            triangles: self.scene.triangles.clone(),
        };
        let mut lines: Vec<V3> = self
            .scene
            .edges
            .iter()
            .flat_map(|e| e.map(|i| positions[i as usize]))
            .collect();
        let mut title = self.quantity.label().to_string();
        if self.deformed {
            let factor = legend::format_value(scale);
            title += &format!("\nDéformée ×{}", factor.trim_start_matches('+'));
        }
        if let (Some((plane, position)), Some(cut)) = (self.plane(), self.cut) {
            // peau coupée, et la coupe elle-même avec le champ à l'intérieur
            let mut kept =
                clip_triangles(&positions, &surface.values, &self.scene.triangles, plane);
            let inside = section(&self.scene.solids, &positions, &surface.values, plane);
            let first = kept.positions.len() as u32;
            kept.positions.extend(inside.positions);
            kept.values.extend(inside.values);
            kept.triangles
                .extend(inside.triangles.iter().map(|t| t.map(|i| i + first)));
            surface = kept;
            lines = clip_segments(&positions, &self.scene.edges, plane)
                .into_iter()
                .flatten()
                .collect();
            let side = if cut.flipped { "≥" } else { "≤" };
            title += &format!(
                "\nCoupe {} {side} {}",
                ["x", "y", "z"][cut.axis],
                legend::format_value(position).trim_start_matches('+')
            );
        }

        let to_f32 = |p: &V3| p.map(|c| c as f32);
        renderer.set_geometry(
            &surface.positions.iter().map(to_f32).collect::<Vec<_>>(),
            &surface.values.iter().map(|&v| v as f32).collect::<Vec<_>>(),
            surface.triangles.as_flattened(),
            &lines.iter().map(to_f32).collect::<Vec<_>>(),
        );
        renderer.set_range((lo, hi));
        renderer.set_camera(&self.camera);
        let layout = legend::layout(
            renderer.size(),
            (lo, hi),
            self.colormap,
            renderer.scale_factor(),
            title.lines().count(),
        );
        renderer.set_legend(&layout, &title);
    }
}
