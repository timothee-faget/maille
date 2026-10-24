//! Cadrage : passage des coordonnées du Mesh aux coordonnées de clip de wgpu ([-1, 1]², y vers le haut).

/// Part de la fenêtre occupée par le Mesh au cadrage initial.
const FILL: f64 = 0.9;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct View {
    /// Point du Mesh affiché au centre de la fenêtre.
    pub center: [f64; 2],
    /// Pixels par unité de longueur du Mesh (identique en x et en y : pas de déformation de l'image).
    pub zoom: f64,
}

impl View {
    /// Cadre la boîte (min, max) dans une fenêtre de `size` pixels dont les `reserved_left` premiers
    /// pixels à gauche sont laissés libres (légende) : centrée dans la partie restante.
    pub fn fit(bounds: ([f64; 2], [f64; 2]), size: [f64; 2], reserved_left: f64) -> Self {
        let (lo, hi) = bounds;
        let width = (size[0] - reserved_left).max(1.0);
        let zoom = FILL * (width / (hi[0] - lo[0])).min(size[1] / (hi[1] - lo[1]));
        Self {
            // le centre de la fenêtre est décalé de reserved_left / 2 à gauche du centre du Mesh
            center: [
                (lo[0] + hi[0]) / 2.0 - reserved_left / 2.0 / zoom,
                (lo[1] + hi[1]) / 2.0,
            ],
            zoom,
        }
    }

    /// Point du Mesh sous le pixel `cursor` (origine en haut à gauche, y vers le bas, comme winit).
    pub fn model_point(&self, cursor: [f64; 2], size: [f64; 2]) -> [f64; 2] {
        [
            self.center[0] + (cursor[0] - size[0] / 2.0) / self.zoom,
            self.center[1] - (cursor[1] - size[1] / 2.0) / self.zoom,
        ]
    }

    /// Zoome d'un facteur `factor` (> 1 : on s'approche) en gardant fixe le point sous le curseur.
    pub fn zoom_at(&mut self, factor: f64, cursor: [f64; 2], size: [f64; 2]) {
        let anchor = self.model_point(cursor, size);
        self.zoom *= factor;
        let moved = self.model_point(cursor, size);
        self.center = [
            self.center[0] + anchor[0] - moved[0],
            self.center[1] + anchor[1] - moved[1],
        ];
    }

    /// Fait suivre la souris au Mesh : `delta` est le déplacement du curseur en pixels.
    pub fn pan(&mut self, delta: [f64; 2]) {
        self.center = [
            self.center[0] - delta[0] / self.zoom,
            self.center[1] + delta[1] / self.zoom,
        ];
    }

    /// (échelle, décalage) tels que clip = échelle · p + décalage, composante par composante.
    pub fn transform(&self, size: [f64; 2]) -> ([f32; 2], [f32; 2]) {
        let scale = [2.0 * self.zoom / size[0], 2.0 * self.zoom / size[1]];
        (
            scale.map(|s| s as f32),
            [
                (-self.center[0] * scale[0]) as f32,
                (-self.center[1] * scale[1]) as f32,
            ],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clip(view: &View, size: [f64; 2], p: [f64; 2]) -> [f64; 2] {
        let (scale, offset) = view.transform(size);
        [
            f64::from(scale[0]) * p[0] + f64::from(offset[0]),
            f64::from(scale[1]) * p[1] + f64::from(offset[1]),
        ]
    }

    #[test]
    fn a_wide_mesh_fills_ninety_percent_of_the_width_and_stays_centered() {
        // Mesh 40 × 10 dans une fenêtre 800 × 600 : la largeur limite
        let bounds = ([0.0, 0.0], [40.0, 10.0]);
        let size = [800.0, 600.0];
        let view = View::fit(bounds, size, 0.0);
        let (lo, hi) = (clip(&view, size, bounds.0), clip(&view, size, bounds.1));
        assert!((lo[0] + 0.9).abs() < 1e-6 && (hi[0] - 0.9).abs() < 1e-6);
        assert!(lo[1] > -0.9 && hi[1] < 0.9 && (lo[1] + hi[1]).abs() < 1e-6);
    }

    #[test]
    fn the_mesh_is_centered_in_the_part_of_the_window_left_free_by_the_legend() {
        // 200 pixels réservés à gauche : le Mesh remplit 90 % des 600 pixels restants, centré dessus
        let bounds = ([0.0, 0.0], [40.0, 10.0]);
        let size = [800.0, 600.0];
        let view = View::fit(bounds, size, 200.0);
        let pixel_x = |p: [f64; 2]| (clip(&view, size, p)[0] + 1.0) / 2.0 * size[0];
        let (left, right) = (pixel_x(bounds.0), pixel_x(bounds.1));
        assert!((right - left - 0.9 * 600.0).abs() < 1e-3, "{left} {right}");
        assert!(((left + right) / 2.0 - 500.0).abs() < 1e-3);
    }

    #[test]
    fn proportions_are_kept_whatever_the_window_shape() {
        // un carré du Mesh reste un carré à l'écran : même nombre de pixels en x et en y
        let size = [1200.0, 300.0];
        let view = View::fit(([0.0, 0.0], [1.0, 1.0]), size, 0.0);
        let (scale, _) = view.transform(size);
        let pixels = [
            f64::from(scale[0]) * size[0] / 2.0,
            f64::from(scale[1]) * size[1] / 2.0,
        ];
        assert!((pixels[0] - pixels[1]).abs() < 1e-3, "{pixels:?}");
        assert!((pixels[1] - 0.9 * 300.0).abs() < 1e-3);
    }

    #[test]
    fn the_window_center_shows_the_view_center_and_y_points_up() {
        let size = [800.0, 600.0];
        let view = View {
            center: [5.0, 2.0],
            zoom: 10.0,
        };
        let near =
            |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() < 1e-12 && (a[1] - b[1]).abs() < 1e-12;
        assert!(near(view.model_point([400.0, 300.0], size), [5.0, 2.0]));
        // 100 pixels plus haut à l'écran : y du Mesh plus grand de 100 / 10
        assert!(near(view.model_point([400.0, 200.0], size), [5.0, 12.0]));
    }

    #[test]
    fn zooming_keeps_the_point_under_the_cursor_fixed() {
        let size = [800.0, 600.0];
        let mut view = View::fit(([0.0, 0.0], [40.0, 40.0]), size, 0.0);
        let cursor = [123.0, 456.0];
        let before = view.model_point(cursor, size);
        view.zoom_at(2.5, cursor, size);
        let after = view.model_point(cursor, size);
        assert!((before[0] - after[0]).abs() < 1e-9 && (before[1] - after[1]).abs() < 1e-9);
        // et le reste s'est bien rapproché
        let (z0, z1) = (
            View::fit(([0.0, 0.0], [40.0, 40.0]), size, 0.0).zoom,
            view.zoom,
        );
        assert!((z1 / z0 - 2.5).abs() < 1e-12);
    }

    #[test]
    fn panning_moves_the_mesh_with_the_cursor() {
        let size = [800.0, 600.0];
        let mut view = View::fit(([0.0, 0.0], [40.0, 40.0]), size, 0.0);
        let grabbed = view.model_point([300.0, 200.0], size);
        view.pan([50.0, -30.0]);
        let under_cursor = view.model_point([350.0, 170.0], size);
        assert!((grabbed[0] - under_cursor[0]).abs() < 1e-9);
        assert!((grabbed[1] - under_cursor[1]).abs() < 1e-9);
    }
}
