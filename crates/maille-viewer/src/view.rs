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
    /// Cadre la boîte (min, max) au centre d'une fenêtre de `size` pixels.
    pub fn fit(bounds: ([f64; 2], [f64; 2]), size: [f64; 2]) -> Self {
        let (lo, hi) = bounds;
        Self {
            center: [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0],
            zoom: FILL * (size[0] / (hi[0] - lo[0])).min(size[1] / (hi[1] - lo[1])),
        }
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
        let view = View::fit(bounds, size);
        let (lo, hi) = (clip(&view, size, bounds.0), clip(&view, size, bounds.1));
        assert!((lo[0] + 0.9).abs() < 1e-6 && (hi[0] - 0.9).abs() < 1e-6);
        assert!(lo[1] > -0.9 && hi[1] < 0.9 && (lo[1] + hi[1]).abs() < 1e-6);
    }

    #[test]
    fn proportions_are_kept_whatever_the_window_shape() {
        // un carré du Mesh reste un carré à l'écran : même nombre de pixels en x et en y
        let size = [1200.0, 300.0];
        let view = View::fit(([0.0, 0.0], [1.0, 1.0]), size);
        let (scale, _) = view.transform(size);
        let pixels = [
            f64::from(scale[0]) * size[0] / 2.0,
            f64::from(scale[1]) * size[1] / 2.0,
        ];
        assert!((pixels[0] - pixels[1]).abs() < 1e-3, "{pixels:?}");
        assert!((pixels[1] - 0.9 * 300.0).abs() < 1e-3);
    }
}
