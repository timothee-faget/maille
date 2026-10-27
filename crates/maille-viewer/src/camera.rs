//! Caméra orbitale orthographique : passage des coordonnées du Mesh aux coordonnées de clip de wgpu
//! (x et y dans [-1, 1], profondeur dans [0, 1], 0 = proche). Vue de face (lacet et tangage nuls) :
//! on regarde vers −z, x à droite, y en haut, ce qui est la vue d'un Mesh 2D du plan z = 0.

use crate::vec3::{add, cross, dot, scale, sub, V3};

/// Part de la fenêtre occupée par le Mesh au cadrage.
const FILL: f64 = 0.9;
/// Tangage limite, en radians : au pôle la vue tournerait sur elle-même.
const MAX_PITCH: f64 = 89.0 * std::f64::consts::PI / 180.0;
/// Radians de rotation par pixel de glissement.
const RADIANS_PER_PIXEL: f64 = 0.01;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Point du Mesh visé : au centre de la fenêtre, décalé de `shift_x`.
    pub target: V3,
    /// Rotation autour de y (lacet) puis de l'axe horizontal de l'écran (tangage), en radians.
    pub yaw: f64,
    pub pitch: f64,
    /// Pixels par unité de longueur du Mesh (identique en x et en y : pas de déformation).
    pub zoom: f64,
    /// Rayon de la sphère qui contient le Mesh : la profondeur [-rayon, rayon] devient [1, 0].
    pub depth_radius: f64,
    /// Décalage de la cible vers la droite, en pixels (place laissée à la légende).
    shift_x: f64,
}

impl Camera {
    /// Vue de face qui cadre la boîte (min, max) dans une fenêtre de `size` pixels dont les
    /// `reserved_left` premiers pixels à gauche sont laissés libres (légende).
    pub fn fit(bounds: (V3, V3), size: [f64; 2], reserved_left: f64) -> Self {
        let mut camera = Self {
            target: [0.0; 3],
            yaw: 0.0,
            pitch: 0.0,
            zoom: 1.0,
            depth_radius: 1.0,
            shift_x: 0.0,
        };
        camera.refit(bounds, size, reserved_left);
        camera
    }

    /// Recadre la boîte (min, max) sans changer l'orientation : elle remplit 90 % de la partie
    /// libre de la fenêtre, centrée dessus, quelle que soit la direction de vue.
    pub fn refit(&mut self, bounds: (V3, V3), size: [f64; 2], reserved_left: f64) {
        let (lo, hi) = bounds;
        let center = scale(add(lo, hi), 0.5);
        let (right, up, _) = self.axes();
        let (mut width, mut height) = (0.0_f64, 0.0_f64);
        for corner in corners(lo, hi) {
            let d = sub(corner, center);
            width = width.max(2.0 * dot(d, right).abs());
            height = height.max(2.0 * dot(d, up).abs());
        }
        let free = (size[0] - reserved_left).max(1.0);
        // boîte réduite à un point ou à un segment : la dimension nulle ne limite pas le zoom
        let fit = |pixels: f64, extent: f64| {
            if extent > 0.0 {
                pixels / extent
            } else {
                f64::INFINITY
            }
        };
        let zoom = FILL * fit(free, width).min(fit(size[1], height));
        self.zoom = if zoom.is_finite() { zoom } else { 1.0 };
        self.target = center;
        self.shift_x = reserved_left / 2.0;
        self.depth_radius = (0.5 * crate::vec3::norm(sub(hi, lo))).max(f64::MIN_POSITIVE);
    }

    /// (droite, haut, œil) : base orthonormée de la caméra ; `œil` va de la cible vers l'observateur.
    pub fn axes(&self) -> (V3, V3, V3) {
        let (sy, cy) = self.yaw.sin_cos();
        let (sp, cp) = self.pitch.sin_cos();
        let eye = [sy * cp, sp, cy * cp];
        let right = [cy, 0.0, -sy];
        (right, cross(eye, right), eye)
    }

    /// Fait tourner le Mesh sous la souris : `delta` est le déplacement du curseur en pixels
    /// (vers la droite : le Mesh tourne vers la droite ; vers le bas : il bascule vers l'avant).
    pub fn rotate(&mut self, delta: [f64; 2]) {
        self.yaw -= delta[0] * RADIANS_PER_PIXEL;
        self.pitch = (self.pitch + delta[1] * RADIANS_PER_PIXEL).clamp(-MAX_PITCH, MAX_PITCH);
    }

    /// Fait suivre la souris au Mesh : `delta` est le déplacement du curseur en pixels.
    pub fn pan(&mut self, delta: [f64; 2]) {
        let (right, up, _) = self.axes();
        self.target = add(
            sub(self.target, scale(right, delta[0] / self.zoom)),
            scale(up, delta[1] / self.zoom),
        );
    }

    /// Point du plan perpendiculaire à la vue qui passe par la cible, sous le pixel `cursor`
    /// (origine en haut à gauche, y vers le bas, comme winit).
    pub fn model_point(&self, cursor: [f64; 2], size: [f64; 2]) -> V3 {
        let (right, up, _) = self.axes();
        let sx = cursor[0] - size[0] / 2.0 - self.shift_x;
        let sy = size[1] / 2.0 - cursor[1];
        add(
            self.target,
            add(scale(right, sx / self.zoom), scale(up, sy / self.zoom)),
        )
    }

    /// Zoome d'un facteur `factor` (> 1 : on s'approche) en gardant fixe le point sous le curseur.
    pub fn zoom_at(&mut self, factor: f64, cursor: [f64; 2], size: [f64; 2]) {
        let anchor = self.model_point(cursor, size);
        self.zoom *= factor;
        let moved = self.model_point(cursor, size);
        self.target = add(self.target, sub(anchor, moved));
    }

    /// Coordonnées de clip (x, y, profondeur) du point `p`.
    #[cfg(test)]
    pub fn project(&self, p: V3, size: [f64; 2]) -> V3 {
        let (right, up, eye) = self.axes();
        let d = sub(p, self.target);
        [
            (dot(d, right) * self.zoom + self.shift_x) * 2.0 / size[0],
            dot(d, up) * self.zoom * 2.0 / size[1],
            0.5 - dot(d, eye) / (2.0 * self.depth_radius),
        ]
    }

    /// Matrice 4 × 4 en colonnes (clip = M · (x, y, z, 1)) pour le GPU.
    pub fn matrix(&self, size: [f64; 2]) -> [[f32; 4]; 4] {
        let (right, up, eye) = self.axes();
        let (kx, ky, kz) = (
            self.zoom * 2.0 / size[0],
            self.zoom * 2.0 / size[1],
            -0.5 / self.depth_radius,
        );
        let offset = |axis: V3, k: f64, shift: f64| shift - dot(self.target, axis) * k;
        let rows = [
            [
                right[0] * kx,
                right[1] * kx,
                right[2] * kx,
                offset(right, kx, self.shift_x * 2.0 / size[0]),
            ],
            [up[0] * ky, up[1] * ky, up[2] * ky, offset(up, ky, 0.0)],
            [
                eye[0] * kz,
                eye[1] * kz,
                eye[2] * kz,
                offset(eye, kz, 0.0) + 0.5,
            ],
            [0.0, 0.0, 0.0, 1.0],
        ];
        std::array::from_fn(|col| std::array::from_fn(|row| rows[row][col] as f32))
    }
}

/// Les huit coins de la boîte (min, max).
pub fn corners(lo: V3, hi: V3) -> [V3; 8] {
    std::array::from_fn(|k| {
        [
            if k & 1 == 0 { lo[0] } else { hi[0] },
            if k & 2 == 0 { lo[1] } else { hi[1] },
            if k & 4 == 0 { lo[2] } else { hi[2] },
        ]
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIZE: [f64; 2] = [800.0, 600.0];

    fn near(a: V3, b: V3, tol: f64) -> bool {
        (0..3).all(|k| (a[k] - b[k]).abs() < tol)
    }

    #[test]
    fn a_wide_mesh_fills_ninety_percent_of_the_width_and_stays_centered() {
        // Mesh 2D de 40 × 10 dans une fenêtre 800 × 600 : la largeur limite
        let bounds = ([0.0; 3], [40.0, 10.0, 0.0]);
        let camera = Camera::fit(bounds, SIZE, 0.0);
        let (lo, hi) = (
            camera.project(bounds.0, SIZE),
            camera.project(bounds.1, SIZE),
        );
        assert!((lo[0] + 0.9).abs() < 1e-9 && (hi[0] - 0.9).abs() < 1e-9);
        assert!(lo[1] > -0.9 && hi[1] < 0.9 && (lo[1] + hi[1]).abs() < 1e-9);
    }

    #[test]
    fn the_mesh_is_centered_in_the_part_of_the_window_left_free_by_the_legend() {
        // 200 pixels réservés à gauche : le Mesh remplit 90 % des 600 pixels restants, centré dessus
        let bounds = ([0.0; 3], [40.0, 10.0, 0.0]);
        let camera = Camera::fit(bounds, SIZE, 200.0);
        let pixel_x = |p: V3| (camera.project(p, SIZE)[0] + 1.0) / 2.0 * SIZE[0];
        let (left, right) = (pixel_x(bounds.0), pixel_x(bounds.1));
        assert!((right - left - 0.9 * 600.0).abs() < 1e-6, "{left} {right}");
        assert!(((left + right) / 2.0 - 500.0).abs() < 1e-6);
    }

    #[test]
    fn proportions_are_kept_whatever_the_window_shape() {
        let size = [1200.0, 300.0];
        let camera = Camera::fit(([0.0; 3], [1.0, 1.0, 0.0]), size, 0.0);
        let o = camera.project([0.0; 3], size);
        let (x, y) = (
            camera.project([1.0, 0.0, 0.0], size),
            camera.project([0.0, 1.0, 0.0], size),
        );
        let (px, py) = ((x[0] - o[0]) * size[0] / 2.0, (y[1] - o[1]) * size[1] / 2.0);
        assert!((px - py).abs() < 1e-9, "{px} {py}");
        assert!((py - 0.9 * 300.0).abs() < 1e-9);
    }

    #[test]
    fn the_front_view_looks_at_a_2d_mesh_with_y_up() {
        let camera = Camera::fit(([0.0; 3], [10.0, 10.0, 0.0]), SIZE, 0.0);
        let (right, up, eye) = camera.axes();
        assert!(near(right, [1.0, 0.0, 0.0], 1e-12));
        assert!(near(up, [0.0, 1.0, 0.0], 1e-12));
        assert!(near(eye, [0.0, 0.0, 1.0], 1e-12));
        // plus haut dans le Mesh, plus haut à l'écran ; plus près de l'œil (z grand), plus petite profondeur
        let (a, b) = (
            camera.project([5.0, 5.0, 0.0], SIZE),
            camera.project([5.0, 6.0, 0.0], SIZE),
        );
        assert!(b[1] > a[1]);
        let c = camera.project([5.0, 5.0, 1.0], SIZE);
        assert!(c[2] < a[2]);
    }

    #[test]
    fn the_axes_stay_orthonormal_and_the_horizon_stays_level() {
        let mut camera = Camera::fit(([0.0; 3], [1.0; 3]), SIZE, 0.0);
        for (dx, dy) in [(37.0, 12.0), (-250.0, 80.0), (400.0, -900.0), (5.0, 5000.0)] {
            camera.rotate([dx, dy]);
            let (r, u, e) = camera.axes();
            assert!((dot(r, r) - 1.0).abs() < 1e-12 && (dot(u, u) - 1.0).abs() < 1e-12);
            assert!(dot(r, u).abs() < 1e-12 && dot(r, e).abs() < 1e-12 && dot(u, e).abs() < 1e-12);
            // droite reste horizontale (y du monde vers le haut) et le haut ne retourne jamais la tête
            assert!(r[1].abs() < 1e-12 && u[1] > 0.0);
        }
        assert!(camera.pitch.abs() <= MAX_PITCH);
    }

    #[test]
    fn dragging_right_turns_the_mesh_toward_the_right() {
        let mut camera = Camera::fit(([-1.0; 3], [1.0; 3]), SIZE, 0.0);
        // le point (0, 0, 1), face avant, passe vers la droite de l'écran
        camera.rotate([30.0, 0.0]);
        assert!(camera.project([0.0, 0.0, 1.0], SIZE)[0] > camera.project([0.0; 3], SIZE)[0]);
        // vers le bas : le haut de la face avant bascule vers l'avant, donc vers le bas de l'écran
        let mut camera = Camera::fit(([-1.0; 3], [1.0; 3]), SIZE, 0.0);
        camera.rotate([0.0, 30.0]);
        assert!(
            camera.project([0.0, 1.0, 0.0], SIZE)[1]
                < camera.project([0.0, 1.0, 0.0], SIZE)[1] + 1.0
        );
        assert!(camera.project([0.0, 0.0, 1.0], SIZE)[1] < camera.project([0.0; 3], SIZE)[1]);
    }

    #[test]
    fn the_whole_box_stays_in_view_and_in_depth_range_from_every_direction() {
        let bounds = ([-3.0, 0.0, 1.0], [5.0, 2.0, 9.0]);
        let mut camera = Camera::fit(bounds, SIZE, 200.0);
        for (dx, dy) in [
            (0.0, 0.0),
            (70.0, 20.0),
            (-140.0, 60.0),
            (300.0, -40.0),
            (90.0, 400.0),
        ] {
            camera.rotate([dx, dy]);
            camera.refit(bounds, SIZE, 200.0);
            for corner in corners(bounds.0, bounds.1) {
                let c = camera.project(corner, SIZE);
                assert!(c[0].abs() <= 1.0 && c[1].abs() <= 0.9 + 1e-9, "{c:?}");
                assert!((0.0..=1.0).contains(&c[2]), "profondeur {c:?}");
            }
        }
    }

    #[test]
    fn the_matrix_gives_the_same_clip_coordinates_as_project() {
        let mut camera = Camera::fit(([0.0; 3], [4.0, 2.0, 3.0]), SIZE, 150.0);
        camera.rotate([55.0, -33.0]);
        let m = camera.matrix(SIZE);
        for p in [[0.0; 3], [4.0, 2.0, 3.0], [1.5, -0.5, 2.2]] {
            let col = |k: usize| m[k].map(f64::from);
            let clip: V3 = std::array::from_fn(|row| {
                col(0)[row] * p[0] + col(1)[row] * p[1] + col(2)[row] * p[2] + col(3)[row]
            });
            assert!(
                near(clip, camera.project(p, SIZE), 1e-5),
                "{clip:?} {:?}",
                camera.project(p, SIZE)
            );
        }
    }

    #[test]
    fn zooming_keeps_the_point_under_the_cursor_fixed() {
        let mut camera = Camera::fit(([0.0; 3], [40.0, 40.0, 0.0]), SIZE, 0.0);
        camera.rotate([40.0, 25.0]);
        let cursor = [123.0, 456.0];
        let before = camera.model_point(cursor, SIZE);
        let zoom = camera.zoom;
        camera.zoom_at(2.5, cursor, SIZE);
        assert!(near(before, camera.model_point(cursor, SIZE), 1e-9));
        assert!((camera.zoom / zoom - 2.5).abs() < 1e-12);
    }

    #[test]
    fn panning_moves_the_mesh_with_the_cursor() {
        let mut camera = Camera::fit(([0.0; 3], [40.0, 40.0, 0.0]), SIZE, 0.0);
        camera.rotate([-60.0, 30.0]);
        let grabbed = camera.model_point([300.0, 200.0], SIZE);
        camera.pan([50.0, -30.0]);
        assert!(near(
            grabbed,
            camera.model_point([350.0, 170.0], SIZE),
            1e-9
        ));
    }

    #[test]
    fn a_degenerate_box_does_not_break_the_fit() {
        let camera = Camera::fit(([2.0; 3], [2.0; 3]), SIZE, 0.0);
        assert!(camera.zoom.is_finite() && camera.zoom > 0.0);
        let flat = Camera::fit(([0.0; 3], [5.0, 0.0, 0.0]), SIZE, 0.0);
        assert!(flat.zoom.is_finite() && flat.zoom > 0.0);
    }
}
