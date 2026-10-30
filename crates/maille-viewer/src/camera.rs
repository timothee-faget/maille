//! Caméra orbitale orthographique : passage des coordonnées du Mesh aux coordonnées de clip de wgpu
//! (x et y dans [-1, 1], profondeur dans [0, 1], 0 = proche). Vue de face (lacet et tangage nuls) :
//! on regarde vers −z, x à droite, y en haut, ce qui est la vue d'un Mesh 2D du plan z = 0.
//! La rotation est un trackball : on tourne autour des axes de l'écran, sans pôle ni borne.

use crate::vec3::{add, cross, dot, norm, scale, sub, V3};

/// Part de la fenêtre occupée par le Mesh au cadrage.
const FILL: f64 = 0.9;
/// Radians de rotation par pixel de glissement.
const RADIANS_PER_PIXEL: f64 = 0.01;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Point du Mesh visé : au centre de la fenêtre, décalé de `shift_x`.
    pub target: V3,
    /// Base orthonormée (droite, haut, œil) de la caméra, dans les coordonnées du Mesh.
    frame: [V3; 3],
    /// Pixels par unité de longueur du Mesh (identique en x et en y : pas de déformation).
    pub zoom: f64,
    /// Centre et rayon de la sphère qui contient le Mesh : la profondeur [-rayon, rayon] autour du
    /// centre devient [1, 0].
    center: V3,
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
            frame: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            zoom: 1.0,
            center: [0.0; 3],
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
        self.enclose(bounds);
    }

    /// Règle la profondeur sur la sphère qui contient la boîte (min, max), sans toucher au cadrage :
    /// la profondeur se mesure depuis le centre de la boîte, pas depuis la cible que le zoom et le
    /// déplacement emmènent ailleurs, donc rien de la boîte n'est coupé, quelle que soit la vue.
    pub fn enclose(&mut self, (lo, hi): (V3, V3)) {
        self.center = scale(add(lo, hi), 0.5);
        self.depth_radius = (0.5 * norm(sub(hi, lo))).max(f64::MIN_POSITIVE);
    }

    /// (droite, haut, œil) : base orthonormée de la caméra ; `œil` va de la cible vers l'observateur.
    pub fn axes(&self) -> (V3, V3, V3) {
        let [right, up, eye] = self.frame;
        (right, up, eye)
    }

    /// Oriente la vue par une rotation `yaw` autour de y (lacet) puis `pitch` autour de l'axe
    /// horizontal de l'écran (tangage), en radians ; (0, 0) est la vue de face.
    pub fn look_from(&mut self, yaw: f64, pitch: f64) {
        let (sy, cy) = yaw.sin_cos();
        let (sp, cp) = pitch.sin_cos();
        let eye = [sy * cp, sp, cy * cp];
        let right = [cy, 0.0, -sy];
        self.frame = [right, cross(eye, right), eye];
    }

    /// Fait tourner le Mesh sous la souris : `delta` est le déplacement du curseur en pixels
    /// (vers la droite : le Mesh tourne vers la droite ; vers le bas : il bascule vers l'avant),
    /// autour des axes de l'écran, quelle que soit l'orientation.
    pub fn rotate(&mut self, delta: [f64; 2]) {
        let [right, up, eye] = self.frame;
        let axis = add(scale(up, delta[0]), scale(right, delta[1]));
        let length = norm(axis);
        if length <= 0.0 {
            return;
        }
        let (axis, angle) = (scale(axis, 1.0 / length), -length * RADIANS_PER_PIXEL);
        // la caméra tourne autour de la cible en sens inverse du Mesh ; on ré-orthonormalise pour
        // que les arrondis ne s'accumulent pas au fil des glissements
        let eye = rotated(eye, axis, angle);
        let eye = scale(eye, 1.0 / norm(eye));
        let right = rotated(right, axis, angle);
        let right = sub(right, scale(eye, dot(right, eye)));
        let right = scale(right, 1.0 / norm(right));
        self.frame = [right, cross(eye, right), eye];
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
    pub fn project(&self, p: V3, size: [f64; 2]) -> V3 {
        let (right, up, eye) = self.axes();
        let d = sub(p, self.target);
        [
            (dot(d, right) * self.zoom + self.shift_x) * 2.0 / size[0],
            dot(d, up) * self.zoom * 2.0 / size[1],
            0.5 - dot(sub(p, self.center), eye) / (2.0 * self.depth_radius),
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
        let depth_offset = 0.5 - dot(self.center, eye) * kz;
        let rows = [
            [
                right[0] * kx,
                right[1] * kx,
                right[2] * kx,
                offset(right, kx, self.shift_x * 2.0 / size[0]),
            ],
            [up[0] * ky, up[1] * ky, up[2] * ky, offset(up, ky, 0.0)],
            [eye[0] * kz, eye[1] * kz, eye[2] * kz, depth_offset],
            [0.0, 0.0, 0.0, 1.0],
        ];
        std::array::from_fn(|col| std::array::from_fn(|row| rows[row][col] as f32))
    }
}

/// `v` tourné de `angle` radians autour de l'axe unitaire `axis` (formule de Rodrigues).
fn rotated(v: V3, axis: V3, angle: f64) -> V3 {
    let (sin, cos) = angle.sin_cos();
    add(
        add(scale(v, cos), scale(cross(axis, v), sin)),
        scale(axis, dot(axis, v) * (1.0 - cos)),
    )
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

    fn assert_orthonormal(camera: &Camera) {
        let (r, u, e) = camera.axes();
        assert!((dot(r, r) - 1.0).abs() < 1e-12 && (dot(u, u) - 1.0).abs() < 1e-12);
        assert!((dot(e, e) - 1.0).abs() < 1e-12);
        assert!(dot(r, u).abs() < 1e-12 && dot(r, e).abs() < 1e-12 && dot(u, e).abs() < 1e-12);
        // base directe : droite × haut = œil
        assert!(near(cross(r, u), e, 1e-12));
    }

    #[test]
    fn the_axes_stay_orthonormal_after_many_drags() {
        let mut camera = Camera::fit(([0.0; 3], [1.0; 3]), SIZE, 0.0);
        for k in 0..10_000 {
            let k = f64::from(k);
            camera.rotate([37.0 * (0.3 * k).sin(), 12.0 * (0.7 * k).cos()]);
        }
        assert_orthonormal(&camera);
    }

    #[test]
    fn full_turns_are_possible_in_every_direction_without_jumps() {
        // un tour = 2π / 0,01 pixels ; glissé par petits pas, vertical, horizontal et en diagonale
        let turn = 2.0 * std::f64::consts::PI / RADIANS_PER_PIXEL;
        let steps = 400;
        for direction in [[0.0, 1.0], [1.0, 0.0], [0.6, -0.8]] {
            let mut camera = Camera::fit(([0.0; 3], [1.0; 3]), SIZE, 0.0);
            camera.rotate([25.0, -40.0]);
            let start = camera.axes();
            let step = turn / f64::from(steps);
            let mut eye = start.2;
            for _ in 0..steps {
                camera.rotate([direction[0] * step, direction[1] * step]);
                let next = camera.axes().2;
                // pas de saut au pôle : l'œil avance d'un pas de rotation, pas plus
                assert!(norm(sub(next, eye)) < 1.01 * step * RADIANS_PER_PIXEL);
                eye = next;
            }
            // un tour complet ramène la vue de départ
            let end = camera.axes();
            assert!(near(end.0, start.0, 1e-9) && near(end.1, start.1, 1e-9));
            assert!(near(end.2, start.2, 1e-9));
            assert_orthonormal(&camera);
        }
    }

    #[test]
    fn dragging_moves_the_front_of_the_mesh_with_the_cursor_whatever_the_view() {
        // même la tête en bas (un demi-tour vertical), glisser vers la droite emmène la face tournée
        // vers l'observateur vers la droite, et vers le bas, vers le bas de l'écran
        let half_turn = std::f64::consts::PI / RADIANS_PER_PIXEL;
        for start in [[0.0, 0.0], [0.0, half_turn], [70.0, 130.0]] {
            for (delta, axis, sign) in [([30.0, 0.0], 0, 1.0), ([0.0, 30.0], 1, -1.0)] {
                let mut camera = Camera::fit(([-1.0; 3], [1.0; 3]), SIZE, 0.0);
                camera.rotate(start);
                let front = camera.axes().2;
                camera.rotate(delta);
                let moved =
                    camera.project(front, SIZE)[axis] - camera.project([0.0; 3], SIZE)[axis];
                assert!(sign * moved > 0.0, "{start:?} {delta:?} : {moved}");
            }
        }
    }

    #[test]
    fn look_from_gives_the_views_along_the_axes() {
        let mut camera = Camera::fit(([0.0; 3], [1.0; 3]), SIZE, 0.0);
        camera.look_from(std::f64::consts::FRAC_PI_2, 0.0);
        assert!(near(camera.axes().2, [1.0, 0.0, 0.0], 1e-12));
        // depuis +y, au pôle : la base reste définie
        camera.look_from(0.0, std::f64::consts::FRAC_PI_2);
        assert!(near(camera.axes().2, [0.0, 1.0, 0.0], 1e-12));
        assert_orthonormal(&camera);
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
    fn zooming_on_a_corner_then_turning_does_not_cut_the_mesh() {
        // le zoom au curseur éloigne la cible du centre ; après un quart de tour, cet écart est
        // dans l'axe de vue : la profondeur doit rester dans [0, 1] pour toute la boîte
        let bounds = ([0.0; 3], [40.0, 40.0, 10.0]);
        let mut camera = Camera::fit(bounds, SIZE, 0.0);
        camera.zoom_at(50.0, [10.0, 590.0], SIZE);
        camera.pan([300.0, -200.0]);
        camera.rotate([157.0, 0.0]);
        for corner in corners(bounds.0, bounds.1) {
            let depth = camera.project(corner, SIZE)[2];
            assert!(
                (0.0..=1.0).contains(&depth),
                "{corner:?} : profondeur {depth}"
            );
        }
    }

    #[test]
    fn enclosing_a_larger_box_keeps_it_in_depth_range() {
        // la déformée déborde de la boîte du cadrage, sans recadrage puisqu'on a zoomé
        let mut camera = Camera::fit(([0.0; 3], [1.0; 3]), SIZE, 0.0);
        camera.rotate([80.0, 40.0]);
        let larger = ([-2.0; 3], [3.0; 3]);
        camera.enclose(larger);
        for corner in corners(larger.0, larger.1) {
            let depth = camera.project(corner, SIZE)[2];
            assert!(
                (0.0..=1.0).contains(&depth),
                "{corner:?} : profondeur {depth}"
            );
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
