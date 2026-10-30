//! Viewer de Maille : affiche un Mesh 2D ou de tétraèdres et ses Fields (statiques ou thermiques)
//! dans une fenêtre (wgpu + winit), ou dans une image PNG (`screenshot_nodal`).
//!
//! Client du cœur : aucun calcul ici, seulement de l'affichage.

mod app;
mod camera;
mod colormap;
mod cut;
mod display;
mod embed;
mod legend;
mod png;
mod render;
mod scene;
mod solid;
#[cfg(test)]
mod testmesh;
mod vec3;

pub use display::{Cut, ShownSet};
pub use embed::ResultsView;
pub use render::Theme;
pub use scene::NodalResults;

use maille_core::mesh::Mesh;
use maille_core::results::Field;

/// Intervalle (min, max) du champ de rang `quantity` (même numérotation que `Shot::quantity`) : pour
/// garder la même légende sur plusieurs Increments.
pub fn quantity_range(
    results: &NodalResults,
    quantity: usize,
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let scene = scene::Scene::from_nodal(results)?;
    let quantities = scene.quantities();
    let q = quantities.get(quantity).ok_or_else(|| {
        format!(
            "champ {quantity} : ce résultat en propose {}, numérotés depuis 0",
            quantities.len()
        )
    })?;
    Ok(scene::range(&scene.values(*q)))
}

/// Ouvre une fenêtre sur le Mesh et son Field ; rend la main quand la fenêtre est fermée.
/// Erreur si le Mesh n'est pas fait de triangles ou de tétraèdres, ou si aucune carte graphique n'est
/// utilisable.
pub fn show(mesh: &Mesh, field: &Field) -> Result<(), Box<dyn std::error::Error>> {
    let title = format!(
        "Maille — {} Elements, {} Nodes, {} DOF",
        mesh.elements.len(),
        mesh.nodes.len(),
        // un seul des deux est rempli, selon l'analyse
        field.displacements.len() + field.temperatures.len()
    );
    run(title, scene::Scene::new(mesh, field)?)
}

/// Ouvre une fenêtre sur un résultat déjà lissé aux Nodes (relu d'un `.vtu`). Mêmes erreurs que
/// `show` : maillage qui n'est ni de triangles ni de tétraèdres, champs de tailles incohérentes.
pub fn show_nodal(results: &NodalResults) -> Result<(), Box<dyn std::error::Error>> {
    let title = format!(
        "Maille — {} Elements, {} Nodes",
        results.cells.len(),
        results.nodes.len()
    );
    run(title, scene::Scene::from_nodal(results)?)
}

fn run(title: String, scene: scene::Scene) -> Result<(), Box<dyn std::error::Error>> {
    println!("{}", app::HELP);
    let event_loop = winit::event_loop::EventLoop::new()?;
    let mut app = app::App::new(title, scene);
    event_loop.run_app(&mut app)?;
    app.into_result()
}

/// Réglages d'une capture d'écran (`screenshot_nodal`).
#[derive(Debug, Clone, PartialEq)]
pub struct Shot {
    /// Taille de l'image, en pixels.
    pub width: u32,
    pub height: u32,
    /// Direction de vue, en degrés : lacet autour de y et tangage ; (0, 0) regarde vers −z.
    pub yaw: f64,
    pub pitch: f64,
    /// Rang du champ dans la liste des touches 1 à 9 (0 : le premier).
    pub quantity: usize,
    pub deformed: bool,
    /// Plan de coupe : axe (0 = x, 1 = y, 2 = z) et position dans la boîte du Mesh, de 0 à 1.
    pub cut: Option<(usize, f64)>,
}

impl Default for Shot {
    fn default() -> Self {
        Self {
            width: 900,
            height: 650,
            yaw: 0.0,
            pitch: 0.0,
            quantity: 0,
            deformed: false,
            cut: None,
        }
    }
}

/// Dessine le résultat hors écran et rend le PNG de l'image, sans ouvrir de fenêtre. Mêmes erreurs
/// que `show_nodal`, plus un rang de champ inexistant et l'absence de carte graphique.
pub fn screenshot_nodal(
    results: &NodalResults,
    shot: &Shot,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let pixels = render_shot(results, shot, Theme::Light)?;
    Ok(png::encode(shot.width.max(1), shot.height.max(1), &pixels))
}

/// Pixels RGBA de la capture, ligne par ligne du haut vers le bas.
fn render_shot(
    results: &NodalResults,
    shot: &Shot,
    theme: Theme,
) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
    let mut display = display::Display::new(scene::Scene::from_nodal(results)?);
    let quantities = display.scene.quantities();
    display.quantity = *quantities.get(shot.quantity).ok_or_else(|| {
        format!(
            "champ {} : ce résultat en propose {}, numérotés depuis 0",
            shot.quantity,
            quantities.len()
        )
    })?;
    display.deformed = shot.deformed && display.scene.has_displacements();
    display.cut = shot
        .cut
        .filter(|_| display.scene.is_solid())
        .map(|(axis, fraction)| display::Cut {
            axis,
            fraction: fraction.clamp(0.0, 1.0),
            flipped: false,
        });
    display
        .camera
        .look_from(shot.yaw.to_radians(), shot.pitch.to_radians());
    let mut renderer = render::Renderer::new_offscreen(shot.width, shot.height, display.colormap)?;
    renderer.set_theme(theme);
    display.refit(&renderer);
    display.refresh(&mut renderer);
    renderer.render_pixels()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::colormap::Colormap;
    use maille_core::results::Stress;

    /// Résultat sans GPU : on saute le test sur une machine qui n'en a pas.
    fn render(results: &NodalResults, shot: &Shot) -> Option<Vec<u8>> {
        render_themed(results, shot, Theme::Light)
    }

    fn render_themed(results: &NodalResults, shot: &Shot, theme: Theme) -> Option<Vec<u8>> {
        match render_shot(results, shot, theme) {
            Ok(pixels) => Some(pixels),
            Err(e) if e.to_string().contains("aucune carte graphique") => {
                eprintln!("test sauté : {e}");
                None
            }
            Err(e) => panic!("{e}"),
        }
    }

    fn pixel(pixels: &[u8], shot: &Shot, x: u32, y: u32) -> [u8; 3] {
        let k = ((y * shot.width + x) * 4) as usize;
        [pixels[k], pixels[k + 1], pixels[k + 2]]
    }

    /// Le pixel est la couleur `expected` de la colormap, à 3 niveaux près.
    fn assert_color(got: [u8; 3], expected: [u8; 4], what: &str) {
        assert!(
            (0..3).all(|k| got[k].abs_diff(expected[k]) <= 3),
            "{what} : {got:?} au lieu de {expected:?}"
        );
    }

    /// Centre de la partie de l'image que la légende laisse libre.
    fn center(shot: &Shot) -> (u32, u32) {
        let reserved = legend::reserved_width(1.0);
        (
            (shot.width as f64 / 2.0 + reserved / 2.0) as u32,
            shot.height / 2,
        )
    }

    fn stress(p: [f64; 3]) -> Stress {
        Stress {
            x: p[0],
            y: p[1],
            z: p[2],
            sigma_xx: 0.0,
            sigma_yy: 0.0,
            sigma_zz: 0.0,
            tau_xy: 0.0,
            tau_yz: 0.0,
            tau_xz: 0.0,
        }
    }

    /// Cube de côté 2 centré, déplacement (0,001 (x + 1), 0, 0) : |u| = 0,001 (x + 1), linéaire,
    /// de 0 sur la face x = −1 à 0,002 sur la face x = +1.
    fn cube_results() -> NodalResults {
        let (nodes, cells) = crate::testmesh::cube();
        NodalResults {
            displacements: nodes
                .iter()
                .map(|p| [0.001 * (p[0] + 1.0), 0.0, 0.0])
                .collect(),
            stresses: nodes.iter().map(|&p| stress(p)).collect(),
            nodes,
            cells,
            ..NodalResults::default()
        }
    }

    #[test]
    fn a_flat_mesh_is_drawn_in_the_exact_color_of_its_field_on_the_background() {
        // un triangle à champ de température constant : toute la colormap tombe sur son bas
        let results = NodalResults {
            nodes: vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            cells: vec![vec![0, 1, 2]],
            temperatures: vec![5.0; 3],
            ..NodalResults::default()
        };
        let shot = Shot::default();
        let Some(pixels) = render(&results, &shot) else {
            return;
        };
        let (cx, cy) = center(&shot);
        // le centre de l'image libre est dans le triangle (son centre de boîte est (1, 1), sur l'hypoténuse :
        // un peu plus bas et à gauche)
        let inside = pixel(&pixels, &shot, cx - 30, cy + 30);
        assert_color(inside, Colormap::Rainbow.lut()[0], "triangle");
        // coin bas droit : fond gris clair, ni blanc ni noir
        assert_color(
            pixel(&pixels, &shot, shot.width - 1, shot.height - 1),
            [204, 204, 204, 255],
            "fond",
        );
    }

    #[test]
    fn the_dark_theme_draws_a_dark_background_and_keeps_the_field_color() {
        let results = NodalResults {
            nodes: vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            cells: vec![vec![0, 1, 2]],
            temperatures: vec![5.0; 3],
            ..NodalResults::default()
        };
        let shot = Shot::default();
        let Some(pixels) = render_themed(&results, &shot, Theme::Dark) else {
            return;
        };
        let (cx, cy) = center(&shot);
        assert_color(
            pixel(&pixels, &shot, cx - 30, cy + 30),
            Colormap::Rainbow.lut()[0],
            "triangle",
        );
        assert_color(
            pixel(&pixels, &shot, shot.width - 1, shot.height - 1),
            [37, 39, 43, 255],
            "fond sombre",
        );
    }

    #[test]
    fn a_cut_shows_the_field_inside_instead_of_the_face_behind_it() {
        let shot = Shot {
            yaw: 90.0, // regard depuis +x, vers le centre
            quantity: 7,
            ..Shot::default()
        };
        let lut = Colormap::Rainbow.lut();
        let (cx, cy) = center(&shot);
        // sans coupe : la face x = +1, où |u| est maximal
        let Some(plain) = render(&cube_results(), &shot) else {
            return;
        };
        assert_color(
            pixel(&plain, &shot, cx, cy),
            *lut.last().unwrap(),
            "face x = 1",
        );
        // coupe au quart de la boîte (x = −0,5) : la face coupée montre |u| = 0,001 · 0,5 à l'intérieur, soit
        // le quart du champ, et passe devant la face x = −1 restée derrière
        let cut = Shot {
            cut: Some((0, 0.25)),
            ..shot.clone()
        };
        let Some(sliced) = render(&cube_results(), &cut) else {
            return;
        };
        let band = (0.25 * lut.len() as f64) as usize;
        assert_color(pixel(&sliced, &cut, cx, cy), lut[band], "coupe x = -0,5");
        assert_ne!(lut[band], lut[0]);
    }

    /// Vue embarquée sur son propre device (comme le ferait l'interface graphique) ; `None` sans carte
    /// graphique.
    fn embedded(results: &NodalResults, size: [u32; 2]) -> Option<ResultsView> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&Default::default())).ok()?;
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default())).ok()?;
        Some(ResultsView::new(device, queue, results, size, 1.0).unwrap())
    }

    #[test]
    fn an_embedded_view_draws_the_same_pixels_as_the_screenshot() {
        let shot = Shot {
            yaw: 90.0,
            pitch: 0.0,
            quantity: 7,
            ..Shot::default()
        };
        let Some(expected) = render(&cube_results(), &shot) else {
            return;
        };
        let Some(mut view) = embedded(&cube_results(), [shot.width, shot.height]) else {
            return;
        };
        view.set_quantity(7);
        view.look_along(0);
        assert_eq!(view.pixels(), expected);
    }

    #[test]
    fn an_embedded_view_follows_its_size_and_its_controls() {
        let Some(mut view) = embedded(&cube_results(), [300, 200]) else {
            return;
        };
        assert_eq!(view.texture().width(), 300);
        view.resize([640, 480], 2.0);
        assert_eq!(
            (view.texture().width(), view.texture().height()),
            (640, 480)
        );
        // le champ choisi, la déformée et la coupe se règlent sans fenêtre
        assert!(view.quantities().len() > 7 && view.has_displacements() && view.is_solid());
        view.set_quantity(2);
        assert_eq!(view.quantity(), 2);
        view.set_quantity(99);
        assert_eq!(view.quantity(), 2, "un rang inexistant est ignoré");
        view.set_deformed(true);
        assert!(view.deformed());
        view.set_cut(Some(Cut {
            axis: 0,
            fraction: 7.0,
            flipped: false,
        }));
        assert_eq!(view.cut().unwrap().fraction, 1.0);
        view.set_cut(None);
        assert!(view.cut().is_none());
        view.draw();
    }

    #[test]
    fn a_highlighted_triangle_is_orange_and_the_rest_is_untouched() {
        let shot = Shot::default();
        let Some(mut view) = embedded(&cube_results(), [shot.width, shot.height]) else {
            return;
        };
        let (cx, cy) = center(&shot);
        let before = view.pixels();
        // moitié basse droite de la face avant (z = +1, vue de face)
        view.set_highlight(vec![[[-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 1.0]]]);
        let after = view.pixels();
        // hors des diagonales de la face, que les arêtes du maillage tracent en gris
        let inside = pixel(&after, &shot, cx + 30, cy + 8);
        assert!(
            inside[0] > 240 && (150..235).contains(&inside[1]) && inside[2] < 20,
            "orange attendu : {inside:?}"
        );
        // moitié haute gauche : inchangée
        let (x, y) = (cx - 30, cy - 8);
        assert_eq!(pixel(&after, &shot, x, y), pixel(&before, &shot, x, y));
        // le retirer rend l'image d'origine
        view.set_highlight(Vec::new());
        assert_eq!(view.pixels(), before);
    }

    #[test]
    fn a_set_is_translucent_and_stays_under_the_highlight() {
        let shot = Shot::default();
        let Some(mut view) = embedded(&cube_results(), [shot.width, shot.height]) else {
            return;
        };
        let (cx, cy) = center(&shot);
        let before = view.pixels();
        let half = vec![[[-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 1.0]]];
        view.set_sets(vec![([255, 0, 0], half.clone())]);
        let after = view.pixels();
        // rouge mélangé au champ (bleu ici) : ni le champ seul, ni le rouge pur
        let (field, set) = (
            pixel(&before, &shot, cx + 30, cy + 8),
            pixel(&after, &shot, cx + 30, cy + 8),
        );
        assert!(set[0] > field[0] && set[0] < 255, "{field:?} → {set:?}");
        assert!(set[2] < field[2] && set[2] > 0, "{field:?} → {set:?}");
        let (x, y) = (cx - 30, cy - 8);
        assert_eq!(pixel(&after, &shot, x, y), pixel(&before, &shot, x, y));
        // le surlignage, opaque, couvre le Set
        view.set_highlight(half.clone());
        let both = view.pixels();
        view.set_sets(Vec::new());
        assert_eq!(view.pixels(), both);
        view.set_highlight(Vec::new());
        assert_eq!(view.pixels(), before);
    }

    #[test]
    fn a_point_is_projected_where_it_is_drawn() {
        let shot = Shot::default();
        let Some(view) = embedded(&cube_results(), [shot.width, shot.height]) else {
            return;
        };
        // le centre du cube est au centre de la partie que la légende laisse libre
        let (cx, cy) = center(&shot);
        let p = view.project([0.0, 0.0, 0.0]);
        assert!(
            (p[0] - f64::from(cx)).abs() < 1.0 && (p[1] - f64::from(cy)).abs() < 1.0,
            "{p:?}"
        );
        // y monte, x va à droite, un point devant est plus proche (profondeur plus petite)
        let right = view.project([1.0, 0.0, 0.0]);
        let up = view.project([0.0, 1.0, 0.0]);
        let near = view.project([0.0, 0.0, 1.0]);
        assert!(right[0] > p[0] && up[1] < p[1]);
        assert!(near[2] < p[2]);
    }

    /// Le cube, avec des déplacements multipliés par `k` : un autre Increment du même maillage.
    fn scaled(k: f64) -> NodalResults {
        let mut r = cube_results();
        for u in &mut r.displacements {
            *u = u.map(|c| c * k);
        }
        r
    }

    #[test]
    fn another_increment_replaces_the_fields_and_keeps_the_view() {
        let shot = Shot::default();
        let Some(mut view) = embedded(&scaled(1.0), [shot.width, shot.height]) else {
            return;
        };
        view.set_quantity(7);
        view.look_along(0);
        view.set_results(&scaled(2.0)).unwrap();
        assert_eq!(view.quantity(), 7, "le champ choisi est gardé");
        // même image qu'une vue neuve sur le second Increment, réglée de la même façon
        let Some(mut fresh) = embedded(&scaled(2.0), [shot.width, shot.height]) else {
            return;
        };
        fresh.set_quantity(7);
        fresh.look_along(0);
        assert_eq!(view.pixels(), fresh.pixels());
    }

    #[test]
    fn a_close_camera_does_not_cut_the_part() {
        // pièce poussée hors du cadre puis zoom au bord de la fenêtre : la cible est loin du cube ;
        // un quart de tour met cet écart dans l'axe de vue, et le cube doit rester dessiné au centre,
        // sans coupe parallèle à l'écran
        let shot = Shot::default();
        let Some(mut view) = embedded(&cube_results(), [shot.width, shot.height]) else {
            return;
        };
        view.set_quantity(7);
        view.drag([-1000.0, 0.0], true);
        view.zoom_at(
            4.0,
            [f64::from(shot.width) - 5.0, f64::from(shot.height) / 2.0],
        );
        view.drag([157.0, 0.0], false);
        // on regarde la face x = −1, où |u| est nul
        let (cx, cy) = center(&shot);
        assert_color(
            pixel(&view.pixels(), &shot, cx, cy),
            Colormap::Rainbow.lut()[0],
            "face x = -1",
        );
    }

    #[test]
    fn a_locked_range_keeps_the_same_colours_between_increments() {
        let shot = Shot::default();
        let (lo, hi) = quantity_range(&scaled(1.0), 7).unwrap();
        let (_, hi2) = quantity_range(&scaled(2.0), 7).unwrap();
        assert!(
            lo < hi && (hi2 - 2.0 * hi).abs() < 1e-9 * hi2,
            "|u| double avec les déplacements"
        );
        let Some(mut small) = embedded(&scaled(1.0), [shot.width, shot.height]) else {
            return;
        };
        small.set_quantity(7);
        small.look_along(0);
        let own = small.pixels();
        small.set_range_lock(Some((lo, hi2)));
        assert_ne!(small.pixels(), own, "la légende suit l'intervalle imposé");
        small.set_range_lock(None);
        assert_eq!(small.pixels(), own);
        assert!(quantity_range(&scaled(1.0), 99).is_err());
    }

    #[test]
    fn the_png_of_a_view_has_the_size_of_its_texture() {
        let Some(mut view) = embedded(&cube_results(), [320, 200]) else {
            return;
        };
        let png = view.png().unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        // IHDR : largeur puis hauteur
        assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 320);
        assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 200);
    }

    #[test]
    fn node_positions_follow_the_deformation() {
        let Some(mut view) = embedded(&scaled(1.0), [300, 200]) else {
            return;
        };
        let rest = view.node_positions();
        view.set_deformed(true);
        let moved = view.node_positions();
        assert_eq!(rest.len(), moved.len());
        assert_ne!(rest, moved);
    }

    #[test]
    fn the_field_rank_is_checked() {
        let shot = Shot {
            quantity: 20,
            ..Shot::default()
        };
        let e = screenshot_nodal(&cube_results(), &shot)
            .unwrap_err()
            .to_string();
        assert!(e.contains("champ 20"), "{e}");
    }
}
