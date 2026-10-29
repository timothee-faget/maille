//! Le viewer embarqué : les résultats dessinés dans une texture d'un device que possède un autre
//! client (l'interface graphique, qui l'affiche comme une image). Même rendu, mêmes champs, même
//! caméra que la fenêtre de `maille view` ; les entrées (souris, clavier) sont à la charge du client,
//! qui appelle ces méthodes.

use std::error::Error;

use crate::display::{Cut, Display};
use crate::render::Renderer;
use crate::scene::{NodalResults, Scene};

/// Lacet et tangage des vues le long des axes x, y, z, en degrés (voir `Camera::axes`) : on regarde
/// depuis +x, depuis +y (tangage limité à 89°, la caméra ne passe pas le pôle) et depuis +z.
const VIEWS: [(f64, f64); 3] = [(90.0, 0.0), (0.0, 89.0), (0.0, 0.0)];

pub struct ResultsView {
    display: Display,
    renderer: Renderer,
}

impl ResultsView {
    /// Prépare la vue sur `results`, dans une texture de `size` pixels (`scale` : pixels physiques par
    /// pixel logique). Erreur si le maillage n'est ni de triangles ni de tétraèdres, ou si les champs
    /// ont des tailles incohérentes.
    pub fn new(
        device: wgpu::Device,
        queue: wgpu::Queue,
        results: &NodalResults,
        size: [u32; 2],
        scale: f64,
    ) -> Result<Self, Box<dyn Error>> {
        let display = Display::new(Scene::from_nodal(results)?);
        let renderer =
            Renderer::new_shared(device, queue, size[0], size[1], scale, display.colormap);
        let mut view = Self { display, renderer };
        view.fit();
        Ok(view)
    }

    /// La texture à afficher (format `Rgba8UnormSrgb`) : à relire après chaque `resize`, elle change.
    pub fn texture(&self) -> &wgpu::Texture {
        self.renderer
            .texture()
            .expect("une vue embarquée dessine toujours dans une texture")
    }

    /// Dessine l'image courante dans la texture.
    pub fn draw(&mut self) {
        self.renderer.draw();
    }

    /// Change la taille de la texture (en pixels) ; la vue se recadre tant qu'on n'a ni zoomé ni
    /// déplacé.
    pub fn resize(&mut self, size: [u32; 2], scale: f64) {
        self.renderer.set_scale_factor(scale);
        self.renderer.resize(size[0], size[1]);
        if self.display.follow_window {
            self.display.refit(&self.renderer);
        }
        self.refresh();
    }

    /// Déplacement du curseur de `delta` pixels, bouton enfoncé : en 3D la pièce tourne, sauf si `pan`
    /// (ou en 2D) : elle suit la souris.
    pub fn drag(&mut self, delta: [f64; 2], pan: bool) {
        if self.display.scene.is_solid() && !pan {
            self.display.camera.rotate(delta);
        } else {
            self.display.camera.pan(delta);
        }
        self.moved();
    }

    /// Zoom de `factor` autour du pixel `cursor` (origine en haut à gauche de la texture).
    pub fn zoom_at(&mut self, factor: f64, cursor: [f64; 2]) {
        self.display
            .camera
            .zoom_at(factor, cursor, self.renderer.size());
        self.moved();
    }

    fn moved(&mut self) {
        self.display.follow_window = false;
        self.renderer.set_camera(&self.display.camera);
    }

    /// Surligne ces triangles (sommets dans les coordonnées de la pièce), par-dessus le maillage ;
    /// une liste vide retire le surlignage.
    pub fn set_highlight(&mut self, triangles: Vec<[[f64; 3]; 3]>) {
        self.display.highlight = triangles;
        self.refresh();
    }

    /// Où le point `p` de la pièce tombe dans la texture : (x, y) en pixels depuis le coin haut gauche,
    /// et la profondeur (0 : devant, 1 : derrière). Sert à désigner ce qu'il y a sous le curseur.
    pub fn project(&self, p: [f64; 3]) -> [f64; 3] {
        let size = self.renderer.size();
        let clip = self.display.camera.project(p, size);
        [
            (clip[0] + 1.0) * 0.5 * size[0],
            (1.0 - clip[1]) * 0.5 * size[1],
            clip[2],
        ]
    }

    /// Base de la caméra dans les coordonnées de la pièce : (droite, haut, œil), où `œil` va de la pièce
    /// vers l'observateur. Sert à dessiner le trièdre d'orientation.
    pub fn orientation(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        self.display.camera.axes()
    }

    /// Recadre la pièce, sans changer la direction de vue.
    pub fn fit(&mut self) {
        self.display.refit(&self.renderer);
        self.refresh();
    }

    /// Regarde la pièce le long de l'axe `axis` (0 = x, 1 = y, 2 = z) et la recadre. Sans effet sur
    /// un maillage 2D, qu'on ne voit que de face.
    pub fn look_along(&mut self, axis: usize) {
        if !self.display.scene.is_solid() {
            return;
        }
        let (yaw, pitch) = VIEWS[axis.min(2)];
        self.display.camera.yaw = yaw.to_radians();
        self.display.camera.pitch = pitch.to_radians();
        self.fit();
    }

    fn refresh(&mut self) {
        self.display.refresh(&mut self.renderer);
    }

    /// Noms des champs affichables, dans l'ordre des touches 1 à 9 de `maille view`.
    pub fn quantities(&self) -> Vec<&'static str> {
        self.display
            .scene
            .quantities()
            .iter()
            .map(|q| q.label())
            .collect()
    }

    /// Rang du champ affiché dans `quantities`.
    pub fn quantity(&self) -> usize {
        self.display
            .scene
            .quantities()
            .iter()
            .position(|&q| q == self.display.quantity)
            .unwrap_or(0)
    }

    /// Affiche le champ de rang `index` de `quantities` ; sans effet s'il n'existe pas.
    pub fn set_quantity(&mut self, index: usize) {
        if let Some(&q) = self.display.scene.quantities().get(index) {
            self.display.quantity = q;
            self.refresh();
        }
    }

    pub fn has_displacements(&self) -> bool {
        self.display.scene.has_displacements()
    }

    pub fn deformed(&self) -> bool {
        self.display.deformed
    }

    /// Montre la déformée (sans effet s'il n'y a pas de déplacements) et recadre sur elle.
    pub fn set_deformed(&mut self, deformed: bool) {
        self.display.deformed = deformed && self.has_displacements();
        if self.display.follow_window {
            self.display.refit(&self.renderer);
        }
        self.refresh();
    }

    pub fn is_solid(&self) -> bool {
        self.display.scene.is_solid()
    }

    pub fn cut(&self) -> Option<Cut> {
        self.display.cut
    }

    /// Plan de coupe (3D seulement) ; `None` le retire.
    pub fn set_cut(&mut self, cut: Option<Cut>) {
        self.display.cut = cut.filter(|_| self.is_solid()).map(|c| Cut {
            fraction: c.fraction.clamp(0.0, 1.0),
            ..c
        });
        self.refresh();
    }

    /// Remplace les champs par ceux d'un autre Increment du même maillage (animation) : le champ
    /// affiché, la déformée, la coupe, la colormap et la caméra sont gardés. Mêmes erreurs que `new`.
    pub fn set_results(&mut self, results: &NodalResults) -> Result<(), Box<dyn Error>> {
        let scene = Scene::from_nodal(results)?;
        self.display.scene = scene;
        let quantities = self.display.scene.quantities();
        if !quantities.contains(&self.display.quantity) {
            self.display.quantity = quantities[0];
        }
        if self.display.follow_window {
            self.display.refit(&self.renderer);
        }
        self.refresh();
        Ok(())
    }

    /// Impose l'intervalle de la légende (`None` : celui des valeurs affichées).
    pub fn set_range_lock(&mut self, range: Option<(f64, f64)>) {
        self.display.range_lock = range;
        self.refresh();
    }

    /// Position affichée de chaque Node (déformée comprise si on la montre), dans les coordonnées de
    /// la pièce : ce que `project` attend.
    pub fn node_positions(&self) -> Vec<[f64; 3]> {
        self.display
            .scene
            .positions(self.display.deformation_scale())
    }

    /// L'image courante en PNG : les mêmes pixels que `maille view --png` à cette taille.
    pub fn png(&mut self) -> Result<Vec<u8>, Box<dyn Error>> {
        let pixels = self.renderer.render_pixels()?;
        let [width, height] = self.renderer.size();
        Ok(crate::png::encode(
            (width.round() as u32).max(1),
            (height.round() as u32).max(1),
            &pixels,
        ))
    }

    /// Passe à la colormap suivante.
    pub fn next_colormap(&mut self) {
        self.display.colormap = self.display.colormap.next();
        self.renderer.set_colormap(self.display.colormap);
        self.refresh();
    }

    /// Pixels RGBA de l'image courante (tests).
    #[cfg(test)]
    pub(crate) fn pixels(&mut self) -> Vec<u8> {
        self.renderer
            .render_pixels()
            .expect("lecture de la texture")
    }
}
