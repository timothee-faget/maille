//! Boucle d'événements winit : fenêtre, clavier (champ, déformée, colormap), souris (zoom, déplacement).

use std::error::Error;
use std::sync::Arc;

use winit::application::ApplicationHandler;
use winit::event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::keyboard::{Key, KeyCode, NamedKey, PhysicalKey};
use winit::window::{Window, WindowId};

use crate::colormap::Colormap;
use crate::legend;
use crate::render::Renderer;
use crate::scene::{range, Quantity, Scene};
use crate::view::View;

/// Touches et souris, rappelées dans le terminal à l'ouverture.
pub const HELP: &str = "Touches : 1 von Mises, 2 σxx, 3 σyy, 4 τxy, 5 |u|, D déformée, \
C colormap, R recadrer, Échap quitter. Souris : molette pour zoomer, glisser pour déplacer.";

/// Facteur de zoom par cran de molette.
const WHEEL_ZOOM: f64 = 1.2;
/// Facteur de zoom par pixel de défilement (pavé tactile).
const PIXEL_ZOOM: f64 = 1.005;

pub struct App {
    title: String,
    scene: Scene,
    quantity: Quantity,
    colormap: Colormap,
    deformed: bool,
    view: View,
    /// Vrai tant que l'utilisateur n'a ni zoomé ni déplacé : la vue se recadre avec la fenêtre.
    follow_window: bool,
    cursor: [f64; 2],
    dragging: bool,
    renderer: Option<Renderer>,
    /// Erreur survenue dans la boucle (qui ne peut pas la renvoyer), rendue par `into_result`.
    error: Option<Box<dyn Error>>,
}

impl App {
    pub fn new(title: String, scene: Scene) -> Self {
        Self {
            title,
            view: View::fit(scene.bounds(0.0), [1.0, 1.0], 0.0),
            scene,
            quantity: Quantity::VonMises,
            colormap: Colormap::Rainbow,
            deformed: false,
            follow_window: true,
            cursor: [0.0, 0.0],
            dragging: false,
            renderer: None,
            error: None,
        }
    }

    pub fn into_result(self) -> Result<(), Box<dyn Error>> {
        self.error.map_or(Ok(()), Err)
    }

    fn create_renderer(&self, event_loop: &ActiveEventLoop) -> Result<Renderer, Box<dyn Error>> {
        let window =
            event_loop.create_window(Window::default_attributes().with_title(&self.title))?;
        Renderer::new(event_loop, Arc::new(window), &self.scene, self.colormap)
    }

    fn deformation_scale(&self) -> f64 {
        if self.deformed {
            self.scene.deformation_scale()
        } else {
            0.0
        }
    }

    /// Envoie au GPU le champ, les positions, la vue et la légende courants.
    fn refresh(&mut self) {
        let values = self.scene.values(self.quantity);
        let (lo, hi) = range(&values);
        let scale = self.deformation_scale();
        let Some(renderer) = self.renderer.as_mut() else {
            return;
        };
        renderer.set_values(&values, (lo, hi));
        renderer.set_positions(&self.scene.positions(scale));
        renderer.set_view(self.view);
        let layout = legend::layout(
            renderer.size(),
            (lo, hi),
            self.colormap,
            renderer.scale_factor(),
        );
        let mut title = self.quantity.label().to_string();
        if self.deformed {
            let factor = legend::format_value(scale);
            title += &format!("\nDéformée ×{}", factor.trim_start_matches('+'));
        }
        renderer.set_legend(&layout, &title);
    }

    fn refit(&mut self) {
        if let Some(renderer) = &self.renderer {
            let reserved = legend::reserved_width(renderer.scale_factor());
            let bounds = self.scene.bounds(self.deformation_scale());
            self.view = View::fit(bounds, renderer.size(), reserved);
            self.follow_window = true;
        }
    }

    fn on_key(&mut self, event: &KeyEvent, event_loop: &ActiveEventLoop) {
        if event.state != ElementState::Pressed || event.repeat {
            return;
        }
        // chiffres par touche physique : sur un clavier AZERTY, « 1 » est sous « & » sans Maj
        let digit = match event.physical_key {
            PhysicalKey::Code(KeyCode::Digit1 | KeyCode::Numpad1) => Some(0),
            PhysicalKey::Code(KeyCode::Digit2 | KeyCode::Numpad2) => Some(1),
            PhysicalKey::Code(KeyCode::Digit3 | KeyCode::Numpad3) => Some(2),
            PhysicalKey::Code(KeyCode::Digit4 | KeyCode::Numpad4) => Some(3),
            PhysicalKey::Code(KeyCode::Digit5 | KeyCode::Numpad5) => Some(4),
            _ => None,
        };
        if let Some(k) = digit {
            self.quantity = Quantity::ALL[k];
        } else {
            match &event.logical_key {
                Key::Named(NamedKey::Escape) => return event_loop.exit(),
                Key::Character(c) => match c.to_lowercase().as_str() {
                    "d" => {
                        self.deformed = !self.deformed;
                        // la déformée déborde du cadrage du Mesh non déformé
                        if self.follow_window {
                            self.refit();
                        }
                    }
                    "c" => {
                        self.colormap = self.colormap.next();
                        if let Some(renderer) = self.renderer.as_mut() {
                            renderer.set_colormap(self.colormap);
                        }
                    }
                    "r" => self.refit(),
                    _ => return,
                },
                _ => return,
            }
        }
        self.refresh();
    }

    /// Seule la vue change : ni le champ ni la légende ne sont renvoyés au GPU.
    fn zoom(&mut self, factor: f64) {
        if let Some(renderer) = self.renderer.as_mut() {
            self.view.zoom_at(factor, self.cursor, renderer.size());
            self.follow_window = false;
            renderer.set_view(self.view);
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.renderer.is_some() {
            return;
        }
        match self.create_renderer(event_loop) {
            Ok(r) => {
                self.renderer = Some(r);
                self.refit();
                self.refresh();
            }
            Err(e) => {
                self.error = Some(e);
                event_loop.exit();
            }
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _: WindowId, event: WindowEvent) {
        let Some(renderer) = self.renderer.as_mut() else {
            return;
        };
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                renderer.resize(size.width, size.height);
                if self.follow_window {
                    self.refit();
                }
                self.refresh();
            }
            WindowEvent::RedrawRequested => renderer.render(),
            WindowEvent::KeyboardInput { event, .. } => self.on_key(&event, event_loop),
            WindowEvent::MouseWheel { delta, .. } => self.zoom(match delta {
                MouseScrollDelta::LineDelta(_, lines) => WHEEL_ZOOM.powf(f64::from(lines)),
                MouseScrollDelta::PixelDelta(p) => PIXEL_ZOOM.powf(p.y),
            }),
            WindowEvent::MouseInput {
                state,
                button: MouseButton::Left,
                ..
            } => self.dragging = state == ElementState::Pressed,
            WindowEvent::CursorMoved { position, .. } => {
                let cursor = [position.x, position.y];
                if self.dragging {
                    self.view
                        .pan([cursor[0] - self.cursor[0], cursor[1] - self.cursor[1]]);
                    self.follow_window = false;
                    renderer.set_view(self.view);
                }
                self.cursor = cursor;
            }
            _ => {}
        }
    }
}
