//! Boucle d'événements winit : fenêtre, clavier (champ, déformée, colormap, coupe), souris
//! (rotation en 3D, zoom, déplacement).

use std::error::Error;
use std::sync::Arc;

use winit::application::ApplicationHandler;
use winit::event::{ElementState, KeyEvent, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::keyboard::{Key, KeyCode, NamedKey, PhysicalKey};
use winit::window::{Window, WindowId};

use crate::display::{Cut, Display};
use crate::render::Renderer;
use crate::scene::Scene;

/// Touches et souris, rappelées dans le terminal à l'ouverture.
pub const HELP: &str =
    "Touches : 1 à 9 champ (von Mises, contraintes, |u|, p si plastique ; thermique : 1 T), \
D déformée, C colormap, R recadrer, Échap quitter. Souris : molette pour zoomer, glisser pour déplacer \
(en 3D : glisser tourne, Maj + glisser ou clic droit déplace). \
Mesh 3D : X, Y, Z coupent selon l'axe (même touche : retire la coupe), ↑ et ↓ déplacent le plan, F change de côté.";

/// Facteur de zoom par cran de molette.
const WHEEL_ZOOM: f64 = 1.2;
/// Facteur de zoom par pixel de défilement (pavé tactile).
const PIXEL_ZOOM: f64 = 1.005;
/// Déplacement du plan de coupe par appui, en fraction de la boîte du Mesh.
const CUT_STEP: f64 = 0.02;

pub struct App {
    title: String,
    display: Display,
    cursor: [f64; 2],
    dragging: Option<MouseButton>,
    shift: bool,
    renderer: Option<Renderer>,
    /// Erreur survenue dans la boucle (qui ne peut pas la renvoyer), rendue par `into_result`.
    error: Option<Box<dyn Error>>,
}

impl App {
    pub fn new(title: String, scene: Scene) -> Self {
        Self {
            title,
            display: Display::new(scene),
            cursor: [0.0, 0.0],
            dragging: None,
            shift: false,
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
        Renderer::new_window(event_loop, Arc::new(window), self.display.colormap)
    }

    fn refresh(&mut self) {
        if let Some(renderer) = self.renderer.as_mut() {
            self.display.refresh(renderer);
        }
    }

    fn refit(&mut self) {
        if let Some(renderer) = &self.renderer {
            self.display.refit(renderer);
        }
    }

    fn on_key(&mut self, event: &KeyEvent, event_loop: &ActiveEventLoop) {
        if event.state != ElementState::Pressed {
            return;
        }
        // flèches : le plan de coupe peut se déplacer en maintenant la touche
        if let PhysicalKey::Code(code @ (KeyCode::ArrowUp | KeyCode::ArrowDown)) =
            event.physical_key
        {
            if let Some(cut) = self.display.cut.as_mut() {
                let step = if code == KeyCode::ArrowUp {
                    CUT_STEP
                } else {
                    -CUT_STEP
                };
                cut.fraction = (cut.fraction + step).clamp(0.0, 1.0);
                self.refresh();
            }
            return;
        }
        if event.repeat {
            return;
        }
        // chiffres par touche physique : sur un clavier AZERTY, « 1 » est sous « & » sans Maj
        let digit = match event.physical_key {
            PhysicalKey::Code(KeyCode::Digit1 | KeyCode::Numpad1) => Some(0),
            PhysicalKey::Code(KeyCode::Digit2 | KeyCode::Numpad2) => Some(1),
            PhysicalKey::Code(KeyCode::Digit3 | KeyCode::Numpad3) => Some(2),
            PhysicalKey::Code(KeyCode::Digit4 | KeyCode::Numpad4) => Some(3),
            PhysicalKey::Code(KeyCode::Digit5 | KeyCode::Numpad5) => Some(4),
            PhysicalKey::Code(KeyCode::Digit6 | KeyCode::Numpad6) => Some(5),
            PhysicalKey::Code(KeyCode::Digit7 | KeyCode::Numpad7) => Some(6),
            PhysicalKey::Code(KeyCode::Digit8 | KeyCode::Numpad8) => Some(7),
            PhysicalKey::Code(KeyCode::Digit9 | KeyCode::Numpad9) => Some(8),
            _ => None,
        };
        if let Some(k) = digit {
            match self.display.scene.quantities().get(k) {
                Some(&q) => self.display.quantity = q,
                None => return,
            }
        } else {
            match &event.logical_key {
                Key::Named(NamedKey::Escape) => return event_loop.exit(),
                Key::Character(c) => match c.to_lowercase().as_str() {
                    "d" if self.display.scene.has_displacements() => {
                        self.display.deformed = !self.display.deformed;
                        // la déformée déborde du cadrage du Mesh non déformé
                        if self.display.follow_window {
                            self.refit();
                        }
                    }
                    "c" => {
                        self.display.colormap = self.display.colormap.next();
                        if let Some(renderer) = self.renderer.as_mut() {
                            renderer.set_colormap(self.display.colormap);
                        }
                    }
                    "r" => self.refit(),
                    axis @ ("x" | "y" | "z") if self.display.scene.is_solid() => {
                        let axis = "xyz".find(axis).unwrap_or(0);
                        self.display.cut = match self.display.cut {
                            Some(cut) if cut.axis == axis => None,
                            _ => Some(Cut {
                                axis,
                                fraction: 0.5,
                                flipped: false,
                            }),
                        };
                    }
                    "f" => match self.display.cut.as_mut() {
                        Some(cut) => cut.flipped = !cut.flipped,
                        None => return,
                    },
                    _ => return,
                },
                _ => return,
            }
        }
        self.refresh();
    }

    fn zoom(&mut self, factor: f64) {
        if let Some(renderer) = self.renderer.as_mut() {
            self.display
                .camera
                .zoom_at(factor, self.cursor, renderer.size());
            self.display.follow_window = false;
            renderer.set_camera(&self.display.camera);
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
                if self.display.follow_window {
                    self.refit();
                }
                self.refresh();
            }
            WindowEvent::RedrawRequested => renderer.render(),
            WindowEvent::KeyboardInput { event, .. } => self.on_key(&event, event_loop),
            WindowEvent::ModifiersChanged(modifiers) => self.shift = modifiers.state().shift_key(),
            WindowEvent::MouseWheel { delta, .. } => self.zoom(match delta {
                MouseScrollDelta::LineDelta(_, lines) => WHEEL_ZOOM.powf(f64::from(lines)),
                MouseScrollDelta::PixelDelta(p) => PIXEL_ZOOM.powf(p.y),
            }),
            WindowEvent::MouseInput {
                state,
                button: button @ (MouseButton::Left | MouseButton::Right),
                ..
            } => self.dragging = (state == ElementState::Pressed).then_some(button),
            WindowEvent::CursorMoved { position, .. } => {
                let cursor = [position.x, position.y];
                if let Some(button) = self.dragging {
                    let delta = [cursor[0] - self.cursor[0], cursor[1] - self.cursor[1]];
                    // en 3D le clic gauche tourne, le droit (ou Maj) déplace ; en 2D il déplace toujours
                    if self.display.scene.is_solid() && button == MouseButton::Left && !self.shift {
                        self.display.camera.rotate(delta);
                    } else {
                        self.display.camera.pan(delta);
                    }
                    self.display.follow_window = false;
                    renderer.set_camera(&self.display.camera);
                }
                self.cursor = cursor;
            }
            _ => {}
        }
    }
}
