//! Boucle d'événements winit : création de la fenêtre, redimensionnement, fermeture.

use std::error::Error;
use std::sync::Arc;

use winit::application::ApplicationHandler;
use winit::event::WindowEvent;
use winit::event_loop::ActiveEventLoop;
use winit::window::{Window, WindowId};

use crate::colormap::Colormap;
use crate::render::Renderer;
use crate::scene::{range, Quantity, Scene};
use crate::view::View;

pub struct App {
    title: String,
    scene: Scene,
    renderer: Option<Renderer>,
    /// Erreur survenue dans la boucle (qui ne peut pas la renvoyer), rendue par `into_result`.
    error: Option<Box<dyn Error>>,
}

impl App {
    pub fn new(title: String, scene: Scene) -> Self {
        Self {
            title,
            scene,
            renderer: None,
            error: None,
        }
    }

    pub fn into_result(self) -> Result<(), Box<dyn Error>> {
        self.error.map_or(Ok(()), Err)
    }

    fn create_renderer(&self, event_loop: &ActiveEventLoop) -> Result<Renderer, Box<dyn Error>> {
        let quantity = Quantity::VonMises;
        let values = self.scene.values(quantity);
        let (lo, hi) = range(&values);
        let title = format!(
            "{} — {} de {lo:.4e} à {hi:.4e}",
            self.title,
            quantity.label()
        );
        let window = event_loop.create_window(Window::default_attributes().with_title(title))?;
        let mut renderer =
            Renderer::new(event_loop, Arc::new(window), &self.scene, Colormap::Rainbow)?;
        renderer.set_positions(&self.scene.positions(0.0));
        renderer.set_values(&values, (lo, hi));
        renderer.set_view(View::fit(self.scene.bounds(), renderer.size()));
        Ok(renderer)
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.renderer.is_some() {
            return;
        }
        match self.create_renderer(event_loop) {
            Ok(r) => self.renderer = Some(r),
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
                renderer.set_view(View::fit(self.scene.bounds(), renderer.size()));
            }
            WindowEvent::RedrawRequested => renderer.render(),
            _ => {}
        }
    }
}
