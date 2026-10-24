//! Boucle d'événements winit : création de la fenêtre, redimensionnement, fermeture.

use std::error::Error;
use std::sync::Arc;

use winit::application::ApplicationHandler;
use winit::event::WindowEvent;
use winit::event_loop::ActiveEventLoop;
use winit::window::{Window, WindowId};

use crate::render::Renderer;

pub struct App {
    title: String,
    renderer: Option<Renderer>,
    /// Erreur survenue dans la boucle (qui ne peut pas la renvoyer), rendue par `into_result`.
    error: Option<Box<dyn Error>>,
}

impl App {
    pub fn new(title: String) -> Self {
        Self {
            title,
            renderer: None,
            error: None,
        }
    }

    pub fn into_result(self) -> Result<(), Box<dyn Error>> {
        self.error.map_or(Ok(()), Err)
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.renderer.is_some() {
            return;
        }
        let attributes = Window::default_attributes().with_title(&self.title);
        let renderer = event_loop
            .create_window(attributes)
            .map_err(Into::into)
            .and_then(|window| Renderer::new(event_loop, Arc::new(window)));
        match renderer {
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
            WindowEvent::Resized(size) => renderer.resize(size.width, size.height),
            WindowEvent::RedrawRequested => renderer.render(),
            _ => {}
        }
    }
}
