//! Viewer de Maille : affiche un Mesh 2D et ses Fields dans une fenêtre (wgpu + winit).
//!
//! Client du cœur : aucun calcul mécanique ici, seulement de l'affichage.

mod app;
mod render;

use maille_core::mesh::Mesh;
use maille_core::results::Field;

/// Ouvre une fenêtre sur le Mesh et son Field ; rend la main quand la fenêtre est fermée.
pub fn show(mesh: &Mesh, field: &Field) -> Result<(), Box<dyn std::error::Error>> {
    let title = format!(
        "Maille — {} Elements, {} Nodes, {} DOF",
        mesh.elements.len(),
        mesh.nodes.len(),
        field.displacements.len()
    );
    let event_loop = winit::event_loop::EventLoop::new()?;
    let mut app = app::App::new(title);
    event_loop.run_app(&mut app)?;
    app.into_result()
}
