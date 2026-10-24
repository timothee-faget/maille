//! Viewer de Maille : affiche un Mesh 2D et ses Fields dans une fenêtre (wgpu + winit).
//!
//! Client du cœur : aucun calcul mécanique ici, seulement de l'affichage.

mod app;
mod colormap;
mod legend;
mod render;
mod scene;
mod view;

use maille_core::mesh::Mesh;
use maille_core::results::Field;

/// Ouvre une fenêtre sur le Mesh et son Field ; rend la main quand la fenêtre est fermée.
/// Erreur si le Mesh n'est pas fait de triangles ou si aucune carte graphique n'est utilisable.
pub fn show(mesh: &Mesh, field: &Field) -> Result<(), Box<dyn std::error::Error>> {
    let title = format!(
        "Maille — {} Elements, {} Nodes, {} DOF",
        mesh.elements.len(),
        mesh.nodes.len(),
        field.displacements.len()
    );
    let scene = scene::Scene::new(mesh, field)?;
    println!("{}", app::HELP);
    let event_loop = winit::event_loop::EventLoop::new()?;
    let mut app = app::App::new(title, scene);
    event_loop.run_app(&mut app)?;
    app.into_result()
}
