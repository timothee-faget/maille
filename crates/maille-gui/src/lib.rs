//! Interface graphique de Maille (M13) : trois espaces (Maillage, Modèle, Résultats), une vue 3D
//! plein écran et des panneaux flottants par-dessus. Client du cœur, de `maille-io` et de
//! `maille-viewer`. L'état de l'interface (`state`) est pur et testé ; `app` ne fait que l'afficher.

mod app;
pub mod compute;
pub mod mesh_space;
pub mod model_space;
mod model_tools;
pub mod quality;
pub mod results;
mod results_ui;
pub mod series;
pub mod sets;
pub mod state;
mod triad;

use std::error::Error;
use std::path::PathBuf;

/// Ouvre la fenêtre de Maille ; rend la main quand elle est fermée. `open` : un `.pvd` ou `.vtu` à
/// montrer d'emblée dans l'espace Résultats.
pub fn run(open: Option<PathBuf>) -> Result<(), Box<dyn Error>> {
    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: eframe::egui::ViewportBuilder::default()
            .with_title("Maille")
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Maille",
        options,
        Box::new(move |cc| {
            // les icônes sont des SVG
            egui_extras::install_image_loaders(&cc.egui_ctx);
            Ok(Box::new(app::App::new(open)))
        }),
    )
    .map_err(|e| format!("l'interface n'a pas pu démarrer : {e}").into())
}
