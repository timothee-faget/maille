//! La légende des champs, dessinée par l'interface dans une carte comme les autres panneaux : le titre
//! (champ, déformée, coupe), la barre de couleurs et ses graduations. Le viewer fournit le contenu
//! (`ResultsView::legend`) et ne la dessine plus dans sa texture.

use eframe::egui::{
    self, Align2, Color32, CornerRadius, FontId, Pos2, Rect, RichText, Sense, Stroke, StrokeKind,
    Vec2,
};
use maille_viewer::Legend;

use crate::theme::{self, Palette};

/// Largeur du contenu de la carte, et de la barre de couleurs.
const CONTENT_WIDTH: f32 = 136.0;
const BAR_WIDTH: f32 = 14.0;
const TICK_LENGTH: f32 = 5.0;
const LABEL_GAP: f32 = 4.0;
/// Hauteur des étiquettes : la barre est décalée d'une demi-ligne de chaque bout pour que
/// la première et la dernière tiennent dans la carte.
const LABEL_HEIGHT: f32 = 14.0;

/// La carte de la légende, collée en haut à gauche ; la barre fait au plus `max_bar` de haut.
pub fn card(ctx: &egui::Context, legend: &Legend, max_bar: f32, offset: [f32; 2]) {
    egui::Area::new(egui::Id::unique("legend"))
        .anchor(Align2::LEFT_TOP, offset)
        .order(egui::Order::Foreground)
        .interactable(false)
        .show(ctx, |ui| {
            theme::panel_frame(ui.visuals()).show(ui, |ui| {
                ui.set_width(CONTENT_WIDTH);
                let palette = Palette::of(ui.visuals());
                for (k, line) in legend.title.lines().enumerate() {
                    let text = if k == 0 {
                        RichText::new(line).color(palette.text)
                    } else {
                        RichText::new(line).small().color(palette.text_weak)
                    };
                    ui.add(egui::Label::new(text).wrap());
                }
                ui.add_space(theme::SPACE_M);
                bar(ui, legend, max_bar, &palette);
            });
        });
}

fn bar(ui: &mut egui::Ui, legend: &Legend, max_bar: f32, palette: &Palette) {
    let height = max_bar.max(60.0);
    let (rect, _) = ui.allocate_exact_size(
        Vec2::new(CONTENT_WIDTH, height + LABEL_HEIGHT),
        Sense::hover(),
    );
    let bar = Rect::from_min_size(
        Pos2::new(rect.left(), rect.top() + LABEL_HEIGHT / 2.0),
        Vec2::new(BAR_WIDTH, height),
    );
    let painter = ui.painter();
    // une bande par couleur de la table, du minimum (en bas) au maximum (en haut) ; la bande
    // déborde d'un demi-pixel sur sa voisine pour ne laisser aucun filet entre elles
    let n = legend.colors.len().max(1) as f32;
    for (k, [r, g, b]) in legend.colors.iter().enumerate() {
        let (low, high) = (k as f32 / n, (k + 1) as f32 / n);
        let band = Rect::from_min_max(
            Pos2::new(bar.left(), bar.bottom() - high * bar.height() - 0.5),
            Pos2::new(bar.right(), bar.bottom() - low * bar.height()),
        );
        painter.rect_filled(band, CornerRadius::ZERO, Color32::from_rgb(*r, *g, *b));
    }
    painter.rect_stroke(
        bar,
        CornerRadius::same(2),
        Stroke::new(1.0, palette.panel_stroke),
        StrokeKind::Inside,
    );
    let label = FontId::monospace(12.0);
    for tick in &legend.ticks {
        let y = bar.bottom() - tick.fraction as f32 * bar.height();
        painter.line_segment(
            [
                Pos2::new(bar.right(), y),
                Pos2::new(bar.right() + TICK_LENGTH, y),
            ],
            Stroke::new(1.0, palette.text_weak),
        );
        painter.text(
            Pos2::new(bar.right() + TICK_LENGTH + LABEL_GAP, y),
            Align2::LEFT_CENTER,
            &tick.label,
            label.clone(),
            palette.text,
        );
    }
}
