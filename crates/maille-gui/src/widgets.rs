//! Les composants partagés de l'interface : ce que plusieurs panneaux dessinent à l'identique. Ils
//! tirent leurs couleurs et leurs mesures de `theme`, comme tout le reste.

use eframe::egui::{
    self, Align, Align2, Button, Color32, CornerRadius, Layout, Response, RichText, Role, Sense,
    TextStyle, Ui, Vec2, WidgetInfo,
};
use egui_lucide::Lucide;

use crate::state::Status;
use crate::theme::{self, Palette};

const ICON_SIZE: f32 = 16.0;
/// Marge du texte d'une ligne à gauche, et diamètre de la pastille d'état.
const ROW_PADDING: f32 = 8.0;
const DOT: f32 = 10.0;

/// L'en-tête d'un panneau : son titre, un bouton qui le ferme ou le masque (rend `true` au clic),
/// puis un filet.
pub fn panel_header(ui: &mut Ui, title: &str, icon: Lucide, tip: &str) -> bool {
    let mut clicked = false;
    ui.horizontal(|ui| {
        ui.heading(title);
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            let close = Button::image(
                icon.size(ICON_SIZE)
                    .color(ui.visuals().weak_text_color())
                    .stroke_width(1.75)
                    .image(),
            )
            .frame(false);
            clicked = ui.add(close).on_hover_text(tip).clicked();
        });
    });
    ui.separator();
    clicked
}

/// Une ligne de choix pleine largeur dans une liste verticale : fond d'accent si elle est choisie,
/// fond discret au survol.
pub fn choice_row(ui: &mut Ui, selected: bool, label: &str) -> Response {
    row(ui, selected, label, None)
}

/// Une étape d'une liste de contrôle : sa pastille d'état puis son libellé, sur une ligne pleine
/// largeur qui se survole.
pub fn step_row(ui: &mut Ui, status: Status, label: &str) -> Response {
    row(ui, false, label, Some(status))
}

fn row(ui: &mut Ui, selected: bool, label: &str, status: Option<Status>) -> Response {
    let p = Palette::of(ui.visuals());
    let size = Vec2::new(ui.available_width(), ui.spacing().interact_size.y);
    let (rect, response) = ui.allocate_exact_size(size, Sense::click());
    response.widget_info(|| WidgetInfo::selected(Role::Button, ui.is_enabled(), selected, label));
    let fill = if selected {
        p.accent
    } else if response.hovered() {
        p.control_hover
    } else {
        Color32::TRANSPARENT
    };
    let painter = ui.painter();
    painter.rect_filled(rect, CornerRadius::same(theme::CONTROL_RADIUS), fill);
    let mut x = rect.left() + ROW_PADDING;
    if let Some(status) = status {
        dot(
            ui,
            egui::Rect::from_center_size(
                egui::pos2(x + DOT / 2.0, rect.center().y),
                Vec2::splat(DOT),
            ),
            status,
        );
        x += DOT + ROW_PADDING;
    }
    let color = if selected { p.on_accent } else { p.text };
    ui.painter().text(
        egui::pos2(x, rect.center().y),
        Align2::LEFT_CENTER,
        label,
        TextStyle::Body.resolve(ui.style()),
        color,
    );
    response.on_hover_cursor(egui::CursorIcon::PointingHand)
}

/// La pastille d'état d'un espace ou d'une étape. « Pas encore fait » est un anneau neutre : le
/// rouge reste pour les erreurs.
pub fn status_dot(ui: &mut Ui, status: Status) {
    let (rect, _) = ui.allocate_exact_size(Vec2::splat(DOT), Sense::hover());
    dot(ui, rect, status);
}

fn dot(ui: &Ui, rect: egui::Rect, status: Status) {
    let p = Palette::of(ui.visuals());
    let painter = ui.painter();
    match status {
        Status::Done => {
            painter.circle_filled(rect.center(), 4.5, p.ok);
        }
        Status::Warning => {
            painter.circle_filled(rect.center(), 4.5, p.warn);
        }
        Status::Missing => {
            painter.circle_stroke(rect.center(), 3.75, egui::Stroke::new(1.5, p.pending));
        }
    }
}

/// Le bouton de l'action principale d'un panneau : plein d'accent.
pub fn primary_button(ui: &mut Ui, text: &str) -> Response {
    let p = Palette::of(ui.visuals());
    ui.add(Button::new(RichText::new(text).color(p.on_accent)).fill(p.accent))
}

/// Comme `primary_button`, désactivable (le calcul est bloqué, ou déjà en cours).
pub fn primary_button_enabled(ui: &mut Ui, enabled: bool, text: &str) -> Response {
    let p = Palette::of(ui.visuals());
    ui.add_enabled(
        enabled,
        Button::new(RichText::new(text).color(p.on_accent)).fill(p.accent),
    )
}
