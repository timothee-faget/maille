//! Le thème de l'interface : les jetons (couleurs, espacements, rayons, tailles de texte) et leur
//! application à egui. Les autres modules n'écrivent aucune couleur ni aucune mesure de style en dur :
//! ils demandent à la `Palette` (`Palette::of`) ou prennent une constante d'ici. Un test le vérifie.
//!
//! Deux palettes, sombre et claire, à accent unique (corail). Le rouge d'`error` ne sert qu'aux vraies
//! erreurs : un état « pas encore fait » est neutre (`Palette::pending`).

use eframe::egui::{
    self, Color32, CornerRadius, FontFamily, FontId, Margin, Shadow, Stroke, Style, TextStyle,
    Vec2, Visuals,
};

/// Rayon des panneaux flottants, des cartes et des contrôles.
pub const PANEL_RADIUS: u8 = 14;
pub const CONTROL_RADIUS: u8 = 8;
/// Marge intérieure d'un panneau flottant.
pub const PANEL_MARGIN: i8 = 12;

/// Échelle d'espacement (px) : 4, 8, 12, 16.
pub const SPACE_S: f32 = 4.0;
pub const SPACE_M: f32 = 8.0;
pub const SPACE_L: f32 = 12.0;

/// Échelle typographique : titre de panneau, corps, détail, valeurs en chasse fixe.
const SIZE_HEADING: f32 = 15.0;
const SIZE_BODY: f32 = 13.5;
const SIZE_SMALL: f32 = 11.5;
const SIZE_MONO: f32 = 12.0;

/// Les couleurs d'un thème. Les fonds translucides sont en alpha non prémultiplié.
#[derive(Clone, Copy, Debug)]
pub struct Palette {
    pub dark: bool,
    /// Fond de la vue, quand il n'y a rien à dessiner (même valeur que le fond du viewer).
    pub canvas: Color32,
    /// Panneau flottant (opaque : la vue ne doit pas se lire à travers), son contour et son ombre.
    pub panel: Color32,
    pub panel_stroke: Color32,
    pub shadow: Color32,
    /// Contrôle (bouton, champ) au repos, survolé, enfoncé ; fond d'un champ de saisie.
    pub control: Color32,
    pub control_hover: Color32,
    pub control_active: Color32,
    pub input: Color32,
    pub text: Color32,
    pub text_weak: Color32,
    /// L'accent : outil, onglet ou ligne active, curseur de slider.
    pub accent: Color32,
    /// Texte posé sur l'accent.
    pub on_accent: Color32,
    pub ok: Color32,
    pub warn: Color32,
    pub error: Color32,
    /// Un état pas encore fait : neutre, ce n'est pas une erreur.
    pub pending: Color32,
    /// Axes x, y, z du trièdre et des boutons de vue.
    pub axes: [Color32; 3],
    /// Courbe d'un graphe.
    pub plot: Color32,
}

impl Palette {
    pub fn of(visuals: &Visuals) -> Self {
        palette(visuals.dark_mode)
    }

    /// Un fond de carte ou de graphe, en retrait du panneau.
    pub fn inset(&self) -> Color32 {
        if self.dark {
            Color32::from_black_alpha(90)
        } else {
            Color32::from_black_alpha(24)
        }
    }
}

pub fn palette(dark: bool) -> Palette {
    if dark {
        Palette {
            dark,
            canvas: Color32::from_rgb(37, 39, 43),
            panel: Color32::from_rgb(30, 31, 36),
            panel_stroke: Color32::from_white_alpha(22),
            shadow: Color32::from_black_alpha(110),
            control: Color32::from_white_alpha(14),
            control_hover: Color32::from_white_alpha(28),
            control_active: Color32::from_white_alpha(40),
            input: Color32::from_black_alpha(80),
            text: Color32::from_rgb(236, 237, 240),
            text_weak: Color32::from_rgb(150, 153, 162),
            accent: Color32::from_rgb(242, 96, 112),
            on_accent: Color32::WHITE,
            ok: Color32::from_rgb(88, 196, 126),
            warn: Color32::from_rgb(236, 178, 70),
            error: Color32::from_rgb(255, 92, 82),
            pending: Color32::from_rgb(110, 113, 122),
            axes: [
                Color32::from_rgb(226, 80, 80),
                Color32::from_rgb(90, 190, 100),
                Color32::from_rgb(80, 130, 235),
            ],
            plot: Color32::from_rgb(120, 190, 255),
        }
    } else {
        Palette {
            dark,
            canvas: Color32::from_rgb(204, 204, 204),
            panel: Color32::from_rgb(247, 247, 249),
            panel_stroke: Color32::from_black_alpha(30),
            shadow: Color32::from_black_alpha(48),
            control: Color32::from_black_alpha(10),
            control_hover: Color32::from_black_alpha(20),
            control_active: Color32::from_black_alpha(32),
            input: Color32::from_white_alpha(210),
            text: Color32::from_rgb(28, 29, 33),
            text_weak: Color32::from_rgb(108, 111, 120),
            accent: Color32::from_rgb(222, 58, 82),
            on_accent: Color32::WHITE,
            ok: Color32::from_rgb(34, 150, 84),
            warn: Color32::from_rgb(190, 120, 10),
            error: Color32::from_rgb(200, 36, 36),
            pending: Color32::from_rgb(160, 163, 172),
            axes: [
                Color32::from_rgb(214, 60, 60),
                Color32::from_rgb(52, 160, 70),
                Color32::from_rgb(50, 100, 220),
            ],
            plot: Color32::from_rgb(20, 90, 190),
        }
    }
}

/// Le cadre commun des panneaux flottants.
pub fn panel_frame(visuals: &Visuals) -> egui::Frame {
    let p = Palette::of(visuals);
    egui::Frame::new()
        .fill(p.panel)
        .corner_radius(CornerRadius::same(PANEL_RADIUS))
        .inner_margin(Margin::same(PANEL_MARGIN))
        .stroke(Stroke::new(1.0, p.panel_stroke))
        .shadow(Shadow {
            offset: [0, 4],
            blur: 16,
            spread: 0,
            color: p.shadow,
        })
}

/// Applique les deux thèmes à egui, une fois, au démarrage.
pub fn install(ctx: &egui::Context) {
    for (theme, dark) in [(egui::Theme::Dark, true), (egui::Theme::Light, false)] {
        ctx.style_mut_of(theme, |style| apply(style, &palette(dark)));
    }
}

fn apply(style: &mut Style, p: &Palette) {
    let radius = CornerRadius::same(CONTROL_RADIUS);
    let v = &mut style.visuals;
    v.override_text_color = Some(p.text);
    v.weak_text_color = Some(p.text_weak);
    v.hyperlink_color = p.text;
    v.warn_fg_color = p.warn;
    v.error_fg_color = p.error;
    v.panel_fill = p.canvas;
    v.window_fill = p.panel;
    v.window_stroke = Stroke::new(1.0, p.panel_stroke);
    v.window_corner_radius = CornerRadius::same(PANEL_RADIUS);
    v.menu_corner_radius = CornerRadius::same(CONTROL_RADIUS);
    v.extreme_bg_color = p.input;
    v.text_edit_bg_color = Some(p.input);
    v.faint_bg_color = p.control;
    v.selection.bg_fill = p.accent;
    v.selection.stroke = Stroke::new(1.0, p.on_accent);
    v.slider_trailing_fill = true;

    let w = &mut v.widgets;
    w.noninteractive.bg_fill = p.panel;
    w.noninteractive.weak_bg_fill = p.panel;
    w.noninteractive.bg_stroke = Stroke::new(1.0, p.panel_stroke);
    w.noninteractive.fg_stroke = Stroke::new(1.0, p.text);
    for (state, fill) in [
        (&mut w.inactive, p.control),
        (&mut w.hovered, p.control_hover),
        (&mut w.active, p.control_active),
        (&mut w.open, p.control_hover),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = Stroke::NONE;
        state.fg_stroke = Stroke::new(1.0, p.text);
        state.corner_radius = radius;
        state.expansion = 0.0;
    }
    w.noninteractive.corner_radius = radius;
    // cases, ronds et rails de curseur (`bg_fill`) plus marqués que les boutons (`weak_bg_fill`), pour
    // se voir sur le fond d'une section, qui a la couleur d'un bouton
    w.inactive.bg_fill = p.control_active;
    w.hovered.bg_fill = p.control_active;

    let s = &mut style.spacing;
    s.item_spacing = Vec2::new(SPACE_M, 6.0);
    s.button_padding = Vec2::new(10.0, 5.0);
    s.interact_size.y = 26.0;
    s.indent = 14.0;
    s.slider_width = 150.0;

    style.text_styles = [
        (
            TextStyle::Heading,
            FontId::new(SIZE_HEADING, FontFamily::Proportional),
        ),
        (
            TextStyle::Body,
            FontId::new(SIZE_BODY, FontFamily::Proportional),
        ),
        (
            TextStyle::Button,
            FontId::new(SIZE_BODY, FontFamily::Proportional),
        ),
        (
            TextStyle::Small,
            FontId::new(SIZE_SMALL, FontFamily::Proportional),
        ),
        (
            TextStyle::Monospace,
            FontId::new(SIZE_MONO, FontFamily::Monospace),
        ),
    ]
    .into();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Contraste de luminance relative (WCAG) entre deux couleurs opaques.
    fn contrast(a: Color32, b: Color32) -> f32 {
        let lum = |c: Color32| {
            let f = |v: u8| {
                let v = f32::from(v) / 255.0;
                if v <= 0.03928 {
                    v / 12.92
                } else {
                    ((v + 0.055) / 1.055).powf(2.4)
                }
            };
            0.2126 * f(c.r()) + 0.7152 * f(c.g()) + 0.0722 * f(c.b())
        };
        let (la, lb) = (lum(a), lum(b));
        (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
    }

    /// Le texte, le texte discret et le texte sur l'accent restent lisibles dans les deux thèmes.
    #[test]
    fn text_contrast_is_readable() {
        for dark in [true, false] {
            let p = palette(dark);
            // le panneau est presque opaque : on lit le texte sur sa teinte pleine
            let panel = Color32::from_rgb(p.panel.r(), p.panel.g(), p.panel.b());
            assert!(contrast(p.text, panel) >= 7.0, "texte, dark = {dark}");
            assert!(
                contrast(p.text_weak, panel) >= 3.0,
                "texte discret, dark = {dark}"
            );
            assert!(
                contrast(p.on_accent, p.accent) >= 3.0,
                "accent, dark = {dark}"
            );
        }
    }

    /// Aucun module hors `theme` n'écrit une couleur littérale : tout passe par la `Palette`.
    /// (Une couleur calculée, comme celle de la colormap, ne commence pas par un chiffre.)
    #[test]
    fn no_hard_coded_color_outside_theme() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/src");
        let mut offenders = Vec::new();
        for entry in fs::read_dir(dir).expect("src lisible") {
            let path = entry.expect("entrée").path();
            if path.extension().is_none_or(|e| e != "rs") || path.ends_with("theme.rs") {
                continue;
            }
            let text = fs::read_to_string(&path).expect("fichier lisible");
            // on s'arrête au module de tests du fichier : ses couleurs sont des données de test
            let code = text.split("#[cfg(test)]").next().unwrap_or(&text);
            for (n, line) in code.lines().enumerate() {
                let mut rest = line;
                while let Some(at) = rest.find("Color32::from_") {
                    let tail = &rest[at..];
                    let after_paren = tail.find('(').map(|i| &tail[i + 1..]);
                    if after_paren
                        .is_some_and(|a| a.trim_start().starts_with(|c: char| c.is_ascii_digit()))
                    {
                        offenders.push(format!("{}:{}", path.display(), n + 1));
                    }
                    rest = &tail["Color32::from_".len()..];
                }
            }
        }
        assert!(
            offenders.is_empty(),
            "couleurs en dur, à prendre dans `theme::Palette` : {offenders:?}"
        );
    }
}
