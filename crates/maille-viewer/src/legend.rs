//! Légende façon Abaqus : barre de couleurs en haut à gauche, graduée, avec le nom du champ.
//! Mise en page en pixels physiques (origine en haut à gauche, y vers le bas).

use crate::colormap::{Colormap, RAINBOW_BANDS};

/// Dimensions en pixels logiques, multipliées par le facteur d'échelle de l'écran (2 en Retina).
const MARGIN: f64 = 16.0;
const TITLE_HEIGHT: f64 = 52.0;
const BAR_WIDTH: f64 = 22.0;
const BAR_MAX_HEIGHT: f64 = 360.0;
const TICK_LENGTH: f64 = 5.0;
const LABEL_GAP: f64 = 8.0;
/// Largeur d'une étiquette « +3.007e+02 » ; avec la barre, fixe la largeur réservée au Mesh.
const LABEL_WIDTH: f64 = 96.0;
pub const FONT_SIZE: f64 = 13.0;
pub const TITLE_FONT_SIZE: f64 = 16.0;

/// Graduations de viridis, réparties régulièrement (l'arc-en-ciel en a une par limite de bande).
const VIRIDIS_TICKS: usize = 6;

#[derive(Debug, Clone, PartialEq)]
pub struct Tick {
    /// Ordonnée de la graduation sur la barre.
    pub y: f64,
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    /// Coin haut gauche du titre.
    pub title: [f64; 2],
    /// Barre : [gauche, haut, droite, bas]. Le minimum est en bas, le maximum en haut.
    pub bar: [f64; 4],
    pub tick_length: f64,
    /// Abscisse du bord gauche des étiquettes.
    pub label_x: f64,
    pub ticks: Vec<Tick>,
}

/// Largeur réservée à la légende sur la gauche de la fenêtre.
pub fn reserved_width(scale: f64) -> f64 {
    scale * (2.0 * MARGIN + BAR_WIDTH + TICK_LENGTH + LABEL_GAP + LABEL_WIDTH)
}

/// Légende d'un champ de (min, max) = `range` dans une fenêtre de `size` pixels.
pub fn layout(size: [f64; 2], range: (f64, f64), colormap: Colormap, scale: f64) -> Layout {
    let left = scale * MARGIN;
    let top = scale * (MARGIN + TITLE_HEIGHT);
    let height = (scale * BAR_MAX_HEIGHT).min(size[1] - top - scale * MARGIN);
    let bar = [left, top, left + scale * BAR_WIDTH, top + height.max(0.0)];
    let intervals = match colormap {
        Colormap::Rainbow => RAINBOW_BANDS,
        Colormap::Viridis => VIRIDIS_TICKS - 1,
    };
    let ticks = (0..=intervals)
        .map(|k| {
            let t = k as f64 / intervals as f64;
            Tick {
                y: bar[3] - t * (bar[3] - bar[1]),
                label: format_value(range.0 + t * (range.1 - range.0)),
            }
        })
        .collect();
    Layout {
        title: [left, scale * MARGIN],
        bar,
        tick_length: scale * TICK_LENGTH,
        label_x: bar[2] + scale * (TICK_LENGTH + LABEL_GAP),
        ticks,
    }
}

/// Nombre au format d'Abaqus : signe toujours écrit, 4 chiffres significatifs, exposant sur
/// 2 chiffres signé (« +3.007e+02 », « -1.250e-03 »).
pub fn format_value(v: f64) -> String {
    let sign = if v < 0.0 { '-' } else { '+' };
    let s = format!("{:.3e}", v.abs());
    let (mantissa, exponent) = s.split_once('e').unwrap_or((&s, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let exp_sign = if exponent < 0 { '-' } else { '+' };
    format!("{sign}{mantissa}e{exp_sign}{:02}", exponent.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_are_formatted_like_abaqus() {
        assert_eq!(format_value(300.69), "+3.007e+02");
        assert_eq!(format_value(-8.0199), "-8.020e+00");
        assert_eq!(format_value(4.7611e-3), "+4.761e-03");
        assert_eq!(format_value(0.0), "+0.000e+00");
        // l'arrondi qui fait passer à la puissance de 10 suivante
        assert_eq!(format_value(9.9996), "+1.000e+01");
        assert_eq!(format_value(1.5e123), "+1.500e+123");
    }

    #[test]
    fn rainbow_has_one_tick_per_band_limit_from_min_at_the_bottom_to_max_at_the_top() {
        let l = layout([800.0, 900.0], (-10.0, 230.0), Colormap::Rainbow, 1.0);
        assert_eq!(l.ticks.len(), RAINBOW_BANDS + 1);
        assert_eq!(l.ticks[0].label, format_value(-10.0));
        assert!((l.ticks[0].y - l.bar[3]).abs() < 1e-9);
        assert_eq!(l.ticks[RAINBOW_BANDS].label, format_value(230.0));
        assert!((l.ticks[RAINBOW_BANDS].y - l.bar[1]).abs() < 1e-9);
        // limite entre la 1re et la 2e bande : min + (max − min) / 12 = 10
        assert_eq!(l.ticks[1].label, "+1.000e+01");
        assert!((l.bar[3] - l.bar[1] - BAR_MAX_HEIGHT).abs() < 1e-9);
    }

    #[test]
    fn viridis_has_six_ticks() {
        let l = layout([800.0, 900.0], (0.0, 1.0), Colormap::Viridis, 1.0);
        let labels: Vec<&str> = l.ticks.iter().map(|t| t.label.as_str()).collect();
        assert_eq!(
            labels,
            [
                "+0.000e+00",
                "+2.000e-01",
                "+4.000e-01",
                "+6.000e-01",
                "+8.000e-01",
                "+1.000e+00"
            ]
        );
    }

    #[test]
    fn the_legend_scales_with_the_screen_and_shrinks_in_a_short_window() {
        let normal = layout([800.0, 900.0], (0.0, 1.0), Colormap::Rainbow, 1.0);
        let retina = layout([1600.0, 1800.0], (0.0, 1.0), Colormap::Rainbow, 2.0);
        for k in 0..4 {
            assert!((retina.bar[k] - 2.0 * normal.bar[k]).abs() < 1e-9);
        }
        assert!((reserved_width(2.0) - 2.0 * reserved_width(1.0)).abs() < 1e-9);
        // fenêtre basse : la barre tient dans la fenêtre, marge du bas comprise
        let short = layout([800.0, 300.0], (0.0, 1.0), Colormap::Rainbow, 1.0);
        assert!((short.bar[3] - (300.0 - MARGIN)).abs() < 1e-9);
        // tout le texte tient dans la largeur réservée
        assert!(normal.label_x + LABEL_WIDTH <= reserved_width(1.0));
    }
}
