//! Colormaps : valeur normalisée t ∈ [0, 1] → couleur sRGB, et tables (LUT) pour une texture 1D.

/// Nombre de bandes de l'arc-en-ciel, comme le réglage par défaut d'Abaqus.
pub const RAINBOW_BANDS: usize = 12;

/// Nombre de texels de la LUT viridis : assez pour qu'aucune bande ne se voie.
const VIRIDIS_TEXELS: usize = 256;

/// 17 points de la table viridis de matplotlib (`_viridis_data`, indices 0, 16, …, 240, 255),
/// données publiées sous CC0 par leurs auteurs (N. Smith, S. van der Walt). Valeurs sRGB.
const VIRIDIS: [(u8, [f64; 3]); 17] = [
    (0, [0.267004, 0.004874, 0.329415]),
    (16, [0.282327, 0.094955, 0.417331]),
    (32, [0.278826, 0.175490, 0.483397]),
    (48, [0.258965, 0.251537, 0.524736]),
    (64, [0.229739, 0.322361, 0.545706]),
    (80, [0.199430, 0.387607, 0.554642]),
    (96, [0.172719, 0.448791, 0.557885]),
    (112, [0.149039, 0.508051, 0.557250]),
    (128, [0.127568, 0.566949, 0.550556]),
    (144, [0.120638, 0.625828, 0.533488]),
    (160, [0.157851, 0.683765, 0.501686]),
    (176, [0.246070, 0.738910, 0.452024]),
    (192, [0.369214, 0.788888, 0.382914]),
    (208, [0.515992, 0.831158, 0.294279]),
    (224, [0.678489, 0.863742, 0.189503]),
    (240, [0.845561, 0.887322, 0.099702]),
    (255, [0.993248, 0.906157, 0.143936]),
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colormap {
    /// Arc-en-ciel en bandes, du bleu (minimum) au rouge (maximum), façon Abaqus.
    Rainbow,
    /// Viridis continue, perceptuellement uniforme et lisible par les daltoniens.
    Viridis,
}

impl Colormap {
    pub fn next(self) -> Self {
        match self {
            Self::Rainbow => Self::Viridis,
            Self::Viridis => Self::Rainbow,
        }
    }

    /// Couleur sRGB (composantes dans [0, 1]) en t, ramené dans [0, 1].
    pub fn color(self, t: f64) -> [f64; 3] {
        let t = t.clamp(0.0, 1.0);
        match self {
            // la bande k couvre [k/12, (k+1)/12] ; ses couleurs vont du bleu (k = 0) au rouge (k = 11)
            Self::Rainbow => {
                let band = ((t * RAINBOW_BANDS as f64) as usize).min(RAINBOW_BANDS - 1);
                hue(240.0 * (1.0 - band as f64 / (RAINBOW_BANDS - 1) as f64))
            }
            Self::Viridis => {
                let x = t * 255.0;
                let k = VIRIDIS
                    .windows(2)
                    .position(|w| x <= f64::from(w[1].0))
                    .unwrap_or(VIRIDIS.len() - 2);
                let ((i0, c0), (i1, c1)) = (VIRIDIS[k], VIRIDIS[k + 1]);
                let s = (x - f64::from(i0)) / f64::from(i1 - i0);
                std::array::from_fn(|c| c0[c] + s * (c1[c] - c0[c]))
            }
        }
    }

    /// Table de couleurs sRGBA 8 bits, à échantillonner sans interpolation (texels de largeur égale) :
    /// une couleur par bande pour l'arc-en-ciel, 256 couleurs pour viridis.
    pub fn lut(self) -> Vec<[u8; 4]> {
        let n = match self {
            Self::Rainbow => RAINBOW_BANDS,
            Self::Viridis => VIRIDIS_TEXELS,
        };
        (0..n)
            .map(|i| {
                // centre du texel i
                let [r, g, b] = self.color((i as f64 + 0.5) / n as f64);
                [to_u8(r), to_u8(g), to_u8(b), 255]
            })
            .collect()
    }
}

/// Teinte h (degrés, 0 = rouge, 120 = vert, 240 = bleu) à saturation et luminosité maximales.
fn hue(h: f64) -> [f64; 3] {
    let channel = |n: f64| {
        let k = (n + h / 60.0) % 6.0;
        1.0 - (k.min(4.0 - k).clamp(0.0, 1.0))
    };
    [channel(5.0), channel(3.0), channel(1.0)]
}

fn to_u8(c: f64) -> u8 {
    (c.clamp(0.0, 1.0) * 255.0).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: [f64; 3], b: [f64; 3]) -> bool {
        a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9)
    }

    #[test]
    fn rainbow_goes_from_blue_to_red_through_green() {
        assert!(close(Colormap::Rainbow.color(0.0), [0.0, 0.0, 1.0]));
        assert!(close(Colormap::Rainbow.color(1.0), [1.0, 0.0, 0.0]));
        // hors bornes : ramené aux extrémités
        assert!(close(Colormap::Rainbow.color(-3.0), [0.0, 0.0, 1.0]));
        assert!(close(Colormap::Rainbow.color(7.0), [1.0, 0.0, 0.0]));
        // milieu : entre cyan et jaune, le vert domine
        let mid = Colormap::Rainbow.color(0.5);
        assert!(mid[1] > mid[0] && mid[1] > mid[2], "{mid:?}");
    }

    #[test]
    fn rainbow_is_constant_within_a_band_and_changes_between_bands() {
        let w = 1.0 / RAINBOW_BANDS as f64;
        for k in 0..RAINBOW_BANDS {
            let (start, end) = (k as f64 * w + 1e-9, (k + 1) as f64 * w - 1e-9);
            assert!(close(
                Colormap::Rainbow.color(start),
                Colormap::Rainbow.color(end)
            ));
            if k + 1 < RAINBOW_BANDS {
                assert!(!close(
                    Colormap::Rainbow.color(end),
                    Colormap::Rainbow.color(end + 2e-9)
                ));
            }
        }
    }

    #[test]
    fn viridis_matches_the_matplotlib_table_at_its_ends_and_anchors() {
        assert!(close(Colormap::Viridis.color(0.0), VIRIDIS[0].1));
        assert!(close(Colormap::Viridis.color(1.0), VIRIDIS[16].1));
        assert!(close(Colormap::Viridis.color(128.0 / 255.0), VIRIDIS[8].1));
        // entre deux points : interpolation linéaire
        let half = Colormap::Viridis.color(8.0 / 255.0);
        let expected = std::array::from_fn(|c| (VIRIDIS[0].1[c] + VIRIDIS[1].1[c]) / 2.0);
        assert!(close(half, expected));
    }

    #[test]
    fn viridis_lightness_increases_monotonically() {
        // propriété clé de viridis : luminance croissante (lisible en niveaux de gris)
        let luminance = |c: [f64; 3]| 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
        let mut previous = -1.0;
        for i in 0..=100 {
            let l = luminance(Colormap::Viridis.color(i as f64 / 100.0));
            assert!(l > previous, "t = {}", i as f64 / 100.0);
            previous = l;
        }
    }

    #[test]
    fn luts_have_one_texel_per_band_or_256_for_viridis() {
        let rainbow = Colormap::Rainbow.lut();
        assert_eq!(rainbow.len(), RAINBOW_BANDS);
        assert_eq!(rainbow[0], [0, 0, 255, 255]);
        assert_eq!(rainbow[RAINBOW_BANDS - 1], [255, 0, 0, 255]);
        let viridis = Colormap::Viridis.lut();
        assert_eq!(viridis.len(), 256);
        assert!(viridis.iter().all(|c| c[3] == 255));
    }
}
