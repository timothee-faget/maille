use crate::error::MailleError;

pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, MailleError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MailleError::InvalidProperty { name, value })
    }
}

/// Matériau isotrope. Chaque propriété n'est exigée que par l'analyse qui s'en sert (vérifié par
/// `Model`) : E (et ν pour les triangles et tétraèdres) en statique, et α pour la dilatation d'un
/// champ de températures ; k en thermique, et en plus ρ et c en thermique transitoire. Avec une
/// `Plasticity`, le Material est élasto-plastique (calcul par `Model::solve_quasi_static`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    /// Module d'Young E.
    pub young: Option<f64>,
    /// Coefficient de Poisson ν, requis par les triangles et les tétraèdres (absent pour Bar et Beam).
    pub poisson: Option<f64>,
    /// Conductivité thermique k : puissance par unité de longueur et par degré (ex. W/(mm·K)).
    pub conductivity: Option<f64>,
    /// Masse volumique ρ (ex. kg/mm³ ; acier : 7,85·10⁻⁶).
    pub density: Option<f64>,
    /// Chaleur massique c : énergie par unité de masse et par degré (ex. J/(kg·K) avec des W, des
    /// mm et des s ; acier : 460). Seul le produit ρc compte en thermique.
    pub specific_heat: Option<f64>,
    /// Coefficient de dilatation thermique α : déformation par degré (ex. 1/K ; acier : 1,2·10⁻⁵).
    pub expansion: Option<f64>,
    /// Plasticité de von Mises, absente pour un Material élastique.
    pub plasticity: Option<Plasticity>,
    /// Loi hyperélastique en grandes déformations, absente : Saint-Venant–Kirchhoff.
    pub hyperelasticity: Option<Hyperelasticity>,
}

/// Loi hyperélastique d'un Material en grandes déformations (`Kinematics::Large`), de paramètres
/// λ et μ tirés de E et ν. Sans elle, la loi est Saint-Venant–Kirchhoff.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hyperelasticity {
    /// Néo-Hookéen compressible : W = μ/2 (I₁ − 3) − μ ln J + λ/2 (ln J)².
    NeoHookean,
}

/// Plasticité de von Mises à écrouissage isotrope linéaire : limite d'élasticité σy(p) = σy0 + H p,
/// p étant la déformation plastique cumulée.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plasticity {
    /// Limite d'élasticité initiale σy0 (une contrainte, ex. MPa).
    pub yield_stress: f64,
    /// Module d'écrouissage H (une contrainte, ex. MPa) : pente de σy en fonction de p. H = 0 :
    /// parfaitement plastique.
    pub hardening: f64,
}

impl Material {
    /// Matériau élastique pour les Bar et Beam.
    pub fn new(young: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: Some(positive("E (module d'Young)", young)?),
            poisson: None,
            conductivity: None,
            density: None,
            specific_heat: None,
            expansion: None,
            plasticity: None,
            hyperelasticity: None,
        })
    }

    /// Matériau conducteur pour une analyse thermique seule : ni E ni ν.
    pub fn conductive(conductivity: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: None,
            poisson: None,
            conductivity: Some(positive("k (conductivité thermique)", conductivity)?),
            density: None,
            specific_heat: None,
            expansion: None,
            plasticity: None,
            hyperelasticity: None,
        })
    }

    /// Ajoute la conductivité k à un matériau élastique.
    pub fn with_conductivity(self, conductivity: f64) -> Result<Self, MailleError> {
        Ok(Self {
            conductivity: Some(positive("k (conductivité thermique)", conductivity)?),
            ..self
        })
    }

    /// Ajoute la masse volumique ρ.
    pub fn with_density(self, density: f64) -> Result<Self, MailleError> {
        Ok(Self {
            density: Some(positive("ρ (masse volumique)", density)?),
            ..self
        })
    }

    /// Ajoute la chaleur massique c.
    pub fn with_specific_heat(self, specific_heat: f64) -> Result<Self, MailleError> {
        Ok(Self {
            specific_heat: Some(positive("c (chaleur massique)", specific_heat)?),
            ..self
        })
    }

    /// Ajoute le coefficient de dilatation α.
    pub fn with_expansion(self, expansion: f64) -> Result<Self, MailleError> {
        Ok(Self {
            expansion: Some(positive("α (coefficient de dilatation)", expansion)?),
            ..self
        })
    }

    /// Rend le Material élasto-plastique (von Mises) : limite d'élasticité initiale σy0 > 0 et module
    /// d'écrouissage H ≥ 0 (0 : parfaitement plastique), voir `Plasticity`.
    pub fn with_plasticity(self, yield_stress: f64, hardening: f64) -> Result<Self, MailleError> {
        let yield_stress = positive("σy0 (limite d'élasticité)", yield_stress)?;
        if !(hardening.is_finite() && hardening >= 0.0) {
            return Err(MailleError::NegativeProperty {
                name: "H (module d'écrouissage)",
                value: hardening,
            });
        }
        Ok(Self {
            plasticity: Some(Plasticity {
                yield_stress,
                hardening,
            }),
            ..self
        })
    }

    /// Loi hyperélastique du Material (E et ν requis) : voir `Hyperelasticity`. Pas de plasticité en
    /// grandes déformations, donc ni les deux ensemble.
    pub fn with_hyperelasticity(self, law: Hyperelasticity) -> Self {
        Self {
            hyperelasticity: Some(law),
            ..self
        }
    }

    /// Matériau élastique isotrope pour les triangles et tétraèdres, avec −1 < ν < 0,5.
    pub fn with_poisson(young: f64, poisson: f64) -> Result<Self, MailleError> {
        if !(poisson > -1.0 && poisson < 0.5) {
            return Err(MailleError::PoissonOutOfRange { poisson });
        }
        Ok(Self {
            poisson: Some(poisson),
            ..Self::new(young)?
        })
    }
}

/// Section d'une barre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Section {
    pub area: f64,
}

impl Section {
    pub fn new(area: f64) -> Result<Self, MailleError> {
        Ok(Self {
            area: positive("A (aire de section)", area)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plasticity_needs_a_positive_yield_stress_and_a_non_negative_hardening() {
        let m = Material::with_poisson(210e3, 0.3).unwrap();
        let perfect = m.with_plasticity(250.0, 0.0).unwrap();
        assert_eq!(
            perfect.plasticity,
            Some(Plasticity {
                yield_stress: 250.0,
                hardening: 0.0
            })
        );
        for sy in [0.0, -1.0, f64::NAN] {
            assert!(matches!(
                m.with_plasticity(sy, 0.0),
                Err(MailleError::InvalidProperty { .. })
            ));
        }
        for h in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                m.with_plasticity(250.0, h),
                Err(MailleError::NegativeProperty { .. })
            ));
        }
    }
}
