use crate::error::MailleError;

pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, MailleError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MailleError::InvalidProperty { name, value })
    }
}

/// Matériau isotrope. Chaque propriété n'est exigée que par l'analyse qui s'en sert (vérifié par
/// `Model`) : E (et ν pour les triangles et tétraèdres) en statique, et α si on impose des
/// températures ; k en thermique, et en plus ρ et c en thermique transitoire.
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

    /// Matériau élastique isotrope pour les triangles et tétraèdres, avec −1 < ν < 0,5.
    pub fn with_poisson(young: f64, poisson: f64) -> Result<Self, MailleError> {
        if !(poisson > -1.0 && poisson < 0.5) {
            return Err(MailleError::InvalidProperty {
                name: "ν (coefficient de Poisson, entre -1 et 0,5 exclus)",
                value: poisson,
            });
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
