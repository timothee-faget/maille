use crate::error::MailleError;

pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, MailleError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MailleError::InvalidProperty { name, value })
    }
}

/// Matériau isotrope. Chaque propriété n'est exigée que par l'analyse qui s'en sert (vérifié par
/// `Model`) : E (et ν pour les triangles et tétraèdres) en statique, k en thermique.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    /// Module d'Young E.
    pub young: Option<f64>,
    /// Coefficient de Poisson ν, requis par les triangles et les tétraèdres (absent pour Bar et Beam).
    pub poisson: Option<f64>,
    /// Conductivité thermique k.
    pub conductivity: Option<f64>,
}

impl Material {
    /// Matériau élastique pour les Bar et Beam.
    pub fn new(young: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: Some(positive("E (module d'Young)", young)?),
            poisson: None,
            conductivity: None,
        })
    }

    /// Matériau conducteur pour une analyse thermique seule : ni E ni ν.
    pub fn conductive(conductivity: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: None,
            poisson: None,
            conductivity: Some(positive("k (conductivité thermique)", conductivity)?),
        })
    }

    /// Ajoute la conductivité k à un matériau élastique.
    pub fn with_conductivity(self, conductivity: f64) -> Result<Self, MailleError> {
        Ok(Self {
            conductivity: Some(positive("k (conductivité thermique)", conductivity)?),
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
