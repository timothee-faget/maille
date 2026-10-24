use crate::error::MailleError;

pub(crate) fn positive(name: &'static str, value: f64) -> Result<f64, MailleError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(MailleError::InvalidProperty { name, value })
    }
}

/// Matériau élastique linéaire.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub young: f64,
    /// Coefficient de Poisson, requis par les éléments 2D (absent pour Bar et Beam).
    pub poisson: Option<f64>,
}

impl Material {
    pub fn new(young: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: positive("E (module d'Young)", young)?,
            poisson: None,
        })
    }

    /// Matériau élastique isotrope pour les éléments 2D, avec −1 < ν < 0,5.
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
