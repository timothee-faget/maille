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
}

impl Material {
    pub fn new(young: f64) -> Result<Self, MailleError> {
        Ok(Self {
            young: positive("E (module d'Young)", young)?,
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
