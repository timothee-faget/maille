//! Bibliothèque de Materials : des valeurs usuelles en SI, converties vers le système d'unités
//! cohérent du modèle. Le cœur ne convertit rien ; la conversion se fait ici, une fois, au moment
//! où un Material de la bibliothèque devient la `MaterialDef` d'un modèle.
//!
//! ```toml
//! [materials.acier]
//! description = "Acier de construction (type S235)"
//! young = 210e9                   # Pa ; les autres : kg/m³, W/(m·K), J/(kg·K), 1/K
//! poisson = 0.3
//! yield_stress = 235e6
//! ```
//!
//! Une bibliothèque est livrée avec Maille (`Library::builtin`) ; `Library::merged` y ajoute celle
//! de l'utilisateur, qui l'emporte sur un nom commun.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::model::MaterialDef;
use crate::IoError;

const BUILTIN: &str = include_str!("../data/materials.toml");

/// Système d'unités cohérent : une longueur, une force et un temps, chacun en unités SI. La masse
/// s'en déduit (m = F T² / L) ; la température reste en kelvin (ou en degrés Celsius, qui ont le
/// même écart).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitSystem {
    length: f64,
    force: f64,
    time: f64,
}

impl UnitSystem {
    /// m, N, s.
    pub const SI: Self = Self {
        length: 1.0,
        force: 1.0,
        time: 1.0,
    };
    /// mm, N, s : contraintes en MPa, masses en tonnes.
    pub const MM_N_S: Self = Self {
        length: 1e-3,
        force: 1.0,
        time: 1.0,
    };

    /// `length`, `force` et `time` valent une unité du système en mètres, newtons et secondes.
    pub fn new(length: f64, force: f64, time: f64) -> Result<Self, IoError> {
        for (name, v) in [("longueur", length), ("force", force), ("temps", time)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(invalid(format!(
                    "unité de {name} : {v} ; elle doit valoir un nombre positif en unités SI"
                )));
            }
        }
        Ok(Self {
            length,
            force,
            time,
        })
    }

    /// Lit `"longueur, force, temps"`, par exemple `"mm, N, s"`. Longueurs : m, cm, mm, in ;
    /// forces : N, kN, lbf ; temps : s.
    pub fn parse(text: &str) -> Result<Self, IoError> {
        let parts: Vec<&str> = text.split(',').map(str::trim).collect();
        let [length, force, time] = parts.as_slice() else {
            return Err(invalid(format!(
                "système d'unités « {text} » : écrire « longueur, force, temps », par exemple « mm, N, s »"
            )));
        };
        let length = match *length {
            "m" => 1.0,
            "cm" => 1e-2,
            "mm" => 1e-3,
            "in" => 0.0254,
            other => return Err(unknown_unit("longueur", other, "m, cm, mm, in")),
        };
        let force = match *force {
            "N" => 1.0,
            "kN" => 1e3,
            "lbf" => 4.448_221_615_260_5,
            other => return Err(unknown_unit("force", other, "N, kN, lbf")),
        };
        let time = match *time {
            "s" => 1.0,
            other => return Err(unknown_unit("temps", other, "s")),
        };
        Self::new(length, force, time)
    }

    /// Masse d'une unité du système, en kilogrammes.
    fn mass(&self) -> f64 {
        self.force * self.time * self.time / self.length
    }

    /// Ce par quoi on multiplie une valeur SI pour l'avoir dans ce système.
    fn factors(&self) -> Factors {
        Factors {
            stress: self.length * self.length / self.force,
            density: self.length.powi(3) / self.mass(),
            specific_heat: self.mass() / (self.force * self.length),
            conductivity: self.time / self.force,
        }
    }
}

struct Factors {
    stress: f64,
    density: f64,
    specific_heat: f64,
    conductivity: f64,
}

fn invalid(problem: String) -> IoError {
    IoError::InvalidModel { problem }
}

fn unknown_unit(what: &str, found: &str, known: &str) -> IoError {
    invalid(format!(
        "unité de {what} « {found} » inconnue ; unités connues : {known}"
    ))
}

/// Un Material de la bibliothèque, en SI. Les propriétés absentes restent absentes.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LibraryMaterial {
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
    /// Module d'Young E (Pa).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub young: Option<f64>,
    /// Coefficient de Poisson ν.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub poisson: Option<f64>,
    /// Conductivité thermique k (W/(m·K)).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conductivity: Option<f64>,
    /// Masse volumique ρ (kg/m³).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub density: Option<f64>,
    /// Chaleur massique c (J/(kg·K)).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub specific_heat: Option<f64>,
    /// Coefficient de dilatation thermique α (1/K).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expansion: Option<f64>,
    /// Limite d'élasticité σy0 (Pa). Elle n'entre pas dans `to_def` : `MaterialDef::yield_stress`
    /// active la plasticité de von Mises, que l'utilisateur doit demander (elle n'existe pas en
    /// contraintes planes).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub yield_stress: Option<f64>,
}

impl LibraryMaterial {
    /// Le Material exprimé dans le système d'unités `units`, prêt pour `[materials.<nom>]`. Il est
    /// élastique et thermique, sans plasticité : voir `yield_stress_in`.
    pub fn to_def(&self, units: UnitSystem) -> MaterialDef {
        let f = units.factors();
        MaterialDef {
            young: self.young.map(|v| v * f.stress),
            poisson: self.poisson,
            conductivity: self.conductivity.map(|v| v * f.conductivity),
            density: self.density.map(|v| v * f.density),
            specific_heat: self.specific_heat.map(|v| v * f.specific_heat),
            expansion: self.expansion,
            yield_stress: None,
            hardening: None,
            hyperelasticity: None,
        }
    }

    /// Le Material d'un modèle écrit dans le système `units`, ramené en SI pour la bibliothèque.
    /// L'écrouissage et l'hyperélasticité ne sont pas gardés.
    pub fn from_def(def: &MaterialDef, units: UnitSystem, description: &str) -> Self {
        let f = units.factors();
        Self {
            description: description.to_string(),
            young: def.young.map(|v| v / f.stress),
            poisson: def.poisson,
            conductivity: def.conductivity.map(|v| v / f.conductivity),
            density: def.density.map(|v| v / f.density),
            specific_heat: def.specific_heat.map(|v| v / f.specific_heat),
            expansion: def.expansion,
            yield_stress: def.yield_stress.map(|v| v / f.stress),
        }
    }

    /// La limite d'élasticité σy0 dans le système `units`, à mettre dans `MaterialDef::yield_stress`
    /// pour une analyse plastique.
    pub fn yield_stress_in(&self, units: UnitSystem) -> Option<f64> {
        self.yield_stress.map(|v| v * units.factors().stress)
    }

    fn check(&self, name: &str) -> Result<(), IoError> {
        let positive = [
            ("young", self.young),
            ("conductivity", self.conductivity),
            ("density", self.density),
            ("specific_heat", self.specific_heat),
            ("yield_stress", self.yield_stress),
        ];
        for (key, value) in positive {
            if let Some(v) = value.filter(|v| !(v.is_finite() && *v > 0.0)) {
                return Err(invalid(format!(
                    "bibliothèque de Materials, « {name} » : {key} = {v} ; une valeur positive est attendue"
                )));
            }
        }
        if let Some(v) = self.poisson.filter(|v| !(*v > -1.0 && *v < 0.5)) {
            return Err(invalid(format!(
                "bibliothèque de Materials, « {name} » : poisson = {v} ; il doit être entre -1 et 0,5"
            )));
        }
        if let Some(v) = self.expansion.filter(|v| !v.is_finite()) {
            return Err(invalid(format!(
                "bibliothèque de Materials, « {name} » : expansion = {v} ; un nombre fini est attendu"
            )));
        }
        Ok(())
    }
}

/// Ensemble de Materials nommés.
#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Library {
    materials: BTreeMap<String, LibraryMaterial>,
}

impl Library {
    /// La bibliothèque livrée avec Maille.
    pub fn builtin() -> Self {
        Self::parse(BUILTIN).expect("la bibliothèque livrée est valide (test `builtin_is_valid`)")
    }

    /// Lit une bibliothèque TOML. Une clé inconnue ou une valeur absurde est une erreur.
    pub fn parse(text: &str) -> Result<Self, IoError> {
        let library: Self = toml::from_str(text).map_err(|e| IoError::Toml(e.to_string()))?;
        for (name, material) in &library.materials {
            material.check(name)?;
        }
        Ok(library)
    }

    /// Écrit la bibliothèque en TOML : `parse(&to_toml()?)` redonne la même.
    pub fn to_toml(&self) -> Result<String, IoError> {
        toml::to_string(self).map_err(|e| IoError::Toml(e.to_string()))
    }

    /// Ajoute `material` sous `name`, ou remplace celui qui porte ce nom. Refuse un nom vide ou des
    /// valeurs absurdes.
    pub fn insert(&mut self, name: &str, material: LibraryMaterial) -> Result<(), IoError> {
        if name.trim().is_empty() {
            return Err(invalid(
                "bibliothèque de Materials : donne un nom au Material".to_string(),
            ));
        }
        material.check(name)?;
        self.materials.insert(name.to_string(), material);
        Ok(())
    }

    /// Ajoute les Materials de `other` ; sur un nom commun, `other` l'emporte.
    pub fn merged(mut self, other: Self) -> Self {
        self.materials.extend(other.materials);
        self
    }

    /// Noms des Materials, dans l'ordre alphabétique.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.materials.keys().map(String::as_str)
    }

    pub fn get(&self, name: &str) -> Option<&LibraryMaterial> {
        self.materials.get(name)
    }

    /// Le Material `name` exprimé dans le système d'unités `units`.
    pub fn material(&self, name: &str, units: UnitSystem) -> Result<MaterialDef, IoError> {
        self.get(name)
            .map(|m| m.to_def(units))
            .ok_or_else(|| IoError::UnknownMaterial {
                key: "bibliothèque de Materials".to_string(),
                name: name.to_string(),
                known: self.names().map(str::to_string).collect(),
            })
    }
}
