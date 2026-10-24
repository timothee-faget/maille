/// Champs de résultats, une valeur par DOF.
#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    pub displacements: Vec<f64>,
    /// Forces de réaction aux DOF contraints, zéro ailleurs.
    pub reactions: Vec<f64>,
}
