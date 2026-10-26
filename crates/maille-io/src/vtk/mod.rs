//! Résultats VTK : un `.vtu` (UnstructuredGrid, ASCII) par Increment et un `.pvd` qui donne le
//! temps de chacun, lisibles par ParaView. Les grandeurs sont aux Nodes seulement.
//!
//! Le `.vtu` et le `.pvd` se relisent (`read_point_data`, `read_pvd`) : c'est ainsi que le
//! couplage thermo-mécanique lit les températures d'un autre calcul. La lecture ne comprend que
//! l'ASCII, ce que ces fonctions écrivent.

mod fields;
mod reader;
mod writer;

pub use fields::result_arrays;
pub use reader::{read_point_data, read_pvd};
pub use writer::{write_pvd, write_vtu};

/// Grandeur aux Nodes : `components` valeurs par Node, Node après Node.
#[derive(Debug, Clone, PartialEq)]
pub struct PointArray {
    pub name: String,
    pub components: usize,
    pub values: Vec<f64>,
}

/// Un `.vtu` d'un `.pvd` : le temps de l'Increment et le nom du fichier.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    pub time: f64,
    pub file: String,
}

fn vtk_error(problem: String) -> crate::IoError {
    crate::IoError::Vtk { problem }
}
