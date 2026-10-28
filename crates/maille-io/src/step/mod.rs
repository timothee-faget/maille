//! Pièce STEP (ISO 10303) : le fichier texte Part 21 (`part21`), puis, à venir (M11), la topologie
//! B-rep et ses surfaces, maillées face par face.
//!
//! Les unités du fichier seront reprises telles quelles, comme pour le STL : le cœur ne convertit
//! rien.

pub(crate) mod part21;

use crate::IoError;

fn error(message: String) -> IoError {
    IoError::Step { message }
}
