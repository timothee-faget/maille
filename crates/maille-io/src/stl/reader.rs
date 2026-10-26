//! Lecture du STL : facettes brutes (trois sommets), avec le `solid` de chacune.

use super::error;
use crate::IoError;

/// Facettes telles que lues, avant soudure. Les normales du fichier sont ignorées : on les
/// recalcule depuis les sommets.
pub(super) struct Facets {
    pub corners: Vec<[[f64; 3]; 3]>,
    /// Indice (dans `solids`) du `solid` de chaque facette.
    pub solid: Vec<usize>,
    /// Noms des `solid`, vides s'ils n'en ont pas (toujours le cas en binaire).
    pub solids: Vec<String>,
}

pub(super) fn read(bytes: &[u8]) -> Result<Facets, IoError> {
    // un STL ASCII peut commencer par « solid » comme l'en-tête d'un binaire : seule la taille tranche
    if bytes.len() >= 84 {
        let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        if bytes.len() == 84 + 50 * count {
            return binary(&bytes[84..], count);
        }
    }
    let text = std::str::from_utf8(bytes).map_err(|_| not_stl())?;
    if !text.trim_start().to_ascii_lowercase().starts_with("solid") {
        return Err(not_stl());
    }
    ascii(text)
}

fn not_stl() -> IoError {
    error(
        "ni STL binaire (taille incohérente avec le nombre de facettes) ni STL ASCII (« solid » attendu en tête) : fichier tronqué ou d'un autre format"
            .to_string(),
    )
}

fn binary(data: &[u8], count: usize) -> Result<Facets, IoError> {
    let float =
        |chunk: &[u8]| f64::from(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
    let mut corners = Vec::with_capacity(count);
    for (i, record) in data.as_chunks::<50>().0.iter().enumerate() {
        let mut c = [[0.0; 3]; 3];
        for (v, corner) in c.iter_mut().enumerate() {
            for (k, x) in corner.iter_mut().enumerate() {
                *x = float(&record[12 + 12 * v + 4 * k..]);
                if !x.is_finite() {
                    return Err(error(format!("facette {} : coordonnée non finie", i + 1)));
                }
            }
        }
        corners.push(c);
    }
    finish(Facets {
        solid: vec![0; corners.len()],
        corners,
        solids: vec![String::new()],
    })
}

type Lines<'a> = std::iter::Peekable<Box<dyn Iterator<Item = (usize, &'a str)> + 'a>>;

fn ascii(text: &str) -> Result<Facets, IoError> {
    let iter: Box<dyn Iterator<Item = (usize, &str)>> = Box::new(
        text.lines()
            .enumerate()
            .map(|(i, l)| (i + 1, l.trim()))
            .filter(|(_, l)| !l.is_empty()),
    );
    let mut lines: Lines = iter.peekable();
    let mut f = Facets {
        corners: Vec::new(),
        solid: Vec::new(),
        solids: Vec::new(),
    };
    let mut open = false;
    while let Some((n, l)) = lines.next() {
        let lower = l.to_ascii_lowercase();
        match lower.split_whitespace().next() {
            Some("solid") => {
                if open {
                    return Err(error(format!(
                        "ligne {n} : « solid » sans « endsolid » avant"
                    )));
                }
                open = true;
                f.solids.push(l[5..].trim().to_string());
            }
            Some("endsolid") if open => open = false,
            Some("facet") if open => {
                expect(&mut lines, "outer loop")?;
                let mut c = [[0.0; 3]; 3];
                for corner in &mut c {
                    *corner = vertex(&mut lines)?;
                }
                expect(&mut lines, "endloop")?;
                expect(&mut lines, "endfacet")?;
                f.corners.push(c);
                f.solid.push(f.solids.len() - 1);
            }
            _ => {
                return Err(error(format!(
                    "ligne {n} : « {l} » inattendu dans un STL ASCII"
                )))
            }
        }
    }
    if open {
        return Err(error("fichier tronqué : « endsolid » manquant".to_string()));
    }
    finish(f)
}

fn finish(f: Facets) -> Result<Facets, IoError> {
    if f.corners.is_empty() {
        return Err(error("le fichier ne contient aucune facette".to_string()));
    }
    Ok(f)
}

fn expect(lines: &mut Lines, keyword: &str) -> Result<(), IoError> {
    match lines.next() {
        Some((_, l))
            if l.to_ascii_lowercase()
                .split_whitespace()
                .eq(keyword.split(' ')) =>
        {
            Ok(())
        }
        Some((n, l)) => Err(error(format!(
            "ligne {n} : « {keyword} » attendu, trouvé « {l} »"
        ))),
        None => Err(error(format!("fichier tronqué : « {keyword} » attendu"))),
    }
}

fn vertex(lines: &mut Lines) -> Result<[f64; 3], IoError> {
    let Some((n, l)) = lines.next() else {
        return Err(error("fichier tronqué : « vertex » attendu".to_string()));
    };
    let lower = l.to_ascii_lowercase();
    let mut words = lower.split_whitespace();
    if words.next() != Some("vertex") {
        return Err(error(format!(
            "ligne {n} : « vertex x y z » attendu, trouvé « {l} »"
        )));
    }
    let mut p = [0.0; 3];
    for x in &mut p {
        *x = words
            .next()
            .and_then(|w| w.parse::<f64>().ok())
            .filter(|v| v.is_finite())
            .ok_or_else(|| {
                error(format!(
                    "ligne {n} : trois nombres attendus après « vertex »"
                ))
            })?;
    }
    Ok(p)
}
