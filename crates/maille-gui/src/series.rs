//! Les résultats d'un calcul vus comme une suite d'Increments : lecture à la demande, valeurs de la
//! sonde en un Node, courbes dans le temps, forces de contact. Pur (hors lecture des fichiers) et
//! testé ; `app` ne fait que l'afficher.

use std::collections::BTreeMap;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};

use maille_core::results::Stress;
use maille_io::msh::{self, MshFile};
use maille_io::vtk::{self, VtuResults};

/// Les Increments d'un `.pvd` (ou l'unique Increment d'un `.vtu`), relus au besoin.
pub struct Series {
    pvd: PathBuf,
    times: Vec<f64>,
    files: Vec<PathBuf>,
    cache: Vec<Option<VtuResults>>,
}

impl Series {
    pub fn open(path: &Path) -> Result<Self, Box<dyn Error>> {
        let read = |p: &Path| {
            fs::read_to_string(p).map_err(|e| format!("impossible de lire {} : {e}", p.display()))
        };
        let (times, files) = match path.extension().and_then(|e| e.to_str()) {
            Some("pvd") => {
                let entries = vtk::read_pvd(&read(path)?)?;
                let dir = path.parent().unwrap_or(Path::new("."));
                (
                    entries.iter().map(|e| e.time).collect(),
                    entries.iter().map(|e| dir.join(&e.file)).collect(),
                )
            }
            Some("vtu") => (vec![0.0], vec![path.to_path_buf()]),
            _ => {
                return Err(format!("{} : attendu un fichier .pvd ou .vtu", path.display()).into())
            }
        };
        if files.is_empty() {
            return Err(format!("{} ne contient aucun incrément", path.display()).into());
        }
        let cache = (0..files.len()).map(|_| None).collect();
        Ok(Self {
            pvd: path.to_path_buf(),
            times,
            files,
            cache,
        })
    }

    pub fn len(&self) -> usize {
        self.files.len()
    }

    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
    }

    pub fn time(&self, k: usize) -> f64 {
        self.times[k]
    }

    pub fn path(&self) -> &Path {
        &self.pvd
    }

    /// L'Increment `k` (depuis 0), lu à la première demande.
    pub fn get(&mut self, k: usize) -> Result<&VtuResults, Box<dyn Error>> {
        if self.cache[k].is_none() {
            let file = &self.files[k];
            let text = fs::read_to_string(file)
                .map_err(|e| format!("impossible de lire {} : {e}", file.display()))?;
            self.cache[k] = Some(vtk::read_results(&text)?);
        }
        Ok(self.cache[k].as_ref().expect("lu ci-dessus"))
    }

    /// Le fichier des forces de contact qui accompagne ce `.pvd` (`<préfixe>_contacts.csv`), s'il existe.
    pub fn contacts_csv(&self) -> Option<PathBuf> {
        let stem = self.pvd.file_stem()?.to_string_lossy().into_owned();
        let path = self.pvd.with_file_name(format!("{stem}_contacts.csv"));
        path.exists().then_some(path)
    }

    /// Le `.msh` du calcul, d'après le modèle à côté des résultats : `maille solve modele.toml` écrit
    /// `modele.pvd` (ou `modele_thermal.pvd`, `modele_0003.vtu`…), et `modele.toml` nomme son
    /// `mesh`. `None` si ce modèle n'est pas là (résultats déplacés, préfixe `-o`).
    pub fn mesh_path(&self) -> Option<PathBuf> {
        let mut stem = self.pvd.file_stem()?.to_string_lossy().into_owned();
        if let Some((head, tail)) = stem.rsplit_once('_') {
            if tail.len() == 4 && tail.bytes().all(|b| b.is_ascii_digit()) {
                stem = head.to_string();
            }
        }
        if let Some(head) = stem.strip_suffix("_thermal") {
            stem = head.to_string();
        }
        let toml = self.pvd.with_file_name(format!("{stem}.toml"));
        let model = maille_io::model::parse(&fs::read_to_string(&toml).ok()?).ok()?;
        Some(toml.with_file_name(model.mesh))
    }

    /// Lit le `.msh` `path` pour en montrer les Sets ; erreur s'il n'a pas les Nodes et les Elements
    /// des résultats (ce n'est pas le maillage du calcul).
    pub fn read_mesh(&mut self, path: &Path) -> Result<MshFile, Box<dyn Error>> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("impossible de lire {} : {e}", path.display()))?;
        let file = msh::read(&text)?;
        let r = self.get(0)?;
        if (file.nodes.len(), file.cells.len(), file.kind) != (r.nodes.len(), r.cells.len(), r.kind)
        {
            return Err(format!(
                "{} n'est pas le maillage de ces résultats : {} Nodes et {} Elements, contre {} et {}",
                path.display(),
                file.nodes.len(),
                file.cells.len(),
                r.nodes.len(),
                r.cells.len()
            )
            .into());
        }
        Ok(file)
    }

    /// Le Node `node` au fil des Increments : une valeur de la grandeur `label` (voir `readings`) par
    /// Increment, avec son temps.
    pub fn curve(&mut self, node: usize, label: &str) -> Result<Vec<[f64; 2]>, Box<dyn Error>> {
        let mut points = Vec::with_capacity(self.len());
        for k in 0..self.len() {
            let time = self.times[k];
            let value = readings(self.get(k)?, node)
                .into_iter()
                .find(|(l, _)| l == label)
                .map(|(_, v)| v)
                .ok_or_else(|| format!("l'incrément {k} n'a pas de « {label} »"))?;
            points.push([time, value]);
        }
        Ok(points)
    }
}

/// Les grandeurs scalaires au Node `node`, étiquetées : chaque composante des champs du fichier, et
/// ce qu'on en déduit (norme du déplacement, von Mises). L'ordre est celui des champs du fichier.
pub fn readings(r: &VtuResults, node: usize) -> Vec<(String, f64)> {
    let mut out = Vec::new();
    for array in &r.arrays {
        let Some(values) = array
            .values
            .get(node * array.components..(node + 1) * array.components)
        else {
            continue;
        };
        match (array.name.as_str(), array.components) {
            ("displacement", 3) => {
                out.extend(
                    ["ux", "uy", "uz"]
                        .iter()
                        .zip(values)
                        .map(|(l, v)| (l.to_string(), *v)),
                );
                out.push((
                    "|u|".to_string(),
                    values.iter().map(|v| v * v).sum::<f64>().sqrt(),
                ));
            }
            ("stress", 6) => {
                let [xx, yy, zz, xy, yz, xz] = [0, 1, 2, 3, 4, 5].map(|i| values[i]);
                for (l, v) in [
                    ("σxx", xx),
                    ("σyy", yy),
                    ("σzz", zz),
                    ("τxy", xy),
                    ("τyz", yz),
                    ("τxz", xz),
                ] {
                    out.push((l.to_string(), v));
                }
                let s = Stress {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    sigma_xx: xx,
                    sigma_yy: yy,
                    sigma_zz: zz,
                    tau_xy: xy,
                    tau_yz: yz,
                    tau_xz: xz,
                };
                out.push(("von Mises".to_string(), s.von_mises()));
            }
            (name, 1) => out.push((name.to_string(), values[0])),
            (name, 3) => {
                for (axis, v) in ["x", "y", "z"].iter().zip(values) {
                    out.push((format!("{name} {axis}"), *v));
                }
            }
            _ => {}
        }
    }
    out
}

/// Le Node le plus proche du curseur (pixels), à `radius` pixels au plus. `projected` : (x, y,
/// profondeur) de chaque Node dans la vue.
pub fn nearest_node(projected: &[[f64; 3]], cursor: [f64; 2], radius: f64) -> Option<usize> {
    projected
        .iter()
        .enumerate()
        .map(|(i, p)| (i, (p[0] - cursor[0]).hypot(p[1] - cursor[1]), p[2]))
        .filter(|(_, d, _)| *d <= radius)
        .min_by(|a, b| a.1.total_cmp(&b.1).then(a.2.total_cmp(&b.2)))
        .map(|(i, ..)| i)
}

/// Courbes nommées : (étiquette, points (temps, valeur)).
pub type NamedCurves = Vec<(String, Vec<[f64; 2]>)>;

/// Les forces sur les surfaces rigides, lues dans le CSV de `maille solve` (`time,contact,fx,fy,fz`) :
/// une courbe (temps, valeur) par contact et par composante, étiquetée « contact c · fy ».
pub fn contact_curves(csv: &str) -> Result<NamedCurves, Box<dyn Error>> {
    let mut by_contact: BTreeMap<usize, Vec<[f64; 4]>> = BTreeMap::new();
    for (n, line) in csv.lines().enumerate().skip(1) {
        let cells: Vec<&str> = line.split(',').collect();
        let bad = || format!("ligne {} du CSV des contacts : « {line} »", n + 1);
        if cells.len() != 5 {
            return Err(bad().into());
        }
        let contact: usize = cells[1].trim().parse().map_err(|_| bad())?;
        let mut row = [0.0; 4];
        for (k, cell) in [0, 2, 3, 4].into_iter().zip(&mut row) {
            *cell = cells[k].trim().parse().map_err(|_| bad())?;
        }
        by_contact.entry(contact).or_default().push(row);
    }
    let mut curves = Vec::new();
    for (contact, rows) in by_contact {
        for (k, component) in ["fx", "fy", "fz"].into_iter().enumerate() {
            curves.push((
                format!("contact {contact} · {component}"),
                rows.iter().map(|r| [r[0], r[k + 1]]).collect(),
            ));
        }
    }
    Ok(curves)
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use maille_core::mesh::Node;
    use maille_io::msh::CellKind;
    use maille_io::vtk::PointArray;

    use super::*;

    fn results(arrays: Vec<PointArray>) -> VtuResults {
        VtuResults {
            nodes: vec![
                Node::new_2d(0.0, 0.0),
                Node::new_2d(1.0, 0.0),
                Node::new_2d(0.0, 1.0),
            ],
            kind: CellKind::Tri3,
            cells: vec![vec![0, 1, 2]],
            arrays,
            units: None,
            large_strain: false,
        }
    }

    fn array(name: &str, components: usize, values: Vec<f64>) -> PointArray {
        PointArray {
            name: name.to_string(),
            components,
            values,
        }
    }

    #[test]
    fn readings_give_each_component_and_the_derived_scalars() {
        let r = results(vec![
            array(
                "displacement",
                3,
                vec![0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0],
            ),
            array(
                "stress",
                6,
                [
                    vec![0.0; 6],
                    vec![100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    vec![0.0; 6],
                ]
                .concat(),
            ),
            array("temperature", 1, vec![10.0, 20.0, 30.0]),
        ]);
        let at = |label: &str| {
            readings(&r, 1)
                .into_iter()
                .find(|(l, _)| l == label)
                .unwrap_or_else(|| panic!("{label}"))
                .1
        };
        assert_eq!((at("ux"), at("uy"), at("|u|")), (3.0, 4.0, 5.0));
        assert_eq!((at("σxx"), at("von Mises")), (100.0, 100.0));
        assert_eq!(at("temperature"), 20.0);
        // un Node hors tableau n'a rien
        assert!(readings(&r, 7).is_empty());
    }

    #[test]
    fn the_nearest_node_is_picked_within_the_radius() {
        let nodes = [[10.0, 10.0, 0.5], [20.0, 10.0, 0.5], [21.0, 10.0, 0.1]];
        assert_eq!(nearest_node(&nodes, [11.0, 10.0], 12.0), Some(0));
        assert_eq!(nearest_node(&nodes, [20.4, 10.0], 12.0), Some(1));
        assert_eq!(nearest_node(&nodes, [100.0, 100.0], 12.0), None);
        // à égalité de distance, le plus proche de l'œil
        let tie = [[0.0, 0.0, 0.9], [0.0, 0.0, 0.2]];
        assert_eq!(nearest_node(&tie, [0.0, 0.0], 5.0), Some(1));
    }

    #[test]
    fn contact_forces_become_one_curve_per_contact_and_component() {
        let csv =
            "time,contact,fx,fy,fz\n0e0,0,0e0,0e0,0e0\n5e-1,0,-1e0,4e2,0e0\n1e0,0,-2e0,9e2,0e0\n";
        let curves = contact_curves(csv).unwrap();
        assert_eq!(curves.len(), 3);
        assert_eq!(curves[1].0, "contact 0 · fy");
        assert_eq!(curves[1].1, vec![[0.0, 0.0], [0.5, 400.0], [1.0, 900.0]]);
        let bad = contact_curves("time,contact,fx,fy,fz\n0,0,oops,0,0\n");
        assert!(bad.unwrap_err().to_string().contains("ligne 2"));
    }

    /// Calcule Hertz dans un dossier neuf, comme le ferait « Calculer » ; rend le `.pvd`.
    fn hertz() -> PathBuf {
        static DONE: std::sync::OnceLock<PathBuf> = std::sync::OnceLock::new();
        DONE.get_or_init(solve_hertz).clone()
    }

    fn solve_hertz() -> PathBuf {
        let dir = std::env::temp_dir().join(format!("maille-series-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for f in ["hertz.toml", "hertz.msh"] {
            let from = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../examples")
                .join(f);
            fs::copy(from, dir.join(f)).unwrap();
        }
        let files = maille_io::solving::solve_file(
            &dir.join("hertz.toml"),
            None,
            &mut |_| maille_core::model::Flow::Continue,
            &mut |_| {},
        )
        .unwrap();
        files
            .into_iter()
            .find(|p| p.extension().is_some_and(|e| e == "pvd"))
            .unwrap()
    }

    #[test]
    fn the_probe_values_are_the_ones_of_the_vtu() {
        let pvd = hertz();
        let mut series = Series::open(&pvd).unwrap();
        assert_eq!(series.len(), 9, "t = 0 puis huit incréments");
        assert!((series.time(8) - 1.0).abs() < 1e-12);
        let last = series.len() - 1;
        // lecture indépendante du dernier .vtu
        let direct =
            vtk::read_results(&fs::read_to_string(pvd.with_file_name("hertz_0008.vtu")).unwrap())
                .unwrap();
        let von_mises: Vec<f64> = direct
            .arrays
            .iter()
            .find(|a| a.name == "stress")
            .unwrap()
            .values
            .chunks(6)
            .map(|c| {
                Stress {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    sigma_xx: c[0],
                    sigma_yy: c[1],
                    sigma_zz: c[2],
                    tau_xy: c[3],
                    tau_yz: c[4],
                    tau_xz: c[5],
                }
                .von_mises()
            })
            .collect();
        // le Node le plus contraint : la sonde y lit exactement le même von Mises
        let node = (0..von_mises.len())
            .max_by(|&a, &b| von_mises[a].total_cmp(&von_mises[b]))
            .unwrap();
        let got = readings(series.get(last).unwrap(), node);
        let vm = got.iter().find(|(l, _)| l == "von Mises").unwrap().1;
        assert_eq!(vm, von_mises[node]);
        // la courbe de la sonde finit sur cette valeur et commence à zéro
        let curve = series.curve(node, "von Mises").unwrap();
        assert_eq!(curve.len(), 9);
        assert_eq!(curve[0], [0.0, 0.0]);
        assert_eq!(curve[8], [1.0, von_mises[node]]);
        assert!(series.curve(node, "inconnu").is_err());
    }

    #[test]
    fn the_contact_curve_ends_on_the_rigid_surface_force_of_the_last_increment() {
        let pvd = hertz();
        let series = Series::open(&pvd).unwrap();
        let csv = fs::read_to_string(series.contacts_csv().expect("hertz_contacts.csv")).unwrap();
        let curves = contact_curves(&csv).unwrap();
        let fy = &curves
            .iter()
            .find(|(l, _)| l == "contact 0 · fy")
            .unwrap()
            .1;
        assert_eq!(fy.len(), 9);
        assert_eq!(fy[0], [0.0, 0.0]);
        // la force croît avec l'enfoncement
        assert!(fy.windows(2).all(|w| w[1][1] > w[0][1]));
    }

    #[test]
    fn the_mesh_of_the_results_is_found_and_checked() {
        let pvd = hertz();
        let mut series = Series::open(&pvd).unwrap();
        let path = series.mesh_path().unwrap();
        assert_eq!(path, pvd.with_file_name("hertz.msh"));
        let file = series.read_mesh(&path).unwrap();
        assert!(!file.sets.is_empty());
        // un Increment seul retrouve aussi son modèle
        let vtu = Series::open(&pvd.with_file_name("hertz_0003.vtu")).unwrap();
        assert_eq!(vtu.mesh_path(), Some(path));
        // un autre maillage est refusé
        let other = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/plaque.msh");
        let err = series.read_mesh(&other).unwrap_err().to_string();
        assert!(err.contains("n'est pas le maillage"), "{err}");
    }

    #[test]
    fn opening_a_vtu_or_a_wrong_file_is_checked() {
        let pvd = hertz();
        let vtu = Series::open(&pvd.with_file_name("hertz_0003.vtu")).unwrap();
        assert_eq!((vtu.len(), vtu.time(0)), (1, 0.0));
        assert!(Series::open(Path::new("/nonexistent/x.pvd")).is_err());
        assert!(Series::open(&pvd.with_file_name("hertz.toml")).is_err());
    }
}
