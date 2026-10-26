//! Triangulation de Delaunay 3D incrémentale (Bowyer-Watson), sur le modèle du 2D : un
//! super-tétraèdre couvre le domaine, chaque point inséré retire la cavité des tétraèdres dont la
//! sphère circonscrite le contient et la retriangule en étoile autour de lui. Pas encore de faces
//! contraintes : la récupération du bord vient ensuite.

use std::collections::{HashMap, HashSet};

use super::predicates::{distance3, insphere, orient3d};
use crate::error::MailleError;

/// Les 4 premiers points sont les sommets du super-tétraèdre.
const SUPER: usize = 4;
/// Taille du super-tétraèdre rapportée à celle du domaine : grand, pour qu'il ne cache pas de
/// face de l'enveloppe convexe entre points presque coplanaires.
const SUPER_SCALE: f64 = 100.0;

/// Face k d'un tétraèdre = celle qui est opposée au sommet k, ordonnée pour que (face, p) soit un
/// tétraèdre positif quand p est du même côté que le sommet k (face tournée vers l'intérieur).
const FACE: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// Tétraèdre positif (`orient3d(v0, v1, v2, v3) > 0`). `n[i]` est le voisin à travers la face
/// opposée au sommet i, `None` sur le bord du super-tétraèdre. Un tétraèdre retiré reste dans le
/// tableau (`alive = false`) pour ne pas renuméroter les autres.
#[derive(Debug, Clone, Copy)]
struct Tet {
    v: [usize; 4],
    n: [Option<usize>; 4],
    alive: bool,
}

/// Cavité d'une insertion : tétraèdres retirés, faces de leur bord avec le tétraèdre d'en face,
/// arêtes orientées de ces faces (arête → face et position).
type Cavity = (
    Vec<usize>,
    Vec<([usize; 3], Option<usize>)>,
    HashMap<(usize, usize), (usize, usize)>,
);

#[derive(Debug, Clone)]
pub struct Delaunay3 {
    points: Vec<[f64; 3]>,
    tets: Vec<Tet>,
    /// Tétraèdre de départ de la recherche du prochain point : le dernier créé.
    hint: usize,
    /// Un tétraèdre vivant qui contient chaque point.
    vertex_tet: Vec<usize>,
}

impl Delaunay3 {
    /// Triangulation vide (le seul super-tétraèdre) d'un domaine contenu dans la boîte [min, max].
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        let c: [f64; 3] = std::array::from_fn(|k| (min[k] + max[k]) / 2.0);
        let size = (0..3).map(|k| max[k] - min[k]).fold(0.0, f64::max);
        let m = SUPER_SCALE * if size > 0.0 { size } else { 1.0 };
        let at = |s: [f64; 3]| [c[0] + m * s[0], c[1] + m * s[1], c[2] + m * s[2]];
        let mut points = vec![
            at([1.0, 1.0, 1.0]),
            at([1.0, -1.0, -1.0]),
            at([-1.0, 1.0, -1.0]),
            at([-1.0, -1.0, 1.0]),
        ];
        if orient3d(points[0], points[1], points[2], points[3]) < 0.0 {
            points.swap(2, 3);
        }
        Self {
            points,
            tets: vec![Tet {
                v: [0, 1, 2, 3],
                n: [None; 4],
                alive: true,
            }],
            hint: 0,
            vertex_tet: vec![0; SUPER],
        }
    }

    /// Insère `p` et renvoie son indice (0 pour le premier point inséré). `MeshingFailed` si `p`
    /// est confondu avec un point existant, hors du super-tétraèdre, ou si l'arrondi rend la
    /// cavité incohérente ; la triangulation n'est alors pas modifiée.
    pub fn insert(&mut self, p: [f64; 3]) -> Result<usize, MailleError> {
        let start = self.locate(p)?;
        let scale = self.edge_scale(start);
        if self.tets[start]
            .v
            .iter()
            .any(|&v| distance3(self.points[v], p) <= 1e-12 * scale)
        {
            return Err(MailleError::MeshingFailed);
        }

        // d'abord sans les tétraèdres dont la sphère passe par p ; si la cavité est incohérente
        // (points cosphériques en nombre), avec eux
        let (cavity, boundary, edges) = [false, true]
            .into_iter()
            .find_map(|ties| self.cavity(start, p, ties))
            .ok_or(MailleError::MeshingFailed)?;
        // position de chaque face du bord dans le tétraèdre d'en face : celle du sommet qui n'en fait pas partie
        let outside: Vec<Option<(usize, usize)>> = boundary
            .iter()
            .map(|&([a, b, c], o)| match o {
                None => Ok(None),
                Some(o) => self.tets[o]
                    .v
                    .iter()
                    .position(|v| ![a, b, c].contains(v))
                    .map(|j| Some((o, j)))
                    .ok_or(MailleError::MeshingFailed),
            })
            .collect::<Result<_, _>>()?;

        let pi = self.points.len();
        self.points.push(p);
        self.vertex_tet.push(0);
        for &t in &cavity {
            self.tets[t].alive = false;
        }
        let first = self.tets.len();
        for (k, (&([a, b, c], o), &link)) in boundary.iter().zip(&outside).enumerate() {
            let mut n = [None, None, None, o];
            for (slot, edge) in n.iter_mut().zip([(c, b), (a, c), (b, a)]) {
                *slot = Some(first + edges[&edge].0);
            }
            self.tets.push(Tet {
                v: [a, b, c, pi],
                n,
                alive: true,
            });
            for v in [a, b, c, pi] {
                self.vertex_tet[v] = first + k;
            }
            if let Some((o, j)) = link {
                self.tets[o].n[j] = Some(first + k);
            }
        }
        self.hint = first;
        Ok(pi - SUPER)
    }

    /// Cavité de `p` à partir du tétraèdre `start` : tétraèdres à retirer, faces (a, b, c) de son
    /// bord et arêtes orientées de ces faces. `ties` : y mettre aussi les tétraèdres dont la
    /// sphère passe par p. `None` si la cavité n'est pas étoilée vue de p.
    fn cavity(&self, start: usize, p: [f64; 3], ties: bool) -> Option<Cavity> {
        // cavité : tétraèdres connexes au tétraèdre de départ dont la sphère circonscrite contient p ;
        // faces (a, b, c) de son bord, tournées vers la cavité, avec le tétraèdre de l'autre côté
        let mut cavity = vec![start];
        let mut inside = HashSet::from([start]);
        let mut boundary = Vec::new();
        let mut k = 0;
        while k < cavity.len() {
            let t = self.tets[cavity[k]];
            k += 1;
            for (i, face) in FACE.iter().enumerate() {
                match t.n[i] {
                    Some(o) if inside.contains(&o) => {}
                    Some(o) if self.in_circumsphere(o, p, ties) => {
                        inside.insert(o);
                        cavity.push(o);
                    }
                    other => boundary.push((face.map(|j| t.v[j]), other)),
                }
            }
        }

        // cohérence avant toute modification : toutes les faces du bord visibles de p (les nouveaux
        // tétraèdres sont positifs), et aucun sommet enfermé dans la cavité (il serait perdu)
        let boundary_vertices: HashSet<usize> = boundary.iter().flat_map(|&(f, _)| f).collect();
        let cavity_vertices: HashSet<usize> = cavity.iter().flat_map(|&t| self.tets[t].v).collect();
        if boundary_vertices != cavity_vertices
            || boundary.iter().any(|&([a, b, c], _)| {
                orient3d(self.points[a], self.points[b], self.points[c], p) <= 0.0
            })
        {
            return None;
        }
        // arêtes orientées du bord : chacune doit avoir son inverse dans une autre face
        let mut edges: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        for (k, &([a, b, c], _)) in boundary.iter().enumerate() {
            // tétraèdre [a, b, c, p] : l'arête b → c est opposée à a, c → a à b, a → b à c
            edges.insert((b, c), (k, 0));
            edges.insert((c, a), (k, 1));
            edges.insert((a, b), (k, 2));
        }
        if edges.len() != 3 * boundary.len()
            || edges.keys().any(|&(x, y)| !edges.contains_key(&(y, x)))
        {
            return None;
        }
        Some((cavity, boundary, edges))
    }

    /// Nombre de points insérés.
    #[cfg(test)]
    pub fn point_count(&self) -> usize {
        self.points.len() - SUPER
    }

    /// Tétraèdres positifs sans sommet du super-tétraèdre (indices de points insérés). Ils
    /// remplissent l'enveloppe convexe des points, à quelques tétraèdres plats près au bord.
    pub fn tetrahedra(&self) -> Vec<[usize; 4]> {
        self.tets
            .iter()
            .filter(|t| t.alive && t.v.iter().all(|&v| v >= SUPER))
            .map(|t| t.v.map(|v| v - SUPER))
            .collect()
    }

    /// Vrai si les points insérés `a` et `b` sont reliés par une arête.
    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        let (a, b) = (a + SUPER, b + SUPER);
        self.tets_around(a).any(|t| self.tets[t].v.contains(&b))
    }

    /// Vrai si les points insérés `a`, `b` et `c` forment une face d'un tétraèdre.
    pub fn has_face(&self, a: usize, b: usize, c: usize) -> bool {
        let (a, b, c) = (a + SUPER, b + SUPER, c + SUPER);
        self.tets_around(a)
            .any(|t| self.tets[t].v.contains(&b) && self.tets[t].v.contains(&c))
    }

    /// Tétraèdres vivants qui contiennent le point `v` (indice brut), de proche en proche par les
    /// faces qui le contiennent.
    fn tets_around(&self, v: usize) -> impl Iterator<Item = usize> {
        let mut seen = HashSet::from([self.vertex_tet[v]]);
        let mut stack = vec![self.vertex_tet[v]];
        let mut out = Vec::new();
        while let Some(t) = stack.pop() {
            out.push(t);
            let tet = self.tets[t];
            for (i, &o) in tet.n.iter().enumerate() {
                if tet.v[i] != v {
                    if let Some(o) = o.filter(|o| seen.insert(*o)) {
                        stack.push(o);
                    }
                }
            }
        }
        out.into_iter()
    }

    fn corners(&self, t: usize) -> [[f64; 3]; 4] {
        self.tets[t].v.map(|v| self.points[v])
    }

    fn in_circumsphere(&self, t: usize, p: [f64; 3], ties: bool) -> bool {
        let [a, b, c, d] = self.corners(t);
        let s = insphere(a, b, c, d, p);
        s > 0.0 || (ties && s == 0.0)
    }

    /// Plus grande arête du tétraèdre `t` : l'échelle du test de point confondu.
    fn edge_scale(&self, t: usize) -> f64 {
        let c = self.corners(t);
        let mut scale = 0.0_f64;
        for i in 0..4 {
            for j in i + 1..4 {
                scale = scale.max(distance3(c[i], c[j]));
            }
        }
        scale
    }

    /// Tétraèdre qui contient `p` (sur une face ou une arête : l'un des tétraèdres voisins), par
    /// marche depuis le dernier créé. `MeshingFailed` si `p` est hors du super-tétraèdre.
    fn locate(&self, p: [f64; 3]) -> Result<usize, MailleError> {
        let mut t = self.hint;
        // la marche peut tourner en rond sur des cas dégénérés : au plus un tour de tous les tétraèdres
        for _ in 0..=self.tets.len() {
            let tet = self.tets[t];
            let c = self.corners(t);
            let exit = FACE
                .iter()
                .position(|f| orient3d(c[f[0]], c[f[1]], c[f[2]], p) < 0.0);
            match exit {
                None => return Ok(t),
                Some(i) => t = tet.n[i].ok_or(MailleError::MeshingFailed)?,
            }
        }
        Err(MailleError::MeshingFailed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Générateur congruentiel déterministe, dans [0, 1[.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn volume(d: &Delaunay3, v: [usize; 4]) -> f64 {
        let c = v.map(|i| d.points[i]);
        orient3d(c[0], c[1], c[2], c[3]) / 6.0
    }

    /// Structure : tétraèdres positifs, voisinage réciproque par une face commune, volumes qui
    /// remplissent exactement le super-tétraèdre, sphères circonscrites vides (Delaunay).
    fn check(d: &Delaunay3, super_volume: f64) {
        let mut total = 0.0;
        for (t, tet) in d.tets.iter().enumerate().filter(|(_, t)| t.alive) {
            let vol = volume(d, tet.v);
            assert!(vol > 0.0, "tétraèdre {t} non positif : {vol}");
            total += vol;
            for (i, face) in FACE.iter().enumerate() {
                let Some(o) = tet.n[i] else { continue };
                assert!(d.tets[o].alive, "voisin {o} de {t} retiré");
                let mut mine = face.map(|j| tet.v[j]);
                mine.sort_unstable();
                let j = d.tets[o].n.iter().position(|&x| x == Some(t));
                let j = j.unwrap_or_else(|| panic!("{o} ne pointe pas vers {t}"));
                let mut theirs = FACE[j].map(|k| d.tets[o].v[k]);
                theirs.sort_unstable();
                assert_eq!(mine, theirs, "face commune de {t} et {o}");
            }
            let [a, b, c, e] = d.corners(t);
            for (i, &q) in d.points.iter().enumerate().skip(SUPER) {
                if !tet.v.contains(&i) {
                    assert!(
                        insphere(a, b, c, e, q) <= 0.0,
                        "point {i} dans la sphère de {t}"
                    );
                }
            }
        }
        assert!(
            (total - super_volume).abs() < 1e-9 * super_volume,
            "{total} contre {super_volume}"
        );
    }

    fn super_volume(d: &Delaunay3) -> f64 {
        volume(d, [0, 1, 2, 3])
    }

    fn real_volume(d: &Delaunay3) -> f64 {
        d.tetrahedra()
            .iter()
            .map(|&t| volume(d, t.map(|v| v + SUPER)))
            .sum()
    }

    #[test]
    fn random_points_give_a_delaunay_tetrahedralization() {
        let mut d = Delaunay3::new([0.0; 3], [1.0; 3]);
        let sv = super_volume(&d);
        let mut rng = Rng(7);
        for _ in 0..200 {
            d.insert([rng.next(), rng.next(), rng.next()]).unwrap();
        }
        assert_eq!(d.point_count(), 200);
        check(&d, sv);
    }

    #[test]
    fn cube_corners_and_random_interior_fill_the_cube() {
        let mut d = Delaunay3::new([0.0; 3], [1.0; 3]);
        let sv = super_volume(&d);
        for i in 0..8 {
            d.insert([(i & 1) as f64, ((i >> 1) & 1) as f64, ((i >> 2) & 1) as f64])
                .unwrap();
        }
        let mut rng = Rng(11);
        for _ in 0..100 {
            d.insert([rng.next(), rng.next(), rng.next()]).unwrap();
        }
        check(&d, sv);
        let v = real_volume(&d);
        assert!((v - 1.0).abs() < 1e-9, "volume du cube : {v}");
    }

    /// Points cosphériques en nombre (grille régulière) : le résultat reste un Delaunay valide.
    #[test]
    fn regular_grid_is_triangulated() {
        let mut d = Delaunay3::new([0.0; 3], [4.0; 3]);
        let sv = super_volume(&d);
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    d.insert([i as f64, j as f64, k as f64]).unwrap();
                }
            }
        }
        check(&d, sv);
        let v = real_volume(&d);
        assert!((v - 64.0).abs() < 1e-9 * 64.0, "volume de la grille : {v}");
    }

    #[test]
    fn duplicate_or_outside_point_fails_without_changing_anything() {
        let mut d = Delaunay3::new([0.0; 3], [1.0; 3]);
        d.insert([0.2, 0.3, 0.4]).unwrap();
        d.insert([0.8, 0.5, 0.1]).unwrap();
        let before = d.tets.len();
        assert!(matches!(
            d.insert([0.2, 0.3, 0.4]),
            Err(MailleError::MeshingFailed)
        ));
        assert!(matches!(
            d.insert([1e4, 0.0, 0.0]),
            Err(MailleError::MeshingFailed)
        ));
        assert_eq!(d.tets.len(), before);
        assert_eq!(d.point_count(), 2);
    }
}
