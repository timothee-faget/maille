//! Triangulation de Delaunay incrémentale (Bowyer-Watson). On part d'un super-triangle qui couvre
//! le domaine, puis chaque point inséré retire la cavité des triangles dont le cercle circonscrit
//! le contient et la retriangule en étoile autour de lui. Les arêtes contraintes (morceaux de bord)
//! arrêtent la cavité : une fois présentes, elles restent (Delaunay contraint).

use std::collections::{HashMap, HashSet};

use super::predicates::{incircle, orient2d};
use crate::error::MailleError;

/// Les 3 premiers points sont les sommets du super-triangle.
const SUPER: usize = 3;
/// Taille du super-triangle rapportée à celle du domaine. Plus il est grand, moins il risque de
/// cacher une arête de l'enveloppe convexe entre points presque alignés.
const SUPER_SCALE: f64 = 100.0;

/// Triangle antihoraire. `n[i]` est le voisin à travers l'arête opposée au sommet i, soit
/// (v[i + 1], v[i + 2]) ; `None` sur le bord du super-triangle. Un triangle retiré reste dans le
/// tableau (`alive = false`) pour ne pas renuméroter les autres.
#[derive(Debug, Clone, Copy)]
struct Tri {
    v: [usize; 3],
    n: [Option<usize>; 3],
    alive: bool,
}

#[derive(Debug, Clone)]
pub struct Delaunay {
    points: Vec<[f64; 2]>,
    tris: Vec<Tri>,
    /// Triangle de départ de la recherche du prochain point : le dernier créé.
    hint: usize,
    /// Un triangle vivant qui contient chaque point.
    vertex_tri: Vec<usize>,
    /// Arêtes contraintes (indices de points bruts, le plus petit d'abord).
    constrained: HashSet<(usize, usize)>,
}

impl Delaunay {
    /// Triangulation vide (le seul super-triangle) d'un domaine contenu dans la boîte [min, max].
    pub fn new(min: [f64; 2], max: [f64; 2]) -> Self {
        let (cx, cy) = ((min[0] + max[0]) / 2.0, (min[1] + max[1]) / 2.0);
        let size = (max[0] - min[0]).max(max[1] - min[1]);
        let m = SUPER_SCALE * if size > 0.0 { size } else { 1.0 };
        Self {
            points: vec![[cx - m, cy - m], [cx + m, cy - m], [cx, cy + m]],
            tris: vec![Tri {
                v: [0, 1, 2],
                n: [None; 3],
                alive: true,
            }],
            hint: 0,
            vertex_tri: vec![0; SUPER],
            constrained: HashSet::new(),
        }
    }

    /// Insère `p` et renvoie son indice (0 pour le premier point inséré). `MeshingFailed` si `p`
    /// est confondu avec un point existant ou si l'arrondi rend la cavité incohérente ; la
    /// triangulation n'est alors pas modifiée.
    pub fn insert(&mut self, p: [f64; 2]) -> Result<usize, MailleError> {
        let start = self.locate(p)?;
        let scale = self.edge_scale(start);
        if self.tris[start]
            .v
            .iter()
            .any(|&v| distance(self.points[v], p) <= 1e-12 * scale)
        {
            return Err(MailleError::MeshingFailed);
        }

        // cavité : triangles connexes au triangle de départ dont le cercle circonscrit contient p
        let mut cavity = vec![start];
        let mut inside = HashSet::from([start]);
        let mut boundary = Vec::new(); // arêtes (a, b) du bord de la cavité et triangle de l'autre côté
        let mut k = 0;
        while k < cavity.len() {
            let t = self.tris[cavity[k]];
            k += 1;
            for i in 0..3 {
                let (a, b) = (t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
                match t.n[i] {
                    Some(o) if inside.contains(&o) => {}
                    Some(o) if !self.is_constrained(a, b) && self.in_circumcircle(o, p) => {
                        inside.insert(o);
                        cavity.push(o);
                    }
                    other => boundary.push((a, b, other)),
                }
            }
        }

        // cohérence avant toute modification : la cavité doit être un polygone simple sans sommet
        // intérieur (T = B − 2), visible depuis p, et sans triangle atteint par un autre chemin
        // derrière une arête contrainte (elle serait détruite)
        let starts: HashSet<usize> = boundary.iter().map(|&(a, _, _)| a).collect();
        if boundary.len() != cavity.len() + 2
            || boundary
                .iter()
                .any(|&(_, _, o)| o.is_some_and(|o| inside.contains(&o)))
            || starts.len() != boundary.len()
            || boundary
                .iter()
                .any(|&(a, b, _)| orient2d(self.points[a], self.points[b], p) <= 0.0)
        {
            return Err(MailleError::MeshingFailed);
        }
        // position de chaque arête du bord dans le triangle d'en face, qui la parcourt en sens inverse
        let outside: Vec<Option<(usize, usize)>> = boundary
            .iter()
            .map(|&(a, b, o)| match o {
                None => Ok(None),
                Some(o) => self
                    .edge_of(o, b, a)
                    .map(|j| Some((o, j)))
                    .ok_or(MailleError::MeshingFailed),
            })
            .collect::<Result<_, _>>()?;

        let pi = self.points.len();
        self.points.push(p);
        self.vertex_tri.push(0);
        for &t in &cavity {
            self.tris[t].alive = false;
        }
        let first = self.tris.len();
        let (mut by_start, mut by_end) = (HashMap::new(), HashMap::new());
        for (k, (&(a, b, o), &link)) in boundary.iter().zip(&outside).enumerate() {
            let id = first + k;
            self.tris.push(Tri {
                v: [a, b, pi],
                n: [None, None, o],
                alive: true,
            });
            by_start.insert(a, id);
            by_end.insert(b, id);
            for v in [a, b, pi] {
                self.vertex_tri[v] = id;
            }
            if let Some((o, j)) = link {
                self.tris[o].n[j] = Some(id);
            }
        }
        // l'arête (b, p) est partagée avec le triangle qui commence en b, (p, a) avec celui qui finit en a
        for k in 0..boundary.len() {
            let [a, b, _] = self.tris[first + k].v;
            self.tris[first + k].n[0] = Some(by_start[&b]);
            self.tris[first + k].n[1] = Some(by_end[&a]);
        }
        self.hint = first;
        Ok(pi - SUPER)
    }

    /// Nombre de points insérés.
    pub fn point_count(&self) -> usize {
        self.points.len() - SUPER
    }

    pub fn point(&self, i: usize) -> [f64; 2] {
        self.points[i + SUPER]
    }

    /// Vrai si les points insérés `a` et `b` sont reliés par une arête.
    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        self.directed_edge(a + SUPER, b + SUPER).is_some()
    }

    /// Marque l'arête (a, b) comme contrainte (`on`) ou la libère.
    pub fn set_constrained(&mut self, a: usize, b: usize, on: bool) {
        let key = (a.min(b) + SUPER, a.max(b) + SUPER);
        if on {
            self.constrained.insert(key);
        } else {
            self.constrained.remove(&key);
        }
    }

    /// Triangles du domaine, en indices de points insérés et dans l'ordre de leurs sommets :
    /// la face k va du sommet k au sommet k + 1. `boundary` donne les arêtes de bord (a, b), le
    /// domaine à leur gauche. On vide l'extérieur de proche en proche à partir du triangle à droite
    /// de chacune, sans traverser d'arête contrainte. `MeshingFailed` si une arête de bord manque
    /// ou si le vidage atteint un triangle du domaine par une fuite (bord non fermé).
    pub fn domain(
        &self,
        boundary: impl Iterator<Item = (usize, usize)>,
    ) -> Result<Vec<[usize; 3]>, MailleError> {
        let mut outside = vec![false; self.tris.len()];
        let mut stack = Vec::new();
        for (a, b) in boundary {
            let (t, _) = self
                .directed_edge(b + SUPER, a + SUPER)
                .ok_or(MailleError::MeshingFailed)?;
            if !outside[t] {
                outside[t] = true;
                stack.push(t);
            }
        }
        while let Some(t) = stack.pop() {
            let tri = self.tris[t];
            for i in 0..3 {
                let (a, b) = (tri.v[(i + 1) % 3], tri.v[(i + 2) % 3]);
                if let Some(o) = tri.n[i] {
                    if !outside[o] && !self.is_constrained(a, b) {
                        outside[o] = true;
                        stack.push(o);
                    }
                }
            }
        }
        let inside: Vec<&Tri> = self
            .tris
            .iter()
            .zip(&outside)
            .filter(|(t, &out)| t.alive && !out)
            .map(|(t, _)| t)
            .collect();
        if inside.iter().any(|t| t.v.iter().any(|&v| v < SUPER)) {
            return Err(MailleError::MeshingFailed);
        }
        Ok(inside.iter().map(|t| t.v.map(|v| v - SUPER)).collect())
    }

    /// Triangles qui ne touchent pas le super-triangle, en indices de points insérés.
    #[cfg(test)]
    pub fn triangles(&self) -> impl Iterator<Item = [usize; 3]> + '_ {
        self.tris
            .iter()
            .filter(|t| t.alive && t.v.iter().all(|&v| v >= SUPER))
            .map(|t| t.v.map(|v| v - SUPER))
    }

    /// Triangle qui contient p (éventuellement sur une arête), par marche : on traverse l'arête dont
    /// p est de l'autre côté jusqu'à ne plus en trouver. Termine dans une triangulation de Delaunay.
    fn locate(&self, p: [f64; 2]) -> Result<usize, MailleError> {
        let mut t = if self.tris[self.hint].alive {
            self.hint
        } else {
            self.tris
                .iter()
                .rposition(|t| t.alive)
                .ok_or(MailleError::MeshingFailed)?
        };
        for _ in 0..self.tris.len() {
            let tri = self.tris[t];
            let crossed = (0..3).find(|&i| {
                let (a, b) = (tri.v[(i + 1) % 3], tri.v[(i + 2) % 3]);
                orient2d(self.points[a], self.points[b], p) < 0.0
            });
            match crossed {
                None => return Ok(t),
                // hors du super-triangle
                Some(i) => t = tri.n[i].ok_or(MailleError::MeshingFailed)?,
            }
        }
        Err(MailleError::MeshingFailed)
    }

    /// Triangle vivant qui parcourt l'arête a → b (indices bruts), et position k de a dans ce
    /// triangle (v[k] = a, v[k + 1] = b). On tourne autour de a de voisin en voisin.
    fn directed_edge(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        let first = self.vertex_tri[a];
        let mut t = first;
        loop {
            let v = self.tris[t].v;
            let k = v.iter().position(|&x| x == a)?;
            if v[(k + 1) % 3] == b {
                return Some((t, k));
            }
            // voisin à travers l'arête (a, v[k + 1]), opposée au sommet k + 2 : il contient aussi a
            t = self.tris[t].n[(k + 2) % 3]?;
            if t == first {
                return None;
            }
        }
    }

    fn is_constrained(&self, a: usize, b: usize) -> bool {
        self.constrained.contains(&(a.min(b), a.max(b)))
    }

    fn in_circumcircle(&self, t: usize, p: [f64; 2]) -> bool {
        let [a, b, c] = self.tris[t].v.map(|v| self.points[v]);
        incircle(a, b, c, p) > 0.0
    }

    /// Indice dans le triangle `t` de l'arête qui va de `a` à `b`.
    fn edge_of(&self, t: usize, a: usize, b: usize) -> Option<usize> {
        let v = self.tris[t].v;
        (0..3).find(|&i| v[(i + 1) % 3] == a && v[(i + 2) % 3] == b)
    }

    /// Plus longue arête du triangle `t`, échelle des tolérances de distance.
    fn edge_scale(&self, t: usize) -> f64 {
        let v = self.tris[t].v;
        (0..3)
            .map(|i| distance(self.points[v[i]], self.points[v[(i + 1) % 3]]))
            .fold(0.0, f64::max)
    }
}

fn distance(p: [f64; 2], q: [f64; 2]) -> f64 {
    (p[0] - q[0]).hypot(p[1] - q[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Points pseudo-aléatoires reproductibles dans [margin, 1 − margin]².
    fn random_points(n: usize, margin: f64) -> Vec<[f64; 2]> {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            margin + (1.0 - 2.0 * margin) * (state >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n).map(|_| [next(), next()]).collect()
    }

    fn triangulate(points: &[[f64; 2]]) -> Vec<[usize; 3]> {
        let mut d = Delaunay::new([0.0, 0.0], [1.0, 1.0]);
        for (i, &p) in points.iter().enumerate() {
            assert_eq!(d.insert(p).unwrap(), i);
        }
        d.triangles().collect()
    }

    /// Triangles antihoraires, cercles circonscrits vides, aire totale du carré unité.
    fn check_delaunay(points: &[[f64; 2]], tris: &[[usize; 3]]) {
        let mut area = 0.0;
        for t in tris {
            let [a, b, c] = t.map(|v| points[v]);
            assert!(orient2d(a, b, c) > 0.0, "triangle {t:?} pas antihoraire");
            area += orient2d(a, b, c) / 2.0;
            for (i, &p) in points.iter().enumerate() {
                assert!(
                    incircle(a, b, c, p) <= 0.0,
                    "point {i} dans le cercle de {t:?}"
                );
            }
        }
        assert!((area - 1.0).abs() < 1e-12, "aire {area}");
    }

    #[test]
    fn random_points_in_square() {
        let mut points = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        points.extend(random_points(300, 0.01));
        let tris = triangulate(&points);
        check_delaunay(&points, &tris);
        // Euler : 2n − 2 − h triangles pour n points dont h sur l'enveloppe convexe
        assert_eq!(tris.len(), 2 * points.len() - 2 - 4);
    }

    /// Grille régulière : 4 points cocirculaires dans chaque case, le cas limite des prédicats.
    #[test]
    fn regular_grid() {
        let n = 10;
        let points: Vec<[f64; 2]> = (0..=n)
            .flat_map(|j| (0..=n).map(move |i| [i as f64 * 0.1, j as f64 * 0.1]))
            .collect();
        let tris = triangulate(&points);
        check_delaunay(&points, &tris);
        assert_eq!(tris.len(), 2 * n * n);
    }

    #[test]
    fn point_on_existing_edge() {
        // (½, ½) tombe sur la diagonale commune des deux premiers triangles
        let points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]];
        let tris = triangulate(&points);
        check_delaunay(&points, &tris);
        assert_eq!(tris.len(), 4);
    }

    #[test]
    fn duplicate_point_is_refused() {
        let mut d = Delaunay::new([0.0, 0.0], [1.0, 1.0]);
        for p in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]] {
            d.insert(p).unwrap();
        }
        assert_eq!(d.insert([1.0, 0.0]), Err(MailleError::MeshingFailed));
        // la triangulation est intacte
        assert_eq!(d.triangles().count(), 1);
        assert_eq!(d.insert([1.0, 1.0]), Ok(3));
    }
}
