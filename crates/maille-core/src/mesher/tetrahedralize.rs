//! Mailleur 3D : tétraèdres conformes à une `Surface`. Avec une taille h (`TetOptions::size`), la
//! surface est d'abord remaillée à la taille h (`remesh`), puis l'intérieur est raffiné une fois le
//! bord récupéré (`refine3`). Récupération du bord : les sommets de la surface
//! entrent dans un Delaunay 3D ; toute arête ou face de la surface qui n'y figure pas est coupée
//! (milieu d'une arête, centre d'une face) jusqu'à ce que chaque morceau soit une arête ou une
//! face du Delaunay. Les points ajoutés sont sur la surface : elle ne change pas de forme. Les
//! tétraèdres dedans sont ceux dont le centre est dans la surface, et leur volume total doit égaler
//! celui de la surface (sinon `MeshingFailed`, jamais un maillage faux). Pas encore de contrôle de
//! la taille ni de la qualité des tétraèdres.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::brep::Brep;
use super::brepmesh::triangulate;
use super::delaunay3::Delaunay3;
use super::predicates::orient3d;
use super::refine3::{self, SizeField3};
use super::remesh::remesh;
use super::surface::Surface;
use super::{BoundarySet, MeshOptions, MAX_NODES};
use crate::error::MailleError;
use crate::material::positive;
use crate::mesh::Node;
use crate::shape::{TET10_EDGES, TET_FACES};

/// Passes de récupération au plus, et points ajoutés au plus par sommet de la surface : au-delà,
/// la récupération diverge (triangles très effilés, angles très aigus).
const MAX_PASSES: usize = 100;
const MAX_GROWTH: usize = 50;
/// Écart relatif toléré entre le volume des tétraèdres et celui de la surface.
const VOLUME_TOLERANCE: f64 = 1e-8;

/// Résultat du mailleur 3D : Nodes, connectivité des tétraèdres (Tet4 pour N = 4) et bords nommés
/// (les patches de la `Surface`). `Mesh::new` en fait un Mesh une fois les Elements créés.
#[derive(Debug, Clone, PartialEq)]
pub struct Tetrahedralization<const N: usize> {
    pub nodes: Vec<Node>,
    /// Connectivité des tétraèdres positifs (comme `Tet4::nodes`).
    pub tetrahedra: Vec<[usize; N]>,
    boundaries: BTreeMap<String, BoundarySet>,
}

impl<const N: usize> Tetrahedralization<N> {
    /// Noms des bords, triés.
    pub fn boundary_names(&self) -> impl Iterator<Item = &str> {
        self.boundaries.keys().map(String::as_str)
    }

    /// Bord nommé `name`. `UnknownBoundary` (avec la liste des noms) si aucun patch ne le porte.
    pub fn boundary(&self, name: &str) -> Result<&BoundarySet, MailleError> {
        self.boundaries
            .get(name)
            .ok_or_else(|| MailleError::UnknownBoundary {
                name: name.to_string(),
                known: self.boundaries.keys().cloned().collect(),
            })
    }
}

impl Tetrahedralization<4> {
    /// Maillage quadratique (Tet10) : un Node au milieu de chaque arête, partagé entre voisins,
    /// dans l'ordre de `Tet10::nodes` (arêtes 0-1, 1-2, 2-0, 0-3, 1-3, 2-3). Les milieux sont au
    /// milieu de la corde : le Tet10 n'accepte que des côtés droits. Les faces de bord gardent
    /// leurs numéros ; leurs Nodes milieux rejoignent les Nodes du bord.
    pub fn to_quadratic(self) -> Tetrahedralization<10> {
        let mut nodes = self.nodes;
        let mut mids: HashMap<(usize, usize), usize> = HashMap::new();
        let mut tetrahedra = Vec::with_capacity(self.tetrahedra.len());
        for t in &self.tetrahedra {
            let mut cell = [0; 10];
            cell[..4].copy_from_slice(t);
            for (k, &(a, b)) in TET10_EDGES.iter().enumerate() {
                let (a, b) = (t[a], t[b]);
                cell[4 + k] = *mids.entry(edge(a, b)).or_insert_with(|| {
                    let (p, q) = (nodes[a], nodes[b]);
                    nodes.push(Node::new(
                        (p.x + q.x) / 2.0,
                        (p.y + q.y) / 2.0,
                        (p.z + q.z) / 2.0,
                    ));
                    nodes.len() - 1
                });
            }
            tetrahedra.push(cell);
        }
        let mut boundaries = self.boundaries;
        for set in boundaries.values_mut() {
            for &(e, f) in &set.faces {
                let [a, b, c] = TET_FACES[f];
                for (x, y) in [(a, b), (b, c), (c, a)] {
                    let k = TET10_EDGES
                        .iter()
                        .position(|&e| e == (x, y) || e == (y, x))
                        .unwrap_or(0);
                    set.nodes.push(tetrahedra[e][4 + k]);
                }
            }
            set.nodes.sort_unstable();
            set.nodes.dedup();
        }
        Tetrahedralization {
            nodes,
            tetrahedra,
            boundaries,
        }
    }
}

/// Morceau de la surface : un triangle de `Surface::triangles`, ou une partie de lui après coupes.
#[derive(Debug, Clone, Copy)]
struct Piece {
    v: [usize; 3],
    triangle: usize,
}

/// Réglages du mailleur 3D.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TetOptions {
    /// Taille h, strictement positive : la surface est remaillée en triangles de rayon
    /// circonscrit au plus h/√3, l'intérieur raffiné jusqu'à un rayon circonscrit d'au plus
    /// h √6 / 4 (celui du tétraèdre régulier de côté h) et un rapport rayon / plus petite arête
    /// d'au plus 2, sauf près de la surface où un point qui l'empiéterait est refusé. `None` : ni
    /// remaillage ni raffinement, les tétraèdres ont la finesse des triangles de la surface.
    pub size: Option<f64>,
    /// Gradation g, sans dimension : de combien la taille peut croître par unité de longueur en
    /// s'éloignant d'un détail fin de la surface, strictement positive.
    pub grading: f64,
}

impl Default for TetOptions {
    fn default() -> Self {
        Self {
            size: None,
            grading: MeshOptions::default().grading,
        }
    }
}

/// Maille l'intérieur de la surface en tétraèdres dont les faces de bord suivent la surface.
/// `InvalidProperty` si la taille h ou la gradation n'est pas strictement positive,
/// `InvalidSurface` si la surface est inutilisable, `MeshingFailed` si la récupération échoue
/// (angles très aigus, arrondi), `TooManyNodes` au-delà de `MAX_NODES`.
pub fn tetrahedralize(
    surface: &Surface,
    options: &TetOptions,
) -> Result<Tetrahedralization<4>, MailleError> {
    let grading = positive("gradation", options.grading)?;
    match options.size {
        Some(h) => fill(&remesh(surface, h, grading)?, Some(h), grading),
        None => fill(surface, None, grading),
    }
}

/// Maille une pièce B-rep (lue d'un STEP) en tétraèdres à la taille h : ses faces sont maillées
/// une à une dans leur plan déroulé (bords nommés : les noms des faces), puis l'intérieur est
/// raffiné comme avec `tetrahedralize`. La taille h est obligatoire (`InvalidBrep` sans elle) ;
/// sur une face ou une arête courbe, elle est en plus bornée par un seizième de tour. Plan,
/// cylindre et cône sont pris en charge. Erreurs de `tetrahedralize`, et `InvalidBrep` si une face
/// ne se laisse pas mailler.
pub fn tetrahedralize_brep(
    brep: &Brep,
    options: &TetOptions,
) -> Result<Tetrahedralization<4>, MailleError> {
    let grading = positive("gradation", options.grading)?;
    let h = options.size.ok_or_else(|| MailleError::InvalidBrep {
        problem: "une pièce B-rep se maille à une taille h : la donner (TetOptions::size, --size)"
            .to_string(),
    })?;
    let h = positive("taille h", h)?;
    fill(&triangulate(brep, h, grading)?, Some(h), grading)
}

/// Essais d'insertion dans un autre ordre quand le Delaunay 3D bute sur un arrondi.
const REORDERINGS: u64 = 3;

/// Tétraèdres conformes à `surface`, intérieur raffiné à la taille h s'il y en a une. Les points
/// coplanaires d'une face plane rendent le Delaunay 3D sensible à l'arrondi, selon l'ordre où ils
/// entrent : sur `MeshingFailed`, la surface est reprise avec ses sommets dans un autre ordre.
fn fill(
    surface: &Surface,
    size: Option<f64>,
    grading: f64,
) -> Result<Tetrahedralization<4>, MailleError> {
    let mut result = fill_once(surface, size, grading);
    for seed in 1..=REORDERINGS {
        if result != Err(MailleError::MeshingFailed) {
            break;
        }
        result = fill_once(&reordered(surface, seed), size, grading);
    }
    result
}

/// La même surface, ses sommets dans un ordre pseudo-aléatoire fixé par `seed`.
fn reordered(surface: &Surface, seed: u64) -> Surface {
    let key = |i: usize| {
        (i as u64 ^ seed.wrapping_mul(0x9E37_79B9_7F4A_7C15))
            .wrapping_mul(0xBF58_476D_1CE4_E5B9)
            .rotate_left(31)
    };
    let mut order: Vec<usize> = (0..surface.vertices.len()).collect();
    order.sort_by_key(|&i| key(i));
    let mut number = vec![0; order.len()];
    for (k, &i) in order.iter().enumerate() {
        number[i] = k;
    }
    Surface {
        vertices: order.iter().map(|&i| surface.vertices[i]).collect(),
        triangles: surface
            .triangles
            .iter()
            .map(|t| t.map(|v| number[v]))
            .collect(),
        patches: surface.patches.clone(),
    }
}

fn fill_once(
    surface: &Surface,
    size: Option<f64>,
    grading: f64,
) -> Result<Tetrahedralization<4>, MailleError> {
    let patch_of = check(surface)?;
    let mut points = surface.vertices.clone();
    let (lo, hi) = bounds(&points);
    let mut del = Delaunay3::new(lo, hi);
    for (i, &p) in points.iter().enumerate() {
        if del.insert(p)? != i {
            return Err(MailleError::MeshingFailed);
        }
    }
    let pieces = conform(surface, &mut points, &mut del).map_err(|e| match e {
        MailleError::MeshingFailed => MailleError::InvalidSurface {
            problem: "la surface ne se laisse pas récupérer (triangles très effilés ou angles très aigus, par exemple un disque en éventail de centaines de triangles) : la remailler avec des triangles plus réguliers".to_string(),
        },
        e => e,
    })?;

    let caster = Caster::new(surface);
    del.classify(|t| caster.contains(centroid(&points, &t)));
    if let Some(h) = size {
        del.set_walls(pieces.iter().map(|p| p.v));
        refine3::refine(&mut del, &mut points, &SizeField3::new(surface, h, grading))?;
    }
    let tetrahedra = del.inside_tetrahedra();
    let volume: f64 = tetrahedra
        .iter()
        .map(|t| {
            let [a, b, c, d] = t.map(|i| points[i]);
            orient3d(a, b, c, d) / 6.0
        })
        .sum();
    if (volume - surface.volume()).abs() > VOLUME_TOLERANCE * surface.volume() {
        return Err(MailleError::MeshingFailed);
    }

    // chaque morceau de la surface est la face d'un seul tétraèdre du domaine
    let mut faces: HashMap<[usize; 3], (usize, usize)> = HashMap::new();
    for (e, t) in tetrahedra.iter().enumerate() {
        for (f, face) in TET_FACES.iter().enumerate() {
            faces.insert(sorted(face.map(|k| t[k])), (e, f));
        }
    }
    let mut boundaries: BTreeMap<String, BoundarySet> = BTreeMap::new();
    for piece in &pieces {
        let face = *faces
            .get(&sorted(piece.v))
            .ok_or(MailleError::MeshingFailed)?;
        let set = boundaries
            .entry(surface.patches[patch_of[piece.triangle]].name.clone())
            .or_default();
        set.nodes.extend(piece.v);
        set.faces.push(face);
    }
    for set in boundaries.values_mut() {
        set.nodes.sort_unstable();
        set.nodes.dedup();
        set.faces.sort_unstable();
    }
    Ok(Tetrahedralization {
        nodes: points.iter().map(|p| Node::new(p[0], p[1], p[2])).collect(),
        tetrahedra,
        boundaries,
    })
}

/// Contrôle minimal de la surface ; renvoie le patch de chaque triangle.
pub(super) fn check(surface: &Surface) -> Result<Vec<usize>, MailleError> {
    let invalid = |problem: &str| MailleError::InvalidSurface {
        problem: problem.to_string(),
    };
    let n = surface.vertices.len();
    if n < 4 || surface.triangles.len() < 4 {
        return Err(invalid("il faut au moins 4 sommets et 4 triangles"));
    }
    if surface.vertices.iter().flatten().any(|x| !x.is_finite()) {
        return Err(invalid("coordonnée non finie"));
    }
    if surface.triangles.iter().flatten().any(|&v| v >= n) {
        return Err(invalid("un triangle cite un sommet qui n'existe pas"));
    }
    let mut patch_of = vec![usize::MAX; surface.triangles.len()];
    for (p, patch) in surface.patches.iter().enumerate() {
        for &t in &patch.triangles {
            *patch_of
                .get_mut(t)
                .ok_or_else(|| invalid("un patch cite un triangle qui n'existe pas"))? = p;
        }
    }
    if patch_of.contains(&usize::MAX) {
        return Err(invalid("un triangle n'est dans aucun patch"));
    }
    Ok(patch_of)
}

fn bounds(points: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

fn sorted(mut v: [usize; 3]) -> [usize; 3] {
    v.sort_unstable();
    v
}

fn centroid(points: &[[f64; 3]], t: &[usize; 4]) -> [f64; 3] {
    std::array::from_fn(|k| t.iter().map(|&i| points[i][k]).sum::<f64>() / 4.0)
}

/// Coupe les arêtes puis les faces de la surface qui manquent au Delaunay, jusqu'à ce qu'il n'en
/// manque plus. Les arêtes d'abord : une face ne peut y figurer que si ses trois arêtes y sont.
fn conform(
    surface: &Surface,
    points: &mut Vec<[f64; 3]>,
    del: &mut Delaunay3,
) -> Result<Vec<Piece>, MailleError> {
    let mut pieces: Vec<Piece> = surface
        .triangles
        .iter()
        .enumerate()
        .map(|(triangle, &v)| Piece { v, triangle })
        .collect();
    for _ in 0..MAX_PASSES {
        if points.len() > MAX_GROWTH * surface.vertices.len() {
            return Err(MailleError::MeshingFailed);
        }
        let missing: BTreeSet<(usize, usize)> = pieces
            .iter()
            .flat_map(|p| (0..3).map(|k| edge(p.v[k], p.v[(k + 1) % 3])))
            .filter(|&(a, b)| !del.has_edge(a, b))
            .collect();
        if !missing.is_empty() {
            let mut mids = HashMap::new();
            for &(a, b) in &missing {
                let m = std::array::from_fn(|k| (points[a][k] + points[b][k]) / 2.0);
                mids.insert((a, b), add(del, points, m)?);
            }
            pieces = pieces.into_iter().flat_map(|p| split(p, &mids)).collect();
            continue;
        }
        let mut next = Vec::with_capacity(pieces.len());
        let mut changed = false;
        for p in pieces {
            let [a, b, c] = p.v;
            if del.has_face(a, b, c) {
                next.push(p);
                continue;
            }
            changed = true;
            let g = std::array::from_fn(|k| (points[a][k] + points[b][k] + points[c][k]) / 3.0);
            let m = add(del, points, g)?;
            for v in [[a, b, m], [b, c, m], [c, a, m]] {
                next.push(Piece { v, ..p });
            }
        }
        pieces = next;
        if !changed {
            return Ok(pieces);
        }
    }
    Err(MailleError::MeshingFailed)
}

/// Insère `p` dans le Delaunay et dans la liste des points.
fn add(del: &mut Delaunay3, points: &mut Vec<[f64; 3]>, p: [f64; 3]) -> Result<usize, MailleError> {
    if points.len() >= MAX_NODES {
        return Err(MailleError::TooManyNodes { max: MAX_NODES });
    }
    let i = del.insert(p)?;
    points.push(p);
    Ok(i)
}

fn edge(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Coupe `piece` à chaque arête de `mids` (arête → point milieu), récursivement.
fn split(piece: Piece, mids: &HashMap<(usize, usize), usize>) -> Vec<Piece> {
    let [a, b, c] = piece.v;
    for v in [[a, b, c], [b, c, a], [c, a, b]] {
        if let Some(&m) = mids.get(&edge(v[0], v[1])) {
            return [[v[0], m, v[2]], [m, v[1], v[2]]]
                .into_iter()
                .flat_map(|v| split(Piece { v, ..piece }, mids))
                .collect();
        }
    }
    vec![piece]
}

/// Appartenance d'un point à la surface par lancer de rayon le long de +x : le nombre de triangles
/// traversés est impair au-dedans. Les triangles sont rangés dans une grille en (y, z). Un rayon
/// qui touche une arête ou un sommet est décalé symboliquement de (ε, ε²) en (y, z) : un seul des
/// triangles voisins compte, sans cas particulier.
struct Caster<'a> {
    surface: &'a Surface,
    lo: [f64; 2],
    cell: [f64; 2],
    n: usize,
    grid: Vec<Vec<usize>>,
}

impl<'a> Caster<'a> {
    fn new(surface: &'a Surface) -> Self {
        let (lo, hi) = bounds(&surface.vertices);
        let n = ((surface.triangles.len() as f64).sqrt().ceil() as usize).clamp(1, 128);
        let cell = [1, 2].map(|k| ((hi[k] - lo[k]) / n as f64).max(f64::MIN_POSITIVE));
        let mut caster = Self {
            surface,
            lo: [lo[1], lo[2]],
            cell,
            n,
            grid: vec![Vec::new(); n * n],
        };
        for (i, t) in surface.triangles.iter().enumerate() {
            let ys = t.map(|v| surface.vertices[v][1]);
            let zs = t.map(|v| surface.vertices[v][2]);
            let (y0, y1) = (caster.index(0, min3(ys)), caster.index(0, max3(ys)));
            let (z0, z1) = (caster.index(1, min3(zs)), caster.index(1, max3(zs)));
            for y in y0..=y1 {
                for z in z0..=z1 {
                    caster.grid[y * n + z].push(i);
                }
            }
        }
        caster
    }

    fn index(&self, axis: usize, x: f64) -> usize {
        (((x - self.lo[axis]) / self.cell[axis]).floor().max(0.0) as usize).min(self.n - 1)
    }

    fn contains(&self, p: [f64; 3]) -> bool {
        let cell = &self.grid[self.index(0, p[1]) * self.n + self.index(1, p[2])];
        let hits = cell.iter().filter(|&&i| self.crosses(i, p)).count();
        hits % 2 == 1
    }

    /// Vrai si le rayon de p vers +x traverse le triangle `i`.
    fn crosses(&self, i: usize, p: [f64; 3]) -> bool {
        let [a, b, c] = self.surface.triangles[i].map(|v| self.surface.vertices[v]);
        let signs =
            [(a, b), (b, c), (c, a)].map(|(s, e)| side([s[1], s[2]], [e[1], e[2]], [p[1], p[2]]));
        if !(signs.iter().all(|&s| s) || signs.iter().all(|&s| !s)) {
            return false;
        }
        // abscisse du plan du triangle sur le rayon
        let n = [
            (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
            (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
            (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]),
        ];
        if n[0] == 0.0 {
            return false;
        }
        let x = a[0] - (n[1] * (p[1] - a[1]) + n[2] * (p[2] - a[2])) / n[0];
        x > p[0]
    }
}

/// Signe de l'aire signée (a, b, p) en (y, z), p décalé de (ε, ε²) si elle est nulle.
fn side(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    let d = [b[0] - a[0], b[1] - a[1]];
    let (l, r) = (d[0] * (p[1] - a[1]), d[1] * (p[0] - a[0]));
    let e = l - r;
    if e.abs() > 1e-12 * (l.abs() + r.abs()) {
        e > 0.0
    } else if d[1] != 0.0 {
        d[1] < 0.0
    } else {
        d[0] > 0.0
    }
}

fn min3(v: [f64; 3]) -> f64 {
    v[0].min(v[1]).min(v[2])
}

fn max3(v: [f64; 3]) -> f64 {
    v[0].max(v[1]).max(v[2])
}
