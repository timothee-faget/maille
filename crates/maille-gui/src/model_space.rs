//! L'espace Modèle sans fenêtre (partie statique) : le brouillon du modèle que l'utilisateur remplit
//! (Material, BC, Load, Pressure), ses diagnostics en langage clair avant le calcul, et son écriture
//! en TOML par `model::write`. `app` ne fait que l'afficher.

use maille_io::model::{
    self, DofDef, MeshingDef, ModelFile, Physics, PlaneStateDef, PressureDef, SectionDef, ValueDef,
};
use maille_io::msh::{CellKind, MshFile, DOMAIN_SET};
use maille_io::IoError;

use crate::state::Tool;

/// Système d'unités cohérent choisi au départ : il s'écrit dans la clé `units` et s'affiche à côté de
/// chaque champ. Le cœur ne convertit rien.
#[derive(Debug, Clone, PartialEq)]
pub enum UnitSystem {
    /// N, m, Pa.
    Si,
    /// N, mm, MPa.
    MmNMpa,
    /// Texte libre : aucune unité affichée.
    Custom(String),
}

impl UnitSystem {
    /// Le texte de la clé `units`.
    pub fn text(&self) -> String {
        match self {
            UnitSystem::Si => "N, m, Pa".to_string(),
            UnitSystem::MmNMpa => "N, mm, MPa".to_string(),
            UnitSystem::Custom(text) => text.clone(),
        }
    }

    /// Retrouve le système d'un texte de `units` ; un texte inconnu reste libre.
    pub fn from_text(text: &str) -> Self {
        match text.trim() {
            "N, m, Pa" => UnitSystem::Si,
            "N, mm, MPa" => UnitSystem::MmNMpa,
            other => UnitSystem::Custom(other.to_string()),
        }
    }

    pub fn length(&self) -> &'static str {
        match self {
            UnitSystem::Si => "m",
            UnitSystem::MmNMpa => "mm",
            UnitSystem::Custom(_) => "",
        }
    }

    pub fn force(&self) -> &'static str {
        match self {
            UnitSystem::Si | UnitSystem::MmNMpa => "N",
            UnitSystem::Custom(_) => "",
        }
    }

    /// Contrainte, module d'Young, pression.
    pub fn stress(&self) -> &'static str {
        match self {
            UnitSystem::Si => "Pa",
            UnitSystem::MmNMpa => "MPa",
            UnitSystem::Custom(_) => "",
        }
    }
}

/// Le Material élastique du modèle (un seul pour l'instant : il couvre tout le domaine).
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialForm {
    pub name: String,
    pub young: f64,
    pub poisson: f64,
}

/// Ce que l'utilisateur a rempli dans l'espace Modèle.
#[derive(Debug, Clone, PartialEq)]
pub struct Draft {
    pub units: UnitSystem,
    pub material: Option<MaterialForm>,
    /// Triangles seulement.
    pub plane_state: PlaneStateDef,
    pub thickness: f64,
    pub bcs: Vec<ValueDef>,
    /// Efforts par Node (`Load`).
    pub loads: Vec<ValueDef>,
    pub pressures: Vec<PressureDef>,
}

impl Default for Draft {
    fn default() -> Self {
        Self {
            units: UnitSystem::MmNMpa,
            material: None,
            plane_state: PlaneStateDef::Stress,
            thickness: 1.0,
            bcs: Vec::new(),
            loads: Vec::new(),
            pressures: Vec::new(),
        }
    }
}

impl Draft {
    /// Un acier, dans le système d'unités choisi.
    pub fn steel(&self) -> MaterialForm {
        MaterialForm {
            name: "acier".to_string(),
            young: if self.units == UnitSystem::Si {
                210e9
            } else {
                210e3
            },
            poisson: 0.3,
        }
    }

    /// Le fichier de modèle correspondant. `meshing` : la recette du maillage si on la connaît.
    pub fn to_file(
        &self,
        mesh_name: &str,
        mesh: &MshFile,
        meshing: Option<MeshingDef>,
    ) -> ModelFile {
        let planar = is_planar(mesh);
        let mut materials = std::collections::BTreeMap::new();
        let mut sections = Vec::new();
        if let Some(m) = &self.material {
            materials.insert(
                m.name.clone(),
                model::MaterialDef {
                    young: Some(m.young),
                    poisson: Some(m.poisson),
                    conductivity: None,
                    density: None,
                    specific_heat: None,
                    expansion: None,
                    yield_stress: None,
                    hardening: None,
                    hyperelasticity: None,
                },
            );
            sections.push(SectionDef {
                set: DOMAIN_SET.to_string(),
                material: m.name.clone(),
                thickness: (planar && (self.thickness - 1.0).abs() > f64::EPSILON)
                    .then_some(self.thickness),
                plane_state: planar.then_some(self.plane_state),
            });
        }
        ModelFile {
            units: self.units.text(),
            mesh: mesh_name.to_string(),
            analysis: None,
            kinematics: None,
            materials,
            sections,
            bcs: self.bcs.clone(),
            loads: self.loads.clone(),
            heat_sources: Vec::new(),
            convections: Vec::new(),
            pressures: self.pressures.clone(),
            contacts: Vec::new(),
            temperatures: None,
            step: None,
            meshing,
            thermal_step: None,
        }
    }

    /// Le TOML du modèle.
    pub fn to_toml(
        &self,
        mesh_name: &str,
        mesh: &MshFile,
        meshing: Option<MeshingDef>,
    ) -> Result<String, IoError> {
        model::write(&self.to_file(mesh_name, mesh, meshing))
    }

    /// Le calcul que le contenu demande, s'il est assez rempli pour qu'on le déduise.
    pub fn physics(&self, mesh: &MshFile) -> Option<Physics> {
        model::physics(&self.to_file("", mesh, None)).ok()
    }

    /// Il y a de quoi déformer la pièce : un effort, une pression ou un déplacement imposé non nul.
    pub fn has_loading(&self) -> bool {
        !self.loads.is_empty()
            || !self.pressures.is_empty()
            || self.bcs.iter().any(|b| b.value != 0.0)
    }

    /// Noms des Sets que le modèle utilise.
    fn used_sets(&self) -> Vec<&str> {
        self.bcs
            .iter()
            .chain(&self.loads)
            .map(|v| v.set.as_str())
            .chain(self.pressures.iter().map(|p| p.set.as_str()))
            .collect()
    }
}

/// Phrase qui décrit le calcul déduit.
pub fn physics_label(physics: Physics) -> &'static str {
    match physics {
        Physics::Static => "statique",
        Physics::Thermal => "thermique",
        Physics::Coupled => "thermo-mécanique",
    }
}

fn is_planar(mesh: &MshFile) -> bool {
    matches!(mesh.kind, CellKind::Tri3 | CellKind::Tri6)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Bloque le calcul.
    Error,
    /// Le calcul se lance, mais le résultat risque d'être décevant.
    Warning,
}

/// Un problème dit en clair : ce qui ne va pas et quoi faire. Un clic mène à l'outil qui le corrige.
#[derive(Debug, Clone, PartialEq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub message: String,
    pub tool: Tool,
}

fn error(tool: Tool, message: String) -> Diagnostic {
    Diagnostic {
        severity: Severity::Error,
        message,
        tool,
    }
}

fn warning(tool: Tool, message: String) -> Diagnostic {
    Diagnostic {
        severity: Severity::Warning,
        message,
        tool,
    }
}

/// Ce qui ne va pas dans le brouillon, avant le calcul. `user_sets` : les Sets créés par l'utilisateur,
/// dont on signale ceux que rien n'utilise.
pub fn diagnose(draft: &Draft, mesh: &MshFile, user_sets: &[&str]) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    match &draft.material {
        None => out.push(error(
            Tool::Material,
            "Aucun Material : définis-en un (module d'Young et coefficient de Poisson).".into(),
        )),
        Some(m) => {
            if m.young <= 0.0 {
                out.push(error(
                    Tool::Material,
                    format!("Le module d'Young de « {} » doit être positif.", m.name),
                ));
            }
            if m.poisson <= -1.0 || m.poisson >= 0.5 {
                out.push(error(
                    Tool::Material,
                    format!(
                        "Le coefficient de Poisson de « {} » doit être entre -1 et 0,5 (exclus).",
                        m.name
                    ),
                ));
            } else if m.poisson > 0.45 && matches!(mesh.kind, CellKind::Tri3 | CellKind::Tet4) {
                out.push(warning(
                    Tool::MeshSize,
                    format!(
                        "ν = {} avec des éléments linéaires : le résultat risque d'être trop raide, coche quadratique.",
                        m.poisson
                    ),
                ));
            }
        }
    }
    if is_planar(mesh) && draft.thickness <= 0.0 {
        out.push(error(
            Tool::Material,
            "L'épaisseur doit être positive.".into(),
        ));
    }

    let mut check_set = |tool: Tool, what: &str, set: &str, needs_faces: bool| {
        match mesh
        .sets
        .get(set)
    {
        None => out.push(error(
            tool,
            format!("{what} renvoie au Set « {set} », qui n'existe pas dans le maillage (remaillage ?) : choisis un autre Set."),
        )),
        Some(s) if needs_faces && s.faces.is_empty() => out.push(error(
            tool,
            format!("{what} demande un bord, et le Set « {set} » n'a pas de faces : choisis un Set de bord."),
        )),
        Some(_) => {}
    }
    };
    for bc in &draft.bcs {
        check_set(
            Tool::Support,
            &format!("La BC « {} »", dof_name(bc.dof)),
            &bc.set,
            false,
        );
    }
    for load in &draft.loads {
        check_set(
            Tool::Load,
            &format!("Le Load « {} »", dof_name(load.dof)),
            &load.set,
            false,
        );
    }
    for p in &draft.pressures {
        check_set(Tool::Load, "La Pressure", &p.set, true);
    }

    for name in user_sets {
        if !draft.used_sets().contains(name) {
            out.push(warning(
                Tool::Support,
                format!("Le Set « {name} » n'est utilisé par aucune BC ni aucun chargement."),
            ));
        }
    }

    if draft.bcs.is_empty() {
        out.push(error(
            Tool::Support,
            "Aucune BC : la pièce flotte. Ajoute un encastrement ou un déplacement imposé.".into(),
        ));
    } else if out.iter().all(|d| d.tool != Tool::Support) {
        let free = free_rigid_modes(mesh, &draft.bcs);
        if !free.is_empty() {
            out.push(error(
                Tool::Support,
                format!(
                    "La pièce peut encore bouger d'un bloc ({}) : ajoute une BC qui bloque ce mouvement.",
                    free.join(", ")
                ),
            ));
        }
    }
    if !draft.bcs.is_empty() && !draft.has_loading() {
        out.push(warning(
            Tool::Load,
            "Aucun chargement : sans effort ni déplacement imposé, la pièce ne se déforme pas."
                .into(),
        ));
    }

    // les règles du cœur, quand rien d'autre n'a été trouvé
    if out.iter().all(|d| d.severity != Severity::Error) {
        if let Err(e) = model::build(&draft.to_file("", mesh, None), mesh, None) {
            out.push(error(Tool::Run, e.to_string()));
        }
    }
    out
}

fn dof_name(dof: DofDef) -> &'static str {
    match dof {
        DofDef::Ux => "ux",
        DofDef::Uy => "uy",
        DofDef::Uz => "uz",
        DofDef::T => "t",
    }
}

/// Les mouvements de corps rigide que les BC ne bloquent pas, en clair ; vide si la pièce est tenue.
///
/// On projette les modes rigides (3 en 2D, 6 en 3D) sur les DOF bloqués : un mouvement rigide reste
/// libre si une combinaison des modes s'annule sur tous les DOF bloqués, soit un rang de la matrice
/// des modes inférieur à leur nombre.
pub fn free_rigid_modes(mesh: &MshFile, bcs: &[ValueDef]) -> Vec<String> {
    let planar = is_planar(mesh);
    let names: &[&str] = if planar {
        &[
            "translation en X",
            "translation en Y",
            "rotation dans le plan",
        ]
    } else {
        &[
            "translation en X",
            "translation en Y",
            "translation en Z",
            "rotation autour de X",
            "rotation autour de Y",
            "rotation autour de Z",
        ]
    };
    let n = names.len();
    let nodes = &mesh.nodes;
    if nodes.is_empty() {
        return Vec::new();
    }
    let count = nodes.len() as f64;
    let c = [
        nodes.iter().map(|p| p.x).sum::<f64>() / count,
        nodes.iter().map(|p| p.y).sum::<f64>() / count,
        nodes.iter().map(|p| p.z).sum::<f64>() / count,
    ];
    let extent = nodes
        .iter()
        .flat_map(|p| [(p.x - c[0]).abs(), (p.y - c[1]).abs(), (p.z - c[2]).abs()])
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);

    // valeur de chaque mode sur le DOF `axis` du Node `i`
    let mode = |m: usize, i: usize, axis: usize| -> f64 {
        let p = &nodes[i];
        let r = [
            (p.x - c[0]) / extent,
            (p.y - c[1]) / extent,
            (p.z - c[2]) / extent,
        ];
        if planar {
            match m {
                0 => f64::from(axis == 0),
                1 => f64::from(axis == 1),
                _ => [-r[1], r[0]][axis],
            }
        } else {
            match m {
                0..=2 => f64::from(axis == m),
                // ω × r pour une rotation unitaire autour de l'axe m - 3
                _ => {
                    let mut w = [0.0; 3];
                    w[m - 3] = 1.0;
                    [
                        w[1] * r[2] - w[2] * r[1],
                        w[2] * r[0] - w[0] * r[2],
                        w[0] * r[1] - w[1] * r[0],
                    ][axis]
                }
            }
        }
    };

    // Gram des modes sur les DOF bloqués
    let mut gram = vec![vec![0.0; n]; n];
    let dim = if planar { 2 } else { 3 };
    for bc in bcs {
        let axis = match bc.dof {
            DofDef::Ux => 0,
            DofDef::Uy => 1,
            DofDef::Uz => 2,
            DofDef::T => continue,
        };
        if axis >= dim {
            continue;
        }
        let Some(set) = mesh.sets.get(&bc.set) else {
            continue;
        };
        for &i in &set.nodes {
            let v: Vec<f64> = (0..n).map(|m| mode(m, i, axis)).collect();
            for a in 0..n {
                for b in 0..n {
                    gram[a][b] += v[a] * v[b];
                }
            }
        }
    }

    let scale = (0..n).map(|a| gram[a][a]).fold(0.0, f64::max);
    let tol = 1e-9 * scale.max(1e-300);
    let mut free: Vec<String> = (0..n)
        .filter(|&a| gram[a][a] <= tol)
        .map(|a| names[a].to_string())
        .collect();
    if free.is_empty() && rank(gram, tol) < n {
        free.push(
            "un mouvement combiné, par exemple une rotation autour d'un point d'appui".into(),
        );
    }
    free
}

/// Rang d'une matrice symétrique positive, par élimination avec pivot sur la diagonale.
fn rank(mut a: Vec<Vec<f64>>, tol: f64) -> usize {
    let n = a.len();
    let mut rank = 0;
    for _ in 0..n {
        let Some((p, &pivot)) = (0..n)
            .map(|i| (i, &a[i][i]))
            .max_by(|x, y| x.1.total_cmp(y.1))
        else {
            break;
        };
        if pivot <= tol {
            break;
        }
        rank += 1;
        let row = a[p].clone();
        for (i, line) in a.iter_mut().enumerate() {
            if i == p {
                continue;
            }
            let f = line[p] / pivot;
            for (value, r) in line.iter_mut().zip(&row) {
                *value -= f * r;
            }
        }
        // le pivot ne sert plus
        a[p].iter_mut().for_each(|v| *v = 0.0);
        a.iter_mut().for_each(|line| line[p] = 0.0);
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::*;
    use maille_core::mesh::Node;
    use maille_io::msh::Set;

    /// Plaque 2 × 1 en deux triangles : Sets `gauche` (nodes 0, 3), `droite` (1, 2), `bas` (0, 1).
    fn plate() -> MshFile {
        let node = |x, y| Node { x, y, z: 0.0 };
        let mut sets = std::collections::BTreeMap::new();
        let nodes_set = |nodes: Vec<usize>, faces: Vec<(usize, usize)>| Set {
            nodes,
            faces,
            elements: Vec::new(),
        };
        sets.insert("gauche".to_string(), nodes_set(vec![0, 3], vec![(1, 2)]));
        sets.insert("droite".to_string(), nodes_set(vec![1, 2], vec![(0, 1)]));
        sets.insert("bas".to_string(), nodes_set(vec![0, 1], vec![(0, 0)]));
        sets.insert(
            DOMAIN_SET.to_string(),
            Set {
                nodes: vec![0, 1, 2, 3],
                faces: Vec::new(),
                elements: vec![0, 1],
            },
        );
        MshFile {
            nodes: vec![
                node(0.0, 0.0),
                node(2.0, 0.0),
                node(2.0, 1.0),
                node(0.0, 1.0),
            ],
            kind: CellKind::Tri3,
            cells: vec![vec![0, 1, 2], vec![0, 2, 3]],
            sets,
        }
    }

    fn bc(set: &str, dof: DofDef) -> ValueDef {
        ValueDef {
            set: set.to_string(),
            dof,
            value: 0.0,
        }
    }

    fn filled() -> Draft {
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.bcs = vec![bc("gauche", DofDef::Ux), bc("gauche", DofDef::Uy)];
        d.pressures = vec![PressureDef {
            set: "droite".into(),
            value: -10.0,
        }];
        d
    }

    fn messages(d: &Draft, mesh: &MshFile, user: &[&str]) -> Vec<String> {
        diagnose(d, mesh, user)
            .into_iter()
            .map(|d| d.message)
            .collect()
    }

    #[test]
    fn a_clamped_plate_with_a_pressure_has_no_diagnostic() {
        assert_eq!(messages(&filled(), &plate(), &[]), Vec::<String>::new());
    }

    #[test]
    fn a_missing_material_is_reported() {
        let mut d = filled();
        d.material = None;
        let m = messages(&d, &plate(), &[]);
        assert!(m.iter().any(|s| s.contains("Aucun Material")), "{m:?}");
    }

    #[test]
    fn a_part_that_can_still_move_as_a_block_is_reported() {
        // blocage de ux seulement sur le bord gauche : la plaque glisse en Y et tourne
        let mut d = filled();
        d.bcs = vec![bc("gauche", DofDef::Ux)];
        let m = messages(&d, &plate(), &[]);
        assert!(
            m.iter()
                .any(|s| s.contains("bouger d'un bloc") && s.contains("translation en Y")),
            "{m:?}"
        );
        // deux DOF sur un seul Node d'un bord de deux Nodes : la rotation reste libre
        d.bcs = vec![bc("gauche", DofDef::Ux), bc("bas", DofDef::Uy)];
        assert!(
            messages(&d, &plate(), &[]).is_empty(),
            "ux gauche + uy bas tient la plaque"
        );
    }

    #[test]
    fn a_single_support_point_leaves_a_combined_rotation() {
        // un seul Node encastré : les deux translations sont bloquées, la rotation reste libre
        let mut mesh = plate();
        mesh.sets.insert(
            "coin".to_string(),
            Set {
                nodes: vec![0],
                faces: Vec::new(),
                elements: Vec::new(),
            },
        );
        let free = free_rigid_modes(&mesh, &[bc("coin", DofDef::Ux), bc("coin", DofDef::Uy)]);
        assert_eq!(free.len(), 1);
        assert!(free[0].contains("mouvement combiné"), "{free:?}");
    }

    #[test]
    fn supports_that_cross_block_everything() {
        // ux sur le bord x = 0 et uy en un point : tout est bloqué ; uy seul sur deux Nodes d'un bord
        // horizontal laisse la translation en X
        let mut mesh = plate();
        mesh.sets.insert(
            "coin".to_string(),
            Set {
                nodes: vec![0],
                faces: Vec::new(),
                elements: Vec::new(),
            },
        );
        let free = free_rigid_modes(&mesh, &[bc("gauche", DofDef::Ux), bc("coin", DofDef::Uy)]);
        assert!(
            free.is_empty(),
            "ux sur x = 0 et uy en un point bloquent tout : {free:?}"
        );
        let free = free_rigid_modes(&mesh, &[bc("bas", DofDef::Uy), bc("coin", DofDef::Uy)]);
        assert!(
            free.iter().any(|s| s.contains("translation en X")),
            "{free:?}"
        );
    }

    #[test]
    fn a_three_dimensional_part_needs_six_blocked_motions() {
        let node = |x, y, z| Node { x, y, z };
        let mut sets = std::collections::BTreeMap::new();
        sets.insert(
            "pied".to_string(),
            Set {
                nodes: vec![0, 1, 2],
                faces: vec![(0, 3)],
                elements: Vec::new(),
            },
        );
        let mesh = MshFile {
            nodes: vec![
                node(0.0, 0.0, 0.0),
                node(1.0, 0.0, 0.0),
                node(0.0, 1.0, 0.0),
                node(0.0, 0.0, 1.0),
            ],
            kind: CellKind::Tet4,
            cells: vec![vec![0, 1, 2, 3]],
            sets,
        };
        let all = [
            bc("pied", DofDef::Ux),
            bc("pied", DofDef::Uy),
            bc("pied", DofDef::Uz),
        ];
        assert!(free_rigid_modes(&mesh, &all).is_empty());
        assert!(
            !free_rigid_modes(&mesh, &all[..2]).is_empty(),
            "sans uz, la pièce bouge encore en Z"
        );
    }

    #[test]
    fn a_set_that_vanished_is_reported_by_name() {
        let mut d = filled();
        d.bcs.push(bc("disparu", DofDef::Ux));
        let m = messages(&d, &plate(), &[]);
        assert!(
            m.iter()
                .any(|s| s.contains("disparu") && s.contains("n'existe pas")),
            "{m:?}"
        );
    }

    #[test]
    fn a_pressure_needs_a_boundary_set() {
        let mut d = filled();
        d.pressures = vec![PressureDef {
            set: DOMAIN_SET.into(),
            value: 1.0,
        }];
        let m = messages(&d, &plate(), &[]);
        assert!(m.iter().any(|s| s.contains("n'a pas de faces")), "{m:?}");
    }

    #[test]
    fn an_unused_user_set_is_a_warning() {
        let diags = diagnose(&filled(), &plate(), &["bas", "gauche"]);
        let unused: Vec<_> = diags
            .iter()
            .filter(|d| d.message.contains("« bas »"))
            .collect();
        assert_eq!(unused.len(), 1);
        assert_eq!(unused[0].severity, Severity::Warning);
        assert!(!diags.iter().any(|d| d.message.contains("« gauche »")));
    }

    #[test]
    fn a_nearly_incompressible_material_with_linear_elements_is_a_warning() {
        let mut d = filled();
        d.material.as_mut().unwrap().poisson = 0.49;
        let diags = diagnose(&d, &plate(), &[]);
        assert!(diags
            .iter()
            .any(|d| d.severity == Severity::Warning && d.message.contains("quadratique")));
        assert!(diags.iter().all(|d| d.severity != Severity::Error));
    }

    #[test]
    fn an_impossible_material_is_an_error() {
        let mut d = filled();
        d.material.as_mut().unwrap().young = -1.0;
        d.material.as_mut().unwrap().poisson = 0.5;
        let diags = diagnose(&d, &plate(), &[]);
        assert_eq!(
            diags
                .iter()
                .filter(|d| d.severity == Severity::Error)
                .count(),
            2
        );
    }

    #[test]
    fn no_load_at_all_is_a_warning_but_an_imposed_displacement_is_a_load() {
        let mut d = filled();
        d.pressures.clear();
        assert!(messages(&d, &plate(), &[])
            .iter()
            .any(|s| s.contains("Aucun chargement")));
        d.bcs.push(ValueDef {
            set: "droite".into(),
            dof: DofDef::Ux,
            value: 0.1,
        });
        assert!(messages(&d, &plate(), &[]).is_empty());
    }

    #[test]
    fn the_calculation_is_deduced_from_the_content() {
        let mesh = plate();
        assert_eq!(Draft::default().physics(&mesh), None);
        assert_eq!(filled().physics(&mesh), Some(Physics::Static));
    }

    #[test]
    fn the_draft_round_trips_through_the_toml() {
        let mesh = plate();
        let d = filled();
        let text = d.to_toml("plaque.msh", &mesh, None).unwrap();
        let file = model::parse(&text).unwrap();
        assert_eq!(file, d.to_file("plaque.msh", &mesh, None));
        assert_eq!(file.units, "N, mm, MPa");
        assert_eq!(file.sections[0].plane_state, Some(PlaneStateDef::Stress));
        // et le modèle se construit
        model::build(&file, &mesh, None).unwrap();
    }

    /// Le parcours de l'interface sur `plaque.dxf` : maillage, Material, BC, calcul, comparé à la
    /// solution exacte de `examples/traction.toml` (ux = 0,001 x, uy = −ν 0,001 y).
    #[test]
    fn the_traction_plate_built_in_the_model_space_gives_the_exact_solution() {
        use crate::mesh_space::{import, mesh, Imported};
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/plaque.dxf");
        let Imported::Cad(cad) = import(&path).unwrap() else {
            panic!("plaque.dxf est un dessin");
        };
        let built = mesh(&cad, 0.5, false, &Default::default()).unwrap();
        let mut d = Draft::default();
        d.material = Some(d.steel());
        d.thickness = 2.0;
        d.bcs = vec![
            bc("gauche", DofDef::Ux),
            ValueDef {
                set: "droite".into(),
                dof: DofDef::Ux,
                value: 0.004,
            },
            bc("bas", DofDef::Uy),
        ];
        assert_eq!(messages(&d, &built.file, &[]), Vec::<String>::new());
        let text = d.to_toml("plaque.msh", &built.file, None).unwrap();
        let loaded = model::build(&model::parse(&text).unwrap(), &built.file, None).unwrap();
        let field = loaded.model.solve().unwrap();
        for (i, p) in built.file.nodes.iter().enumerate() {
            let (ux, uy) = (field.displacements[2 * i], field.displacements[2 * i + 1]);
            assert!((ux - 0.001 * p.x).abs() < 1e-9, "ux du Node {i}");
            assert!((uy + 0.3 * 0.001 * p.y).abs() < 1e-9, "uy du Node {i}");
        }
    }

    #[test]
    fn units_follow_the_chosen_system() {
        assert_eq!(UnitSystem::from_text("N, m, Pa"), UnitSystem::Si);
        assert_eq!(UnitSystem::from_text("N, mm, MPa").stress(), "MPa");
        let custom = UnitSystem::from_text("lbf, in, psi");
        assert_eq!(custom.stress(), "");
        assert_eq!(custom.text(), "lbf, in, psi");
        let d = Draft {
            units: UnitSystem::Si,
            ..Draft::default()
        };
        assert_eq!(d.steel().young, 210e9);
    }
}
