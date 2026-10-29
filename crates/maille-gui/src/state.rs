//! État de l'interface, sans fenêtre ni GPU : l'espace actif, l'outil de chacun, le contenu des
//! barres d'outils et des cartes du parcours. Les widgets (`app`) ne font que l'afficher.

/// Les trois espaces de travail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Space {
    Mesh,
    Model,
    Results,
}

impl Space {
    pub const ALL: [Space; 3] = [Space::Mesh, Space::Model, Space::Results];

    pub fn label(self) -> &'static str {
        match self {
            Space::Mesh => "Maillage",
            Space::Model => "Modèle",
            Space::Results => "Résultats",
        }
    }

    fn index(self) -> usize {
        self as usize
    }
}

/// Un outil de la barre d'outils : son inspecteur flotte à sa gauche.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    /// Navigation : tourner, déplacer, zoomer (commun aux trois espaces).
    Navigate,
    // Maillage
    Import,
    MeshSize,
    FaceSets,
    VertexSets,
    ManageSets,
    Quality,
    // Modèle
    Material,
    Support,
    Load,
    Thermal,
    Contact,
    StepAmplitude,
    Run,
    // Résultats
    Open,
    Field,
    Deformation,
    Cut,
    Probe,
    Increments,
    Curve,
    Capture,
}

impl Tool {
    pub fn label(self) -> &'static str {
        match self {
            Tool::Navigate => "Flèche",
            Tool::Import => "Importer",
            Tool::MeshSize => "Taille h",
            Tool::FaceSets => "Faces → Set",
            Tool::VertexSets => "Sommets → Set",
            Tool::ManageSets => "Sets",
            Tool::Quality => "Qualité",
            Tool::Material => "Material",
            Tool::Support => "Encastrement",
            Tool::Load => "Load, Pressure",
            Tool::Thermal => "Thermique",
            Tool::Contact => "Contact",
            Tool::StepAmplitude => "Step, Amplitude",
            Tool::Run => "Calculer",
            Tool::Open => "Ouvrir",
            Tool::Field => "Champ",
            Tool::Deformation => "Déformée",
            Tool::Cut => "Coupe",
            Tool::Probe => "Sonde",
            Tool::Increments => "Increments",
            Tool::Curve => "Courbe XY",
            Tool::Capture => "Capture PNG",
        }
    }

    /// Faux tant que l'outil n'a pas son espace : l'inspecteur dit alors quand il arrive.
    pub fn available(self) -> bool {
        matches!(
            self,
            Tool::Navigate
                | Tool::Import
                | Tool::MeshSize
                | Tool::FaceSets
                | Tool::ManageSets
                | Tool::Quality
                | Tool::Material
                | Tool::Support
                | Tool::Load
                | Tool::Thermal
                | Tool::Contact
                | Tool::StepAmplitude
                | Tool::Run
                | Tool::Open
                | Tool::Field
                | Tool::Deformation
                | Tool::Cut
        )
    }
}

/// Action immédiate de la barre d'outils (pas d'inspecteur).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// Cadrer la pièce (touche F).
    Fit,
    /// Regarder le long de l'axe 0 = x, 1 = y, 2 = z.
    View(usize),
}

impl Action {
    pub fn label(self) -> &'static str {
        match self {
            Action::Fit => "Cadrer",
            Action::View(0) => "Vue X",
            Action::View(1) => "Vue Y",
            Action::View(_) => "Vue Z",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Item {
    Tool(Tool),
    Action(Action),
}

/// Une famille d'outils de la barre : un séparateur la distingue des autres.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group {
    pub items: Vec<Item>,
}

/// Contenu de la barre d'outils d'un espace : les outils communs, puis ceux de l'espace.
pub fn toolbar(space: Space) -> Vec<Group> {
    use Item::{Action as A, Tool as T};
    let common = Group {
        items: vec![
            T(Tool::Navigate),
            A(Action::Fit),
            A(Action::View(0)),
            A(Action::View(1)),
            A(Action::View(2)),
        ],
    };
    let group = |tools: &[Tool]| Group {
        items: tools.iter().map(|&t| T(t)).collect(),
    };
    let own = match space {
        Space::Mesh => vec![
            group(&[Tool::Import, Tool::MeshSize]),
            group(&[Tool::FaceSets, Tool::VertexSets, Tool::ManageSets]),
            group(&[Tool::Quality]),
        ],
        Space::Model => vec![
            group(&[Tool::Material]),
            group(&[Tool::Support, Tool::Load, Tool::Thermal, Tool::Contact]),
            group(&[Tool::StepAmplitude, Tool::Run]),
        ],
        Space::Results => vec![
            group(&[Tool::Open]),
            group(&[Tool::Field, Tool::Deformation, Tool::Cut]),
            group(&[Tool::Probe, Tool::Increments, Tool::Curve, Tool::Capture]),
        ],
    };
    std::iter::once(common).chain(own).collect()
}

/// État d'une étape du parcours, et d'un espace (✅ / ⚠️ / ❌ du sélecteur).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Done,
    Warning,
    Missing,
}

/// Une ligne de la carte du parcours : un clic active l'outil correspondant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Step {
    pub label: &'static str,
    pub status: Status,
    pub tool: Tool,
    /// Étape facultative : elle ne compte pas dans l'état de l'espace.
    pub optional: bool,
}

/// Où en est le travail de l'utilisateur. L'espace Modèle y ajoutera ses champs avec ses outils.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Progress {
    /// Une géométrie ou un maillage est importé.
    pub geometry: bool,
    /// Il y a un maillage (importé tel quel, ou fait par le mailleur).
    pub meshed: bool,
    /// Nombre de Sets de bord du maillage.
    pub sets: usize,
    pub results_open: bool,
    /// Espace Modèle : un Material, des BC, un chargement, et le nombre de diagnostics bloquants.
    pub material: bool,
    pub bcs: bool,
    pub loaded: bool,
    pub errors: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct State {
    pub space: Space,
    /// Outil actif de chaque espace, dans l'ordre de `Space::ALL`.
    tools: [Tool; 3],
    /// Cartes du parcours affichées (désactivables).
    pub cards: bool,
    /// Inspecteur de l'outil actif affiché (masquable, pour voir la pièce en entier).
    pub inspector: bool,
    pub progress: Progress,
}

impl Default for State {
    fn default() -> Self {
        Self {
            space: Space::Results,
            tools: [Tool::Navigate; 3],
            cards: true,
            inspector: true,
            progress: Progress::default(),
        }
    }
}

impl State {
    pub fn tool(&self) -> Tool {
        self.tools[self.space.index()]
    }

    /// Active un outil et montre son inspecteur.
    pub fn set_tool(&mut self, tool: Tool) {
        self.tools[self.space.index()] = tool;
        self.inspector = true;
    }

    /// Clic sur un outil de la barre : un autre outil s'active ; l'outil déjà actif masque ou remontre
    /// son inspecteur.
    pub fn pick(&mut self, tool: Tool) {
        if tool == self.tool() {
            self.inspector = !self.inspector;
        } else {
            self.set_tool(tool);
        }
    }

    /// Change d'espace en gardant l'outil que chacun avait.
    pub fn go(&mut self, space: Space) {
        self.space = space;
    }

    /// Les étapes de la carte d'un espace.
    pub fn checklist(&self, space: Space) -> Vec<Step> {
        let step = |label, status, tool, optional| Step {
            label,
            status,
            tool,
            optional,
        };
        use Status::{Done, Missing};
        match space {
            Space::Mesh => {
                let p = &self.progress;
                let done = |ok| if ok { Done } else { Missing };
                vec![
                    step(
                        "Importer une géométrie ou un maillage",
                        done(p.geometry),
                        Tool::Import,
                        false,
                    ),
                    step("Mailler", done(p.meshed), Tool::MeshSize, false),
                    step(
                        "Avoir au moins un Set",
                        done(p.sets > 0),
                        Tool::FaceSets,
                        false,
                    ),
                ]
            }
            Space::Model => {
                let p = &self.progress;
                let done = |ok| if ok { Done } else { Missing };
                let filled = p.material && p.bcs && p.loaded;
                let ready = if !filled {
                    Missing
                } else if p.errors > 0 {
                    Status::Warning
                } else {
                    Done
                };
                vec![
                    step(
                        "Choisir un Material",
                        done(p.material),
                        Tool::Material,
                        false,
                    ),
                    step("Poser les BC", done(p.bcs), Tool::Support, false),
                    step("Poser un chargement", done(p.loaded), Tool::Load, false),
                    step("Prêt à calculer", ready, Tool::Run, false),
                ]
            }
            Space::Results => {
                let open = if self.progress.results_open {
                    Done
                } else {
                    Missing
                };
                vec![
                    step("Ouvrir un calcul", open, Tool::Open, false),
                    step("Sonde, coupe ou courbe", Missing, Tool::Cut, true),
                ]
            }
        }
    }

    /// État d'un espace : fait quand toutes ses étapes obligatoires le sont, manquant sinon.
    pub fn space_status(&self, space: Space) -> Status {
        let steps = self.checklist(space);
        let required: Vec<_> = steps.iter().filter(|s| !s.optional).collect();
        if required.iter().all(|s| s.status == Status::Done) {
            Status::Done
        } else if required.iter().any(|s| s.status == Status::Done) {
            Status::Warning
        } else {
            Status::Missing
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_space_remembers_its_own_tool() {
        let mut s = State::default();
        s.go(Space::Results);
        s.set_tool(Tool::Cut);
        s.go(Space::Mesh);
        assert_eq!(s.tool(), Tool::Navigate);
        s.set_tool(Tool::Import);
        s.go(Space::Results);
        assert_eq!(s.tool(), Tool::Cut);
        s.go(Space::Mesh);
        assert_eq!(s.tool(), Tool::Import);
    }

    #[test]
    fn the_active_tool_hides_and_shows_its_inspector() {
        let mut s = State::default();
        s.set_tool(Tool::Cut);
        assert!(s.inspector);
        s.pick(Tool::Cut);
        assert!(
            !s.inspector && s.tool() == Tool::Cut,
            "le même outil masque"
        );
        s.pick(Tool::Cut);
        assert!(s.inspector, "le même outil remontre");
        s.inspector = false;
        s.pick(Tool::Field);
        assert!(
            s.inspector && s.tool() == Tool::Field,
            "un autre outil remontre"
        );
        s.inspector = false;
        s.set_tool(Tool::Open);
        assert!(s.inspector, "une carte qui active un outil remontre");
    }

    #[test]
    fn every_toolbar_starts_with_the_common_group_and_lists_each_tool_once() {
        for space in Space::ALL {
            let bar = toolbar(space);
            assert_eq!(bar[0].items[0], Item::Tool(Tool::Navigate), "{space:?}");
            assert!(bar[0].items.contains(&Item::Action(Action::Fit)));
            let all: Vec<Item> = bar.iter().flat_map(|g| g.items.clone()).collect();
            for item in &all {
                assert_eq!(all.iter().filter(|i| *i == item).count(), 1, "{item:?}");
            }
        }
    }

    #[test]
    fn a_space_shows_only_its_own_tools() {
        let tools = |space| -> Vec<Tool> {
            toolbar(space)
                .iter()
                .flat_map(|g| g.items.iter())
                .filter_map(|i| match i {
                    Item::Tool(t) => Some(*t),
                    Item::Action(_) => None,
                })
                .collect()
        };
        assert!(tools(Space::Results).contains(&Tool::Cut));
        assert!(!tools(Space::Results).contains(&Tool::Import));
        assert!(tools(Space::Mesh).contains(&Tool::Import));
        assert!(!tools(Space::Model).contains(&Tool::Field));
    }

    #[test]
    fn a_checklist_step_leads_to_a_tool_of_its_space() {
        let s = State::default();
        for space in Space::ALL {
            let bar: Vec<Item> = toolbar(space).into_iter().flat_map(|g| g.items).collect();
            for step in s.checklist(space) {
                assert!(
                    bar.contains(&Item::Tool(step.tool)),
                    "{space:?} : « {} » mène à un outil absent de la barre",
                    step.label
                );
            }
        }
    }

    #[test]
    fn opening_results_turns_the_results_space_done() {
        let mut s = State::default();
        assert_eq!(s.space_status(Space::Results), Status::Missing);
        s.progress.results_open = true;
        assert_eq!(s.space_status(Space::Results), Status::Done);
        // les autres espaces n'ont encore rien de fait
        assert_eq!(s.space_status(Space::Mesh), Status::Missing);
        assert_eq!(s.space_status(Space::Model), Status::Missing);
    }

    #[test]
    fn the_mesh_space_follows_the_import_the_meshing_and_the_sets() {
        let mut s = State::default();
        assert_eq!(s.space_status(Space::Mesh), Status::Missing);
        s.progress.geometry = true;
        assert_eq!(s.space_status(Space::Mesh), Status::Warning);
        s.progress.meshed = true;
        s.progress.sets = 4;
        assert_eq!(s.space_status(Space::Mesh), Status::Done);
        assert!(s
            .checklist(Space::Mesh)
            .iter()
            .all(|st| st.status == Status::Done));
    }

    #[test]
    fn the_model_space_is_ready_when_filled_without_blocking_diagnostics() {
        let mut s = State::default();
        assert_eq!(s.space_status(Space::Model), Status::Missing);
        s.progress.material = true;
        assert_eq!(s.space_status(Space::Model), Status::Warning);
        s.progress.bcs = true;
        s.progress.loaded = true;
        s.progress.errors = 1;
        let ready = |s: &State| s.checklist(Space::Model)[3].status;
        assert_eq!(ready(&s), Status::Warning, "rempli mais bloqué");
        s.progress.errors = 0;
        assert_eq!(ready(&s), Status::Done);
        assert_eq!(s.space_status(Space::Model), Status::Done);
    }

    #[test]
    fn only_the_tools_of_the_finished_spaces_are_available() {
        for tool in [
            Tool::Cut,
            Tool::Navigate,
            Tool::Import,
            Tool::MeshSize,
            Tool::FaceSets,
            Tool::ManageSets,
            Tool::Quality,
        ] {
            assert!(tool.available(), "{tool:?}");
        }
        for tool in [Tool::Probe, Tool::VertexSets, Tool::Curve, Tool::Capture] {
            assert!(!tool.available(), "{tool:?}");
        }
    }
}
