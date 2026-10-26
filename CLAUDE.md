# Maille — calcul par éléments finis, mécanique et thermique

Logiciel libre de simulation FEM en Rust, façon Abaqus mais plus simple. Projet d'apprentissage
qui doit pouvoir être utilisé par d'autres, y compris des non-experts. Projet perso public,
quelques heures par semaine, sans échéance.

Les notes de cadrage vivent dans le vault, catégorie 28 (`20-29 Mes Loisirs/28 Maille/`) :
état de l'art (28.11), architecture (28.21), stack (28.31), roadmap (28.41), journal (28.51).
Les lire avant toute décision d'architecture. Ce fichier ne les duplique pas.

## État du projet

M0 à M6 terminés : FEM 1D barre + poutre ; élasticité plane Tri3/Tri6 (dont Kirsch) ; viewer wgpu avec contraintes lissées aux Nodes ; 3D Tet4/Tet10 ; thermique stationnaire (conduction, source, convection) sur triangles et tétraèdres ; solveur creux maison (ligne de ciel + Cuthill-McKee inverse) ; mailleur 2D maison (`maille_core::mesher` : Delaunay contraint, raffinement de Ruppert, segments et arcs, tailles h par bord avec gradation, bords nommés pour les BC) ; thermique transitoire (`Step` à pas fixe, θ-méthode de Euler implicite à Crank-Nicolson, capacité cohérente avec ρ et c du Material, un `Increment` par pas) ; couplage thermo-mécanique séquentiel (champ de températures d'un calcul thermique donné au calcul statique, dilatation α sur Bar, triangles et tétraèdres) ; élasto-plasticité de von Mises en petites déformations (écrouissage isotrope linéaire, retour radial et tangente cohérente, Newton-Raphson par incrément, triangles en déformations planes et tétraèdres). Tout est validé contre l'analytique. Suite (M7 à M11 : formats de fichiers, performance, import CAO, mailleur 3D) dans la roadmap (28.41). Workspace Cargo : `crates/maille-core` (lib pure, zéro dépendance), `crates/maille-io` (fichiers : `.msh` Gmsh 4.1 ASCII et modèle TOML, M7), `crates/maille-viewer` (client 2D : wgpu, winit, glyphon) et `crates/maille-cli` (client : `maille-cli kirsch`, `maille-cli cylindre`, `maille-cli trempe`, `maille-cli cylindre-thermo` et `maille-cli cylindre-plastique` maillent leur géométrie puis ouvrent le viewer).

Commandes :
- `cargo build`
- `cargo test`
- `cargo clippy --all-targets -- -D warnings`
- `cargo fmt`

## Périmètre

- Dedans : mécanique statique linéaire, thermique (stationnaire puis transitoire), couplage thermo-mécanique, plasticité
  de von Mises en petites déformations (M6), fichiers de maillage, de modèle et de résultats (M7), performance du solveur
  (M8), import CAO en 2D (DXF, M9) et en 3D (STL, mailleur 3D, M10 ; STEP, M11).
- Dehors, sauf demande explicite : fluides, électromagnétisme, contact général, grandes déformations,
  calcul parallèle/GPU, interface graphique (viendra une fois les formats stables).
- Ne pas élargir le périmètre de ta propre initiative. En cas de doute, demander.

## Ordre de travail

Le plan suit les jalons M0 à M11 de la roadmap (28.41) : FEM 1D, puis 2D, viewer minimal, 3D,
mailleur 2D, thermique transitoire, non-linéarités (M0 à M6, terminés), puis flux de fichiers (M7),
performance (M8), import DXF 2D (M9), mailleur 3D et viewer 3D (M10), import STEP (M11). Les autres
non-linéarités (contact, grandes déformations) passent en M12+. Ne pas sauter de jalon : le viewer et
le mailleur arrivent après un solveur validé, pas avant.

Flux de fichiers cible : CAO (DXF, STL, STEP) → mailleur → fichier de maillage (Gmsh `.msh` 4.1, ensembles nommés =
physical groups) ; fichier de modèle (TOML) qui y renvoie par son nom ; `maille solve` écrit le fichier de résultats
(VTK `.vtu` + `.pvd`). Choix du STEP (crate Rust, OpenCASCADE ou outil externe) : spike avant M11, et décision de Tim
(licence : attention au LGPL).

## Règles transverses

- **Un solveur faux n'a aucune valeur.** Tout élément et toute physique se valident contre une solution
  analytique (poutre en flexion, plaque trouée de Kirsch, barre en régime thermique…) avant d'avancer.
  Écrire le test de validation avant ou avec le code, pas après.
- Le cœur de calcul (maillage, éléments, assemblage, solveur) est pur : pas d'I/O, pas de dépendance
  au viewer ni à l'interface. Le CLI et le viewer sont des clients du cœur.
- Unités : le cœur ne convertit rien. L'utilisateur fournit un système d'unités cohérent (comme Abaqus).
  Documenter ce choix dans le fichier d'entrée.
- Pas de `unwrap()` ni de panique dans le cœur sur une entrée utilisateur. Erreurs typées, avec un message
  qui dit quoi corriger : le public visé n'est pas expert en FEM.
- Flottants : `f64` partout dans le solveur. Comparer avec une tolérance, jamais avec `==`.
- Pas de dépendance lourde ajoutée sans la noter dans 28.31 Stack (choix + raison + licence).
  La licence du projet n'est pas encore choisie : vérifier la compatibilité avant d'ajouter une crate (attention au LGPL/GPL).

## Vocabulaire

Contraignant, ne pas inventer de synonymes : **Mesh** (maillage), **Node**, **Element**, **DOF** (degré de liberté),
**Material**, **BC** (condition aux limites : Dirichlet / Neumann), **Load**, **Step**, **Field** (champ de résultats).
Thermique (M3) : **Analysis** (`Static` / `Thermal`, la physique résolue), **HeatSource** (source volumique),
**Convection** (échange h (T − T∞) sur une face, champs `film` et `sink` comme Abaqus). **Step** reste réservé à
l'étape de calcul avec son temps (M5).
Transitoire (M5) : **Step** porte le temps (durée, nombre d'incréments, θ) et s'applique à un Model dont l'Analysis est
`Thermal` ; l'Analysis reste sur le Model. **Increment** : résultat d'un pas (temps + Field).
Couplage (M5) : séquentiel, comme Abaqus : le **champ de températures** d'un Field thermique est donné au Model statique
(`with_temperatures`, dilatation α (T − T_ref)). Ce n'est pas une BC : « BC de température » reste le Dirichlet sur `T`.
Mailleur (M4) : **Geometry** (domaine : contour extérieur + trous), **Curve** (morceau de bord, segment ou arc, avec sa
taille h et son nom), **Triangulation** (sortie du mailleur : Nodes, triangles, bords nommés ; ce n'est pas un Mesh),
**BoundarySet** (Nodes et faces d'un bord nommé).
Plasticité (M6) : **Plasticity** (σy0 et H du Material, `with_plasticity`), **VonMises** (la loi au point de Gauss),
**PlasticState** (εp et p d'un point de Gauss), **p** (déformation plastique cumulée, `equivalent_plastic_strains` du
Field ; « εp » désigne le tenseur), **Amplitude** (λ(t) linéaire par morceaux qui multiplie les Load et les valeurs des BC).
Le calcul **quasi-statique** (`Step::quasi_static`, `Model::solve_quasi_static`) résout un équilibre par Newton-Raphson
à chaque Increment ; ses `residuals` sont les résidus relatifs des itérations. La Step reste le temps seul.
Fichiers (M7) : **Set** (ensemble nommé de Nodes, de faces ou d'Elements : un groupe physique de Gmsh), **MshFile**
(contenu d'un `.msh` : Nodes, connectivité, Sets ; ce n'est pas un Mesh, il n'a pas de Material). Dans le futur fichier de
modèle, « sections » reprend « section assignment » d'Abaqus ; dans le cœur, `Section` reste l'aire d'une Bar.

## Workflow

Rôles : Claude écrit tout le code. Tim dirige, relit et valide ; il ne code pas. Les règles de dépôt (branches,
commits, PR) sont dans `CONTRIBUTING.md` et s'appliquent à Claude sans exception.

- Dépôt : `github.com/timothee-faget/maille` (public, MIT OR Apache-2.0). `main` est protégée : jamais de commit
  ni de push direct dessus, tout passe par une PR avec la CI verte (fmt, clippy, test, cargo-deny).
- Une branche par sujet, depuis `main` à jour : `<type>/<sujet>` (`feat/`, `fix/`, `refactor/`, `docs/`, `ci/`,
  `chore/`). Un jalon peut se découper en plusieurs PR, une par étape relisible.
- Commits en Conventional Commits, en français (`feat(core): …`, portées `core`, `mesher`, `viewer`, `cli`).
  Petits et cohérents ; chaque commit compile et passe les tests.
- Avant de pousser : `cargo fmt --all`, `cargo clippy --all-targets -- -D warnings`, `cargo test`,
  `cargo deny check licenses sources`. Montrer la sortie plutôt que d'affirmer que ça passe.
- PR : description selon `.github/pull_request_template.md` (quoi, pourquoi, validation analytique avec les écarts).
  Ensuite, s'arrêter et attendre l'OK explicite de Tim dans le chat. Après l'OK : `gh pr merge --squash
  --delete-branch`, puis revenir sur `main` à jour. Ne jamais fusionner sans cet OK, même si la CI est verte.
- Fin de jalon : tag annoté `vX.Y.0` sur `main`, proposé à Tim avant d'être poussé.
- Décider avant de coder : cadrage d'un jalon, choix d'architecture, nouvelle dépendance, changement d'API publique
  → proposer à Tim et attendre sa décision. Les détails d'implémentation restent à Claude.
- Surgical changes : ne toucher que ce que la demande exige, ne pas refactorer au passage.
- Quand un jalon avance ou qu'une décision est prise, la noter dans le vault (28.51 Journal ou 28.31 Stack)
  en respectant les règles du vault du CLAUDE.md global.
