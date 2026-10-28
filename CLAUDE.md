# Maille — calcul par éléments finis, mécanique et thermique

Logiciel libre de simulation FEM en Rust, façon Abaqus mais plus simple. Projet d'apprentissage
qui doit pouvoir être utilisé par d'autres, y compris des non-experts. Projet perso public,
quelques heures par semaine, sans échéance.

Les notes de cadrage vivent dans le vault, catégorie 28 (`20-29 Mes Loisirs/28 Maille/`) :
état de l'art (28.11), architecture (28.21), stack (28.31), roadmap (28.41), journal (28.51).
Les lire avant toute décision d'architecture. Ce fichier ne les duplique pas.

## État du projet

M0 à M11 (dont M7b et M10b) terminés : FEM 1D barre + poutre ; élasticité plane Tri3/Tri6 (dont Kirsch) ; viewer wgpu avec contraintes lissées aux Nodes ; 3D Tet4/Tet10 ; thermique stationnaire (conduction, source, convection) sur triangles et tétraèdres ; solveur creux maison (ligne de ciel + Cuthill-McKee inverse) ; mailleur 2D maison (`maille_core::mesher` : Delaunay contraint, raffinement de Ruppert, segments et arcs, tailles h par bord avec gradation, bords nommés pour les BC) ; thermique transitoire (`Step` à pas fixe, θ-méthode de Euler implicite à Crank-Nicolson, capacité cohérente avec ρ et c du Material, un `Increment` par pas) ; couplage thermo-mécanique séquentiel (champ de températures d'un calcul thermique donné au calcul statique, dilatation α sur Bar, triangles et tétraèdres) ; élasto-plasticité de von Mises en petites déformations (écrouissage isotrope linéaire, retour radial et tangente cohérente, Newton-Raphson par incrément, triangles en déformations planes et tétraèdres) ; performance (M8) : banc criterion (`cargo bench -p maille-core`), Cholesky deux fois plus rapide, solveur itératif (gradients conjugués préconditionnés par IC(0), `SolverKind`, `Auto` à 100 000 DOF libres ; une console Tet10 de 10⁶ DOF se résout en ~4 min) ; import DXF 2D (M9) : `maille_io::dxf::read_geometry` lit LINE, ARC, CIRCLE et LWPOLYLINE (arrondis compris) et le calque de chaque entité nomme le bord ; mailleur 3D et viewer 3D (M10) : `maille_io::stl::read_surface` lit un STL ASCII ou binaire en `Surface` fermée (un Patch par `solid`, sinon découpage aux arêtes vives), `maille_core::mesher::tetrahedralize` en fait des tétraèdres conformes à la surface (Delaunay 3D, récupération du bord, volume contrôlé), le viewer montre la peau des tétraèdres avec un plan de coupe et `maille view --png` capture hors écran. Taille h en 3D (M10b) : la surface est regroupée en facettes planes, remaillée facette par facette par le mailleur 2D, puis l'intérieur est raffiné (rapport rayon / arête visé à 2, surface figée) ; validé par le tube de Lamé 3D en Tet10. Import STEP (M11) : `maille_io::step::read_brep` lit un solide (faces planes,
cylindriques, coniques, sphériques et toriques, arêtes en segments, arcs et B-splines) et
`maille_core::mesher::tetrahedralize_brep` le maille à la taille h, chaque face dans une carte plane (déroulée, ou
conforme pour la sphère et le tore) ; validé par le tube de Lamé lu en B-rep exact et la sphère creuse de Lamé
(`examples/sphere.step`). Tout est validé contre l'analytique. Suite (M12+) dans la roadmap (28.41). Workspace Cargo : `crates/maille-core` (lib pure, zéro dépendance), `crates/maille-io` (fichiers : `.msh` Gmsh 4.1 ASCII, modèle TOML et résultats VTK, M7 ; dessin DXF 2D, M9 ; surface STL, M10 ; pièce STEP, M11), `crates/maille-viewer` (client 2D et 3D : wgpu, winit, glyphon) et `crates/maille-cli` (le binaire `maille` : `maille solve modele.toml` écrit les résultats `.vtu` et `.pvd`, `maille view resultats.pvd` ouvre le viewer, `maille mesh dessin.dxf --size h [--quadratic]` écrit un `.msh` depuis un dessin DXF, `maille mesh piece.stl [--size h] [--feature-angle a] [--quadratic]` depuis un STL, `maille mesh piece.step --size h
[--quadratic]` depuis un STEP ; les cas `maille kirsch`, `cylindre`, `trempe`, `cylindre-thermo` et `cylindre-plastique` maillent leur géométrie puis ouvrent le viewer). Exemples et guide des fichiers dans `examples/`.

Commandes :
- `cargo build`
- `cargo test`
- `cargo clippy --all-targets -- -D warnings`
- `cargo fmt`

## Périmètre

- Dedans : mécanique statique linéaire, thermique (stationnaire puis transitoire), couplage thermo-mécanique, plasticité
  de von Mises en petites déformations (M6), fichiers de maillage, de modèle et de résultats (M7), performance du solveur
  (M8, fait), import CAO en 2D (DXF, M9, fait) et en 3D (STL, mailleur 3D, M10 et M10b, faits ; STEP, M11, fait).
  Grandes déformations (M12, en cours) : Lagrangien total sur triangles en déformations planes et tétraèdres,
  Saint-Venant–Kirchhoff et néo-Hookéen, pression suiveuse, découpe automatique de l'incrément ; cadrage dans
  28.23 Non-linéarités.
- Dehors, sauf demande explicite : fluides, électromagnétisme, contact général,
  calcul parallèle/GPU, interface graphique (viendra une fois les formats stables).
- Ne pas élargir le périmètre de ta propre initiative. En cas de doute, demander.

## Ordre de travail

Le plan suit les jalons M0 à M11 de la roadmap (28.41) : FEM 1D, puis 2D, viewer minimal, 3D,
mailleur 2D, thermique transitoire, non-linéarités (M0 à M6), puis flux de fichiers (M7),
performance (M8), import DXF 2D (M9), mailleur 3D et viewer 3D (M10), import STEP (M11), tous terminés. Les autres
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
Dessin (M9) : **Drawing** n'existe pas dans le code, on dit « dessin DXF » ; le **calque** DXF d'une entité devient le nom du bord (une **Curve** de la **Geometry**), donc le **Set** où poser les BC. `read_geometry(texte, h)` donne une `Geometry` avec la même taille h partout.
Mailleur 3D (M10) : **Surface** (surface fermée d'une pièce : sommets soudés, triangles à normale sortante, entrée du
mailleur 3D ; elle vit dans `maille_core::mesher` et `maille_io::stl` la réexporte), **Patch** (morceau nommé de la Surface :
un `solid` du STL ou une zone limitée par des arêtes vives, `surface_1`, `surface_2`… ; son nom devient le Set des BC),
**Tetrahedralization** (sortie du mailleur 3D : Nodes, tétraèdres, bords nommés ; ce n'est pas un Mesh, comme la
Triangulation en 2D). **TetOptions** (M10b : `size` = taille h optionnelle, `grading` ; `tetrahedralize(&surface, &options)`) :
avec une taille h, la surface est remaillée facette par facette puis l'intérieur raffiné ; `None` garde la finesse du STL.
Une **facette** est un morceau plan et connexe d'un Patch (interne au mailleur, ce n'est pas un Set).
Import STEP (M11) : **Brep** (pièce d'un seul solide lue d'un STEP : sommets, arêtes, faces ; `maille_io::step::read_brep`
la donne avec l'unité de longueur du fichier, sans rien convertir), **BrepFace** (face de la pièce, à ne pas confondre avec
la face d'un Element ; son nom devient un Patch), **BrepEdge** (arête entre deux sommets), **EdgeCurve** (courbe d'une arête :
segment, arc de cercle ou d'ellipse, morceau de **BSpline** ; ce n'est pas la Curve 2D du mailleur), **FaceSurface**
(surface d'une BrepFace : plan, cylindre, cône, sphère, tore ; ce n'est pas la `Surface` triangulée d'un STL ; la face
est orientée par ses contours, et une sphère dit en plus si elle est bombée, `convex`). `tetrahedralize_brep(&brep, &options)` exige une taille h.
Solveur (M8) : **SolverKind** (`Auto` / `Direct` / `Iterative`, `Model::with_solver`) : le Cholesky en ligne de ciel (direct) ou les gradients conjugués préconditionnés par un Cholesky incomplet IC(0) (itératif). Ce n'est pas un trait : deux solveurs, un `enum`.
Grandes déformations (M12) : **Kinematics** (`Small` / `Large`, `Model::with_kinematics` ; `Large` = **grandes déformations**, par `solve_quasi_static` seulement), **raideur géométrique** (le terme de la tangente en contrainte S, en plus de la raideur du Material), **contrainte de Cauchy** (σ = J⁻¹ F S Fᵀ sur la configuration déformée : ce que le Field porte en `Large`). Loi par défaut : Saint-Venant–Kirchhoff, avec E et ν du Material. **Hyperelasticity** (`Material::with_hyperelasticity`, `NeoHookean` : loi du Material en `Large` seulement, à la place de SVK). **Pression suiveuse** : en `Large`, une Pressure reste normale à la face déformée (pas une force morte). **Coupe de l'incrément** : si Newton ne converge pas, `solve_quasi_static` coupe l'incrément en deux, jusqu'à 5 fois (1/32) ; seuls les Increments de la Step sont gardés (vaut aussi pour la plasticité).
Pression (M7b) : **Pressure** (pression uniforme sur une face d'Element, `with_pressures`, statique ; positive, elle comprime). Distincte de la **Convection** (thermique) et du **Load** (force par Node).
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
- Commits en Conventional Commits, en français (`feat(core): …`, portées `core`, `mesher`, `io`, `viewer`, `cli`).
  Petits et cohérents ; chaque commit compile et passe les tests.
- Avant de pousser : `cargo fmt --all`, `cargo clippy --all-targets -- -D warnings`, `cargo test`,
  `cargo deny check licenses sources`. Montrer la sortie plutôt que d'affirmer que ça passe.
- PR : description selon `.github/pull_request_template.md` (quoi, pourquoi, validation analytique avec les écarts).
  Ensuite, s'arrêter et attendre l'OK explicite de Tim dans le chat. Après l'OK : `gh pr merge --squash
  --delete-branch`, puis revenir sur `main` à jour. Ne jamais fusionner sans cet OK, même si la CI est verte.
- Fin de jalon : la PR de clôture passe la version du workspace (`Cargo.toml`) à `X.Y.0`, puis tag annoté `vX.Y.0` sur `main`,
  proposé à Tim avant d'être poussé.
- Décider avant de coder : cadrage d'un jalon, choix d'architecture, nouvelle dépendance, changement d'API publique
  → proposer à Tim et attendre sa décision. Les détails d'implémentation restent à Claude.
- Surgical changes : ne toucher que ce que la demande exige, ne pas refactorer au passage.
- Quand un jalon avance ou qu'une décision est prise, la noter dans le vault (28.51 Journal ou 28.31 Stack)
  en respectant les règles du vault du CLAUDE.md global.
