# Maille — calcul par éléments finis, mécanique et thermique

Logiciel libre de simulation FEM en Rust, façon Abaqus mais plus simple. Projet d'apprentissage
qui doit pouvoir être utilisé par d'autres, y compris des non-experts. Projet perso public,
quelques heures par semaine, sans échéance.

Les notes de cadrage vivent dans le vault, catégorie 28 (`20-29 Mes Loisirs/28 Maille/`) :
état de l'art (28.11), architecture (28.21), stack (28.31), roadmap (28.41), journal (28.51).
Les lire avant toute décision d'architecture. Ce fichier ne les duplique pas.

## État du projet

M0 à M4 terminés : FEM 1D barre + poutre ; élasticité plane Tri3/Tri6 (dont Kirsch) ; viewer wgpu avec contraintes lissées aux Nodes ; 3D Tet4/Tet10 ; thermique stationnaire (conduction, source, convection) sur triangles et tétraèdres ; solveur creux maison (ligne de ciel + Cuthill-McKee inverse) ; mailleur 2D maison (`maille_core::mesher` : Delaunay contraint, raffinement de Ruppert, segments et arcs, tailles h par bord avec gradation, bords nommés pour les BC). Tout est validé contre l'analytique. Prochain jalon : M5, thermique transitoire. Workspace Cargo : `crates/maille-core` (lib pure, zéro dépendance), `crates/maille-viewer` (client 2D : wgpu, winit, glyphon) et `crates/maille-cli` (client : `maille-cli kirsch` et `maille-cli cylindre` maillent leur géométrie puis ouvrent le viewer).

Commandes :
- `cargo build`
- `cargo test`
- `cargo clippy --all-targets -- -D warnings`
- `cargo fmt`

## Périmètre

- Dedans : mécanique statique linéaire, thermique (stationnaire puis transitoire), puis couplage thermo-mécanique.
- Dehors, sauf demande explicite : fluides, électromagnétisme, contact général, grandes déformations,
  calcul parallèle/GPU, import CAO (STEP).
- Ne pas élargir le périmètre de ta propre initiative. En cas de doute, demander.

## Ordre de travail

Le plan suit les jalons M0 à M6 de la roadmap (28.41) : FEM 1D, puis 2D, viewer minimal, 3D,
mailleur 2D, thermique transitoire, non-linéarités. Ne pas sauter de jalon : le viewer et
le mailleur arrivent après un solveur validé, pas avant.

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
l'étape de calcul avec son temps (M5) : une Step portera une Analysis.
Mailleur (M4) : **Geometry** (domaine : contour extérieur + trous), **Curve** (morceau de bord, segment ou arc, avec sa
taille h et son nom), **Triangulation** (sortie du mailleur : Nodes, triangles, bords nommés ; ce n'est pas un Mesh),
**BoundarySet** (Nodes et faces d'un bord nommé).

## Workflow

- Surgical changes : ne toucher que ce que la demande exige, ne pas refactorer au passage.
- Avant de proposer un commit : `cargo fmt`, `cargo clippy --all-targets -- -D warnings`, `cargo test`.
  Montrer la sortie plutôt que d'affirmer que ça passe.
- Quand un jalon avance ou qu'une décision est prise, la noter dans le vault (28.51 Journal ou 28.31 Stack)
  en respectant les règles du vault du CLAUDE.md global.
