# Maille — calcul par éléments finis, mécanique et thermique

Logiciel libre de simulation FEM en Rust, façon Abaqus mais plus simple. Projet d'apprentissage
qui doit pouvoir être utilisé par d'autres, y compris des non-experts. Projet perso public,
quelques heures par semaine, sans échéance.

Les notes de cadrage vivent dans le vault, catégorie 28 (`20-29 Mes Loisirs/28 Maille/`) :
état de l'art (28.11), architecture (28.21), stack (28.31), roadmap (28.41), journal (28.51).
Les lire avant toute décision d'architecture. Ce fichier ne les duplique pas.

## État du projet

M0 terminé (FEM 1D barre + poutre, validé contre l'analytique). M1 à venir. Workspace Cargo : `crates/maille-core` (lib pure, zéro dépendance) et `crates/maille-cli` (client).

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

## Workflow

- Surgical changes : ne toucher que ce que la demande exige, ne pas refactorer au passage.
- Avant de proposer un commit : `cargo fmt`, `cargo clippy --all-targets -- -D warnings`, `cargo test`.
  Montrer la sortie plutôt que d'affirmer que ça passe.
- Quand un jalon avance ou qu'une décision est prise, la noter dans le vault (28.51 Journal ou 28.31 Stack)
  en respectant les règles du vault du CLAUDE.md global.
