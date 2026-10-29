# Maille

Logiciel libre de calcul par éléments finis, en mécanique et en thermique, écrit en Rust. L'esprit est celui
d'Abaqus, en plus simple : moins de capacités et une interface plus claire, pour un public qui ne connaît pas
forcément les éléments finis.

Chaque élément et chaque physique sont validés contre une solution analytique avant d'être intégrés.

## Ce que Maille sait faire

- Barres et poutres 1D, élasticité plane (Tri3, Tri6), 3D (Tet4, Tet10)
- Thermique stationnaire et transitoire (conduction, source, convection)
- Couplage thermo-mécanique séquentiel (dilatation thermique)
- Plasticité de von Mises en petites déformations (écrouissage isotrope linéaire)
- Grandes déformations (Lagrangien total, Saint-Venant–Kirchhoff et néo-Hookéen), pression suiveuse
- Contact sans frottement contre une surface rigide (plan, cylindre, sphère), en petites et grandes déformations
- Pression sur une face, en 2D et en 3D
- Solveur creux direct (Cholesky en ligne de ciel) et itératif (gradients conjugués préconditionnés), choisi
  selon la taille du problème
- Mailleur 2D (Delaunay contraint, raffinement de Ruppert) et mailleur 3D (tétraèdres, à la taille h)
- Import CAO : dessin DXF 2D, surface STL 3D, pièce STEP 3D (faces planes, cylindriques, coniques,
  sphériques et toriques)
- Viewer natif 2D et 3D (wgpu), avec plan de coupe et capture en image PNG
- Interface graphique (egui) : importer une pièce, la mailler, créer des Sets à la souris, remplir le modèle
  avec des diagnostics en clair, lancer le calcul en tâche de fond, parcourir les résultats (incréments,
  sonde, courbe, capture) ; le modèle reste un fichier TOML rejouable par `maille solve`

Jalons M0 à M13 terminés (version 0.13.0). Le projet est en cours de développement : l'API n'est pas stable.

## Démarrer

```sh
cargo build --release
./target/release/maille                              # ouvre l'interface graphique
cd examples
../target/release/maille solve traction.toml     # maillage .msh + modèle .toml → résultats .vtu/.pvd
../target/release/maille view traction.pvd       # ouvre le viewer
../target/release/maille mesh sphere.step --size 2 --quadratic   # maille une pièce STEP → sphere.msh
```

Le flux de travail (CAO DXF, STL ou STEP, maillage Gmsh `.msh`, modèle `.toml`, résultats VTK lisibles par
ParaView) et les exemples (dont `console_gd`, `hertz` et `balle` pour les grandes déformations et le contact)
sont décrits dans [`examples/README.md`](examples/README.md).

Des cas d'exemple, qui maillent leur géométrie, comparent le résultat à la solution analytique, puis ouvrent le
viewer : `maille console`, `kirsch`, `cylindre`, `trempe`, `cylindre-thermo` et `cylindre-plastique`
(par exemple `cargo run --release -p maille-cli -- kirsch`).

## Organisation

| Crate | Rôle |
|---|---|
| `maille-core` | Cœur de calcul pur (maillage, éléments, assemblage, solveur, mailleurs 2D et 3D), sans dépendance ni I/O |
| `maille-io` | Fichiers : dessin DXF, surface STL, pièce STEP, maillage Gmsh `.msh`, modèle TOML, résultats VTK |
| `maille-viewer` | Viewer 2D et 3D (wgpu, winit, glyphon), embarquable dans l'interface |
| `maille-gui` | Interface graphique (egui) : espaces Maillage, Modèle et Résultats |
| `maille-cli` | Le binaire `maille` : l'interface (sans argument), `solve`, `view`, `mesh` et les cas d'exemple |

## Unités

Le cœur ne convertit aucune unité. C'est à vous de fournir un système cohérent, par exemple N, mm, MPa
(voir la documentation de `maille-core`).

## Contribuer

Voir [CONTRIBUTING.md](CONTRIBUTING.md).

## Licence

Au choix, sous l'une ou l'autre de ces licences :

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT ([LICENSE-MIT](LICENSE-MIT))

Sauf mention contraire de votre part, toute contribution que vous soumettez volontairement pour inclusion dans
Maille, au sens de la licence Apache-2.0, est sous cette double licence, sans condition supplémentaire.
