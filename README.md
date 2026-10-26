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
- Mailleur 2D (Delaunay contraint, raffinement de Ruppert) et viewer natif (wgpu)

Le projet est en cours de développement : l'API n'est pas stable.

## Démarrer

```sh
cargo build --release
cd examples
../target/release/maille solve traction.toml     # maillage .msh + modèle .toml → résultats .vtu/.pvd
../target/release/maille view traction.pvd       # ouvre le viewer 2D
```

Le flux de travail (maillage Gmsh `.msh`, modèle `.toml`, résultats VTK lisibles par ParaView) et les exemples
sont décrits dans [`examples/README.md`](examples/README.md).

Des cas d'exemple, qui maillent leur géométrie, comparent le résultat à la solution analytique, puis ouvrent le
viewer : `maille console`, `kirsch`, `cylindre`, `trempe`, `cylindre-thermo` et `cylindre-plastique`
(par exemple `cargo run --release -p maille-cli -- kirsch`).

## Organisation

| Crate | Rôle |
|---|---|
| `maille-core` | Cœur de calcul pur (maillage, éléments, assemblage, solveur), sans dépendance ni I/O |
| `maille-io` | Fichiers : maillage Gmsh `.msh`, modèle TOML, résultats VTK |
| `maille-viewer` | Viewer 2D (wgpu, winit, glyphon) |
| `maille-cli` | Le binaire `maille` : `solve`, `view` et les cas d'exemple |

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
