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

Le projet est en cours de développement : l'API n'est pas stable et il n'y a pas encore de format de fichier
d'entrée.

## Démarrer

```sh
cargo build --release
cargo run --release -p maille-cli -- kirsch
```

Les cas disponibles sont `console`, `kirsch`, `cylindre`, `trempe`, `cylindre-thermo` et
`cylindre-plastique`. Chacun maille sa géométrie, compare le résultat à la solution analytique, puis ouvre le
viewer.

## Organisation

| Crate | Rôle |
|---|---|
| `maille-core` | Cœur de calcul pur (maillage, éléments, assemblage, solveur), sans dépendance ni I/O |
| `maille-viewer` | Viewer 2D (wgpu, winit, glyphon) |
| `maille-cli` | Cas d'exemple en ligne de commande |

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
