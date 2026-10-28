# Contribuer à Maille

## Workflow

1. Créer une branche depuis `main` : `<type>/<sujet-court>` (par exemple `feat/hexa8` ou `fix/convection-tet10`).
2. Faire des commits petits et cohérents, au format [Conventional Commits](https://www.conventionalcommits.org/fr/).
3. Ouvrir une Pull Request vers `main`, en suivant le modèle de description. La CI doit être verte.
4. La PR est fusionnée en squash après relecture. `main` est protégée : aucun push direct.

## Messages de commit

```
<type>(<portée>): <résumé à l'impératif, en minuscules, sans point final>

<corps facultatif : le pourquoi, pas le comment>
```

- Types : `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `build`, `ci`, `chore`.
- Portées : `core`, `mesher`, `io`, `viewer`, `cli`. On l'omet pour un changement transverse.
- Un changement incompatible de l'API se signale par `!` (`feat(core)!: …`) et une ligne `BREAKING CHANGE:`.

Exemples : `feat(core): ajouter l'élément Hexa8`, `fix(viewer): recadrer la déformée au zoom`.

## Avant d'ouvrir une PR

```sh
cargo fmt --all
cargo clippy --all-targets -- -D warnings
cargo test
cargo deny check licenses sources
```

## Règles du code

- **Un solveur faux n'a aucune valeur** : tout nouvel élément ou toute nouvelle physique arrive avec un test qui
  le valide contre une solution analytique.
- `maille-core` reste pur : pas d'I/O, pas de dépendance vers le viewer ni vers le CLI.
- Pas de `unwrap()` ni de panique dans le cœur sur une entrée utilisateur. Les erreurs sont typées, avec un
  message qui dit quoi corriger.
- `f64` partout dans le solveur. Les flottants se comparent avec une tolérance, jamais avec `==`.
- Une nouvelle dépendance doit avoir une licence permissive compatible avec MIT et Apache-2.0 (voir
  `deny.toml`).
