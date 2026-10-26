# Exemples et guide des fichiers

Maille résout un modèle en trois fichiers : **maillage → modèle → résultats**. En 2D, le maillage
peut venir d'un dessin DXF (`maille mesh`).

```
dessin (.dxf) ── maille mesh ──▶
maillage (.msh) ──┐
                  ├─ maille solve ──▶ résultats (.pvd + .vtu) ──▶ maille view / ParaView
modèle (.toml) ───┘                                  │
        ▲                                            │
        └──── [temperatures] (un calcul lit les résultats thermiques d'un autre)
```

Depuis ce dossier :

```sh
maille solve traction.toml        # écrit traction.pvd et traction.vtu
maille view traction.pvd          # ouvre le viewer (2D) ; dernier incrément par défaut
maille view trempe.pvd 50         # un incrément précis, depuis 0
```

| Exemple | Ce qu'il montre | Solution exacte |
|---|---|---|
| `traction` | plaque en contraintes planes, Tri6 | champ linéaire |
| `trempe` | thermique transitoire, 200 incréments | série de Fourier |
| `tube_thermique` puis `tube_contraintes` | couplage thermo-mécanique (lancer dans cet ordre) | Timoshenko et Goodier |
| `tube_pression` | pression sur un bord, déformations planes | Lamé |
| `cube` | 3D, Tet4 (ouvrir dans ParaView) | champ linéaire |

Les tests lancent chaque exemple et le comparent à sa solution exacte : ils restent justes.

## Unités

Maille ne convertit rien : choisissez un système d'unités cohérent et écrivez-le dans `units`
(texte libre, recopié dans les résultats). Par exemple N, mm, MPa en mécanique ; W, mm, K en
thermique, plus s et kg en transitoire.

## Le dessin : DXF 2D ASCII

```sh
maille mesh plaque.dxf --size 0.5 --quadratic   # écrit plaque.msh (Tri6 ; sans --quadratic, Tri3)
```

`--size` est la taille h des Elements, dans les unités du dessin (Maille ne convertit rien).
`-o` choisit le fichier de sortie. Entités lues : `LINE`, `ARC`, `CIRCLE` et `LWPOLYLINE` (arrondis
compris), dans l'ordre et le sens qu'on veut. Les contours doivent être **fermés** : les extrémités
se touchent (à un millionième de la diagonale du dessin) et deux entités seulement se rejoignent en
un point. Le contour de plus grande aire est le contour extérieur, les autres sont des trous, et un
angle de moins de 60° entre deux entités est refusé. Le **calque** de chaque entité est le nom du
bord, donc le Set où poser les BC (`0` donne `bord`). Le Set `domaine` contient tous les triangles.
Splines, ellipses, blocs et DXF binaire ne sont pas lus. `plaque.dxf`, `trempe.dxf` et `tube.dxf`
sont les dessins des exemples.

## Le maillage : Gmsh `.msh` 4.1 ASCII

Produit par [Gmsh](https://gmsh.info) (outil externe) avec `Mesh.MshFileVersion = 4.1` et
`Mesh.Binary = 0`. Éléments lus : triangles à 3 ou 6 Nodes (2D), tétraèdres à 4 ou 10 Nodes (3D),
dans le plan z = 0 pour le 2D. Les **groupes physiques** de Gmsh deviennent des **Sets** nommés :

- un groupe de surface (2D) ou de volume (3D) est un Set d'Elements : `domaine`, par exemple ;
- un groupe de lignes (2D) ou de surfaces (3D) est un bord : ses Nodes reçoivent les BC, ses faces
  la convection ;
- un groupe de points est un Set de Nodes.

Les maillages d'ordre 2 à côtés courbes sont refusés (Maille veut des côtés droits).

## Le modèle : TOML

Une clé inconnue ou mal écrite est une erreur, avec la ligne et la colonne. Les chemins (`mesh`,
`from`) sont relatifs au dossier du modèle.

| Clé | Rôle |
|---|---|
| `units`, `mesh`, `analysis` | obligatoires ; `analysis` vaut `"static"` ou `"thermal"` |
| `[materials.nom]` | `young`, `poisson`, `conductivity`, `density`, `specific_heat`, `expansion`, et `yield_stress` + `hardening` pour la plasticité de von Mises |
| `[[sections]]` | `set`, `material` ; triangles : `thickness` (1 par défaut) et `plane_state` (`"stress"` ou `"strain"`, obligatoire en statique). Tout Element a exactement une section |
| `[[bcs]]` | `set` de Nodes, `dof` (`ux`, `uy`, `uz` en statique ; `t` en thermique), `value` |
| `[[loads]]` | comme `[[bcs]]` : une force **par Node** (pas la force totale) |
| `[[heat_sources]]` | `set` d'Elements, `value` (puissance par volume) |
| `[[convections]]` | `set` de bord, `film` (h) et `sink` (T∞) |
| `[[pressures]]` | `set` de bord, `value` : une pression positive comprime (contre la normale sortante), négative elle tire. Statique, triangles et tétraèdres |
| `[temperatures]` | couplage : `from` (un `.pvd` thermique sur le même maillage), `increment`, `reference` (T_ref) |
| `[step]` | absent : calcul linéaire. Sinon `kind`, `duration`, `increments` |

`[step]` :
- `kind = "transient"` (analyse thermique) : `initial_temperature` obligatoire, `theta` optionnel
  (1 = Euler implicite, 0,5 = Crank-Nicolson) ;
- `kind = "quasi_static"` (analyse statique, plasticité) : `amplitude = [[0, 0], [1, 1]]` multiplie
  les charges et les valeurs des BC (rampe de 0 à 1 par défaut). Le premier point est (0, 0).

Pas encore possible : une valeur de BC qui varie le long d'un bord.

## Les résultats : VTK `.vtu` + `.pvd`

Un `.vtu` par incrément, plus un `.pvd` qui donne le temps de chacun (un calcul linéaire écrit un
seul `.vtu` et un `.pvd`). Format ASCII, lisible par ParaView. Grandeurs aux Nodes :

- statique : `displacement`, `reaction`, `stress` (xx, yy, zz, xy, yz, xz, lissées), `von_mises`,
  et `plastic_strain` pour un calcul plastique ;
- thermique : `temperature`, `reaction` (puissance aux Nodes imposés).

## Régénérer les maillages

Les `.msh` 2D d'ici viennent des `.dxf` voisins, par `maille mesh` (le cube, 3D, est écrit à la main) :

```sh
cargo test -p maille-cli --test examples -- --ignored
```
