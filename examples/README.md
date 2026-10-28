# Exemples et guide des fichiers

Maille résout un modèle en trois fichiers : **maillage → modèle → résultats**. En 2D, le maillage
peut venir d'un dessin DXF 2D ; en 3D, d'une surface STL ou d'une pièce STEP (`maille mesh`).

```
dessin (.dxf) ──┐
surface (.stl) ─┤
pièce (.step) ──┴ maille mesh ──▶
maillage (.msh) ──┐
                  ├─ maille solve ──▶ résultats (.pvd + .vtu) ──▶ maille view / ParaView
modèle (.toml) ───┘                                  │
        ▲                                            │
        └──── [temperatures] (un calcul lit les résultats thermiques d'un autre)
```

Depuis ce dossier :

```sh
maille solve traction.toml        # écrit traction.pvd et traction.vtu
maille view traction.pvd          # ouvre le viewer (2D ou 3D) ; dernier incrément par défaut
maille view trempe.pvd 50         # un incrément précis, depuis 0
```

| Exemple | Ce qu'il montre | Solution exacte |
|---|---|---|
| `traction` | plaque en contraintes planes, Tri6 | champ linéaire |
| `trempe` | thermique transitoire, 200 incréments | série de Fourier |
| `tube_thermique` puis `tube_contraintes` | couplage thermo-mécanique (lancer dans cet ordre) | Timoshenko et Goodier |
| `tube_pression` | pression sur un bord, déformations planes | Lamé |
| `cube` | 3D, Tet4 (viewer 3D ou ParaView) | champ linéaire |
| `tube3d` | 3D depuis un STL : `maille mesh tube3d.stl --size 2 --quadratic` d'abord, Tet10 | Lamé |
| `sphere` | 3D depuis un STEP (faces sphériques exactes) : `maille mesh sphere.step --size 2 --quadratic` d'abord, Tet10 | Lamé (sphère) |

**Viewer 3D.** Un maillage de tétraèdres s'affiche par sa peau, colorée par le champ (touches 1 à 8 :
von Mises, σxx, σyy, σzz, τxy, τyz, τxz, |u| ; 9 pour p en plasticité). Glisser tourne, Maj + glisser (ou
clic droit) déplace, la molette zoome, `D` montre la déformée, `R` recadre. **Plan de coupe** : `X`, `Y`
ou `Z` coupent selon l'axe (la même touche retire la coupe), `↑` et `↓` déplacent le plan, `F` change
de côté ; la face coupée montre le champ à l'intérieur (interpolé sur les sommets des tétraèdres).
Les mêmes réglages donnent une image sans ouvrir de fenêtre :

```sh
maille view cube.pvd --png cube.png --field 8 --deformed --cut x=0.5 --angles 35,25
```

(`--size LxH`, `--angles lacet,tangage` en degrés, `--field k` pour la touche k, `--cut axe=position` avec
la position de 0 à 1 dans la boîte du maillage.) Il faut une carte graphique, même sans fenêtre.

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

## La surface : STL 3D, ASCII ou binaire

```sh
maille mesh cube.stl --quadratic             # écrit cube.msh (Tet10 ; sans --quadratic, Tet4)
maille mesh tube3d.stl --size 2 --quadratic  # tétraèdres de taille h = 2 mm
```

Le STL doit être une surface **fermée d'une seule pièce** (pas de trou, pas de cavité, chaque arête
partagée par deux triangles). Maille soude les sommets, vérifie l'orientation et corrige des normales
rentrantes. Le **nom du bord** est celui du `solid` (STL ASCII à plusieurs `solid`, comme
`cube.stl`) ; sinon la surface est découpée aux arêtes vives (`--feature-angle`, 30° par défaut) en
bords `surface_1`, `surface_2`… Le Set `domaine` contient tous les tétraèdres.

**Taille h (`--size`).** Maille regroupe les triangles du STL en faces planes, remaille chacune à la
taille h (triangles de rayon circonscrit au plus h/√3, comme en 2D), puis remplit l'intérieur de
tétraèdres de taille au plus h (rayon circonscrit au plus h √6/4, rapport rayon / plus petite arête
visé à 2), plus fins près des petits détails (gradation 0,3). Limites : une surface courbe d'un STL est
faite de bandes planes, et une bande plus étroite que h garde sa largeur (la finesse est alors celle
du facettage) ; près de la surface, la qualité est visée mais pas garantie, et des tétraèdres plats
(slivers) peuvent rester ; un coin ou un dièdre de moins de 6° environ est refusé, avec un message.

**Sans `--size`**, la finesse est celle des triangles du STL (Maille ajoute des points sur la
surface, mais n'affine pas) : un STL grossier donne un maillage grossier, avec des tétraèdres parfois
aplatis, et des triangles très effilés (un disque en éventail de centaines de triangles) sont refusés.
`--size` lève ce refus, puisque la surface est remaillée.

## La pièce : STEP 3D

```sh
maille mesh sphere.step --size 2 --quadratic  # tétraèdres de taille h = 2 mm, en Tet10
```

Un fichier STEP texte (AP203, AP214, AP242) d'**un seul solide** : plusieurs solides (un assemblage),
une surface ouverte ou un solide à cavités sont refusés, avec un message. Maille lit la vraie forme
de la pièce : faces planes, cylindriques, coniques, sphériques et toriques ; arêtes en segments,
arcs de cercle ou d'ellipse et B-splines. Une autre surface (B-spline, extrusion…) est refusée, avec
sa ligne dans le fichier. `--size` est obligatoire. Il est dans l'**unité de longueur du fichier**,
que `maille mesh` affiche (`mm`, `m`, `inch`…) : Maille ne convertit rien, le modèle TOML doit être
dans la même unité. Le **nom de chaque face** dans le fichier est le nom du bord (`face_1`,
`face_2`… si elle n'en a pas), et le Set `domaine` contient tous les tétraèdres.

Chaque face est maillée sur sa surface exacte, en triangles de taille au plus h, puis l'intérieur
est rempli comme pour un STL (`--size`). Sur une face ou une arête courbe, la taille est en plus
bornée par un seizième de tour : un trou de rayon r a au moins 16 morceaux, quelle que soit h. Les
sommets sont sur la vraie surface, les triangles en dessous : l'erreur de géométrie décroît comme
h². Limites : une pointe où deux faces se rejoignent tangentes (le chanfrein d'un écrou contre un de
ses pans) peut refuser les petites tailles h, et une taille h grossière devant les détails de la
pièce peut faire échouer le remplissage (le message propose une autre taille h). `sphere.step` est
écrit à la main : un huitième de sphère creuse de rayons 10 et 20 mm, faces `interieur`,
`exterieur`, `sym_x`, `sym_y` et `sym_z`.

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

Les `.msh` 2D d'ici viennent des `.dxf` voisins, par `maille mesh` (le cube, 3D, est écrit à la main ; `cube.stl` en est la surface, pour `maille mesh`). `tube3d.msh`
n'est pas livré : `maille mesh tube3d.stl --size 2 --quadratic` l'écrit en une seconde, de même que
`sphere.msh` (`maille mesh sphere.step --size 2 --quadratic`) :

```sh
cargo test -p maille-cli --test examples -- --ignored
```
