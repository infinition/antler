# Phase 8B - fermeture d'un patch fini de controles mixtes

## Resultat

Les briques de marcheur derivees localement (`X`, `Y`, `Z`) se completent en
un patch stabilisateur fini a sept liens. Ce patch est une reference algebrique
compacte, pas encore une surface de dislocation complete ni une realisation du
ladder ANTLER.

Les six controles, tous de poids quatre, sont :

```text
YXZXIII    XZIIXZI    XZXYIII
YXIIYXI    XIXIXIX    YIZIYIY
```

Les deux premiers sont le couple mixte deja controle par Schur exact. Les
quatre autres sont des mots de meme longueur dans le vocabulaire de marcheur
phase-controle `X/Y/Z`.

## Gate algebrique du code

| propriete | valeur |
|---|---:|
| nombre de liens | `7` |
| rang stabilisateur | `6` |
| qubits encodes | `1` |
| degenerescence du fondamental | `2` |
| distance minimale | `3` |
| gap de syndrome de `-J sum S_a` | `2J` |
| Paulis de poids `<3` non scalaires sur le code | `0 / 210` |

Les representants logiques minimaux trouves sont `XZIIIIX` et `ZYIIIIZ`; ils
anticommutent. Le resultat est deduit du rang et du centralisateur exacts.

## Fermeture microscopique par paires

Les 15 paires de controles sont chacune downfoldees dans le Hamiltonien de
deux marcheurs explicites a `g/Delta=0.05`. Chaque bloc a 2 048 etats, dont
128 etats de jauge bas. Toutes les paires se ferment sur
`{I,S_i,S_j,S_iS_j}` :

| quantite pire cas | valeur |
|---|---:|
| coefficient hors algebre generee | `1.66e-16` |
| residu d'hermiticite | `9.20e-19` |
| plus petit coefficient de controle | `1.268997e-4` |
| plus grand produit de deux controles | `9.694863e-9` |

Les produits de deux controles sont donc plus petits que les controles
principaux d'environ quatre ordres de grandeur au point enregistre.

## Lemme de fermeture multi-marcheur

Chaque marcheur vit sur un cycle `C4`. Le support modulo deux de tout chemin
qui revient au vide est un cycle de `C4`; il n'existe que deux possibilites :

```text
0000  (identite)       1111  (tour complet)
```

Un tour complet produit le mot `S_i` du marcheur. Les operateurs de Pauli ne
peuvent differer par reordonnancement que par une phase scalaire. Dans la
serie de Schrieffer--Wolff formelle du Hamiltonien a six marcheurs, tout chemin
bas-vers-bas est donc un produit de sous-ensemble des six `S_i`. Comme les
`S_i` commutent, l'algebre effective est le groupe stabilisateur de taille
`2^6=64`, et agit scalairement sur le qubit encode.

Ce n'est pas une extrapolation numerique : le code enumere independamment les
parites de tous les chemins retournes jusqu'a l'ordre 12, et la preuve utilise
la dimension un du cycle-espace de `C4`.

Au niveau tronque aux ordres quatre et huit, la somme des 15 produits est
`7.60e-8`; elle laisse une borne de gap de syndrome positive `2.53647e-4`
dans les unites du Hamiltonien effectif. Les ordres superieurs ne sont pas
bornes numeriquement par cette ligne.

## Position scientifique

Le patch montre que le probleme ne se bloque plus sur la selectivite locale ni
sur le crosstalk ideal des marcheurs : dans cette grammaire declaree, le
crosstalk reste dans l'algebre des stabilisateurs et ne produit pas un Pauli
logique nouveau. C'est une condition necessaire importante pour une memoire
codee.

Il ne prouve pas une topologie intrinseque thermodynamique ni une statistique
non abelienne. Un seul qubit code fini n'est pas un espace de fusion de twists.
La prochaine construction devra ajouter une vraie geometrie de coupure,
plusieurs defauts et des deformations locales derivees avant tout test de braid.

## Claim boundary

Le lemme tous-ordres vaut dans la serie perturbative de la grammaire ideale de
marcheurs `C4` independants; il ne constitue ni diagonalisation exacte du bloc
a six marcheurs, ni preuve de convergence a couplage fini. Les marcheurs
neutres, la phase conditionnelle `pi/2`, le patch complet et son embedding
`U(1)` n'ont pas ete derives du Hamiltonien ANTLER fige. Aucun twist complet,
fusion, mouvement de defaut, braid non abelien, universalite ou seuil de faute
n'est revendique.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_mixed_check_patch_algebra.py
python experiments/phase7/run_phase8b_mixed_check_patch_pairwise_overlap.py
python experiments/phase7/run_phase8b_mixed_patch_closed_walk_closure.py
```
