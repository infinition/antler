# Phase 8B  -  préflight d’extrémité de twist

## Positionnement

Un twist `e <-> m` est l’extrémité d’une coupure de branche qui échange les
charges électrique et magnétique. Dans la construction de Bombin, la coupure
modifie la géométrie/stabilisateurs et son extrémité porte un stabilisateur
pentagonal non-CSS. Les twists peuvent reproduire la fusion et le braid des
anyons d’Ising, mais les seuls braids Ising restent de type Clifford, donc non
universels. [Bombin (2010)](https://arxiv.org/abs/1004.1838)

Cette littérature motive la cible. Elle ne rend pas l’implémentation ANTLER
automatique.

## Pentagon local dérivé

Un marcheur neutre à cinq positions est couplé successivement à cinq liens par
`X Z X Z X`. Le Schur complément exact (`160` états) donne le mot
non-CSS `XZXZX` :

| quantité profonde (`g/Delta <= 0.075`) | résultat |
|---|---:|
| puissance du pentagone | `5.0136` |
| coefficient non scalaire hors `{XZXZX}` | `< 4.2e-17` |

La boucle à cinq pas est importante : elle évite la règle de sélection paire
qui avait réfuté le gadget médiateur simple de Phase 8B.

## Premier chevauchement de coupure

Le pentagone `XZXZXII` est ensuite joint à une plaquette voisine
`ZXIIIZX`. Ils partagent deux liens sur lesquels les Pauli sont échangés, donc
commutent. Le bloc exact a `2 560` états.

| quantité profonde | puissance |
|---|---:|
| pentagone | `5.0502` |
| plaquette | `4.0365` |
| produit commutant | `9.0613` |
| mot hors `{I,pentagone,plaquette,produit}` | `< 7.7e-17` |

Le commutateur des stabilisateurs cibles est nul. C’est un test de
compatibilité locale, pas encore une coupure de branche complète.

## Ce qui manque avant de parler de twist

1. Construire un patch 2D complet avec **deux** extrémités pentagonales et
   démontrer que tous ses stabilisateurs microscopiquement dérivés commutent.
2. Démontrer le gap et l’indistinguabilité locale du sous-espace de fusion.
3. Construire des opérateurs de chaînes qui échangent `e` et `m` à travers la
   coupure et vérifier l’algèbre de Majorana/fusion sans l’insérer à la main.
4. Seulement alors, déplacer deux défauts par déformations locales et mesurer
   deux échanges adjacents non commutatifs.

Les approches Floquet à twists connues emploient précisément une dynamique
non commutante de déformations/mesures ; la simple modulation de notre famille
de stabilisateurs commutants est insuffisante. [Ellison, Sullivan & Dua
(2023)](https://arxiv.org/abs/2306.08027)

## Claim boundary

Le mécanisme de pentagone et son premier chevauchement sont des contrôles
locaux exacts de nouvelles primitives de marcheur phase-contrôlées. Ils ne
constituent ni une géométrie de twist complète, ni un espace de fusion, ni une
tresse, ni une réalisation ANTLER native. Ils ne démontrent donc pas encore de
statistique non abélienne.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_mixed_pentagon_walker_audit.py
python experiments/phase7/run_phase8b_pentagon_plaquette_overlap_audit.py
```
