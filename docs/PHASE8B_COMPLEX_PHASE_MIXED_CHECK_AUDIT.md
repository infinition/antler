# Phase 8B - lien Y et controle mixte central

## Resultat

Le vocabulaire de marcheurs Phase 8B peut produire un lien de Pauli `Y` sans
ajouter `Y` comme terme effectif postule. Dans la base `Z` des liens de jauge,
le saut hermitien vers l'etat de marcheur suivant est choisi comme

```
<flipped gauge | H | gauge> = i g z(gauge).
```

Il applique d'abord `Z` a l'amplitude-source puis retourne le bit, donc
`i X Z = Y` exactement. La phase coherente `pi/2` est l'ingredient de
controle supplementaire; elle n'est pas une matrice logique ni une matrice de
tresse inseree.

Cette extension est necessaire pour les geometries de twist usuelles : les
plaquettes le long de la coupure sont mixtes et le controle central contient
une action `Y` dans la representation triangulaire. Voir
[Yoder & Kim (2017)](https://arxiv.org/abs/1612.04795), en particulier la
description des plaquettes mixtes et du `Y` central, et
[Bombin (2010)](https://arxiv.org/abs/1004.1838) pour le role des extremites
de coupure e/m.

## Audit du lien Y isole

Un marcheur neutre a quatre positions execute explicitement `Y X Z X` sur
quatre liens. Le Schur complement exact a energie nulle agit sur 16 etats de
jauge, dans un Hamiltonien total de dimension 64.

| quantite, regime profond `g/Delta <= 0.075` | resultat |
|---|---:|
| puissance de `YXZX` | `4.0091` |
| coefficient non scalaire hors cible | `< 2.8e-17` |
| residu d'hermiticite | `0` |

Le signe du coefficient depend de l'orientation de la boucle et n'est pas
une phase logique. Le fait utile est l'isolement tomographique du mot `YXZX`.

## Premier chevauchement de controle mixte

Le controle central `YXZXII` est joint a la plaquette `XZIIXZ`. Ils partagent
deux liens : `Y`/`X` et `X`/`Z` anticommutent chacun, donc le commutateur total
est nul. Les deux marcheurs explicites donnent un bloc de dimension 1 024,
avec 64 etats de jauge bas.

| quantite profonde | puissance |
|---|---:|
| controle mixte `YXZXII` | `4.0365` |
| plaquette `XZIIXZ` | `4.0365` |
| produit commutant | `8.0514` |
| mot hors `{I,mixte,plaquette,produit}` | `< 1.6e-16` |
| `||[mixte, plaquette]||_F` | `0` |

Le resultat teste les crosstalks du bloc joint : le produit d'ordre huit est
le seul terme croise resolu. Ce n'est pas une hypothese de stabilisateurs
independants.

## Ce que cela debloque, et ce que cela ne prouve pas

La grammaire locale sait maintenant compiler des mots sur liens `X`, `Y` et
`Z` a partir de boucles de marcheur phase-controlees. Elle peut donc exprimer
les controles mixtes requis par une geometrie complete de coupure de branche.

Il reste, dans cet ordre :

1. generer un patch de dislocation complet avec deux extremites et verifier
   **tous** les commutateurs de chevauchement ;
2. calculer rang stabilisateur, dimension de fusion, distance et
   indistinguabilite locale ;
3. deriver les deformations locales qui deplacent les extremites ;
4. seulement si les deux holonomies adjacentes ont un commutateur non nul,
   tester Yang--Baxter.

## Claim boundary

Les deux resultats sont des controles locaux exacts d'une ressource de
marcheur neutre phase-controlee. Le marcheur, sa phase conditionnelle et son
integration dans le Hamiltonien ANTLER fige restent a deriver. Aucun patch de
twist complet, espace de fusion, mouvement de defaut, braid non abelien,
universalite ou seuil de faute n'est etabli ici.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_y_link_walker_audit.py
python experiments/phase7/run_phase8b_y_mixed_plaquette_overlap_audit.py
```
