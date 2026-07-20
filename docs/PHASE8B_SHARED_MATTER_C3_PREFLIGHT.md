# Phase 8B - preflight de composition C3 a marcheur partage

## Objet

La bibliotheque locale a matiere partagee derive statiquement les trois liens
`X/Y/Z_rail tensor X_walker`. Avant de proposer une boucle C4, ce preflight
demande si trois de ces liens peuvent etre assembles dans le plus petit anneau
ferme explicite.

Le bloc a charge totale cinq contient trois qubits de rail, trois sites de
marcheur charge-deux `d0,d1,d2`, et un reservoir de paire distinct pour chaque
arete orientee. Le code bas a le marcheur sur `d0`; `d1,d2` sont virtuels avec
un detuning `0.5`. Les canaux de conversion incident a un meme site de
marcheur sont declares distincts car ils adressent des reservoirs distincts.

## Resultat

Le Schur exact sur le secteur a 1488 etats produit bien le mot ferme `XXX`.
Il n'est cependant pas selectif : le parasite principal est `XIX` et il est
plus grand dans tout le regime profond teste.

| g/Delta | `c_XXX` | plus grand parasite | rapport parasite/cible |
|---:|---:|---:|---:|
| 0.10 | `-4.24e-4` | `XIX = -1.67e-3` | 3.93 |
| 0.075 | `-4.74e-5` | `XIX = -3.95e-4` | 8.33 |
| 0.05 | `-3.14e-6` | `XIX = -6.98e-5` | 22.2 |

Le fit profond est `XXX ~ (g/Delta)^7.05`, contre
`parasite ~ (g/Delta)^4.56`. Aller plus profond aggrave donc la selectivite,
au lieu de la proteger.

Le scan independant de detuning du marcheur `0.5,1,2,5,10` a
`g/Delta=0.10,0.05` ne trouve aucune ligne avec
`parasite/cible < 0.25`; la meilleure vaut encore `3.93`.

## Verdict

La composition directe C3 des reservoirs/canaux actuels est rejetee. Le
resultat est utile : la bibliotheque X/Y/Z locale reste validee, mais le
premier mecanisme de mise en commun introduit un terme a deux rails d'ordre
inferieur qui noie le retour ferme multi-corps. Il est donc interdit de
promouvoir ces liens a un marcheur C4, a un patch ou a une tresse dans cette
grammaire directe.

La suite exige un mecanisme derive qui supprime `XIX` (contre-terme physique,
echo de plus haut ordre, encodage de marcheur different ou selection rule),
puis la repetition des memes audits de selectivite avant toute geometrie de
defaut.

## Claim boundary

Ce no-go est limite a l'anneau C3 explicite, au contrat de canaux et a la
fenetre de detunings testes. Il ne refute ni les liens X/Y/Z a un seul lien,
ni un contre-terme effectivement derive, ni un encodage de marcheur different,
ni toute extension ANTLER multi-liens.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_shared_matter_c3_walker_preflight.py
python experiments/phase7/run_phase8b_shared_matter_c3_detuning_scan.py
```
