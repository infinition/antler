# Audit de tresse U(1) : benchmark Iemini

## Objet exact

Ce document audite les opérateurs de tresse à support fini du Hamiltonien à
deux fils de F. Iemini *et al.*, PRL 115, 156402 (2015), sur sa ligne
exactement soluble `lambda=1`. Il calibre l'infrastructure ANTLER sur une
référence publiée à conservation de nombre ; il ne dérive pas cette référence
du Hamiltonien ANTLER gelé et ne décrit pas un échange adiabatique physique.

La convention de Fock du projet est rung-major. Le passage depuis la convention
fil-par-fil de la référence ajoute la parité d'inversion
`sum_i n_b,i sum_(k>i) n_a,k`. Le frame analytique ainsi corrigé est annulé
par le parent Hamiltonien à `L=N=6` avec un résidu maximal
`1.72e-15`. C'est un garde-fou obligatoire avant toute extrapolation sans ED.

## Opérateurs et garde-fous

Les deux générateurs finis sont :

- `R_aR,bR=(I+Z_aR,bR)/sqrt(2)`, supporté sur les `j` derniers barreaux ;
- `R_aR,aL=(I+Z_aR,aL)/sqrt(2)`, bilocal sur les `j` premiers et derniers
  barreaux, conformément à la conservation de nombre.

Pour un remplissage `nu=1/2`, la troncature analytique publiée est
`nu^(2j)+(1-nu)^(2j)`. Le résidu Yang--Baxter brut n'est rapporté que parce
que la norme du commutateur projeté dépasse explicitement `1e-3` à chaque
point. La polarisation des matrices projetées n'est qu'un diagnostic : elle
ne remplace jamais l'opérateur brut.

## Résultats au support maximal calculé

| L | support j | queue analytique | fuite d'amplitude | défaut d'unitarité projeté | ||[R1,R2]|| | résidu YB brut |
|---|---:|---:|---:|---:|---:|---:|
| 8 | 3 | 3.125e-2 | 3.113e-1 | 1.223e-1 | 1.170 | 1.114e-1 |
| 10 | 4 | 7.8125e-3 | 2.232e-1 | 6.683e-2 | 1.281 | 6.364e-2 |

Les opérateurs sont bien supportés aux bords prescrits, le commutateur est
robustement non nul, et les défauts de troncature diminuent lorsque le support
augmente. C'est exactement le comportement attendu pour la construction
finie publiée.

La conclusion correcte est néanmoins **préliminaire** : à `L=10`, la fuite
d'amplitude `0.223` et le résidu YB brut `0.0636` ne permettent pas de dire
qu'une tresse finie est déjà exacte ou sans excitation hors code. Il faut un
extrapolation contrôlée `L,j -> infini` avec une fraction de support restant
loin de `L/2`, puis une vraie dynamique d'échange. L'algèbre est une cible
externe reconnue par la soufflerie, pas encore une porte ANTLER.

## Reproductibilité

- `experiments/phase5/run_phase5_iemini_braid_audit.py` : ED à `L=8` ;
- `experiments/phase5/run_phase5_iemini_braid_scaling.py` : frames exacts,
  `L=8,10`, et contrôle de convention ;
- `results/phase5/iemini_braid_audit.json` ;
- `results/phase5/iemini_braid_scaling.json`.
