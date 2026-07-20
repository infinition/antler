# Phase 6  -  statut de la calibration externe après les contrôles 6L–6M

## Verdict

Le benchmark externe d'Iemini reste utile pour calibrer les diagnostics, mais
il ne franchit pas les portes nécessaires à un qubit topologique protégé ni à
une tresse physique dans cette archive. Aucun calcul `L=12` d'opérateur de
bord n'est lancé à partir de cette série : il ne constituerait pas une
nouvelle preuve de protection.

Ce verdict concerne le parent externe à `lambda=1`. Il ne réfute pas une
future construction native ANTLER conçue avec une limite soluble et un état
parent distinct.

## 1. Contrôle positif du harnais d'opérateur de bord

Le critère utilisé est

\[
\epsilon_{\rm edge} = \frac{\|(1-P)[H,O]P\|_F}{\|OP\|_F}.
\]

Sur une chaîne de Kitaev ouverte au point exact `mu=0, t=Delta`, le Majorana
gauche `a_0` donne exactement `epsilon_edge = 0`. Le Majorana de bulk `a_1`
donne `epsilon_edge = 2.0`. Le harnais distingue donc effectivement un
zéro-mode exact d'un opérateur local non protégé.

Cette calibration ne rend pas le parent Iemini équivalent à la chaîne de
Kitaev : elle valide uniquement l'implémentation numérique du test projeté.

## 2. Axes taille et support : correction d'interprétation

Le parcours maximal disponible changeait simultanément la taille et le
support : `(L,j)=(6,2),(8,3),(10,4)`. La baisse du résidu normalisé
`7.22 -> 6.32 -> 5.37` est donc principalement une amélioration de support,
pas un scaling en taille.

Au support fixé `j=3`, le résidu vaut `6.3184` à `L=8` et `6.3057` à `L=10`.
Il est donc essentiellement plat et reste très loin de zéro, alors même que
la queue analytique de troncature décroît géométriquement. Le générateur
tronqué testé n'est pas qualifié comme mode de bord quasi-conservé, ni comme
primitive de braid.

## 3. Gap neutre : nouveau point L=10

Une action matricielle compilée, recoupée avec l'ED creuse à `L=6` avec une
erreur absolue `1.33e-14`, donne le premier gap neutre à charge fixée :

| L | gap neutre |
|---:|---:|
| 4 | 1.4621251 |
| 6 | 0.7174655 |
| 8 | 0.4181120 |
| 10 | 0.2719927 |

La suite décroît strictement. Un ajustement descriptif pur en loi de puissance
donne `gap ~ 18.82 L^-1.834` (`R²=0.9995` en espace logarithmique), mais quatre
tailles ne démontrent ni un exposant thermodynamique ni une fermeture rigoureuse.
Elles ne démontrent surtout pas de saturation vers un gap neutre non nul.

Le parent est de plus compressible dans les secteurs de charge : les énergies
d'ajout/retrait sont nulles à la précision numérique aux tailles auditées. Une
charge préparée/fixée ou une énergie de charging externe serait donc une
ressource expérimentale supplémentaire, pas une propriété à inférer du
doublet de charge fixe.

## Conséquence méthodologique

Le prochain travail natif ne doit pas être un scan de paramètres ni une
extension de support présentée comme une tresse. Il doit commencer par une
construction état-vers-Hamiltonien avec :

- les symétries de parité et la charge pondérée exactes ;
- une contrainte locale ou un invariant analytique ;
- un opérateur de bord candidat dont `epsilon_edge(L,j)` peut être étudié à
  support fixé puis à support croissant ;
- des gaps neutre et de charge audités séparément avant toute dynamique.

Une éventuelle reprise du benchmark externe devra d'abord montrer, sur les
deux axes indépendants, un `epsilon_edge` qui décroît vers zéro et un gap
neutre dont la saturation est établie. Sans ces prérequis, ni Yang--Baxter ni
une évolution adiabatique ne seront promus comme preuve de calcul topologique.

## Sources reproductibles

- `experiments/phase6/run_phase6l_kitaev_edge_harness_control.py` ;
- `results/phase6/kitaev_edge_harness_control.json` ;
- `experiments/phase6/run_phase6m_neutral_gap_scaling.py` ;
- `results/phase6/neutral_gap_scaling.json` ;
- `experiments/phase6/run_phase6k_edge_support_scaling.py` ;
- `results/phase6/edge_support_scaling.json`.
