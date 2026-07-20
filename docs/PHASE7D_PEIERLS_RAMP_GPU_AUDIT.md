# Phase 7D  -  audit GPU des rampes finies de phase de Peierls

## Méthode

L'audit CPU de rampes n'avait pas pu sérialiser ses métriques avant la limite
locale. Le même Hamiltonien 472D est ici propagé sur la RTX 4070 Ti par
décomposition spectrale exacte de chaque segment hermitien, en `complex128`.
Les 16 colonnes du frame monomère sont propagées simultanément. Aucun
Hamiltonien effectif supplémentaire n'est introduit.

Le contrôle garde `t_leg=1`, crosstalk `0.01g` et 16 sous-cycles. Toutes les
jambes passent linéairement de phase Peierls `0` à `pi`, puis reviennent à
`0`, tandis que les conversions charge-2 restent actives.

## Résultats et convergence de segmentation

Le switch instantané donne une fuite de `6.97e-8`, un résidu de parité de
`1.72e-5`, une distance logique de `3.12e-5` et une valeur singulière minimale
de `0.9999999652`. Les rampes linéaires sont ensuite réévaluées avec 2, 4 et 8
segments de phase de même durée totale :

| fraction de sous-cycle par rampe | segments | fuite | résidu de parité | distance logique | valeur singulière min. |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.005 | 2 | `1.0628e-6` | `1.6078e-5` | `2.9222e-5` | `0.9999994686` |
| 0.005 | 4 | `9.1478e-7` | `1.6237e-5` | `2.9506e-5` | `0.9999995426` |
| 0.005 | 8 | `8.8194e-7` | `1.6272e-5` | `2.9569e-5` | `0.9999995590` |
| 0.020 | 2 | `1.6235e-5` | `1.1582e-5` | `1.2719e-5` | `0.9999918823` |
| 0.020 | 4 | `1.3865e-5` | `1.0058e-5` | `1.1953e-5` | `0.9999930673` |
| 0.020 | 8 | `1.3340e-5` | `9.8224e-6` | `1.1962e-5` | `0.9999933302` |

La rampe courte introduit une perte contrôlée plutôt qu'un effondrement de la
fermeture Rabi. C'est qualitativement différent des kicks de potentiel avec
Hamiltonien actif, qui produisaient des fuites de l'ordre de `1e-2`.

Les différences 2-vers-4 segments décroissent d'un facteur `4.508` aux deux
fractions. Les différences 4-vers-8 donnent un ordre observé de `2.173` pour
la fuite aux deux points, et de `2.168` / `2.695` pour la parité. La
représentation numérique par segments de la rampe est donc convergée dans cette
fenêtre. Le calcul emploie une propagation spectrale exacte de chaque segment
sur RTX 4070 Ti, en `complex128`.

## Limite de promotion

Cette convergence ne certifie que la discrétisation du profil de phase dans le
bloc microscopique exact 472D. Elle ne valide pas la bande passante physique,
la calibration, le jitter temporel, le bruit de phase, ni l'indépendance du
contrôle de phase et de la conversion charge-2. Ces erreurs, puis l'intégration
au ladder complet, restent des tests requis avant toute promotion matérielle.

## Décision

Le pont Peierls possède maintenant un contrôle numérique local convergé : des
rampes finies courtes n'écrasent pas la primitive dans le bloc exact. La preuve
matérielle reste ouverte ; ce résultat ne démontre ni phase topologique, ni
mode de bord, code 2D, braid, non-abélianité, universalité ou tolérance aux
fautes.

Résultats machine : `results/phase7/peierls_phase_ramp_gpu_audit.json`,
`results/phase7/peierls_phase_ramp_gpu_4step_refinement.json`,
`results/phase7/peierls_phase_ramp_gpu_8step_refinement.json` et
`results/phase7/peierls_phase_ramp_segmentation_convergence.json`.
