# Phase 7D  -  refocalisation continue CDT

## Construction

Le générateur continu est `K_cont=sum_j j(n_aj+n_bj)+sum_(j,chi)(2j+1)n_d(j,chi)`.
Le poids d'un médiateur charge-2 sur le lien `j,j+1` vaut la somme des poids de
sa paire. Par conséquent, `[K_cont,H_pair]=0` pour les canaux same-rail et
opposite-rail, tandis que le hopping intrajambe change ce poids de un.

Un drive sinusoïdal à `xi=2.4048255577`, premier zéro de `J0`, annule donc le
hopping dans la moyenne haute fréquence sans annuler la conversion de paire.
Le frame tournant est fermé par un nombre entier de cycles par pulse Rabi.
Dans le bloc complet 472D, le commutateur de conversion de paire est exactement
nul pour les deux familles de canaux ; le commutateur de hopping a norme de
Frobenius `33.2265`.

## Audit dynamique

Le cas sévère reste `t_leg=1`, crosstalk `0.01g`. Sans refocalisation, la fuite
est `1.98e-2` et le résidu de parité `2.62e-2`. Après convergence interne de
l'intégration (`16` contre `32` pas par cycle), la modulation CDT donne :

| cycles par pulse Rabi | fuite | résidu de parité | distance logique |
| ---: | ---: | ---: | ---: |
| 1 | `1.12e-2` | `8.81e-3` | `8.85e-3` |
| 2 | `1.89e-3` | `1.59e-3` | `1.74e-3` |
| 4 | `4.35e-4` | `3.70e-4` | `4.14e-4` |

La suppression est réelle et convergente mais ne satisfait pas encore une
porte de haute fidélité. À quatre cycles, le potentiel maximal de rail est
déjà `2334.46` dans les unités du modèle ; une extrapolation exige une borne
matérielle de bande passante.

## Décision

Cette voie est algébriquement qualifiée et réduit fortement le problème des
kicks : la conversion de paire n'est pas réduite avec le hopping. Dans la
fenêtre dynamique testée, elle n'est pas qualifiée comme primitive de porte.
Il faut une analyse haute fréquence, une amplitude optimisée à fréquence finie
ou une ressource de modulation dont la bande passante est justifiée.

Cela ne démontre ni un Hamiltonien ANTLER natif complet, ni gap, mode de bord,
phase topologique, code 2D, braid, non-abélianité, universalité ou tolérance
aux fautes.

Le scan local d'amplitude à fréquence finie est clos dans
`PHASE7D_CDT_FINITE_FREQUENCY_OPTIMIZATION.md` : le premier zéro de Bessel
reste optimal pour la fuite dans la fenêtre `xi=2.10..2.75`; aucun réglage
voisin ne transforme le point à quatre cycles en porte.

Résultat machine : `results/phase7/continuous_cdt_refocusing_audit.json`.
