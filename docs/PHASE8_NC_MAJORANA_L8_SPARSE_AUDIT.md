# Phase 8  -  candidat Floquet à nombre conservé, audit exact jusqu'à L=8

## Résultat

Un candidat externe à densité fixe a franchi le premier test de taille non
trivial. Le Hamiltonien est le modèle Floquet effectif à nombre conservé de
Defossez *et al.* :

\[
H_{\rm eff}=\alpha H_0+(1-\alpha)P^\dagger H_0P,
\qquad P=e^{-i\pi J_x/2}.
\]

Le point `U0=-2`, `alpha=0.5` était déjà préenregistré : c'était le meilleur
ratio split/gap du benchmark Phase 7D à `L=6,N=4`. Il est maintenant réévalué
avec charge paire et densité constante `N/(2L)=1/4`, aux tailles `L=4,N=2` et
`L=8,N=4`.

| taille | dimension Fock fixée-N | split de parité | gap neutre min. | split/gap |
| ---: | ---: | ---: | ---: | ---: |
| `L=4, N=2` | 28 | `5.87e-2` | `5.59e-1` | `1.049e-1` |
| `L=8, N=4` | 1820 | `9.10e-3` | `3.02e-1` | `3.012e-2` |

La représentation sparse du cycle est validée avant l'audit : à `L=6,N=4`,
elle coïncide avec l'implémentation dense antérieure à une norme de Frobenius
`1.62e-14`; l'unitarité de la rotation et le commutateur de parité valent
respectivement `1.37e-14` et `2.85e-15` à `L=8`.

Le spectre de Schmidt de l'état du secteur de parité opposé est pairé, à la
précision numérique, sur les niveaux dominants. Le secteur fondamental montre
des multiplets mais pas une double dégénérescence uniforme ; cette structure
est donc rapportée comme un signal qualitatif, pas comme un invariant certifié.

## Contrôles négatifs enregistrés à L=8,N=4

| contrôle | gap neutre min. | split/gap | conclusion |
| --- | ---: | ---: | --- |
| retirer le mélange Floquet (`alpha=1`) | `8.88e-15` | divergent | gapless dans ce secteur fini |
| retirer l'interaction (`U0=0`) | `3.82e-14` | divergent | gapless dans ce secteur fini |

Le gap du candidat exige donc les deux ingrédients : interaction et mélange
Floquet. Cela réfute l'explication triviale « simple effet de taille du
hopping libre » dans les deux contrôles testés.

## Interprétation exacte

Il s'agit d'une **cible externe finie encourageante** : le ratio split/gap se
réduit fortement de `L=4` à `L=8` sans fermeture observée du gap, et les
contrôles correspondants sont gapless. Deux tailles ne permettent ni un fit
asymptotique crédible ni la démonstration d'une phase. La prochaine porte est
une réplication à densité fixe et taille longue par MPS/DMRG à conservation
`U(1)`, suivie seulement ensuite de la dérivation des ressources `U0` et `P`
depuis le Hamiltonien microscopique ANTLER.

Ce résultat ne démontre pas encore une phase ANTLER native, une dégénérescence
de bord protégée, un braid, une non-abélianité, l'universalité ou la tolérance
aux fautes.

Résultat machine : `results/phase7/nc_majorana_l8_sparse_audit.json`.
