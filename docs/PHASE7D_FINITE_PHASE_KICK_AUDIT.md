# Phase 7D  -  audit des kicks de phase à durée finie

## Question

L'écho alterné idéal utilisait le kick instantané

\[
Q=e^{-i\pi[\sum_j j(n_{aj}+n_{bj})+N_{\rm med}]}.
\]

Les potentiels de site sont une famille de contrôle déjà présente dans les
protocoles digitaux ANTLER. Ce test demande cependant si un pulse de potentiel
**fini**, avec conversion médiateur et hopping laissés actifs, suffit à
approcher ce kick.

La simulation utilise `t_leg=1`, crosstalk inactif `0.01g` et 16 sous-cycles
d'écho. Pendant chaque kick de durée `pi/kappa`, le Hamiltonien est

\[
H_{\rm kick}=H_{\rm pulse}+\kappa K,
\qquad K=\sum_j j(n_{aj}+n_{bj})+N_{\rm med}.
\]

Conserver `H_pulse` actif est volontairement conservateur : il évite de
supposer que les couplages peuvent être coupés instantanément.

## Résultat : non qualifié dans la fenêtre testée

Pour `kappa=100,200,400,800`, les valeurs maximales de potentiel de rail sont
respectivement `300,600,1200,2400`. Aucune ne reproduit le kick idéal :

| `kappa` | fuite | résidu de parité | distance logique |
| ---: | ---: | ---: | ---: |
| 100 | `1.66e-2` | `6.87e-2` | `7.37e-1` |
| 200 | `4.44e-2` | `8.91e-3` | `3.06e-1` |
| 400 | `2.03e-2` | `1.55e-3` | `1.50e-1` |
| 800 | `4.59e-2` | `8.29e-4` | `8.96e-2` |

La non-monotonie de la fuite est une interférence cohérente entre le kick et
la dynamique active, non un bruit statistique. Les kicks deviennent courts en
`1/kappa`, mais leur très grand nombre dans une séquence à 16 sous-cycles rend
la perturbation cumulée encore dominante dans cette fenêtre.

## Décision

L'écho à kick instantané reste un contrôle mathématique exact, mais son
implémentation « potentiel fort avec tous les couplages encore actifs » est
**rejetée dans la fenêtre enregistrée**. On ne peut pas présenter les kicks
idéaux comme une primitive matérielle ANTLER acquise.

Les alternatives à dériver avant toute phase étendue sont :

1. une fenêtre de contrôle qui éteint/refocalise les conversions et le hopping
   pendant le kick ;
2. une séquence composite qui réduit fortement le nombre de kicks ;
3. une estimation matérielle de bande passante et de rapport
   `kappa/(U,Delta,g)` justifiant un régime de kick rapide.

Une première variante continue de l'alternative 2 est documentée dans
`PHASE7D_CONTINUOUS_CDT_REFOCUSING_AUDIT.md`. Elle préserve exactement la
conversion de paire et supprime partiellement le hopping, mais ne franchit pas
encore le seuil de porte dans la fenêtre de fréquence testée.

Ce résultat n'affecte pas les no-go ou les contrôles antérieurs ; il borne
seulement la disponibilité physique de la correction d'écho. Aucun claim de
phase topologique, code 2D, tresse, non-abélianité, universalité ou tolérance
aux fautes n'est établi.

Résultat machine : `results/phase7/finite_phase_kick_audit.json`.
