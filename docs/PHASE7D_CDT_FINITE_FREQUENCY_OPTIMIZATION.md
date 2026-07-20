# Phase 7D  -  optimisation d'amplitude CDT à fréquence finie

## Question préenregistrée

À fréquence finie, les corrections de Magnus peuvent déplacer le meilleur
rapport `xi=A/omega` par rapport au premier zéro asymptotique de Bessel,
`xi0=2.4048255577`. Le scan ne varie donc qu'une seule variable, dans la
fenêtre `xi=2.10..2.75`, à quatre cycles par pulse Rabi et avec la même
dynamique 472D, `t_leg=1` et crosstalk `0.01g`.

## Résultat

Le zéro de Bessel reste le meilleur point suivant le critère préenregistré
« fuite puis parité » :

- `xi=2.40482556` : fuite `4.439e-4`, parité `3.769e-4` ;
- voisin gauche `xi=2.35` : fuite `4.682e-4` ;
- voisin droit `xi=2.46` : fuite `4.514e-4` et parité légèrement plus basse,
  `3.680e-4` ;
- aux bords `xi=2.10` et `2.75` : fuite `1.031e-3` et `9.177e-4`.

Il n'y a donc pas de point finement accordé qui transforme le protocole
quatre-cycles en porte à haute fidélité. L'optimum de fuite est stable dans la
fenêtre et le faible gain de parité à `xi=2.46` ne compense pas la hausse de
fuite.

## Décision

La correction d'amplitude seule est épuisée dans cette fenêtre. La suite ne
doit pas être un scan plus large : elle nécessite soit une construction
Floquet qui annule les corrections d'ordre supérieur, soit une bande passante
physiquement justifiée permettant d'augmenter la fréquence, soit un nouveau
mécanisme de contrôle. Cette décision ne change pas les limites de phase :
aucun qubit topologique, gap protecteur, code 2D ou braid n'est établi.

Résultat machine : `results/phase7/cdt_finite_frequency_optimization.json`.
