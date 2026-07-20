# Phase 4.7 : fermeture numérique de la quantification digitale

Les cinq profondeurs `D=4,6,8,10,12` sont maintenant toutes évaluées au pas
fin `dt=0.125`, avec la durée adiabatique `T(D) proportional D^2`. Le fit

`odd_slope(D) = -1 + a D^-p`

donne `a=0.66493`, `p=2.28624 ± 0.00764`, `R²=0.999983` et une RMSE de
`3.91e-5`. Les cinq points ont leakage `<5e-5`, `sigma_min>0.999975` et un
gap d'isolation de handoff positif.

Les contrôles indépendants en pas de temps à `D=6,8`, les déformations de
chemin et la composition `n=1,2,4,8` sont déjà archivés. La clôture qualifie
la **quantification numérique de la primitive Z abélienne** dans le modèle
gelé ; elle ne qualifie ni une tresse non abélienne, ni l'universalité, ni la
tolérance aux fautes.

Script d'assemblage : `experiments/phase4_7/run_phase47_publication_closure.py`.
Résultat : `results/phase4_7/publication_closure.json`.
