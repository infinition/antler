# Phase 4.7  -  état des campagnes de fermeture digitale

Ce document empêche de promouvoir des tendances intermédiaires en preuve de
quantification. Les scripts et les logs sont conservés ; seuls les points
terminés, sérialisés et indépendamment convergés pourront alimenter une
extrapolation finale.

## Livrables

- `experiments/phase4_7/run_phase47_deep_limit.py` : profondeur, durée et
  convergence en pas de temps ;
- `experiments/phase4_7/run_phase47_path_invariance.py` : déformations de
  rampes et de séquence à topologie fixe ;
- `experiments/phase4_7/run_phase47_composition.py` : cycles `n=1,2,4,8` ;
- `docs/DIGITAL_TRANSFER_LEMMA.md` : transfert à deux niveaux et convention
  de signe.

## État à la date de cette archive de travail

Les campagnes sont maintenant clôturées numériquement. La promotion reste
strictement limitée à la primitive logique abélienne du modèle gelé :

- le scan primaire à `dt=0.25` est terminé : les pentes progressent de
  `-0.9735` (`D=4`) à `-0.99799` (`D=12`) et son fit donne
  `p=2.484`, `R²=0.999775` ;
- pour `D=12`, le leakage primaire (`≈3.43e-2`) est incompatible avec les
  points moins profonds et doit être considéré comme potentiellement dominé
  par l’intégrateur ;
- `D=6, dt=0.5` échoue nettement, ce qui démontre pourquoi la convergence
  numérique doit être séparée de la limite physique ;
- les contrôles fins `D=6,8, dt=0.125` sont terminés. Ils montrent que
  `dt=0.5` échoue fortement et que `dt=0.25` garde un biais de pente d’environ
  `1.1e-3` à `D=8` et un leakage plus grand que le calcul fin ;
- le raffinement indépendant `D=12, dt=0.125` est terminé :
  `slope=-0.997730`, leakage `4.97e-5`, `sigma_min=0.999975` et mélange
  hors diagonale `2.68e-9`. Il confirme que le leakage primaire à `D=12`
  était dominé par l'intégrateur, mais ne suffit pas seul à promouvoir le fit
  profond mélangé en résolution ;
- toutes les déformations de chemin sont terminées à `D=8`. Leur décalage de
  phase maximal par rapport au baseline vaut `7.65e-5` et leur distance
  unitaire maximale `5.41e-5`, avec leakage entre `3.55e-4` et `5.02e-4`.
  C'est une invariance numérique robuste dans cette famille de
  paramétrisations à profondeur finie, pas une preuve autonome de protection
  topologique ;
- la composition `n=1,2,4,8` est terminée au point de fonctionnement testé.
  Le test `n=8` a nécessité une reconstruction de branche à partir des
  transports bruts `+/-theta` : la racine carrée principale seule aliasait la
  phase. Après correction, l'erreur d'additivité par rapport à une exécution
  directe vaut au plus `2.57e-2 rad` (`n=8`), le leakage pire cas reste au
  plus `2.46e-4` et aucune croissance superlinéaire n'est observée dans ce
  test. Cela qualifie la composition de la primitive `Z` abélienne au point
  testé, pas une porte universelle ni une garantie sous bruit.

Les sorties brutes se trouvent dans `results/raw/phase4_7/`; les JSON finaux,
quand les campagnes auront complété leurs contrôles, doivent être sous
`results/phase4_7/`.

## Conditions de promotion

La mention « quantification digitale fermée » exige simultanément :

1. un fit profond utilisant uniquement des points à pas convergé ;
2. une table explicite `dt=0.5,0.25,0.125` à `D=6` et `D=8` ;
3. leakage, valeurs singulières, mélange hors diagonale et gap de handoff
   pour chaque point retenu ;
4. déformations de chemin et composition complètes ; **satisfait au point de
   fonctionnement testé** ;
5. absence de contradiction entre le fit et le contrôle numérique fin.

La grille fine est maintenant complète avec `D=4,6,8,10,12` à `dt=0.125`.
Le fit uniquement fin donne `p=2.28624 ± 0.00764`, `R²=0.999983`, et toutes
les portes numériques sont satisfaites. La bonne qualification est donc :
**quantification numérique Phase 4.7 fermée pour la primitive Z abélienne**.
Voir `docs/PHASE4_7_PUBLICATION_CLOSURE.md`.
