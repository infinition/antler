# Phase 8C-T4  -  Calibration d’espace de code à quatre twists de bord

## Verdict

**PASS comme référence de codage/fusion ; pas de fusion mesurée ni de tresse.** Une famille de patches planaires CSS de distances `3` et `5` possède quatre coins où une frontière `X` rencontre une frontière `Z`. Ces coins sont la réalisation de référence de twists de bord. Le doublet de code reste localement indistinguable, tandis que son support logique minimal croît avec la séparation des coins.

Le script reproductible est `experiments/phase8c/run_phase8c_corner_twist_fusion_calibration.py`; la sortie est `results/phase8c/corner_twist_fusion_calibration.json`.

## Parent de référence

Pour une grille carrée impaire de taille `d x d`, les qubits neutres portent :

- des plaquettes internes de poids quatre, alternativement `X` et `Z` ;
- des checks de bord de poids deux, `Z` en haut/bas et `X` à gauche/droite, dans le motif alterné qui ferme le patch.

Les quatre coins joignent donc des frontières de type opposé. Cette association entre coins de surface code et twists est une construction de code connue ; elle est ici employée seulement comme référence stabilisateur, sans Hamiltonien ANTLER microscopique. [Yoder & Kim (2016)](https://arxiv.org/abs/1612.04795)

## Tests exacts

| Quantité | `d=3` | `d=5` |
|---|---:|---:|
| Qubits | 9 | 25 |
| Rang stabilisateur | 8 | 24 |
| Dimension du code | 2 | 2 |
| Séparation minimale des coins | 2 | 4 |
| Distance logique | 3 | 5 |
| Gap de syndrome | `2J` | `2J` |
| Sondes de Pauli sous la distance | 351 | 1 089 525 |
| Action logique locale sous la distance | 0 | 0 |

Les cordes logiques `X` et `Z` anticommuttent dans les deux tailles. Le résultat important est le scaling : en éloignant les coins, la dimension du doublet reste deux mais toute action logique exige un support plus grand. C’est le test de protection locale requis pour interpréter le doublet comme la calibration d’un secteur de codage à quatre twists, et non comme une dégénérescence accidentelle.

## Portée physique

Cette étape donne une cible honnête pour la fusion : avec plusieurs twists, le secteur protégé doit être une dimension de code localement illisible qui survit à l’augmentation de leur séparation. Elle ne mesure pas encore une issue de fusion `1` contre `psi`, parce qu’aucun défaut intérieur n’est créé ni mesuré. Elle n’autorise donc pas encore les mots « anyon non abélien réalisé » ou « tresse ».

## Porte T5

T5 doit remplacer les coins fixes par des twists **intérieurs et déformables**, définir un Hamiltonien de déformation locale, puis propager deux échanges adjacents. Les sorties minimales sont `U1`, `U2`, leakage, gaps, projection locale et `||[U1,U2]||`. L’équation de Yang--Baxter n’est interprétable que si ce commutateur est nettement non nul.

## Claim boundary

Affirmé : la référence de code à quatre coins de type twist possède un doublet de distance `3` puis `5`, et passe l’audit exhaustif de localité sous la distance.

Non affirmé : des twists mobiles, une mesure de fusion, une statistique non abélienne, une tresse, une universalité, une dérivation ANTLER, ou une réalisation expérimentale.
