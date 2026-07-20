# Phase 8C-T5i : compilation ancilla des checks de déformation

Les boucles T5f/T5g nécessitent dix-neuf checks de graphe distincts de poids
`4` ou `6`. Ils ne sont plus laissés comme projecteurs abstraits : chacun est
compilé dans un circuit de mesure standard avec un ancilla de lecture, des
rotations de base `H/S` et `4` ou `6` CNOT data-vers-ancilla.

Sur trois vecteurs d'état reproductibles de dimension `2^9`, les deux branches
de l'instrument sont vérifiées directement : `M_+ = (I + P)/2` et
`M_- = (I - P)/2`. Le pire résidu est `1.44e-15`.

Ce résultat établit seulement un contrat de circuit de référence. L'ancilla,
les CNOT, leur connectivité, les hook errors, la fidélité et la dérivation
depuis ANTLER restent des ressources ou problèmes externes.

Script : `experiments/phase8c/run_phase8c_ancilla_check_circuit_audit.py`.
Résultat : `results/phase8c/ancilla_check_circuit_audit.json`.
