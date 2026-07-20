# Phase 8C-T5h : contrôle anti-fausse holonomie locale

T5f/T5g emploient des checks de graphe de poids élevé, plus trois mesures
ponctuelles externes. Avant de leur attribuer même une signification de
déformation de défaut, il faut exclure le mécanisme plus banal suivant : une
suite courte de mesures à un seul qubit fabriquerait seule une boucle logique
non commutative.

Le contrôle part du même code de graphe périodique `3 x 3` (rang `6`, GSD
`8`) et explore exhaustivement toutes les mesures de Pauli non identitaires à
un sommet (`27` choix), dans la branche d'outcome `+1`, jusqu'à une profondeur
de six. Chaque étape doit préserver le rang et le gate local : les `27` Paulis
à un qubit doivent rester scalaires ou nuls après projection.

Résultat : aucune boucle non vide ne revient au span stabilisateur initial aux
profondeurs `1` à `6`. Les nombres d'états protégés nouveaux explorés sont
`20, 284, 1704, 5260, 8364, 5416`. Dans cette grammaire et cette borne, les
checks de déformation de poids élevé sont donc opérationnellement nécessaires
aux boucles T5f/T5g.

Ce contrôle renforce l'interprétation très limitée de T5g, mais ne la
transforme pas en tresse : il ne traite ni les circuits locaux plus longs, ni
un appareil physique pour les checks, ni les lignes d'univers de défauts.

Script : `experiments/phase8c/run_phase8c_local_measurement_only_control.py`.
Résultat : `results/phase8c/local_measurement_only_control.json`.
