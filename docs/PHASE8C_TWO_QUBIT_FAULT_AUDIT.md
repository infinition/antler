# Phase 8C-T5k : fautes de CNOT à deux qubits

Le contrat de décision lance son premier gate complet : les 15 Paulis non
identitaires injectés après chaque CNOT data--ancilla sont propagés exactement.
Tous les ordres `4!`/`6!` sont comparés aux stabilisateurs courants pour les
24 étapes de checks T5f/T5g.

Chaque étape admet un ordre sans faute CNOT unique induisant un Pauli logique
sur les données. Les flips de bit de lecture subsistent (`912` dans les
horaires optimaux) : ce PASS ne remplace donc ni extraction répétée, ni
décodage, ni étude de bruit.

Script : `experiments/phase8c/run_phase8c_two_qubit_fault_audit.py`.
Résultat : `results/phase8c/two_qubit_fault_audit.json`.
