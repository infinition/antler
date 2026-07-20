# Phase 8C-T5j : audit des hook errors ancilla

Le circuit T5i est soumis au défaut dangereux le plus simple : un `Z` sur
l'ancilla de lecture entre deux CNOT. Il se propage en produit suffixe des
Paulis du check sur les données. Tous les ordres CNOT (`4!` ou `6!`) ont été
testés pour les `24` étapes de check de poids élevé des boucles T5f/T5g.

Résultat : chaque étape possède au moins un ordre sans hook logique. Les
suffixes propagés sont tous détectés par le groupe stabilisateur courant ou
scalaires ; le pire nombre minimal de hooks logiques est `0`.

C'est un PASS important mais étroit : il ne couvre ni fautes CNOT arbitraires,
ni préparation/lecture, ni décodage répété, ni bruit stochastique. Il ne rend
pas le circuit tolérant aux fautes, ni natif ANTLER.

Script : `experiments/phase8c/run_phase8c_ancilla_hook_error_audit.py`.
Résultat : `results/phase8c/ancilla_hook_error_audit.json`.
