# Phase 8C-T5f  -  Holonomie fermée avec check auxiliaire déclaré

Le cycle T5e se ferme si, et seulement si, on ajoute la mesure locale externe
`YIIIIIIII` sur sa jambe de retour. Les quatre jambes reviennent alors au span
stabilisateur initial, passent toutes le test de localité à un qubit et induisent
une transformation symplectique logique non identitaire.

Cette holonomie est **conditionnelle** à une nouvelle ressource de mesure ; elle
n'est ni dérivée d'ANTLER ni une tresse non abélienne. Une seule boucle ne peut
pas établir une non-commutativité. Il faut une seconde boucle indépendante et
un commutateur non nul.

Script : `experiments/phase8c/run_phase8c_auxiliary_closure_holonomy.py`.
Résultat : `results/phase8c/auxiliary_closure_holonomy.json`.
