# Phase 8C-T5e  -  Obstruction de clôture du premier cycle

Trois mutations de graphe ont une séquence de mesures qui passe le test de
localité à un qubit. La quatrième, nécessaire pour revenir au graphe initial,
ne passe dans aucun des `13 699` ordres non répétés des sept checks finaux.

Le cycle n'est donc pas fermé : aucune holonomie, aucun échange et aucun
commutateur ne sont rapportés. Le résultat est un rejet de cette grammaire de
clôture, pas un no-go général des déformations à twists.

Script : `experiments/phase8c/run_phase8c_interior_twist_cycle_closure_audit.py`.
Résultat : `results/phase8c/interior_twist_cycle_closure_audit.json`.

La prochaine ressource à dériver est un check auxiliaire ou une cellulation
intérieure différente, soumis au même audit de rang, GSD et localité avant de
retenter un cycle fermé.
