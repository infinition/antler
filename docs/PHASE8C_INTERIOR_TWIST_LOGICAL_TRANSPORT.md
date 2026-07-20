# Phase 8C-T5d  -  Transport logique et frames de Pauli

## Verdict

**PASS comme transport logique de référence pour une déformation mesurée.** Les
quatre mesures T5c transportent une base complète de trois paires de Pauli
logiques. À chaque étape, les six représentants restent hors du span
stabilisateur, commutent avec les checks courants et conservent leur forme
symplectique canonique.

Le script est `experiments/phase8c/run_phase8c_interior_twist_logical_transport.py`;
la sortie est `results/phase8c/interior_twist_logical_transport.json`.

## Contrôles passés

- La séquence est celle de T5c : quatre mesures de checks finaux à issue `+1`.
- Les six Paulis logiques sont suivis après chaque mise à jour ; leur matrice
  symplectique reste exactement trois blocs `[[0,1],[1,0]]`.
- À l'arrivée, ils sont exprimés dans une base logique finale explicitement
  déclarée dans le JSON, sans imposer une porte logique.
- Les `2^4=16` patterns d'issues ont chacun un frame de Pauli chronologique,
  modulo phase globale, formé par les pivots pré-mesure.

Cette étape résout le problème « quelle information logique a été transportée ? »
pour **un** déplacement abstrait. Elle ne donne pas encore une tresse : il faut
construire un second déplacement adjacent sur le même espace logique, composer
les deux transports et seulement alors évaluer leur commutateur.

## Claim boundary

Affirmé : transport logique stabilisateur à issue `+1` et table de frames pour
toutes les outcomes sur le graphe de référence `3 x 3`.

Non affirmé : mouvement physique, dispositif de mesure, Hamiltonien ANTLER,
fusion lue, échange, tresse non abélienne, universalité ou tolérance aux fautes.
