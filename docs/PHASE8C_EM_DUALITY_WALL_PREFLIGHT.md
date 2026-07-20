# Phase 8C-T3a  -  Calibrage de dualité `e<->m` et rejet des faux murs

## Verdict

**PASS de calibrage, pas de twist.** Le parent `Z2` pur de T2 possède une dualité globale exacte qui transforme étoiles et plaquettes, ainsi qu'une corde logique `Z` en corde logique `X`. Ce calcul fixe une condition nécessaire, mesurable et non ambiguë pour un mur `e<->m` futur. Il montre aussi que deux propositions plus simples ne satisfont pas cette condition.

Le script est `experiments/phase8c/run_phase8c_em_duality_wall_preflight.py`; le résultat est `results/phase8c/em_duality_wall_preflight.json`.

## Dualité calibrée

Sur le tore `3x3`, poser la permutation duale des liens

\[
h(x,y)\mapsto v(x,y-1),\qquad
v(x,y)\mapsto h(x-1,y),
\]

puis appliquer Hadamard sur tous les liens. L'opération globale `D` vérifie

\[
D A_v D^\dagger=B_{p(v)},\qquad D B_p D^\dagger=A_{v(p)}.
\]

Les neuf étoiles deviennent exactement les neuf plaquettes et réciproquement. En particulier, la corde `Z[h(0,0)] Z[h(1,0)] Z[h(2,0)]` devient la corde logique non triviale `X[v(0,2)] X[v(1,2)] X[v(2,2)]`.

Cette transformation établit l'algèbre de permutation d'anyons que le futur mur doit reproduire localement. Elle est **globale**, non bornée spatialement, et ne possède donc aucune extrémité de twist.

## Contrôles négatifs

| Candidat | Résultat | Pourquoi il est rejeté |
|---|---|---|
| Cocycle/signe statique `pi` | Échec | Une phase modifie des coefficients, pas le type de Pauli : aucune corde `Z` ne devient une corde `X`. |
| Hadamard sur une colonne bornée | Échec | Les checks restent commutants et deviennent parfois mixtes, mais aucune étoile ne devient une plaquette ni inversement. C'est un changement de base local, pas une dualité de support. |

Le second contrôle est particulièrement important : des checks non-CSS ou un pentagone local peuvent être nécessaires dans une dislocation, mais leur seule présence ne prouve jamais une permutation `e<->m`.

## Porte T3b

T3b doit spécifier une **cellulation locale complète** avec un segment fini de coupure, ses checks d'extrémité et les cordes physiques de part et d'autre. Il devra établir simultanément :

1. commutation, rang et gap du parent modifié ;
2. cordes de syndromes dont le type est échangé au franchissement ;
3. localité des extrémités ;
4. indistinguabilité locale du sous-espace pertinent.

Un tel résultat serait un défaut de twist de référence. Il faudrait encore deux paires de défauts, une fusion et des transports dérivés pour atteindre le test de tresse non commutative.

## Claim boundary

Affirmé : la dualité globale `e<->m` et les échecs des deux faux murs sont établis exactement sur le modèle de référence à liens neutres.

Non affirmé : un mur local, une dislocation, des twists, une fusion, une tresse non abélienne, une dérivation ANTLER ou une implémentation matérielle.
