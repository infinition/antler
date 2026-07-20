# Phase 8C-T5c  -  Déformation de code résolue par mesures

## Verdict

**PASS comme protocole de référence à mesures de stabilisateurs.** T5b avait
montré qu'un chemin hamiltonien gappé peut néanmoins rendre le code lisible.
Ici, quatre mesures de checks de face finaux remplacent les générateurs qui
anticommutent, sans réduire la dimension du code et sans faire apparaître de
Pauli logique de poids un.

Le script est
`experiments/phase8c/run_phase8c_interior_twist_measurement_deformation.py`;
la sortie est `results/phase8c/interior_twist_measurement_deformation.json`.

## Séquence obtenue

Sur le contrôle périodique `3 x 3`, on fixe l'issue `+1` de chaque mesure. Les
issues `-1` ne changent pas l'algèbre des stabilisateurs mais imposent un frame
de Pauli qui n'est pas encore suivi dans ce livrable.

```text
1. IIIYYIXZI
2. IIIZIYZIX
3. IZXIXZIXZ
4. XZIIIIZXI
```

Chaque mesure anticommute avec exactement deux générateurs courants. La règle
de mise à jour remplace l'un d'eux par le check mesuré et multiplie l'autre par
le pivot ; elle préserve ainsi le rang sans ajouter de terme logique choisi à
la main.

| Étape | Rang | GSD | Checks qui anticommuttent avant mesure | Paulis à un qubit non scalaires |
|---|---:|---:|---:|---:|
| Initiale | 6 | 8 |  -  | 0 / 27 |
| 1 | 6 | 8 | 2 | 0 / 27 |
| 2 | 6 | 8 | 2 | 0 / 27 |
| 3 | 6 | 8 | 2 | 0 / 27 |
| Finale | 6 | 8 | 2 | 0 / 27 |

Le span stabilisateur final est exactement celui du second graphe à quatre
twists. Contrairement à T5b, aucun état intermédiaire de ce protocole
stabilisateur ne laisse un Pauli à un qubit agir de manière logique.

## Ce que ce résultat établit

Il existe maintenant une déformation de code de référence qui respecte les
deux conditions minimales manquantes : dimension de fusion fixe et
indistinguabilité locale à chaque étape enregistrée. Le chemin linéaire
continu est rejeté ; la voie à mesures est retenue comme cible pour le transport
logique.

## Porte T5d

T5d doit suivre une base complète de six Paulis logiques à travers la séquence,
pour toutes les issues de mesure et leurs frames. Ensuite, deux mutations
adjacentes doivent être composées : le commutateur ne sera calculé que si les
deux transports sont réellement définis sur le même espace logique.

## Claim boundary

Affirmé : une séquence abstraite de mesures de stabilisateurs maintient sur le
petit graphe la commutation, le rang, GSD `8` et la scalarité des Paulis à un
qubit.

Non affirmé : un appareil de mesure, une dérivation ANTLER, le transport
physique `e<->m`, la fusion lue expérimentalement, une holonomie, une tresse
non abélienne, l'universalité ou la tolérance aux fautes.
