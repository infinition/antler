# Phase 8C-T5b  -  Rejet du déplacement linéaire naïf

## Verdict

**Rejeté comme déplacement protégé.** Entre deux graphes T5a qui diffèrent par
une mutation d'arête, l'interpolation directe

`H(s)=-(1-s) sum S_A-s sum S_B`

ne ferme pas son gap sur le contrôle exact `3 x 3`. C'est tentant, mais
insuffisant : pendant le trajet, une sonde sur un seul qubit agit de façon non
scalaire dans le sous-espace bas. Le chemin est gappé mais il n'est pas
localement protégé.

Le calcul reproductible est
`experiments/phase8c/run_phase8c_interior_twist_linear_deformation_audit.py`;
sa sortie est `results/phase8c/interior_twist_linear_deformation_audit.json`.

## Contrôle exact

Le graphe périodique `3 x 3` a neuf qubits de sommet (`512` états) et quatre
twists intérieurs. Les deux extrémités du chemin encodent `k=3` qubits, ont
GSD `8` et distance `2`. Ce petit système ne sert pas à démontrer un scaling :
il permet de diagonaliser exactement tout le chemin et de séparer clairement
gap et protection.

| Point | Dimension de bande basse | Gap au-dessus | pire lecture locale |
|---|---:|---:|---:|
| `s=0` | 8 | `4.0000` | `< 1e-8` |
| `s=0.25` | 8 | `2.2519` | `0.1627` |
| `s=0.50` | 8 | `1.1632` | `0.2937` |
| `s=0.75` | 8 | `2.2519` | `0.1627` |
| `s=1` | 8 | `4.0000` | `< 1e-8` |

La bande de huit états reste dégénérée à précision numérique sur tous les
points. Le contre-exemple est donc propre : le gap seul n'autorise pas le mot
« déplacement topologique ». À `s=0.5`, le Pauli local enregistré est de norme
non scalaire `0.2937` après projection sur la bande basse.

## Conséquence

Le mécanisme à bannir est maintenant précis : **allumer et éteindre
linéairement des checks non commutants**. Il produit un tunnel local dans le
sous-espace censé porter l'information de fusion.

T5c doit trouver une déformation à issues de mesure ou un gadget microscopique
résolu qui maintienne simultanément :

1. une bande de dimension fixe ;
2. un gap contrôlé ;
3. `P(s) O_local P(s)` scalaire pour chaque point du chemin ;
4. un transport calculé des opérateurs logiques.

Sans ces quatre conditions, calculer une holonomie serait seulement calculer
une porte d'un code localement lisible.

## Claim boundary

Affirmé : l'interpolation linéaire directe est un contre-exemple exact,
gappé mais non protégé, dans le contrôle imposé à neuf qubits.

Non affirmé : un no-go pour toute déformation de code, une tresse, une fusion,
un Hamiltonien ANTLER, une réalisation matérielle ou un calcul universel.
