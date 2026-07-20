# Phase 8C-T2  -  Patch de code de jauge `Z2` pur

## Verdict

**PASS comme référence de code abélien imposée ; aucune promotion microscopique ANTLER ni non abélienne.** Le tore carré `3 x 3` construit uniquement avec les 18 qubits de lien neutres `tau_e` passe les portes qui n'étaient pas accessibles sur l'étoile T0 ni la plaquette T1 : espace fondamental à quatre dimensions, distance trois, gap de syndrome et indistinguabilité locale complète sous la distance.

Le calcul reproductible est `experiments/phase8c/run_phase8c_z2_code_patch_preflight.py`; son résultat est `results/phase8c/z2_code_patch_preflight.json`.

## Contrat de référence

La matière est absente ou strictement gelée. Le générateur de Gauss devient

\[
A_v=\prod_{e\ni v}\tau^x_e,
\]

et le flux est

\[
B_p=\prod_{e\in\partial p}\tau^z_e.
\]

Le parent de référence est le code torique à projecteurs commutants

\[
H_{\rm ref}=-J_s\sum_v A_v-J_p\sum_p B_p,
\qquad J_s=J_p=1.
\]

Les qubits `tau_e`, les checks étoile et plaquette, et ce Hamiltonien sont des **ressources de référence explicitement insérées**. Ils ne sont pas dérivés du ladder ANTLER gelé ni de ses médiateurs charge-2.

## Résultats exacts

| Quantité | Valeur |
|---|---:|
| Liens / étoiles / plaquettes | `18 / 9 / 9` |
| Hilbert des liens | `2^18 = 262144` |
| Rang Gauss indépendant | `8` |
| Dimension du secteur `A_v=+1` | `1024` |
| Rang stabilisateur total | `16` |
| Qubits encodés / GSD | `2 / 4` |
| Distance minimale | `3` |
| Gap de syndrome exact | `4` |
| Pauli locaux contrôlés (poids 1--2) | `1431` |
| Actions logiques non scalaires sous distance | `0` |

Les relations globales `prod_v A_v=I` et `prod_p B_p=I` expliquent les deux contraintes dépendantes. Deux boucles logiques de poids trois ont été extraites et anticommuttent, comme attendu pour un couple de Paulis logiques sur le tore.

Pour chaque mot de Pauli non identitaire supporté sur au plus deux liens, le test exact donne soit une anticommutation avec un check (donc une projection nulle dans le code), soit un élément du stabilisateur (donc une projection scalaire). Aucun mot local ne réalise une action logique. C'est le test de protection locale requis avant tout défaut.

## Ce que T2 rend possible  -  et ce qu'il ne rend pas possible

T2 qualifie le substrat **abélien** dans lequel un défaut peut être défini. La prochaine porte T3 est un mur de domaine qui transforme réellement les cordes électriques et magnétiques l'une dans l'autre : `e <-> m`. Une phase de hopping `pi`, une chaîne de signes, ou un flux `B_p=-1` ne remplit pas ce contrat : ce sont respectivement une convention de jauge locale ou un vison, pas une permutation d'anyons.

T2 ne démontre ni défaut de twist, ni espace de fusion, ni braid non commutatif ni universalité. Il ne dit pas non plus que le code reste stable après ajout de matière mobile : les charges de matière seront des défauts contrôlés à auditer après le mur, et non le support initial de la mémoire.

## Claim boundary

Affirmé : un parent de jauge `D(Z2)` pur, à liens neutres déclarés, satisfait les critères algébriques et locaux d'un code torique `3 x 3`.

Non affirmé : une dérivation ANTLER, une phase thermodynamique, un mur `e<->m`, des twists, une fusion non abélienne, une tresse, une porte universelle, du bruit matériel ou une tolérance aux fautes.
