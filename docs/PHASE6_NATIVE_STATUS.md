# Phase 6  -  statut de découverte native ANTLER

## Règle de méthode

La Phase 6 ne promeut aucune quasi-dégénérescence de petit système en qubit
topologique. Les candidats ANTLER sont testés dans l'ordre : doublet isolé,
indistinguabilité locale, scaling en taille, puis seulement algèbre de tresse.
Iemini reste une calibration externe gelée.

## 6A  -  Trois jambes scalaire, hopping corrélé et densités locales

La première famille conserve la convention Jordan--Wigner rung-major ANTLER,
la charge U(1), les hoppings scalaires SSH et de rung, plus `U_rung` et
`V_leg`. Elle exclut liens SU(2), pièges de bord et termes du parent Iemini.

Sur le scan préenregistré `L=4,N=2` et `L=5,N=3`, aucun candidat ne passe le
filtre. Le meilleur point a `split/gap=0.895` et une distinguabilité locale de
densité `0.216`, très loin des seuils `1e-3` et `0.1`.

Verdict : **rejet de cette famille minimale**, sans calcul de tresse.

## 6B  -  Ring-exchange scalaire chiral sur plaquettes trois-jambes

Un échange local orienté de plaquette a ensuite été ajouté. Il est propre à la
géométrie à trois jambes et différent du transfert de paires du parent Iemini.
Le meilleur point spectral atteint `split/gap=1.74e-2` avec un gap `0.622`,
mais est fortement lisible par une densité locale (`0.331`).

Verdict : **rejet**. Un bon spectre seul ne suffit pas : ce doublet n'est pas
un espace de fusion protégé.

## Décision de conception suivante

Les deux rejets montrent qu'empiler des termes scalaires locaux sans symétrie
de fusion dédiée ne construit pas un code. La prochaine famille ne doit pas
être un balayage arbitraire supplémentaire. Elle doit partir d'une symétrie
relative explicite  -  par exemple une parité de branches préservée par des
processus locaux pairs  -  puis dériver ou réfuter son origine microscopique
dans les degrés de liberté ANTLER.

Cette contrainte est un résultat de design : sans elle, les quasi-doublets
trouvés restent localement lisibles, donc inutilisables pour une tresse.

## 6C  -  mécanisme microscopique de paire par médiateurs à flux

Un bloc local à deux rails de médiateurs détunés a été testé à flux `pi`.
L'interférence annule exactement le cotunnelling à une particule : la norme
du bloc croisé de Schur vaut `3.48e-19` à `t=0.1`. En présence de
`U_mediator=4`, le transfert entre les deux paires basses est non nul, isolé
et suit `t^3.9938` (`R^2=0.9999996`). Le contrôle `U_mediator=0` le réduit à
`2.66e-15`; l'effet n'est donc pas un terme de paire dissimulé.

Verdict : **mécanisme local qualifié** comme ingrédient de synthèse. Ce n'est
pas encore un Hamiltonien ANTLER en échelle, un code protégé ni une tresse.
La prochaine tâche est une extension tuilée, à médiateurs explicites et corde
JW explicitement redérivée, suivie d'un préflight de protection avant tout
calcul de braid.

## 6D–6E  -  médiateur de charge 2 à parité exacte, puis ladder

Le test 6C ne préserve pas la parité de branche exactement dans son Hamiltonien
microscopique. Une alternative à médiateur `d` de charge deux conserve, elle,
les deux parités de branche exactement et engendre localement une séparation
de paire proportionnelle à `g^1.9964`. Le transfert à une particule est nul.

Le premier ladder explicite de ces médiateurs a ensuite été scanné à `L=4,5`.
Les parités restent exactes et le meilleur point est localement très peu
lisible (`1.89e-3`), mais ne satisfait pas le critère spectral :
`split/gap=1.96e-3` à `L=4`, et le même point se dégrade à `3.53e-2` à `L=5`.

Verdict : **rejet de ce scan de ladder** avant tout calcul de tresse. Le
médiateur de charge 2 reste un ingrédient local valable, mais le design de
phase many-body doit désormais être fixé par un invariant ou une limite
soluble, pas par un scan aveugle supplémentaire.

## 6F–6I  -  construction parent, audit de protection et gaps de charge

La matrice d'interaction de liaison du parent externe Iemini est factorisée
exactement en trois canaux de médiateurs de charge 2. Le Hamiltonien étendu
fixé par cette factorisation converge avec le détuning vers le frame parent
externe (`overlap=0.99536` à `Delta=1280`), avec une séparation logique qui
s'extrapole vers zéro. C'est un pont microscopique vers un benchmark, **pas**
un Hamiltonien topologique ANTLER distinct.

Le protocole est maintenant « protection avant gap ». Un opérateur de bord
fini doit satisfaire `||(1-P)[H,O]P|| -> 0` sous un scaling contrôlé. Le
générateur Iemini tronqué reste non conservé aux tailles `L=6,8`; aucune tresse
physique ne lui est attribuée. Enfin, le parent externe a un gap neutre à
charge fixée mais un coût d'ajout/retrait nul : la préparation à charge fixée
ou un charging energy reste nécessaire expérimentalement.

L'habillage de Schrieffer--Wolff de premier ordre du générateur de bord enlève
la divergence due aux médiateurs virtuels, mais laisse un résidu projeté fini
proche de `9.5` sur `Delta=80..640`. Le manque vient donc du support fini de
l'opérateur, pas d'un faux succès de spectre. Le pont n'est pas promu vers une
tresse.

L'audit matriciel compilé de ce support atteint `L=10` et est vérifié contre
l'ED creuse à `L=6`. À support maximal, fuite et résidu projeté diminuent
(`0.618 -> 0.316` et `7.22 -> 5.37` de `L=6` à `L=10`). C'est un signal de
convergence à poursuivre, mais le résidu reste grand : **aucun mode de bord
protégé ni braid physique n'est revendiqué**.

## 6L-6M -- external-calibration correction

The previously quoted maximal-support trajectory mixed two variables: both
the system size and the edge-generator support changed. It is therefore not
size scaling. At fixed support `j=3`, the normalized projected residual is
`6.3184` at `L=8` and `6.3057` at `L=10`: effectively flat and still large.
The finite-support external generator is not qualified as a protected edge
mode or physical braid.

The same matrix-free route reproduces the `L=6` neutral gap to `1.33e-14` and
finds `0.2719927` at `L=10`, continuing the sequence
`1.4621, 0.7175, 0.4181, 0.2720` for `L=4,6,8,10`. This does not prove a
thermodynamic gap closure, but it establishes no finite-size saturation and
blocks an `L=12` edge-support calculation as evidence for protection. The
Kitaev exact-zero-mode control separately validates the projected-commutator
harness (`epsilon_edge=0` for the edge, `2.0` for a bulk control).

The current native next step remains a state-to-Hamiltonian construction with
an analytic parent constraint, exact symmetries, and separate edge/gap gates;
it is not a larger blind parameter scan. See
`docs/PHASE6_EXTERNAL_CALIBRATION_STATUS.md` for the full boundary.

## Sources reproductibles

- `antler/native_threeleg.py` ;
- `experiments/phase6/run_phase6_native_threeleg_preflight.py` ;
- `experiments/phase6/run_phase6_native_ringexchange_preflight.py` ;
- `results/phase6/native_threeleg_preflight.json` ;
- `results/phase6/native_ringexchange_preflight.json` ;
- `antler/native_fusion.py` ;
- `experiments/phase6/run_phase6c_flux_pair_derivation.py` ;
- `results/phase6/flux_pair_mediator_local_preflight.json` ;
- `docs/PHASE6C_FUSION_SYMMETRY_SPEC.md` ;
- `antler/native_charge2.py` ;
- `experiments/phase6/run_phase6d_charge2_mediator_audit.py` ;
- `results/phase6/charge2_mediator_local_audit.json` ;
- `antler/native_charge2_ladder.py` ;
- `experiments/phase6/run_phase6e_charge2_ladder_preflight.py` ;
- `results/phase6/charge2_mediator_ladder_preflight.json` ;
- `docs/PHASE6D_CHARGE2_MEDIATOR_SPEC.md`.
- `antler/multiplet_mediator_parent.py` ;
- `experiments/phase6/run_phase6f_iemini_mediator_factorization.py` ;
- `results/phase6/iemini_bond_mediator_factorization.json` ;
- `experiments/phase6/run_phase6g_multiplet_parent_convergence.py` ;
- `results/phase6/multiplet_mediator_parent_convergence.json` ;
- `experiments/phase6/run_phase6h_edge_operator_preflight.py` ;
- `results/phase6/edge_operator_protection_preflight.json` ;
- `experiments/phase6/run_phase6i_charge_sector_audit.py` ;
- `results/phase6/charge_sector_audit.json` ;
- `docs/PHASE6_PROTECTION_FIRST_PROTOCOL.md`.
- `experiments/phase6/run_phase6j_sw_dressed_edge_audit.py` ;
- `results/phase6/sw_dressed_edge_operator_audit.json`.
- `experiments/phase6/run_phase6k_edge_support_scaling.py` ;
- `results/phase6/edge_support_scaling.json`.
