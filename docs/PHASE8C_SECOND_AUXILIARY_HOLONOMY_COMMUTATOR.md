# Phase 8C-T5g : deux holonomies conditionnelles et commutateur logique

Ce test ferme une seconde boucle de déformation sur le même code de graphe
`3 x 3`, puis compare sa transformation logique à la boucle conditionnelle
T5f. C'est un étalon algébrique du banc de mesures, **pas** une tresse
d'anyons et pas une réalisation microscopique ANTLER.

La deuxième boucle tourne la seconde arête retirée autour du sommet `(1,1)`
selon `E -> N -> W -> E`. Les trois jambes conservent rang `6`, GSD `8` et le
test local enregistré : les 27 Paulis à un qubit se projettent toujours sur
zéro ou un scalaire. Deux checks ponctuels supplémentaires sont toutefois
nécessaires : `IIIIIIIZI` et `IIIYIIIII`. Avec le check `YIIIIIIII` déjà
déclaré pour T5f, ce sont des ressources externes, non dérivées d'ANTLER.

Les deux boucles reviennent exactement au même groupe stabilisateur. Dans la
base logique `X1,Z1,X2,Z2,X3,Z3`, leurs matrices symplectiques GF(2) ne
commutent pas : `A B != B A`, et le commutateur de groupe `A B A^-1 B^-1`
est non identitaire. C'est la première calibration de non-commutativité du
pipeline de transport logique.

Un correctif de robustesse accompagne ce résultat. Les représentations GF(2)
en forme échelonnée simple dépendent de l'ordre d'insertion ; elles ne doivent
jamais être comparées directement pour conclure que deux spans diffèrent. Les
scripts T5c/T5d utilisent désormais une forme échelonnée réduite canonique.
T5c, T5d, T5e et T5f ont été rejoués après ce correctif : leurs PASS restent
des PASS, et l'échec de fermeture T5e dans la grammaire *sans* auxiliaire
reste un échec.

Script : `experiments/phase8c/run_phase8c_second_auxiliary_holonomy_commutator.py`.
Résultat : `results/phase8c/second_auxiliary_holonomy_commutator.json`.

## Limite des claims

La non-commutativité obtenue est conditionnée par trois mesures locales
imposées. Elle n'établit ni lignes d'univers de défauts, ni secteur de fusion
à deux dimensions, ni appareil de mesure, ni Hamiltonien microscopique
ANTLER, ni holonomie adiabatique, ni relation de Yang--Baxter. Elle ne permet
donc aucune revendication de statistique non abélienne, d'universalité ou
d'ordinateur quantique topologique.

La prochaine porte est physique : spécifier une mesure/commande qui préserve
la symétrie, la dériver d'un Hamiltonien microscopique, puis certifier gap,
leakage, localité de défaut et frames d'outcome avant de reparler d'échange.
