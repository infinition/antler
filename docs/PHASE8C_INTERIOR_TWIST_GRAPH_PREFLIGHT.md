# Phase 8C-T5a  -  Premier graphe à twists intérieurs

## Verdict

**PASS comme référence stabilisateur de twists intérieurs ; aucun mouvement ni
tresse n'est déduit.** Le premier modèle à quatre défauts intérieurs est
maintenant constructible de façon reproductible, sans importer de matrice de
braid. Il utilise le formalisme de graphe de Sarkar et Yoder : les qubits sont
posés sur les sommets, les checks sur les faces, et les sommets de degré impair
identifient les twists. [Sarkar & Yoder (2024)](https://arxiv.org/abs/2101.09349)

Le script est
`experiments/phase8c/run_phase8c_interior_twist_graph_preflight.py`; la sortie
est `results/phase8c/interior_twist_graph_preflight.json`.

## Construction exacte

On part d'un tore carré périodique de taille `L x L`, avec degré quatre partout.
On retire deux arêtes disjointes. Chaque arête retirée laisse ses deux
extrémités de degré trois, et fusionne deux faces carrées en une face de poids
six. Il y a donc exactement quatre sommets impairs, tous intérieurs.

Autour d'un sommet, les opérateurs de secteur sont choisis dans l'ordre
cyclique :

```text
degré 4 : X, Z, X, Z
degré 3 : X, Y, Z
```

Deux secteurs adjacents anticommuttent. Chaque check de face est le produit des
opérateurs de secteur rencontrés en parcourant son bord. Deux faces voisines
ont alors deux anticommutions locales et commutent globalement. C'est la
construction de graphe, pas un ajustement numérique de mots de Pauli.

## Audits statiques

| Taille | Qubits | Twists intérieurs | Rang | Qubits encodés | GSD | Distance mesurée | Sondes sous la distance non scalaires |
|---|---:|---:|---:|---:|---:|---:|---:|
| `L=4` | 16 | 4 | 13 | 3 | 8 | 3 | 0 |
| `L=6` | 36 | 4 | 33 | 3 | 8 | 4 | 0 |

Le tore sans dislocation encode deux qubits dans les deux tailles. La
configuration à quatre sommets impairs en encode trois : le qubit additionnel
est cohérent avec le comptage attendu pour quatre twists, mais ce comptage ne
constitue pas à lui seul une mesure de fusion. Le parent
`-J sum_f S_f` a un gap de syndrome exact `2J` dans cette référence.

## Premier candidat de déformation

Sur `L=6`, remplacer localement une arête retirée horizontale par une arête
retirée verticale conserve les paramètres statiques (`k=3`, GSD `8`, distance
`4`) et déplace l'ensemble des sommets impairs. Neuf checks restent identiques,
cinq sont retirés et cinq sont introduits. Les checks changés entre les deux
configurations anticommuttent douze fois : la transition n'est donc pas une
simple interpolation de Hamiltoniens commutants.

C'est précisément le bon résultat de préflight : l'existence de deux codes
statiques voisins ne donne ni unitaire, ni holonomie, ni tresse. T5b doit
dériver une séquence de mesure/check-deformation ou un chemin hamiltonien gappé
qui résout ces anticommutions, et suivre les opérateurs logiques pendant cette
séquence.

## Portée et prochaine porte

Cette référence dépasse les coins fixes de T4 : les twists sont maintenant des
sommets impairs intérieurs d'un graphe explicitement fourni. En revanche les
qubits de sommet sont une **nouvelle référence externe** ; ils ne sont pas les
liens neutres de Phase 8C et ne sont pas dérivés des médiateurs charge-2
ANTLER.

La baseline T5b d'interpolation linéaire a maintenant été rejetée : elle reste
gappée mais rend le sous-espace localement lisible. T5c doit donc fournir, avant
toute affirmation de tresse :

1. une déformation locale dérivée entre deux graphes à défauts ;
2. le spectre/gap ou les issues de mesure à chaque étape ;
3. le transport explicite des trois paires de Pauli logiques ;
4. deux échanges adjacents avec `||[U1,U2]|| > 0` ;
5. Yang--Baxter seulement après ce quatrième test.

## Claim boundary

Affirmé : une famille de parents stabilisateurs à quatre twists intérieurs
commutants, avec code `k=3`, gap de syndrome et localité sous la distance,
ainsi qu'un candidat de mutation de graphe statique non trivial.

Non affirmé : transport physique `e<->m`, fusion mesurée, mouvement adiabatique,
holonomie, tresse non abélienne, dérivation ANTLER, faisabilité matérielle,
universalité ou tolérance aux fautes.
