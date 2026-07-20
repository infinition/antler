# Phase 6  -  protocole « protection avant gap »

## Principe

Un petit `split/gap` est nécessaire pour un code statique, mais ne démontre ni
la localité de l'information logique, ni l'existence d'un opérateur de bord,
ni une tresse. Toute suite Phase 6 applique maintenant cet ordre :

1. symétries exactes et charge pondérée ;
2. construction par état connu / contraintes locales, pas par scan aveugle ;
3. opérateur de bord quasi-conservé dans le sous-espace de code ;
4. gap neutre du code et gaps entre secteurs de charge, rapportés séparément ;
5. seulement alors dynamique, Berry/Wilson loop et braid.

Le bon test d'opérateur de bord à taille finie est

\[
\epsilon_{\rm edge}=\|(1-P)[H,O_{\rm edge}]P\|.
\]

Pour un parent dont `H P=0`, cela devient `||H O_edge P||`. La norme complète
de `[H,O_edge]` sur tout l'espace contient des états excités arbitraires : elle
est un diagnostic secondaire, pas un certificat logique. L'opérateur doit de
plus avoir une action logique non triviale et un support dont la troncature
converge avec la taille.

## Résultats appliqués

### Brique à flux

Le médiateur à flux `pi` annule le transfert à une particule, mais les hoppings
microscopiques ne préservent pas une parité de branche exacte. Il ne franchit
donc pas la première porte pour un code protégé.

### Médiateur simple de charge 2

Le médiateur de charge 2 conserve les parités exactes et dérive localement un
transfert de paire. Son premier ladder explicite est néanmoins rejeté : le
near-miss `L=4` se dégrade à `L=5`. La symétrie seule ne sélectionne pas la
bonne phase many-body.

### Construction parent par médiateurs multiples

L'interaction de liaison du parent externe d'Iemini est négative semi-définie
de rang trois : ses valeurs propres sont `(-8,-4,-4,0,0,0)`. Trois médiateurs
de charge 2 par liaison la factorisent exactement à l'ordre deux. Deux des
médiateurs portent une charge de parité généralisée `(-1,-1)` ; les parités
exactes sont donc celles du système étendu, et coïncident avec les parités de
rails dans le sous-espace sans médiateur.

Sans scan de couplages, le Hamiltonien microscopique ainsi fixé converge vers
le frame parent externe : recouvrement minimal `0.99536` à `Delta=1280`; la
séparation logique suit `a/Delta+b`, avec `b=-5.98e-4` et `R^2=0.99987` sur la
queue contrôlée. C'est un **pont microscopique vers un benchmark externe**,
pas un nouveau Hamiltonien topologique ANTLER.

## Les deux garde-fous qui restent ouverts

Le générateur de bord Iemini à support fini ne passe pas encore le test de
protection à `L=6,8` : même si sa fuite diminue lorsque le support augmente,
`||(1-P)[H,O]P||` reste grand aux tailles calculées. Il est conservé comme
diagnostic de convergence, jamais comme tresse physique validée.

À remplissage `N/(2L)=1/2`, le parent externe a un gap neutre fini à taille
finie, mais un coût d'ajout/retrait de particule nul à la précision numérique
pour `L=4,6,8`. Il ne constitue donc pas une mémoire incompressible autonome
en contact avec un réservoir. Une préparation à charge fixée ou une énergie de
charging doit faire partie d'une proposition expérimentale.

## Porte obligatoire pour tout nouveau modèle

Avant de lancer une nouvelle ED de Phase 6, le modèle doit fournir :

- une symétrie de parité exacte, y compris les charges de tout médiateur ;
- l'absence exacte de transfert à une particule dans le sous-espace bas, si
  l'architecture repose sur une sélection de paire ;
- un état/une contrainte parent ou un invariant analytique précis ;
- une proposition d'opérateur de bord et une loi de convergence de
  `epsilon_edge` ;
- les gaps neutre et de charge séparés.

L'étape suivante n'est donc pas une nouvelle grille `L=4,5,6`. C'est la
dérivation de l'opérateur de bord **habillé par Schrieffer--Wolff** pour le
pont à médiateurs, puis son audit `epsilon_edge(L,support,Delta)`. Si cet
opérateur ne converge pas, le pont est seulement un reproducer de spectre et
ne doit pas être promu vers une tresse.

Le premier habillage a été effectué. Avec
`S_high,low=-H_high,low/Delta`, l'opérateur physique est
`O_phys=O+[S,O]+O(Delta^-1)`. Il supprime bien la croissance virtuelle du
commutateur brut : à `Delta=640`, la norme projetée normalisée passe de
`82.7` à `9.55`. Elle plafonne toutefois autour de `9.5` sur
`Delta=80..640`, loin de zéro. Ce résultat isole le verrou : le dressing de
médiateur est nécessaire mais le générateur de bord à support fini lui-même
reste insuffisant. Aucun braid n'est promu.

Une application matricielle compilée du parent, vérifiée exactement contre
l'ED creuse à `L=6`, permet ensuite un premier scaling à `L=10` sans
matérialiser la matrice `184756 x 184756`. Au support maximal disponible
`j=2,3,4` pour `L=6,8,10`, la fuite passe de `0.618` à `0.440` puis `0.316` et
le résidu projeté normalisé de `7.22` à `6.32` puis `5.37`. La direction est
encourageante, mais aucune borne asymptotique n'est établie et le résidu reste
très supérieur à zéro. C'est une tendance de convergence, pas un opérateur de
bord qualifié ni une tresse.

## Correction 6L-6M: prerequisite order enforced

The projected-commutator implementation now has an exact positive control:
an ideal Kitaev edge Majorana gives normalized `epsilon_edge=0`, whereas a
bulk Majorana gives `2.0`. The large external-parent residual is therefore a
negative result of the diagnostic, not an uncalibrated convention artifact.

The prior `L=6,8,10` maximal-support path also changed support from `j=2` to
`j=4`; it cannot be interpreted as size scaling. At fixed `j=3`, the residual
is `6.3184` at `L=8` and `6.3057` at `L=10`. Separately, the neutral gap falls
to `0.2719927` at `L=10` after an independent matrix-free `L=6` cross-check.
No `L=12` support calculation is launched as evidence for a protected braid.

Sources : `experiments/phase6/run_phase6f_iemini_mediator_factorization.py`,
`experiments/phase6/run_phase6g_multiplet_parent_convergence.py`,
`experiments/phase6/run_phase6h_edge_operator_preflight.py`,
`experiments/phase6/run_phase6i_charge_sector_audit.py`,
`experiments/phase6/run_phase6j_sw_dressed_edge_audit.py`,
`experiments/phase6/run_phase6k_edge_support_scaling.py` et les JSON associés
sous `results/phase6/`.
