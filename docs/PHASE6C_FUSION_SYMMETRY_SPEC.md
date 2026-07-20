# Phase 6C  -  médiateurs à flux et mécanisme de transfert de paire

## Objet exact

Les préflights 6A et 6B ont rejeté deux familles ANTLER à trois jambes : leurs
quasi-doublets restaient lisibles par une densité locale. La réponse n'est pas
d'ajouter un autre terme de paire effectif arbitraire. Cette étape vérifie un
mécanisme microscopique local qui peut, dans une extension future, préserver
la parité de chaque branche.

Le bloc contient huit modes hard-core : deux rails bas `a0,a1` et `b0,b1`,
deux rails médiateurs détunés `p0,p1` et `m0,m1`, et les paramètres
`Delta=5`, `U_mediator=4`, `E_bind=2`. Chaque `a_r -> b_r` a deux chemins. Le
chemin `m` porte une phase de Peierls `phi`; les paires sur les rails bas sont
liées par `E_bind`, et une double occupation d'un même rail médiateur coûte
`U_mediator`.

Cette construction est un **bloc local de dérivation**, pas une modification
du Hamiltonien ANTLER gelé. Elle emploie des hoppings hard-core ordinaires.
Une insertion dans une échelle ANTLER devra redériver explicitement les
cordes Jordan--Wigner rung-major de l'échelle entière.

## Sélection par interférence

Dans le secteur à une particule, l'élimination de Schur des médiateurs donne,
au second ordre,

\[
H^{(2)}_{a_r b_r}=-\frac{t^2}{\Delta}(1+e^{i\phi})+O(t^4).
\]

Ainsi, à `phi=pi`, le transfert à une particule s'annule. Le calcul exact du
bloc confirme une norme croisée de Schur de `3.48e-19` pour `t=0.1`, contre
`5.66e-3` à flux nul.

Dans le secteur à deux particules, le transfert `|a0 a1> <-> |b0 b1>` conserve
`(-1)^N_a` et `(-1)^N_b`, car il modifie chaque occupation de branche de deux.
Un processus virtuel passant par les médiateurs apparaît à l'ordre quatre :

\[
\kappa_{\rm pair}=c(\Delta,U_{\rm mediator},E_{\rm bind})t^4+O(t^6).
\]

La dépendance en `U_mediator` est une partie indispensable du test. À flux
`pi`, le contrôle avec `U_mediator=0` donne une séparation de paires maximale
de seulement `2.66e-15`; avec `U_mediator=4`, elle est non nulle et le fit
sur `t=0.05..0.20` donne une puissance `3.9938` (`R^2=0.9999996`). Le sous-
espace de paires garde une capture minimale de `0.99675` et un gap minimal de
`1.9909` sur ce scan.

## Ce qui est établi, et ce qui ne l'est pas

Établi localement : une interaction de médiateur et une interférence de flux
peuvent sélectionner un transfert de paire de quatrième ordre tout en
annulant le transfert de particule unique au second ordre. C'est une piste de
dérivation native, distincte du terme parent Iemini imposé.

Non établi : une phase topologique, un doublet localement indistinguable, une
protection en taille, une dynamique adiabatique, une tresse, une algèbre non
abélienne ou une universalité. Les nombres locaux ne doivent jamais être
présentés comme une démonstration de calcul quantique topologique.

## Prochaine porte expérimentale et numérique

Avant tout audit de code ou de tresse, construire un **Hamiltonien de ladder
tuilé avec les médiateurs explicites**  -  sans ajouter `kappa_pair` à la main  - 
et préenregistrer :

1. secteur de charge, ordre global des modes et convention de corde ;
2. symétries exactes de parité de branche ;
3. audit de doublet, indistinguabilité par une base d'opérateurs locaux et
   scaling en longueur ;
4. seulement si ces trois conditions passent, gap dynamique, braids et test
   de Yang--Baxter accompagné d'une norme de commutateur non nulle.

Sources reproductibles : `antler/native_fusion.py`,
`experiments/phase6/run_phase6c_flux_pair_derivation.py` et
`results/phase6/flux_pair_mediator_local_preflight.json`.
