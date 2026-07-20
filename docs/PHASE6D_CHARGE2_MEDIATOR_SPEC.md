# Phase 6D–6E  -  médiateur de charge 2 et ladder à parité exacte

## Pourquoi cette bifurcation

Le bloc 6C à flux `pi` sélectionne un transfert de paire à basse énergie, mais
ses hoppings élémentaires à une particule ne conservent pas exactement les
parités de branche. Il ne peut donc pas, seul, fonder une protection de parité
microscopique.

La Phase 6D introduit explicitement un degré de liberté hard-core `d` de
charge U(1) égale à deux. Sur une liaison, l'interaction est

\[
-g\left[d^\dagger(a_0a_1+b_0b_1)+\mathrm{h.c.}\right].
\]

Elle conserve exactement la charge pondérée `N_a+N_b+2N_d`, ainsi que
`(-1)^N_a` et `(-1)^N_b`. Elle ne connecte pas le secteur à une particule
entre les rails. L'élimination du médiateur détuné donne un transfert de paire
à l'ordre `g^2/Delta`, sans ajouter directement un opérateur de transfert de
paire `a†a†bb`.

## Audit local 6D

Le bloc à cinq modes satisfait les trois tests préenregistrés :

- commutateurs avec chacune des deux parités : exactement `0` ;
- norme de transfert inter-rail dans le secteur de charge un : exactement
  `0` ;
- séparation de paires non nulle et isolée, avec puissance mesurée
  `g^1.9964` (`R^2=0.9999992`).

Le médiateur de charge 2 est donc un **ingrédient local à symétrie exacte**.
Cette assertion dépend d'un nouveau degré de liberté et d'une conversion
atomes–molécule locale ; elle n'est pas encore une dérivation expérimentale à
partir des seuls paramètres de l'échelle ANTLER gelée.

## Préflight du ladder 6E

Le mécanisme a été tuilé sur les liaisons d'une échelle à deux rails, en
conservant les médiateurs explicitement. Le scan préenregistré comporte 144
points (`L=4,5`, charge totale 4), avec hopping intra-rail, conversion `g`,
détuning et interactions de densité. Aucun terme du parent Iemini n'est inséré.

Résultat : **aucun candidat ne passe simultanément**

\[
\frac{\delta E_{\rm logique}}{\Delta_{\rm iso}}<10^{-3},\qquad
d_{\rm local}<0.1,\qquad \Delta_{\rm iso}>10^{-2}.
\]

Le meilleur point `L=4` a une distinguabilité de densité excellente
(`1.89e-3`) et une symétrie de parité exacte, mais
`split/gap=1.96e-3`, donc manque le seuil spectral. Surtout, ce même point
passe à `3.53e-2` à `L=5`; il ne montre pas une tendance protectrice.

## Conséquence de design

La condition de symétrie exacte est désormais satisfaite par une extension
concrète, mais elle ne suffit pas. Le prochain Hamiltonien ne devra pas être
un balayage aveugle autour de ce near-miss. Il devra partir d'un invariant ou
d'une limite soluble qui impose à la fois la structure de parité, les termes
diagonaux induits par l'élimination de `d`, et la bonne phase à grand volume.
Puis seulement : scaling, base complète d'opérateurs locaux, dynamique et
tresse.

Ce résultat ne démontre ni phase topologique, ni code protégé, ni braid, ni
ordinateur quantique topologique. Il identifie précisément la différence entre
un mécanisme local correct et une phase de many-body protégée.

Sources : `antler/native_charge2.py`,
`experiments/phase6/run_phase6d_charge2_mediator_audit.py`,
`results/phase6/charge2_mediator_local_audit.json`,
`antler/native_charge2_ladder.py`,
`experiments/phase6/run_phase6e_charge2_ladder_preflight.py` et
`results/phase6/charge2_mediator_ladder_preflight.json`.
