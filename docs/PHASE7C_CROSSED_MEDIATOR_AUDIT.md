# Phase 7C  -  médiateur charge-2 croisé à support plaquette

## Primitive testée

Deux médiateurs charge-2 sont chacun couplés de manière cohérente à deux
paires disjointes dans une même plaquette de quatre barreaux :

\[
d_1^\dagger(a_0a_1+b_2b_3)+\mathrm{h.c.},\qquad
d_2^\dagger(b_0b_1+a_2a_3)+\mathrm{h.c.}
\]

La signature de parité de chaque médiateur est `(0,0)`, donc la primitive
conserve les parités de rail nues. Il s'agit cependant d'un nouveau support
microscopique de plaquette : le catalogue précédent ne couvrait que les
médiateurs attachés à un unique lien.

## Résultat local exact

À `U=20`, `Delta=10`, `g=0.5`, le bloc exact donne
`c_XXXX=-2.03e-7`, contre environ `-1.37e-10` pour le meilleur canal
lien-par-lien. Le sweep `g=0.30..1.50` mesure

\[
c_{XXXX}\propto g^{3.999},\qquad \|H_{\rm unwanted}\|\propto g^{1.999}
\]

avant compensation. C'est la première confirmation numérique d'un mécanisme
quatre-corps d'ordre quatre dans la grammaire charge-2 élargie.

Les termes d'ordre deux dominants sont diagonaux (`Z`, `ZZ`). Pour chaque
point de couplage, l'audit les extrait du même bloc non compensé et ajoute les
contre-termes statiques opposés autorisés (`rail_biases`, `zz_couplings`).
Après cette calibration à `g=0.5`, la norme de tous les autres mots de Pauli,
normalisée par l'échelle cible fixée `|c_target|=0.5`, est `1.14e-6` et suit
elle-même une loi `g^4`. Cette normalisation ne mesure **pas** la sélectivité
par rapport au très faible `c_XXXX` effectivement observé.

## Mur restant

La compensation ne produit pas un stabilisateur pur. Les opérateurs quatre
corps `YYYY`, `XXYY`, `XYXY`, etc. restent de même ordre que `XXXX`; l'alignement
sur le seul `XXXX` n'est que `0.1135`, soit une norme des autres mots de Pauli
d'environ `2.80` fois la norme du `XXXX` observé. Son coefficient absolu demeure
`-2.03e-7`, très loin de `-0.5`, et le résidu à échelle cible fixée reste
presque unitaire.

La calibration est obtenue à partir du même bloc fini qui est évalué. Elle ne
constitue donc pas encore une dérivation indépendante de Schrieffer--Wolff, ni
une recette expérimentale de contre-termes, ni une preuve que les motifs
peuvent être tuilés sans conflit entre plaquettes.

## Décision

La primitive croisée est **qualifiée comme direction locale**, contrairement
aux canaux lien-par-lien, mais n'est pas promue en parent ANTLER 2D. Avant une
optimisation continue, il faut :

1. dériver analytiquement le coefficient d'ordre quatre et les contre-termes ;
2. identifier la combinaison de stabilisateurs réellement générée ;
3. vérifier la compatibilité des plaquettes qui se chevauchent ;
4. trouver un mécanisme qui sépare `XXXX` des autres termes quatre-corps ou
   modifier honnêtement la cible de code.

Les primitives the submitted candidate de charge 4 et de contrainte auxiliaire ne sont pas
implémentées ici : la première introduit un nouveau opérateur de conversion à
quatre fermions, la seconde impose déjà un terme quatre-corps. Elles demandent
une spécification microscopique séparée avant tout audit de promotion.
