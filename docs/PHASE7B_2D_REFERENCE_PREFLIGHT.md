# Phase 7B  -  référence stabilisateur 2D et frontière microscopique

## Verdict

Le premier contrôle exact d'une extension **genuinement bidimensionnelle**
passe : un parent de type code torique, encodé dans des barreaux ANTLER gelés à
une particule, possède bien une algèbre de stabilisateurs commutante et une
indistinguabilité locale complète jusqu'à la distance du code. Ce contrôle ne
réalise toutefois pas encore un Hamiltonien ANTLER natif à faible corps. Il est
une cible de calibration, pas une découverte de phase du ladder gelé.

Le calcul reproductible est
`experiments/phase7/run_phase7_2d_surface_code_preflight.py`. Pour le tore
`3 x 3`, il donne : 18 qubits de barreau, rang stabilisateur 16, 4 états
fondamentaux, distance 3, et 1 431 Paulis non triviaux de poids 1 ou 2 tous
scalaires (ou nuls) après projection dans l'espace du code. Voir
`results/phase7/2d_surface_code_preflight.json`.

## Parent de référence exact

Un bord `e` de la maille carrée est un double puits avec deux modes de rail
`a_e,b_e`. On impose le sous-espace monomère par

\[
C_e=(n^a_e+n^b_e-1)^2,\qquad P_e=1-C_e.
\]

`C_e` est déjà un projecteur exact. Dans son noyau, les opérateurs physiques
neutres

\[
X_e=a_e^\dagger b_e+b_e^\dagger a_e,\qquad
Z_e=n^a_e-n^b_e
\]

forment les Paulis du qubit de barreau ; en espace complet,
`X_e^2=Z_e^2=P_e`. Pour une étoile `s` et une plaquette `p`, poser

\[
A_s=\prod_{e\ni s}X_e,\quad B_p=\prod_{e\in\partial p}Z_e,\quad
Q_s={P_s-A_s\over2},\quad Q_p={P_p-B_p\over2},
\]

où `P_s` et `P_p` sont les produits des `P_e` sur les mêmes supports. Le
parent de référence est

\[
H_{\rm ref}=U\sum_e C_e+J_s\sum_sQ_s+J_p\sum_pQ_p.
\]

Chaque terme est hermitien et idempotent. Deux termes de types `A` et `B`
partagent zéro ou deux barreaux, donc les deux anticommutations locales se
compensent. Les contraintes Mott commutent également avec tout terme qui les
contient. C'est ainsi un parent frustration-free à projecteurs commutants sur
le **Hilbert complet**, pas seulement après projection monomère.

Il conserve exactement la charge totale `U(1)`. Une étoile contient quatre
transferts de rail, donc elle conserve aussi les parités globales de rail
`P_a` et `P_b`; les plaquettes et contraintes de Mott sont diagonales.

Sur un tore, le parent est le code torique standard : les relations produit de
toutes les étoiles et produit de toutes les plaquettes réduisent le rang de
deux, donnant quatre états fondamentaux. Une boucle `Z` non contractile de
longueur trois est le premier opérateur logique pour le contrôle `3 x 3`.

## Ce que le test établit réellement

Toute projection `P O_S P` d'un opérateur local physique, conservant la
charge, dans le secteur à un atome par barreau est une combinaison de Paulis
sur les barreaux de `S`. Sous la distance du code, chaque Pauli soit
anticommute avec au moins un stabilisateur (projection nulle), soit appartient
au groupe stabilisateur (projection scalaire). Le contrôle énumère exactement
tous les Paulis jusqu'au poids deux et vérifie cette assertion.

Il s'agit d'un contrôle positif pour l'audit d'algèbre locale Phase 7 : il se
distingue de l'Ising/cat archivé, que l'audit rejette déjà sur un seul barreau
à cause de `X_j`.

## Ce qui reste ouvert  -  le mur pertinent

`A_s` et `B_p` sont des termes à quatre barreaux : `A_s` est un produit de
quatre bilinéaires de rail, donc un opérateur à huit opérateurs de création ou
annihilation dans l'écriture microscopique. Le dépôt ne contient **aucune**
dérivation contrôlée montrant que les hoppings corrélés, flux et médiateurs de
charge 2 d'ANTLER produisent ces termes, avec les mauvais termes bornés et une
fenêtre de gap explicite. Les spéculations de type honeycomb ou « deux canaux
médiateurs donnent XX/YY/ZZ » ne passent pas ce seuil sans un calcul de
Schrieffer--Wolff complet et des ingrédients microscopiques explicitement
listés.

Le texte de travail externe reçu le 19 juillet esquisse un no-go pour tout
parent commutant à largeur bornée et une nécessité de largeur croissante. Sa
preuve formelle, avec hypothèses et références vérifiables, n'était pas encore
déposée au moment de ce préflight. Il est donc enregistré comme direction
théorique à auditer, non comme théorème ANTLER déjà validé.

## Porte obligatoire avant toute promotion

Une proposition « ANTLER 2D native » doit fournir, avant un scan numérique :

1. les degrés de liberté supplémentaires et la géométrie 2D qui grandit dans
   les deux directions ;
2. un Hamiltonien microscopique à corps bas et la convention Jordan--Wigner
   complète ;
3. une dérivation perturbative avec petite quantité contrôlée, ordre des
   processus, coefficients, termes parasites et bornes d'erreur ;
4. les tests de symétrie `U(1)`/parités, de gap de charge et de gap neutre ;
5. l'audit de l'algèbre locale exhaustive sur des familles de taille et de
   distance croissantes.

Cette cible est **abélienne**. Elle ne démontre ni tresse non abélienne, ni
universalité, ni tolérance aux fautes, ni même une réalisation expérimentale.
