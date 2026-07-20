# Phase 7C  -  stress-test de recouvrement du médiateur croisé

## Question testée

Une interaction locale ne peut devenir la brique d'un code 2D que si elle
survit lorsqu'elle est répétée. Ce contrôle construit deux motifs croisés sur
les supports `(0,1,2,3)` et `(2,3,4,5)`. Ils partagent les barreaux `2,3`.
Les contre-termes `Z/ZZ` de la calibration à une plaquette sont additionnés,
sans nouveau fit global.

Ce support à six barreaux est un stress-test compact. Il ne représente ni une
maille complète de surface code, ni une géométrie de tresse.

## Contrôles numériques

À `U=20`, `Delta=10`, `g=0.5`, la diagonalisation exacte est effectuée dans
les secteurs de parité de rail, afin de ne pas perdre de vecteurs dans une
dégénérescence à énergie nulle. Le bloc charge-fixée a dimension `3304` et le
sous-espace monomère a dimension `64`.

Après addition des contre-termes :

- capture monomère minimale : `0.99979999` ;
- gap vers le premier état hors sous-espace : `39.94246` ;
- coefficients des deux composantes `XXXXII` et `IIXXXX` :
  `-2.03055e-7` chacune ;
- parité de rail nue : résidu de commutateur nul à la précision numérique.

Le recouvrement ne détruit donc pas le mécanisme local ni l'isolation Mott.

## La combinaison effectivement générée

Les composantes `XXXX`, `YYYY` et les six mots avec deux `Y` forment, avec
leurs signes observés, l'opérateur

\[
F_p=\prod_{j\in p}S_j^+ + \prod_{j\in p}S_j^- .
\]

Après compensation, la projection de l'Hamiltonien effectif sur `F_p` explique
`0.94878` de sa norme sur ce stress-test, avec coefficient
`-1.62444e-6` pour chacune des deux plaquettes. Le médiateur croisé fabrique
donc principalement ce *all-flip* collectif, et non le stabilisateur Pauli
pur `XXXX`.

## Mur de commutativité

Pour les deux supports qui se chevauchent, la norme spectrale exacte vaut

\[
\lVert[F_{(0,1,2,3)},F_{(2,3,4,5)}]\rVert_2=1.
\]

La non-commutativité est d'ordre un. Ce n'est pas une tresse ni une ressource
topologique : elle empêche au contraire d'interpréter ces deux opérateurs comme
des stabilisateurs simultanément imposables dans cette géométrie.

## Décision

Le recouvrement valide seulement la compatibilité structurelle minimale du
motif et identifie sa vraie interaction dominante. La géométrie à recouvrement
deux-barreaux est **rejetée comme parent stabilisateur commutant**. Les étapes
suivantes sont de dériver le coefficient Schrieffer--Wolff de `F_p`, puis de
chercher une autre géométrie ou un autre degré de liberté où les termes
effectifs pertinents commutent et restent sélectifs.
