# Phase 8B  -  lemme de boucle de marcheur programmable

## Énoncé

Soit un marcheur neutre avec une position basse `0` et `ell-1` positions de
coût `Delta`, arrangées en cycle. Sur le lien `j` du cycle, le saut est
conditionné par un Pauli hermitien `O_j` sur un lien de jauge distinct. Si
`S=O_1...O_ell` est un mot Pauli hermitien, alors la branche basse contient :

\[
H_{\rm eff}=c_I I+
2(-1)^{\ell-1}\frac{\prod_j\lambda_j}{\Delta^{\ell-1}}S+
O(\lambda^{\ell+2}/\Delta^{\ell+1}).
\]

Le premier terme non scalaire est donc le tour complet du cycle, à l’ordre
`ell`. Tous les chemins fermés plus courts rebroussent nécessairement chemin
sur au moins une arête et se réduisent à des facteurs `O_j^2=I`.

## Preuve perturbative

Dans `P=|0><0|`, un terme de Schrieffer--Wolff doit commencer et finir à la
position `0`. Pour une marche de longueur strictement inférieure à `ell`, le
cycle ne peut être parcouru entièrement ; tout chemin fermé contient donc une
paire de sauts inverses. Les opérateurs associés s’annulent deux à deux et le
mot de jauge est scalaire. À la longueur `ell`, les deux orientations du cycle
donnent le mot `S` et son adjoint, identiques puisque `S` est hermitien. Les
`ell-1` dénominateurs virtuels donnent le coefficient ci-dessus.

Cette preuve suppose les liens du mot sur des qubits distincts. Elle ne
prétend pas contrôler le crosstalk entre plusieurs marcheurs se chevauchant.

## Contrôles exacts archivés

| boucle | mot obtenu | dimension | puissance mesurée |
|---|---|---:|---:|
| 4 pas | étoile/plaquette | 1 024 (chevauchement joint) | `4.0365` |
| 5 pas | pentagone `XZXZX` | 160 | `5.0136` |
| 5+4 pas | pentagone + plaquette adjacente | 2 560 | `5.0502`, `4.0365`, produit `9.0613` |

Les audits de chevauchement correspondants ne résolvent aucun mot Pauli hors
de l’algèbre générée par les stabilisateurs cibles.

## Contrat de compilation pour un twist

Pour une géométrie de coupure de branche explicitement définie :

1. associer un marcheur à chaque stabilisateur étoile, plaquette, mixte ou
   pentagonal ;
2. vérifier le critère symplectique pair sur chaque support partagé ;
3. auditer tous les motifs de chevauchement locaux, puis un patch fini avec
   tous les marcheurs ;
4. seulement après ces portes, mesurer le sous-espace de fusion et les
   holonomies de déformation.

Le lemme fournit la brique de compilation locale. Il ne fournit ni la
géométrie complète de twist, ni l’origine matérielle des sauts
phase-contrôlés, ni une tresse non abélienne.

## Claim boundary

Les marcheurs `X/Z`-conditionnés sont de nouvelles primitives déclarées,
extérieures au ladder ANTLER gelé. La validité de l’expansion est locale et en
régime profond. Une architecture complète doit encore démontrer le pavage
microscopique, le gap, la localité du code, la fusion et le braid.
