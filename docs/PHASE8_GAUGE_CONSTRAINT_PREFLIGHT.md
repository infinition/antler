# Phase 8  -  préflight de contrainte de Gauss locale

## But

La route vers une tresse non abélienne exige d'abord que le tunneling de
barreau nu, qui est `U(1)`-conservant mais fatal pour la parité de rail, soit
physiquement interdit plutôt que seulement petit. Ce préflight fixe le test
algébrique minimal que devra satisfaire toute extension ANTLER avec ancilla de
jauge.

Il s'agit d'un contrôle abstrait : l'ancilla neutre est insérée, elle n'est
pas dérivée des médiateurs de charge 2 existants.

## Bloc minimal

Dans le secteur de charge de rail fixée `N=1`, de base `|a>, |b>`, ajouter un
ancilla neutre à deux niveaux. Définir

\[
G=(-1)^{N_a}X_g,\qquad H_G=\frac{\lambda}{2}(1-G).
\]

La contrainte est locale et conserve la charge. Elle admet la décomposition

\[
H_G=\frac{\lambda}{2}(1-X_g)+\lambda n_aX_g,
\]

qui identifie la nouvelle primitive matérielle à dériver : un couplage local
densité–ancilla neutre, plus le contrôle de l'ancilla.

## Résultat exact

Pour `lambda=3` :

| Test | Valeur |
|---|---:|
| dimension du secteur physique `G=+1` | 2 |
| gap de contrainte | 3 |
| `||G²-1||` | 0 |
| `||[H_G,G]||` | 0 |
| `||{c_a†c_b+h.c.,G}||` | 0 |
| norme projetée du tunnel nu | `2.24e-17` |
| `||[Z_g(c_a†c_b+h.c.),G]||` | 0 |
| norme non scalaire du tunnel habillé projeté | 1 |

La première ligne utile est positive : le tunnel nu est non invariant de
jauge et sa projection sur le secteur physique est nulle. La dernière ligne
est le garde-fou : le bloc unique n'est **pas** un qubit topologique, car un
opérateur local invariant de jauge peut encore agir logiquement.

## Porte d'implémentation ANTLER

Avant toute jonction T, toute proposition doit :

1. déclarer l'ancilla neutre ou son équivalent microscopique ;
2. dériver `n_a X_g` (ou une contrainte équivalente) du Hamiltonien natif ;
3. reproduire l'algèbre ci-dessus sur son bloc de Fock à charge pondérée ;
4. construire plusieurs sommets puis faire l'audit exhaustif de
   l'indistinguabilité locale sur le sous-espace de code.

## Claim boundary

Ce contrôle ne dérive ni l'ancilla, ni une phase, ni un code multi-sommets,
ni une fusion, ni une tresse. Il établit seulement la différence opérationnelle
entre une contrainte de jauge capable de bannir le tunnel nu et une protection
topologique effectivement obtenue.

## Reproduction

```powershell
python experiments/phase7/run_phase8_local_gauge_constraint_preflight.py
```
