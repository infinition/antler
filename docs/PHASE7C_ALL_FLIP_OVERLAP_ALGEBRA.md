# Phase 7C  -  algèbre de recouvrement du terme all-flip

L'audit du médiateur croisé a identifié l'interaction dominante

\[
F_p=\prod_{j\in p}S_j^+ + \prod_{j\in p}S_j^- .
\]

Avant tout scan de paramètres, son algèbre de répétition doit être connue. Un
contrôle exact sur huit qubits logiques compare un support de quatre barreaux
à un second support de quatre barreaux ayant `0` à `4` barreaux communs.

| Barreaux communs | Norme spectrale de \( [F_p,F_q] \) |
| ---: | ---: |
| 0 | 0 |
| 1 | 1 |
| 2 | 1 |
| 3 | 1 |
| 4 (même support) | 0 |

Ainsi, deux termes `F_p` distincts avec un recouvrement non vide ne commutent
pas. Le résultat à deux plaquettes n'était donc pas un accident de géométrie :
une famille répétée du même terme ne peut pas être un parent stabilisateur
commutant, sauf supports disjoints.

Ce rejet est étroit : il concerne la combinaison *all-flip* réellement
engendrée par le médiateur croisé, pas toute interaction ANTLER ni tout parent
non commutant. Il interdit néanmoins de consacrer un RL à ajuster les seuls
couplages de ce motif dans l'espoir d'obtenir directement le stabilisateur 2D
commutant de référence.
