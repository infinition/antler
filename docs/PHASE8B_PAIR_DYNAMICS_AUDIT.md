# Phase 8B  -  dynamique de navette de paire sous Gauss

## Verdict

Le transfert de paire entre les deux blocs est dynamique, cohérent et reste
dans le secteur physique. C'est la première démonstration de propagation,
plutôt qu'un simple test de commutateur, pour la route Phase 8B.

## Protocole

Deux particules du rail `a` commencent dans le premier bloc. Le seul terme de
frontière est le transfert de paire; les marcheurs sont initialement vides et
les trois qubits de frontière sont dans le secteur `X=+1`. L'évolution exacte
est suivie jusqu'à `pi/J_pair`.

## Résultats

| `g/Delta` | population maximale de la paire cible | population virtuelle max du marcheur | fuite max de Gauss |
|---:|---:|---:|---:|
| 0.10 | 0.96112 | 0.14633 | `4.44e-13` |
| 0.05 | 0.99007 | 0.03143 | `7.63e-13` |
| 0.025 | 0.99755 | 0.00932 | `3.49e-13` |

Le maximum de transfert survient au temps `pi/(2 J_pair)`, comme pour une
navette cohérente à deux niveaux. La réduction de population de marcheur avec
la profondeur SW est cohérente avec son rôle virtuel.

## Claim boundary

Une unique paire est propagée sous des termes de marcheur et de transfert
insérés. Cela ne démontre ni une phase à plusieurs paires, ni un gap
thermodynamique, ni l'indistinguabilité locale, ni une implantation native,
ni une jonction T ou une tresse non abélienne.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_two_block_pair_dynamics.py
```
