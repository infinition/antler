# Phase 8B  -  audit de sélection d'ordre impair du gadget de Gauss

## Verdict

La construction T3 de `PHASE8B_Z2_GAUSS_RESOURCE_THEOREM.md` est **réfutée
telle qu'écrite**. Cela ne réfute ni T1 (ressource de jauge indispensable),
ni T2 (no-go du gauging fin), ni toute construction à gros grain modifiée.

Le Hamiltonien déclaré est

\[
H=\Delta_G n_\mu+\lambda(\mu^\dagger A+A\mu),
\qquad A=X_L+X_R+\eta P_B,
\]

et la cible est (G_B=X_LP_BX_R). Comme chaque vertex change l'occupation de
`mu` d'une unité, une amplitude qui part et revient dans le secteur
`n_mu=0` contient un nombre **pair** de vertices. La branche basse exacte est

\[
H_{\mathrm{eff}}=\frac{\Delta_G-
\sqrt{\Delta_G^2+4\lambda^2 A^2}}{2}=f(A^2).
\]

Elle appartient donc à l'algèbre paire

\[
\operatorname{span}\{1, X_LX_R, X_LP_B, P_BX_R\},
\]

qui ne contient pas (X_LP_BX_R). En particulier, un terme
`K ~ lambda^3 eta / Delta_G^2` est interdit par cette règle de sélection.

## Calcul exact

Sur les trois qubits commutants `(X_L, P_B, X_R)`, pour
`eta=0.37`, `Delta_G=10` et `lambda/Delta_G=0.20..0.025` :

| Quantité | Résultat |
|---|---:|
| max `|c_{X_L P_B X_R}|` | `2.78e-16` |
| max résidu de reconstruction par l'algèbre paire | `2.35e-15` |
| min `|c_{X_LX_R}|`, parasite d'ordre 2 | `1.25e-2` |

Les termes d'ordre pair sont bien dans le commutant abélien, mais ils ne
fabriquent pas la contrainte de Gauss à trois facteurs. Les compenser ne peut
pas créer le terme absent.

## Conséquence constructive

Une révision viable devra déclarer une nouvelle primitive qui brise explicitement
la sélection « nombre pair de conversions de `mu` » : par exemple un couplage
microscopique à trois facteurs, une structure de médiateurs avec mélange
supplémentaire démontré, ou une autre architecture. Cette primitive doit être
écrite avant de lui attribuer un ordre SW, puis repasser les audits de
sélectivité à bas ordre.

## Claim boundary

Cet audit porte uniquement sur le gadget linéaire à un médiateur affiché dans
le théorème Phase 8B. Il ne conclut rien sur une construction T3 modifiée, une
phase appariée, une protection émergente, une fusion, une jonction T ou une
tresse non abélienne.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_gauss_odd_order_audit.py
```
