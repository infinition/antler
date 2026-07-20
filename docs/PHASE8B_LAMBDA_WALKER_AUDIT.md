# Phase 8B  -  audit exact du marcheur Lambda réparé

## Verdict

La correction de la règle de sélection est cohérente avec le gadget du gadget à un seul
médiateur : un marcheur neutre à états internes crée une **boucle** virtuelle,
et non un aller-retour sur un vertex unique. Au niveau algébrique du marcheur,
la boucle engendre exactement le mot de Gauss grossier

\[
G_B=X_L(1-2n_{a,1})(1-2n_{a,2})X_R.
\]

Ce résultat qualifie seulement une nouvelle primitive déclarée  -  le saut de
marcheur conditionné par la densité. Il ne dérive pas cette primitive de la
grammaire ANTLER existante et n'établit aucune phase ni tresse.

## Gadget audité

Le secteur comporte `|vac>, |mu0>, |mu1>, |mu2>` avec au plus un marcheur
neutre. Les vertices sont :

\[
\lambda X_L:|vac\rangle\leftrightarrow|\mu_0\rangle,
\quad w(1-2n_{a,1}):|\mu_0\rangle\leftrightarrow|\mu_1\rangle,
\]
\[
w(1-2n_{a,2}):|\mu_1\rangle\leftrightarrow|\mu_2\rangle,
\quad \lambda X_R:|\mu_2\rangle\leftrightarrow|vac\rangle.
\]

Le parcours de la boucle complète a quatre vertices et porte le produit de
leurs signes : (X_Lp_1p_2X_R). Les parcours rebroussés sont des carrés
d'involutions et donnent des scalaires.

## Résultat exact

Pour `Delta=10`, couplages égaux et `lambda/Delta=0.20..0.025` :

| Mesure | Valeur |
|---|---:|
| coefficient cible à `lambda/Delta=0.025` | `-7.7833e-6` |
| pente log-log sur tout le scan | `3.903` |
| pente profonde (`<=0.075`) | `3.973` |
| maximum de tout coefficient Walsh non scalaire hors cible | `2.09e-17` |

La séparation entre les secteurs `G=+1` et `G=-1` vaut `2|K|`, comme attendu.
L'absence de parasite non scalaire est exacte dans ce modèle à un marcheur,
car le seul invariant de signe de la boucle est le produit des quatre
vertices.

## Ce qui reste bloqué

1. Le vertex `mu_j† mu_(j-1)(1-2n_a)` est une **nouvelle primitive**. Sa
   classe de localité est déclarée mais sa réalisation par ANTLER n'est pas
   dérivée.
2. La formulation « une primitive » ne vaut pas encore preuve du budget
   matériel initial « au plus un nouveau DOF local » : le marcheur doit être
   spécifié comme un unique ancilla multi-niveaux, distinct des qubits de
   jauge de frontière, ou le budget doit être révisé explicitement.
3. Le bloc n'inclut pas encore le transport de paires inter-blocs, les
   canaux Floquet, la dynamique de matière complète, ni l'audit local du
   sous-espace physique.

## État de la théorie Phase 8B

- T1/T2 : théoriques, non ré-auditées ici.
- Gadget à un médiateur : réfuté exactement.
- Couche interne du marcheur Lambda : validée algébriquement et numériquement.
- Construction microscopique native, phase protégée, fusion et tresse : non
  établies.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_lambda_walker_gadget_audit.py
```
