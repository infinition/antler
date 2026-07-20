# Phase 7  -  Audit adversarial de la dérivation état-vers-Hamiltonien

Audit à charge de `docs/PHASE7_STATE_TO_HAMILTONIAN_DERIVATION.md`, mené par
calcul explicite. Chaque point est établi ou réfuté ; les statuts corrigés
sont reportés dans le document principal (marqués « corrigé par audit »).

## Verdict global (avant le détail)

1. **Le code n'est pas localement indistinguable dans l'algèbre physique
   complète** : `P X_j P = X̄` (non scalaire) pour *tout* barreau `j`, bord ou
   bulk. Seule l'indistinguabilité **sous l'algèbre symétrique** (opérateurs
   préservant `P_a, P_b`) est exacte. Cas A : **faux**. Cas B : **vrai**.
2. **« Mode de bord » retiré au point fixe** : `[H_fix, X_j] = 0` pour tout
   `j`. `X_j` est le **paramètre d'ordre local** d'une brisure spontanée Z₂ ;
   la GSD 2 en OBC **et** PBC est la signature SSB, pas une dégénérescence de
   bord. Une localisation de bord n'existe que sous perturbation (`h > 0`),
   avec un diagnostic quantitatif bord/bulk donné au §2.
3. **`C_j` n'est pas un projecteur** (`C_j² ≠ C_j`). « Commuting-projector au
   sens strict » était faux tel qu'énoncé. Réparation exacte : réécriture
   spectrale en projecteurs commutants explicites, ou projecteur minimal
   `K_j = 1_{q_cell ≠ 1}` (formule locale au §3). Les preuves de gap
   survivent inchangées.
4. **Origine microscopique requalifiée** : un seul canal médiateur de charge
   2 (6D) ne produit que la moitié de `X_jX_{j+1}` ; l'autre moitié est un
   superéchange par doublon (natif mais non médié). Les **trois canaux 6F ne
   sont pas utilisés**. À `θ = π`, le signe du superéchange s'inverse
   (calcul de chemins au §4) : le parent sur la ligne de matching est
   `+J̃ Σ Y_jY_{j+1}` (Ising en base Y, staggered), pas `−J̃ Σ XX`. Le
   document principal affirmait le contraire : **corrigé**.
5. Le théorème U(1) tient, avec ses hypothèses précisées (§5) ; il interdit
   les modes de bord *chargés* quasi conservés, et rien de plus ; il ne
   démontre pas l'unicité du choix « degré relatif des rails ».
6. **Reclassement final : mémoire Ising/Z₂ sous symétrie exacte (SSB de la
   parité relative, paramètre d'ordre local P-impair).** Rejet comme « code
   localement indistinguable » et comme candidat topologique. Conservé comme
   benchmark de protection sous symétrie et contrôle positif du harnais
   `ε_edge`.

---

## 1. Localité logique réelle  -  calcul exact de P O P

`P = |Ω_+⟩⟨Ω_+| + |Ω_−⟩⟨Ω_−|`. Actions par barreau (`|→⟩, |←⟩` propres de
`X` ; `Z|→⟩ = |←⟩`, `Y|→⟩ = −i|←⟩`, `n^a|→⟩ = (|→⟩+|←⟩)/2`) :

`X_j|Ω_±⟩ = ±|Ω_±⟩` ; `Z_j|Ω_±⟩, Y_j|Ω_±⟩ ∝ |Ω_∓ sur un seul barreau⟩ ⊥ 𝒞`
(un seul barreau retourné, orthogonal aux deux états pour `L ≥ 2`) ;
`n^{a}_j|Ω_±⟩ = ½|Ω_±⟩ + ½(barreau retourné) `.

D'où, exactement, pour tout `j` et tout `L ≥ 2` :

| opérateur | base `{Ω_+, Ω_−}` | base parité `{+, −}` |
|---|---|---|
| `P X_j P` | `diag(+1, −1)` | `[[0,1],[1,0]] = X̄` |
| `P Y_j P` | `0` | `0` |
| `P Z_j P` | `0` | `0` |
| `P n^a_j P` | `½ I` | `½ I` |
| `P n^b_j P` | `½ I` | `½ I` |

**Question décisive : `P X_j P` est non scalaire  -  oui, de norme 1, pour tout
`j`.** `X_j = a†_j b_j + b†_j a_j` est un opérateur physique local à un
barreau (c'est le hopping de barreau `J_perp` du modèle gelé). Il implémente
`X̄` exactement, depuis n'importe quel barreau.

Distinction sans ambiguïté :

- **A. Indistinguabilité face à tous les opérateurs locaux : FAUX.**
  Réfutée par le calcul ci-dessus.
- **B. Indistinguabilité face aux opérateurs locaux préservant `P_a, P_b` :
  VRAI.** (Lemme 2 du doc principal : pour `[O, P_b] = 0` de support strict,
  `⟨+|O|+⟩ = ⟨−|O|−⟩` et `⟨+|O|−⟩ = 0` ; donc `P O P` scalaire, exactement.)

**Reclassement obligatoire** : mémoire Ising/Z₂ protégée par symétrie, PAS
code topologique localement indistinguable. Toute affirmation
« indistinguabilité locale » du doc principal est restreinte au cas B.
Conséquence physique : un bain qui se couple à `X_j` (transfert croisé
résiduel `ε_⊥`) produit une erreur logique `X̄` **au premier ordre, depuis
n'importe quel site**  -  taux ∝ `ε_⊥` par site, ∝ `ε_⊥ L` au total. La
protection est intégralement conditionnelle à l'exactitude de la symétrie
dans le Hamiltonien **et dans le couplage au bain**.

## 2. « Mode de bord » ou ordre de bulk ?

### 2.1 Au point fixe : ordre de bulk, sans ambiguïté

Pour tout `j` (pas seulement `j = 0`) :
`[X_j, Π^B_{j−1}] = [X_j, Π^B_j] = 0` (algèbre `X N = X`, `[X_j, X_jX_{j±1}] = 0`),
`[X_j, C_k] = [X_j, n^d] = 0` (X conserve toutes les charges de cellule). Donc

```
[H_fix, X_j] = 0    pour j = 0, 1, …, L−1,
```

et tous les `X_j` ont la même action logique `X̄` (§1). **`X_0` n'a aucun
statut de bord** : `X_j` est le paramètre d'ordre local d'un point d'Ising
classique (dans la base x), avec `L` opérateurs P-impairs localement
conservés  -  structure de point fin-réglé, pas de mode localisé. Les
expressions « mode de bord » et « quasi-conservation de bord » sont
**retirées au point fixe**.

### 2.2 GSD 2 en OBC et PBC = signature SSB

Une dégénérescence de bord véritable (chaîne de Kitaev fermionique) donne
OBC 2 / PBC 1 : les modes de bord disparaissent sans bord. Ici la
dégénérescence est **insensible aux conditions aux bords** parce qu'elle
vient de deux états de brisure `|Ω_±⟩` reliés par la symétrie globale
(`P_b|Ω_+⟩ = |Ω_−⟩`), avec paramètre d'ordre local `⟨X_j⟩ = ±1`. C'est la
définition opérationnelle de la brisure spontanée Z₂. Le doc principal
tirait la bonne conclusion (« pas d'ordre topologique intrinsèque ») sans en
tirer la conséquence terminologique : c'est fait ici.

### 2.3 Ce qui reste vrai : localisation de bord uniquement sous perturbation

Sous `h Σ Z_j` (`r = 2h/J < 1`), le point fin-réglé se résout :

- **non habillé** : `[H, X_j] = 2ih Y_j` ⇒ `ε(X_j) ≈ 2h` pour tout `j`,
  bord comme bulk  -  aucune distinction au support 1 ;
- **habillé, bord** : la récurrence de mode zéro fort converge (un seul lien
  adjacent) : `ε(Γ^{(w)}_{bord}) = 2J_I r^w √((1−r²)/(1−r^{2w}))`  -  décroît
  exponentiellement avec le support `w` ;
- **habillé, bulk** : la même récurrence lancée sur un site de bulk branche
  des deux côtés : le lien gauche génère au premier pas le résidu
  `2i(h) X_{j−1}Y_jX_{j+1}` (coefficient `∝ h`, non réductible par
  géométrie) ; en langage fermion libre, le spectre BdG n'a **aucun mode
  sous-gap de bulk**  -  prédiction : `ε_min(bulk, w) = Θ(h)`, **plat en `w`**.

**Diagnostic de localisation quantitatif** (à auditer numériquement) :

```
ε_min(site 0, support w) / ε_min(site L/2, support w)  ≈  (r^w · c₁) / (c₂ h/J)  → 0
```

C'est le seul sens dans lequel « bord » est défendable, et il n'existe qu'à
`h ≠ 0`. Par ailleurs même à `h > 0` le paramètre d'ordre de bulk persiste :
`P X_j P ≈ ± m X̄` avec `m = (1−r²)^{1/8} ≈ 0.995` à `r = 0.2`  -  l'action
logique locale de bulk ne disparaît jamais dans la phase ordonnée.

## 3. Correction algébrique des « projecteurs »

### 3.1 Réfutation

`q_cell(j) ∈ {0,1,2,3,4}` (bulk ; `{0..3}` au bord). Valeurs propres de
`C_j = (q_cell−1)²` : `{1, 0, 1, 4, 9}`. Donc

```
C_j² = (q_cell−1)⁴  a pour valeurs  {1, 0, 1, 16, 81}  ≠  C_j .
```

**`C_j` n'est pas idempotent. « H_fix est un parent commuting-projector au
sens strict » était faux tel qu'énoncé.** Ce qui reste vrai sans
modification : tous les termes sont PSD, commutent deux à deux, et
annihilent les cibles (commuting frustration-free PSD parent).

### 3.2 Réparations exactes (deux, au choix)

**(a) Réécriture spectrale.** Avec `Π^{S}_j = 1_{q_cell(j) ∈ S}` (diagonaux,
idempotents, commutant entre eux et avec tout `H_fix`) :

```
C_j = Π^{{0,2}}_j + 4 Π^{{3}}_j + 9 Π^{{4}}_j ,
```

donc `H_fix` est une combinaison **positive** de projecteurs locaux
mutuellement commutants `{Π^{{0,2}}_j, Π^{{3}}_j, Π^{{4}}_j, n^d_j, Π^B_j}`  - 
la forme « somme de projecteurs commutants » est restaurée avec des
projecteurs explicites, au prix de coefficients `{U, 4U, 9U, Δ, J}`.

**(b) Projecteur minimal.** `K_j = 1_{q_cell(j) ≠ 1}`. Formule locale exacte
(fonctions symétriques élémentaires `e_k` des quatre occupations
`S_j = {n^a_j, n^b_j, n^d_{j−1}, n^d_j}`, idempotentes) :

```
1_{q_cell=1} = e₁ − 2e₂ + 3e₃ − 4e₄ ,        K_j = 1 − e₁ + 2e₂ − 3e₃ + 4e₄ .
```

Vérification exhaustive (`k` modes occupés) : `k−2C(k,2)+3C(k,3)−4C(k,4) =
0,1,0,0,0` pour `k = 0..4`. ✓ `K_j² = K_j`, diagonal, support = barreau `j`
+ ses deux liaisons (termes de densité jusqu'à 4 corps), `ker K_j = ker C_j`,
commute avec tous les autres termes (`X` et le vertex `g` conservent chaque
`q_cell`).

### 3.3 Preuves de gap refaites avec les bons opérateurs

Les deux variantes (`U ΣC_j` spectralement décomposé, ou `U ΣK_j`) donnent
un Hamiltonien = somme de projecteurs locaux commutants à coefficients
positifs ⇒ spectre additif sur les secteurs joints. Alors, exactement, ∀L :

- **noyau** : inchangé (`ker C_j = ker K_j`) ⇒ Lemme 3 intact, GSD = 2 ;
- **gap neutre** (`Q = L`) : `min(J, Δ)`  -  mur de domaine (`Π^B = 1`, coût
  `J`) ou médiateur réel (coût `Δ`) ; inchangé ;
- **gap de charge** : `Q = L ± 1` ⇒ `Σ_j (q_cell,j − 1) = ±1` ⇒ au moins une
  cellule violée ⇒ `E ≥ U` ; atteint (une particule ajoutée/retirée sur un
  barreau) ⇒ `E₀(L±1) − E₀(L) = U` exactement, **dans les deux variantes** ;
- **déviations multiples** : `E ≥ U·m` pour la variante `C` (convexité de
  `(q−1)²`) ; `E ≥ U·⌈m/3⌉` pour la variante `K` (une cellule absorbe au
  plus 3 unités  -  et en pratique les cellules voisines sont violées aussi
  par les médiateurs partagés). Fenêtre grand-canonique : `|μ| < U` (variante
  C), `|μ| < U/3` garanti (variante K).

Statut corrigé : **parent frustration-free à termes commutants ; forme
commuting-projector stricte disponible via (a) ou (b).**

## 4. Origine microscopique native  -  décomposition stricte

### 4.1 Qui produit quoi

```
Π^B_j = ½( N_jN_{j+1} − X_jX_{j+1} ) ,
X_jX_{j+1} = [a†_jb_j a†_{j+1}b_{j+1} + h.c.]   (canal croisé, ΔN_a = ±2)
           + [a†_jb_j b†_{j+1}a_{j+1} + h.c.]   (canal échange, ΔN_a = 0).
```

| élément | statut |
|---|---|
| `Π^B_j` (et l'état cible) | **parent cible abstrait**, dérivé du sous-espace cible  -  sa réalisation microscopique est une obligation séparée |
| canal croisé + diagonal `P^{aa}+P^{bb}` | obtenu par **un** médiateur de charge 2 par liaison (vertex 6D), ordre `g²/Δ` |
| canal échange + diagonal `P^{ab}+P^{ba}` | superéchange par **doublon virtuel** (hoppings intra-rail natifs, coût `E_v ≈ 2U`)  -  natif, mais **pas un processus de médiateur de charge 2** |
| les trois canaux 6F | **non utilisés**. La factorisation 6F concernait la liaison du parent Iemini externe ; Phase 7 n'emploie qu'un canal médiateur (combinaison `aa+bb`) plus un canal doublon neutre. Aucune filiation 6F ne doit être revendiquée. |

**Un médiateur par lien ne suffit pas** : à `t = 0`, le modèle effectif est
le point de Heisenberg caché (gapless, réfutation analytique du voisinage
6E). Le parent n'est atteint que sur la ligne de matching
`2t²/E_v = g²/Δ`, `V_leg = 0`  -  deux ingrédients, une condition de réglage.
La formule « aucun terme inséré à la main » est requalifiée : *l'état et les
projecteurs sont dérivés ; la réalisation microscopique exige un réglage de
couplages (matching), et le canal d'échange n'est pas médié par charge 2.*

Contrainte supplémentaire vs modèle gelé : la ligne de matching exige un
superéchange uniforme, donc soit `J1 = J2 = t` (hopping non dimérisé  -  écart
déclaré à l'échelle SSH gelée), soit un médiateur dimérisé
`g_i² = 2 J_i² Δ / E_v` liaison par liaison.

### 4.2 Signes rung-major à θ = π : correction de fond

Calcul de chemins pour le superéchange sur `|a@j, b@j+1⟩ → |b@j, a@j+1⟩`
(cordes JW évaluées sur les états intermédiaires) :

- chemin 1 : `a†_{j+1}a_j` (corde sur `2j+1`, vide : signe +) puis
  `b†_j b_{j+1}` (corde sur `2j+2`, **occupée** : signe −) ;
- chemin 2 : `b†_j b_{j+1}` (corde vide : +) puis `a†_{j+1}a_j`
  (corde sur `2j+1`, **occupée** : −).

Produit des éléments `⟨f|V|m⟩⟨m|V|i⟩ = −t²` par chemin, somme `−2t²`,
`H_eff = −(1/E_v)(−2t²) = **+2t²/E_v**` : le superéchange **change de signe**
à `θ = π` (contre `−2t²/E_v` à `θ = 0`). Le canal médiateur, lui, garde son
signe (les cordes interviennent au carré : `−g²/Δ` aux deux θ).

Conséquence sur la ligne de matching à `θ = π` :

```
H_eff = −J̃(σ⁺σ⁺ + h.c.) + J̃(σ⁺σ⁻ + h.c.) + diag = + J̃ Σ Y_jY_{j+1} + const .
```

Le parent au point fixe fermionique est un **Ising en base Y,
antiferromagnétique** ; états cibles = produits staggered d'états propres de
`Y_j` ; `Π^B(θ=π)_j = ½(N_jN_{j+1} + Y_jY_{j+1})`. Équivalent au cas `θ = 0`
par unitaire on-site (`Π_{j impair} Z_j` puis rotation `e^{−iπ/4 Σ Z_j}`),
mais **le doc principal affirmait le même parent XX aux deux θ via
`d → −d` : faux, corrigé.** La famille Phase 7 est définie canoniquement à
`θ = 0` ; la version `θ = π` requiert la redéfinition ci-dessus.

## 5. Théorème U(1)  -  hypothèses et portée exactes

**Hypothèses spectrales** : (i) `[H, Q] = 0` ; (ii) `P` = projecteur
spectral sur le multiplet fondamental, entièrement contenu dans un secteur
`Q₀` ; (iii) multiplet exactement dégénéré d'énergie `E₀` ; (iv) `O` de
charge définie `q ≠ 0` (`[Q,O] = qO`) ; (v)
`Δ_c(q) = min spec(H|_{Q₀+q}) − E₀ > 0`. Conclusion : `ε_edge(O) ≥ Δ_c(q)`.

**Cas quasi-dégénéré** : si le multiplet a un étalement `δ`
(`‖(H−Ē)P‖ ≤ δ`), la décomposition
`(1−P)[H,O]P = (1−P)(H−Ē)OP − (1−P)O(H−Ē)P` donne

```
ε_edge(O) ≥ Δ_c(q) − δ .
```

**Opérateurs de charge mixte** : `O = Σ_q O_q` (harmoniques
`O_q = (2π)^{−1}∫dθ e^{−iqθ} e^{iθQ} O e^{−iθQ}`). La borne s'applique à
chaque harmonique chargée ; **l'harmonique neutre `O₀` n'est contrainte en
rien**.

**Ce que le théorème interdit exactement** : un opérateur (ou une composante
d'opérateur) *chargé* qui soit à la fois quasi conservé et d'action non
triviale sur un multiplet fondamental de charge définie, dans un système à
gap de charge. **Rien de plus.** Il n'interdit ni action logique neutre, ni
opérateurs chargés en système compressible, ni quoi que ce soit sur la
dynamique.

**Pourquoi il ne sélectionne pas le degré relatif des rails** : le théorème
élimine les candidats chargés ; le choix de `P_rel` comme algèbre logique
est ensuite un **choix de construction** guidé par l'inventaire des
opérateurs neutres disponibles (neutres P-impairs : `X_j, Y_j` et habillages ;
neutres P-pairs : inertes sur un doublet de parité par le Lemme 2). D'autres
structures neutres (nombres quantiques spatiaux, symétries émergentes) ne
sont pas exclues par le théorème. L'unicité n'est **pas démontrée**  -  le doc
principal la présentait comme « tout en découle » : requalifié en heuristique
de conception.

## 6. Tableau corrigé (L = 4, 6, 8, 10) et verdict

Point fixe `H_fix` (`U = 4, Δ = 2, J = 1`), secteur `Q = L`. Valeurs exactes
(tolérance numérique `1e−10`) ; lignes `h = 0.1` = prédictions dérivées
(§2.3), `r = 2h/J = 0.2`, `J_I = J/2`.

| observable | L=4 | L=6 | L=8 | L=10 |
|---|---|---|---|---|
| GSD OBC | 2 | 2 | 2 | 2 |
| GSD PBC | 2 | 2 | 2 | 2  -  **signature SSB, pas bord** |
| gap neutre | 1.0 | 1.0 | 1.0 | 1.0 |
| gap de charge | 4.0 | 4.0 | 4.0 | 4.0 |
| `max_j ‖P X_j P − ½tr(PX_jP)·I‖` | 1.0 | 1.0 | 1.0 | 1.0  -  **∀j, non scalaire** |
| même test, opérateurs locaux préservant `P_a,P_b` (support ≤ 2 barreaux) | 0 | 0 | 0 | 0 |
| `ε(X_0)` (bord, h=0) | 0 | 0 | 0 | 0 |
| `ε(X_{L/2})` (bulk, h=0) | 0 | 0 | 0 | 0  -  **aucune localisation au point fixe** |
| `ε(Γ^{(w)}_{bord})`, h=0.1, w=1..4 | 0.200 / 0.0392 / 0.00784 / 0.00157 | idem | idem | idem (plat en L) |
| `ε_min(bulk, w)`, h=0.1 | ≈ 0.2 plat en w | idem | idem | idem (prédiction Θ(h)) |
| `P X_j P` sous h=0.1 | `± m X̄`, `m ≈ 0.995` | idem | idem | idem  -  action logique locale de bulk persistante |

**Verdict : « SPT/Ising sous symétrie »  -  plus précisément mémoire Ising/Z₂
par brisure spontanée de la parité relative, protégée uniquement par
l'exactitude microscopique de `P_a, P_b` (Hamiltonien ET bain). « Code
localement indistinguable » : REJET. Candidat topologique : REJET.**

Usages légitimes restants : (i) benchmark de protection sous symétrie avec
gates falsifiables (ce tableau) ; (ii) contrôle positif/négatif du harnais
`ε_edge` (bord vs bulk sous `h`) ; (iii) support du théorème U(1), qui
survit à l'audit et contraint toute construction Phase ≥ 7 : un futur
candidat devra exhiber `P O P` non scalaire **avec** `max_j ‖P X_j P −
scalaire‖ → 0` pour tous les opérateurs locaux physiques  -  ce que la
présente construction ne fait pas.
