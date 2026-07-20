# Phase 7  -  Construction analytique état-vers-Hamiltonien (parent natif à médiateurs de charge 2)

Namespace : **Phase 7 uniquement**. Aucune modification, aucune revendication
rétroactive sur le modèle gelé (`antler/model.py`, convention JW rung-major).
Le présent document est une dérivation algébrique autoportante, directement
traduisible en tests matriciels via `antler/phase7_parent_audit.py`.

> **STATUT CORRIGÉ PAR AUDIT ADVERSARIAL** (`docs/PHASE7_ADVERSARIAL_AUDIT.md`) :
> le résultat est une **mémoire Ising/Z₂ par brisure spontanée de la parité
> relative, protégée uniquement par symétrie exacte**  -  PAS un code
> localement indistinguable, PAS un candidat topologique. `P X_j P = X̄` est
> non scalaire pour **tout** barreau `j` (opérateur logique local partout) ;
> au point fixe `[H_fix, X_j] = 0` ∀j : `X_j` est un **paramètre d'ordre de
> bulk**, pas un mode de bord ; `C_j` n'est pas idempotent (« commuting-
> projector strict » réparé par réécriture spectrale, audit §3) ; à `θ = π`
> le parent de la ligne de matching est `+J̃ΣY_jY_{j+1}` (staggered), pas
> `−J̃ΣXX`. Les sections ci-dessous portant ces claims sont annotées.

Résumé en une phrase : on démontre d'abord un **no-go** (aucun opérateur de
bord *chargé* ne peut être quasi-conservé dans un bulk incompressible), ce qui
réoriente la construction vers des **opérateurs logiques neutres** agissant
sur le degré de liberté *relatif* des deux rails ; on dérive, à partir d'un
sous-espace cible explicite, un parent **frustration-free à termes
commutants, incompressible**, dont la famille microscopique native
(médiateur de charge 2 + superéchange doublon, sur une ligne de matching)
converge vers lui au second ordre. La dérivation **explique analytiquement
l'échec du scan 6E** : le point « médiateur pur » est un point de Heisenberg
caché, donc sans gap. L'audit reclasse l'objet : mémoire Ising sous symétrie
(SSB), utile comme benchmark et contrôle de harnais, rejetée comme candidat
topologique.

---

## 1. Espace de Hilbert local, charges, parités, opérateurs élémentaires

### 1.1 Modes

Échelle à deux rails, `L` barreaux (rungs) `j = 0..L-1`, liaisons (bonds)
`j = 0..L-2`.

- Rail `a` : modes hard-core `a_j`, indice linéaire `2j` ;
- rail `b` : modes hard-core `b_j`, indice linéaire `2j+1` ;
- médiateur de liaison : mode hard-core `d_j` sur la liaison `(j, j+1)`,
  indice `2L + j`, **charge U(1) égale à 2**.

C'est exactement l'indexation de `antler/native_charge2_ladder.py`
(`mode_a`, `mode_b`, `mode_d`). Nombre de modes : `3L - 1`.

### 1.2 Statistique et cordes de Jordan–Wigner

Les modes `a, b` suivent la convention ANTLER rung-major
`a_k = exp(i θ Σ_{m<k} n_m) b^{(B)}_k` sur l'ordre linéaire `k = 2i + σ`.
Le médiateur `d` est **défini** comme un mode hard-core de charge paire,
JW-neutre : il commute avec toutes les cordes (objet de charge 2, statistique
d'échange triviale avec les modes de rail). C'est une définition du modèle
Phase 7, pas une dérivation.

**Lemme 1 (localité des opérateurs de paire  -  restriction sur θ).**
Pour `k < l`, en notant `N_{<k}(s) = Σ_{m<k} n_m(s)` et
`N_mid(s) = Σ_{k<m<l} n_m(s)` :

```
a_k a_l = exp( i θ [ 2 N_{<k} + 1 + N_mid ] ) b^(B)_k b^(B)_l .
```

*Preuve.* `a_l|s⟩ = e^{iθ N_{<l}(s')} b^(B)_l|s⟩` avec `s' = s` privé du site
`l` ; puis `a_k` donne `e^{iθ N_{<k}(s'')}`. Or
`N_{<l} + N_{<k} = 2N_{<k} + n_k + N_mid` et `n_k = 1` sur le support
non nul. ∎

Le facteur `e^{2iθ N_{<k}}` est une **corde non locale** attachée à
l'opérateur de paire, sauf si `2θ ∈ 2πZ`. Donc :

> **Le vertex de médiation `d†_j a_j a_{j+1}` est un opérateur local
> uniquement pour `θ ∈ {0, π}`.** La famille Phase 7 est définie à `θ = 0`
> (bosons hard-core, représentation de `native_charge2_ladder.py`) et à
> `θ = π` (fermions). Pour tout autre θ, il y a non-localité cachée : no-go
> déclaré (§13).

À `θ = π`, règles de signe explicites (dérivées du Lemme 1, `N_mid` réduit au
seul site intermédiaire) :

```
a_j a_{j+1} = −(−1)^{ n^b_j }     b^(B)_{2j}   b^(B)_{2j+2}
b_j b_{j+1} = −(−1)^{ n^a_{j+1} } b^(B)_{2j+1} b^(B)_{2j+3}
a†_j b_j     =                    b^(B)†_{2j}  b^(B)_{2j+1}      (aucune corde)
```

Sur le secteur pertinent (barreaux simplement occupés, §3), `n^b_j = 0` dans
le canal `aa` et `n^a_{j+1} = 0` dans le canal `bb` : les deux canaux portent
le **même signe global −1**, absorbable par `d → −d`. La cohérence de signe
entre canaux est un test obligatoire du contrat (§12).

### 1.3 Charges et parités

```
Q   = Σ_j (n^a_j + n^b_j) + 2 Σ_j n^d_j          (charge U(1) pondérée)
P_a = (−1)^{N_a},   P_b = (−1)^{N_b}             (parités de branche)
```

Contrainte cinématique : `(−1)^{N_a+N_b} = (−1)^Q` dans tout secteur de
charge ; à `Q` fixé, la seule symétrie discrète indépendante est la **parité
relative** `P_rel ≡ P_a`.

### 1.4 Opérateurs élémentaires de barreau

```
N_j = n^a_j + n^b_j − 2 n^a_j n^b_j              (projecteur : barreau simplement occupé)
X_j = a†_j b_j + b†_j a_j                        (échange de rail sur le barreau j)
Z_j = n^a_j − n^b_j                              (déséquilibre de rail)
Y_j = i( b†_j a_j − a†_j b_j )
```

Algèbre exacte (vérifiable en dimension 4 par barreau) :

```
X_j² = Y_j² = Z_j² = N_j ,   X_j N_j = N_j X_j = X_j ,   X_j³ = X_j ,
[Z_j, X_j] = 2i Y_j ,        X_j Y_j = i Z_j N_j  (cycliques),
X_j |0⟩_j = X_j |ab⟩_j = 0 .
```

Sur le sous-espace `N_j = 1` (états `|a⟩_j`, `|b⟩_j`), `(X, Y, Z)` est une
algèbre de Pauli avec `Z|a⟩ = +|a⟩`, `Z|b⟩ = −|b⟩`. `X_j` **conserve la
charge du barreau** et est **impair** sous `P_a` et sous `P_b`.

### 1.5 Charge de cellule

Chaque médiateur (charge 2) est partagé pour moitié entre ses deux barreaux :

```
q_cell(j) = n^a_j + n^b_j + n^d_{j−1} + n^d_j ,      (n^d_{−1} ≡ n^d_{L−1} ≡ 0 en OBC)
Σ_j q_cell(j) = Q       (identité exacte, toutes conditions aux bords ouvertes).
```

---

## 2. Théorème du mur U(1)  -  pas d'opérateur de bord chargé dans un bulk incompressible

C'est le pivot de conception. Il transforme le mur Phase 6 en contrainte
structurelle.

**Théorème 1.** Soit `H` conservant `Q`, `P` le projecteur sur le multiplet
fondamental, supposé entièrement dans le secteur `Q = Q₀` et exactement
dégénéré d'énergie `E₀`. Soit `O` un opérateur de charge définie `q ≠ 0`
(`[Q, O] = q O`). Soit
`Δ_c(q) = min spec( H |_{Q₀+q} ) − E₀ > 0` (gap de charge). Alors

```
ε_edge(O) = ||(1−P)[H,O]P|| / ||O P|| ≥ Δ_c(q) .
```

*Preuve.* `O P` est entièrement dans le secteur `Q₀ + q`, donc orthogonal à
l'espace fondamental : `(1−P) O P = O P`. Alors
`(1−P)[H,O]P = (1−P)(H − E₀) O P = (H − E₀) O P`, et
`||(H−E₀) O P|| ≥ Δ_c(q) ||O P||` puisque `H − E₀ ≥ Δ_c(q)` sur le secteur
`Q₀+q`. ∎

(Si le multiplet est quasi dégénéré avec étalement `δ`, la borne devient
`Δ_c(q) − δ`.)

**Conséquences.**

1. Dans un bulk **incompressible** (`Δ_c > 0`, objectif A), tout opérateur de
   bord **chargé** a un `ε_edge` borné inférieurement par une constante
   indépendante de `L` et du support : la loi `ε_edge → 0` est impossible.
   L'objectif B force donc `[Q, O_edge] = 0` : **l'opérateur de bord doit
   être neutre.**
2. Le parent externe U(1) de la Phase 6 n'échappait à ce théorème que parce
   qu'il est compressible (`Δ_c ≈ 0`)  -  et son générateur de bord échoue
   quand même (`ε_edge ≈ 6.31`, plat en `L`). La route « générateur de bord
   chargé de type Iemini » est fermée **des deux côtés** : interdite par le
   théorème dans le régime incompressible, réfutée numériquement dans le
   régime compressible.
3. Le seul degré de liberté restant qui puisse porter une action logique non
   triviale tout en étant neutre est le **degré relatif des rails**, chargé
   sous `P_a` mais pas sous `U(1)` : les opérateurs `X_j` (et leurs
   habillages). Toute la suite en découle.

---

## 3. État cible et sous-espace cible exact

### 3.1 Définition

États de barreau `|→⟩_j = (|a⟩_j + |b⟩_j)/√2`, `|←⟩_j = (|a⟩_j − |b⟩_j)/√2`
(états propres de `X_j`, charge 1 exacte). Sous-espace cible (secteur
`Q = L`, aucun médiateur réel) :

```
|Ω_+⟩ = ⊗_{j=0}^{L−1} |→⟩_j ,      |Ω_−⟩ = ⊗_{j=0}^{L−1} |←⟩_j ,
𝒞 = span{ |Ω_+⟩, |Ω_−⟩ } .
```

Base logique symétrique (états propres de `P_a`) :

```
|±⟩ = ( |Ω_+⟩ ± (−1)^L |Ω_−⟩ ) / √2 ,        P_a |±⟩ = ± |±⟩ .
```

Actions exactes : `P_b|Ω_±⟩ = |Ω_∓⟩`, `P_a|Ω_±⟩ = (−1)^L|Ω_∓⟩`,
`X_j|Ω_±⟩ = ±|Ω_±⟩`, `Q|Ω_±⟩ = L|Ω_±⟩`.

### 3.2 Indistinguabilité locale  -  exacte, pas asymptotique

**Lemme 2.** Pour tout opérateur `O` de support strictement inclus dans les
`L` barreaux :
(i) `⟨Ω_+| O |Ω_−⟩ = 0` ;
(ii) si de plus `[O, P_b] = 0` (opérateur symétrique), alors
`⟨Ω_+| O |Ω_+⟩ = ⟨Ω_−| O |Ω_−⟩`, donc `⟨+|O|+⟩ = ⟨−|O|−⟩` et
`⟨+|O|−⟩ = 0` dès qu'en outre `[O,P_a]=0`.

*Preuve.* (i) L'élément de matrice factorise sur tout barreau hors support en
`⟨→|←⟩ = 0` ; un seul barreau hors support suffit. (ii)
`⟨Ω_+|O|Ω_+⟩ = ⟨Ω_+|P_b† O P_b|Ω_+⟩ = ⟨Ω_−|O|Ω_−⟩`. ∎

Donc `d_local = 0` exactement, pour tout `L ≥` support + 1, **mais
uniquement sur l'algèbre symétrique** (opérateurs commutant avec
`P_a, P_b`).

**Correction d'audit (décisive)** : dans l'algèbre physique complète, le
code n'est PAS localement indistinguable. Calcul exact (audit §1) :
`P X_j P = diag(+1,−1)` dans `{Ω_+, Ω_−}`, soit `X̄` dans la base de
parité  -  **non scalaire, de norme 1, pour tout barreau `j`**. `X_j` est un
opérateur physique local (le hopping de barreau `J_perp` du modèle gelé) qui
implémente l'opération logique depuis n'importe quel site. Par ailleurs
`P Y_j P = P Z_j P = 0` et `P n^{a,b}_j P = ½ P` (scalaires). Distinction
obligatoire :
- A. indistinguabilité face à tous les opérateurs locaux : **fausse** ;
- B. indistinguabilité face aux opérateurs locaux préservant `P_a, P_b` :
  **vraie** (Lemme 2).

Le résultat est donc une **mémoire Ising/Z₂ protégée par symétrie** (le
paramètre d'ordre `X_j` est P-impair), pas un code topologique. Un bain
couplé à `X_j` produit une erreur logique au premier ordre depuis n'importe
quel site (taux ∝ `ε_⊥ L`).

---

## 4. Dérivation des projecteurs à partir du sous-espace cible

Méthode parent standard : pour chaque support local `S`, le terme parent est
le projecteur sur l'orthogonal du sous-espace engendré par les matrices
densité réduites `Tr_{S^c} |Ω_±⟩⟨Ω_±|`.

### 4.1 Support = 1 barreau

`ρ_j(Ω_±)` est portée par les états de charge 1 du barreau. Le parent local en
découle : pénaliser toute déviation de la charge de cellule,

```
C_j = ( q_cell(j) − 1 )² .
```

Développement exact (tous les `n` sont idempotents)  -  uniquement des termes
de densité à 1 et 2 corps, tous natifs :

```
C_j = 1 − (n^a_j + n^b_j + n^d_{j−1} + n^d_j)
      + 2[ n^a_j n^b_j  +  (n^a_j + n^b_j)(n^d_{j−1} + n^d_j)  +  n^d_{j−1} n^d_j ] .
```

(`U_rung` = `n^a n^b` ; répulsion médiateur–rail ; répulsion
médiateur–médiateur ; potentiels chimiques.) `C_j` est diagonal, PSD, à
valeurs entières, et `C_j|Ω_±⟩ = 0`.

### 4.2 Support = 1 liaison (médiateur)

`ρ(Ω_±)` ne contient aucun médiateur : le parent est le projecteur `n^d_j`
lui-même (idempotent). `n^d_j|Ω_±⟩ = 0`.

### 4.3 Support = 2 barreaux adjacents

Sur `span{|aa⟩,|ab⟩,|ba⟩,|bb⟩}` (les deux barreaux simplement occupés), les
réduites des cibles engendrent `span{|→→⟩, |←←⟩}`. Or, identité exacte en
dimension 4 :

```
|→→⟩⟨→→| + |←←⟩⟨←←| = (1 + X_j X_{j+1}) / 2 .
```

Le projecteur parent sur l'orthogonal, prolongé par 0 hors du secteur
simplement occupé (rappel `X_j = N_j X_j N_j`), est donc

```
Π^B_j = ( N_j N_{j+1} − X_j X_{j+1} ) / 2 .
```

Vérifications : `Π^B_j` est hermitien ; `(Π^B_j)² = Π^B_j` (car
`(X_jX_{j+1})² = N_jN_{j+1}` et `X N = X`) ; `Π^B_j |Ω_±⟩ = 0` ;
`Π^B_j = 0` sur tout état où l'un des deux barreaux n'est pas simplement
occupé.

Développement en opérateurs physiques (les « canaux » de la liaison) :

```
X_j X_{j+1} = [ a†_j b_j a†_{j+1} b_{j+1} + h.c. ]   (transfert de paire croisé, ΔN_a = ±2)
            + [ a†_j b_j b†_{j+1} a_{j+1} + h.c. ]   (échange de rail, ΔN_a = 0).
```

Le premier canal est exactement celui que le médiateur de charge 2 engendre
(`bb → d → aa`) ; le second est un **superéchange par doublon virtuel**
(hoppings intra-rail au second ordre, §10)  -  ce n'est **pas** un processus
de médiateur de charge 2. Statut précis (audit §4) : `Π^B_j` est *dérivé* du
sous-espace cible (aucun terme deviné), mais sa réalisation microscopique
exige **deux** ingrédients (médiateur 6D + doublon natif) et une **condition
de matching** `2t²/E_v = g²/Δ` ; un médiateur seul ne suffit pas (point de
Heisenberg caché, §10.3). Les **trois canaux 6F ne sont pas utilisés** ici  - 
la factorisation 6F appartenait au parent externe Iemini ; aucune filiation
n'est revendiquée.

---

## 5. Candidate Hamiltonian

### 5.1 Parent au point fixe (objet principal, exactement soluble)

```
H_fix = Σ_{j=0}^{L−1} U C_j  +  Σ_{j=0}^{L−2} Δ n^d_j  +  Σ_{j=0}^{L−2} J Π^B_j  +  H_bord ,

U > 0,  Δ > 0,  J > 0 ;      H_bord = 0  (canonique ; ancres facultatives §9.4).
```

Terme par terme :

1. `U C_j` : contrainte de charge de cellule (Mott locale), §4.1  -  **coût
   local vrai** d'ajout/retrait de charge ; ce n'est ni une contrainte de
   charge fixée à la main, ni un terme global `Uc(N−N₀)²` : c'est une somme
   de projecteurs strictement locaux.
2. `Δ n^d_j` : détuning du médiateur ; `Δ > 0` est **obligatoire** (sinon
   dégénérescence extensive, §13.4).
3. `J Π^B_j` : stabiliseur de liaison dérivé, §4.3.

### 5.2 Famille microscopique native (mêmes symétries, mêmes ingrédients que 6D)

```
H_micro = Σ_j U C_j + Σ_j Δ n^d_j
        − g Σ_j [ d†_j ( a_j a_{j+1} + b_j b_{j+1} ) + h.c. ]        (vertex 6D)
        − t Σ_j [ a†_j a_{j+1} + b†_j b_{j+1} + h.c. ]               (hopping intra-rail)
        + V Σ_j [ n^a_j n^a_{j+1} + n^b_j n^b_{j+1} ]                (V_leg)
        + h Σ_j Z_j                                                   (détuning de rail)
        + H_bord .
```

**`J_perp = 0` exactement** : le hopping de barreau `a†_j b_j` est impair
sous les deux parités de branche ; son annulation microscopique est
précisément le rôle de l'interféromètre à flux π audité en Phase 5 (norme de
Schur croisée `3.48e−19`). La Phase 7 consomme ce résultat comme ressource ;
elle ne modifie pas le modèle gelé.

La relation `H_micro → H_fix` (élimination de Schrieffer–Wolff au second
ordre, condition de ligne d'Ising) est dérivée au §10.

---

## 6. Exact symmetries  -  démonstrations

Bilan par terme (`ΔN_a`, `ΔN_b`, `ΔN_d`, `ΔQ`) :

| terme | ΔN_a | ΔN_b | ΔN_d | ΔQ | P_a | P_b |
|---|---|---|---|---|---|---|
| `C_j`, `n^d_j`, `V`, `h Z_j` (diagonaux) | 0 | 0 | 0 | 0 | ✓ | ✓ |
| `X_j X_{j+1}` canal croisé | ±2 | ∓2 | 0 | 0 | ✓ | ✓ |
| `X_j X_{j+1}` canal échange | 0 | 0 | 0 | 0 | ✓ | ✓ |
| `d†_j a_j a_{j+1} + h.c.` | ∓2 | 0 | ±1 | 0 | ✓ | ✓ |
| `d†_j b_j b_{j+1} + h.c.` | 0 | ∓2 | ±1 | 0 | ✓ | ✓ |
| `t` intra-rail | 0 | 0 | 0 | 0 | ✓ | ✓ |

Donc, terme à terme et exactement :

```
[H_fix, Q] = [H_micro, Q] = 0 ,
[H_fix, P_a] = [H_micro, P_a] = 0 ,
[H_fix, P_b] = [H_micro, P_b] = 0 .
```

Ce sont des identités d'algèbre (chaque terme change `N_a`, `N_b` par des
quantités paires et `Q` par zéro), pas des résultats numériques. Le test
numérique correspondant du contrat doit rendre `0` à la précision machine.

Remarque : `H_fix` conserve en outre chaque `q_cell(j)` **au point fixe**
(tous ses termes commutent avec tous les `q_cell(j)`), tandis que `H_micro`
ne conserve que `Q` (le terme `t` déplace la charge entre cellules  -  c'est
lui qui rend la contrainte Mott dynamique et non triviale).

---

## 7. Local-projector algebra  -  noyaux, commutateurs, frustration-freeness

### 7.1 Noyaux locaux

- `ker C_j` = états avec `q_cell(j) = 1` ;
- `ker n^d_j` = médiateur `j` vide ;
- `ker Π^B_j` = `span{|→→⟩, |←←⟩}` ⊕ (tout état où `N_j N_{j+1} = 0`).

### 7.2 Commutateurs du parent  -  H_fix est un parent commuting-projector

Liens disjoints : trivialement nuls (supports disjoints, opérateurs pairs).

Liens chevauchants, la seule vérification non triviale :

```
[Π^B_j, Π^B_{j+1}] = ¼ [ N_jN_{j+1} − X_jX_{j+1} ,  N_{j+1}N_{j+2} − X_{j+1}X_{j+2} ] = 0 ,
```

car (a) les parties `N` sont diagonales et commutent entre elles et avec les
`X_kX_{k+1}` de même support partiel via `X N = X` ; (b)

```
[X_jX_{j+1}, X_{j+1}X_{j+2}] = X_j X_{j+1}² X_{j+2} − X_{j+1}X_{j+2}X_jX_{j+1}
                             = X_j N_{j+1} X_{j+2} − X_j N_{j+1} X_{j+2} = 0 ,
```

en utilisant `[X_j, X_{j+2}] = 0` (supports disjoints) et
`X_{j+1}² = N_{j+1}`, `X N = N X`. De plus `[Π^B_j, C_k] = 0` (les `X`
conservent toutes les charges de cellule) et `[Π^B_j, n^d_k] = 0`. Donc :

> **Tous les termes de `H_fix` commutent deux à deux**  -  mais **correction
> d'audit** : `C_j` n'est PAS un projecteur (`spec(C_j) = {0,1,4,9}`,
> `C_j² ≠ C_j`). Tel qu'écrit, `H_fix` est un parent **frustration-free à
> termes PSD commutants**, pas un « commuting-projector » au sens strict.
> La forme stricte est restaurée exactement par l'une des deux réécritures
> (audit §3) :
> (a) spectrale : `C_j = Π^{{0,2}}_j + 4Π^{{3}}_j + 9Π^{{4}}_j` avec
> `Π^{S}_j = 1_{q_cell(j)∈S}` (projecteurs locaux diagonaux commutants) ;
> (b) projecteur minimal `K_j = 1_{q_cell(j)≠1} = 1 − e₁ + 2e₂ − 3e₃ + 4e₄`
> (fonctions symétriques élémentaires des quatre occupations de la cellule),
> idempotent, même noyau, mêmes commutations. Les preuves de gap (§8)
> survivent inchangées dans les deux variantes.

Pour `H_micro`, en revanche (à déclarer honnêtement) :

```
[Π^B_j, d†_j a_j a_{j+1} + d†_j b_j b_{j+1} + h.c.] ≠ 0 ,
```

valeur explicite au secteur pertinent : le vertex annihile la combinaison
antisymétrique `(|aa⟩−|bb⟩)/√2` (qui est l'image de `Π^B_j` dans le plan
aligné) mais sa partie création `a†a†d` recrée cette combinaison depuis
`|d_j⟩` :

```
[Π^B_j, V_g,j] = (g/√2) ( |φ_j⟩⟨d_j| − |d_j⟩⟨φ_j| ) ⊗ 1_reste ,
φ_j = (|a_j a_{j+1}⟩ − |b_j b_{j+1}⟩)/√2 ,        ||[Π^B_j, V_g,j]|| = g/√2 · √2 = g .
```

`H_micro` n'est donc **pas** un parent commutant ; aucune preuve de gap n'est
revendiquée pour lui au niveau projecteur. Son statut est perturbatif (§10).

### 7.3 Frustration-freeness et intersection globale des noyaux

Chaque terme de `H_fix` est PSD et annihile `|Ω_±⟩` (§4). Donc `H_fix ≥ 0`
et `H_fix |Ω_±⟩ = 0` : **frustration-free, démontré**.

**Lemme 3 (intersection globale).** Pour `L ≥ 2`, en OBC comme en PBC :

```
ker H_fix = 𝒞 = span{ |Ω_+⟩, |Ω_−⟩ } ,      dim ker H_fix = 2 .
```

*Preuve.* `E = 0` exige terme à terme (PSD) : `n^d_j = 0` ∀j (aucun
médiateur) ; alors `q_cell(j) = n^a_j + n^b_j = 1` ∀j : chaque barreau est
exactement simplement occupé  -  le secteur est isomorphe à `(C²)^{⊗L}`,
`Q = L` gelée. Sur ce secteur, `Π^B_j = 0` ⇔ `S_j ≡ X_j X_{j+1} = +1`. Les
`S_j` sont `L−1` involutions commutantes indépendantes (OBC) sur `(C²)^{⊗L}`
 -  un groupe stabilisateur abélien  -  donc leur espace propre commun `+1` a
dimension `2^L / 2^{L−1} = 2`, et il contient `𝒞` : c'est `𝒞`. En PBC les
`L` contraintes satisfont la relation `Π_j S_j = 1`, soit `L−1`
indépendantes : dimension 2 également. ∎

Conséquence diagnostique : **dégénérescence 2 en OBC et en PBC**  -  signature
d'un ordre protégé par symétrie de classe Kitaev/Ising, **pas** d'un ordre
topologique intrinsèque (qui donnerait une GSD dépendant de la topologie).
Prédiction falsifiable, à auditer telle quelle.

---

## 8. Bulk gap and charge gap  -  démonstrations exactes

Tous les termes de `H_fix` commutent : le spectre est **additif**,

```
spec(H_fix) ⊆ { U m_C + Δ m_d + J m_B :  m_C, m_d, m_B ∈ N } ,
```

et chaque combinaison atteinte correspond à des secteurs simultanés des
projecteurs. D'où, exactement et pour tout `L` :

1. **Gap neutre** (secteur `Q = L`) :
   `Δ_neutre = min(J, Δ)`  -  atteint par un mur de domaine (`S_j = −1` sur une
   liaison, coût `J`) ou par la substitution de deux monomères alignés par un
   médiateur (coût `Δ`). Avec le point recommandé `Δ > J` : `Δ_neutre = J`,
   **indépendant de `L`**  -  macroscopique, à comparer à la séquence externe
   `1.4621 → 0.2720` qui décroît.

2. **Gap de charge** :
   pour tout état de charge `Q = L ± m`,
   `Σ_j (q_cell(j) − 1) = Q − L = ±m` ⇒ `Σ_j (q_cell(j) − 1)² ≥ m`, donc
   `E ≥ U m` ; la borne est atteinte (ajouter une particule sur un barreau,
   une seule cellule violée) :

   ```
   E₀(Q = L ± 1) − E₀(Q = L) = U       exactement, pour tout L.
   ```

   Compressibilité nulle, stabilité grand-canonique pour `|μ| < U`. Le coût
   est **local** (une cellule), pas global : c'est le vrai coût d'ajout
   exigé par l'objectif A. Aucune énergie de charging externe n'est
   nécessaire au point fixe. (Pour `H_micro`, le coût devient
   `U − O(t²/U, g²/Δ)`, uniforme dans le bulk  -  pas d'état de bord de charge
   sous le gap, cf. §9.4 ancre 1.)

3. **Dégénérescence** : `dim = 2` exactement (Lemme 3), splitting `0` exact
   au point fixe pour tout `L`.

---

## 9. Edge operator and boundary anchor

### 9.1 Opérateurs logiques  -  correction d'audit : ordre de bulk, pas mode de bord

**Au point fixe, il n'existe aucun mode de bord.** Le calcul complet (audit
§2) donne `[H_fix, X_j] = 0` pour **tout** `j = 0..L−1` (chaque `X_j`
commute avec `X_{j−1}X_j`, `X_jX_{j+1}`, tous les `N`, `C`, `n^d`), et tous
les `X_j` ont la même action logique `X_j|Ω_±⟩ = ±|Ω_±⟩`, soit `X̄` dans la
base de parité. `X_0` n'a donc aucun statut distinctif : **`X_j` est le
paramètre d'ordre local d'une brisure spontanée Z₂** (point d'Ising
classique dans la base x, avec `L` opérateurs P-impairs conservés  -  point
fin-réglé). Les expressions « mode de bord » et « quasi-conservation de
bord » sont **retirées au point fixe**. Restent exacts : `[Q, X_j] = 0`
(neutralité, condition du Théorème 1), et les opérateurs logiques
`X̄ = X_j` (n'importe quel `j`), `Z̄ = P_a`.

La cohérence GSD(OBC) = GSD(PBC) = 2 (Lemme 3) est la **signature SSB** :
deux états de brisure reliés par la symétrie globale, insensibles aux
conditions aux bords  -  par opposition à une dégénérescence de bord
(Kitaev fermionique : OBC 2 / PBC 1).

Une localisation de bord n'apparaît que **sous perturbation** (`h ≠ 0`,
§9.2) : au support 1, `ε(X_j) ≈ 2h` pour tout `j` (bord = bulk) ; c'est
seulement l'habillage qui distingue le bord (récurrence convergente,
`ε ∝ r^w`) du bulk (branchement des deux côtés, plateau prédit `Θ(h)`
indépendant du support  -  aucun mode sous-gap de bulk dans le spectre
fermion-libre). Diagnostic quantitatif à auditer :
`ε_min(bord, w)/ε_min(bulk, w) ≈ c·r^w·J/h → 0`.

### 9.2 Loi ε_edge sous perturbation symétrique  -  mode zéro fort dérivé

Perturbation native dominante : `h Σ_j Z_j` (détuning de rail, symétrique).
Sur le secteur monomère, `H_fix + hΣZ` se réduit à la chaîne d'Ising en champ
transverse `H = −J_I Σ X_jX_{j+1} + h Σ Z_j + const`, `J_I = J/2`. Posant
`r = h / J_I = 2h/J < 1`, l'habillage s'obtient ordre par ordre par la
récurrence de mode zéro fort : le résidu de `[H, ·]` au rang `n` est annulé
par le terme de rang `n+1`,

```
Γ^{(j)} = 𝒩_j Σ_{n=1}^{j} r^{n−1} ( Π_{m<n} Z_m ) X_n            (support = j barreaux de bord),
[H, Γ^{(j)}] = 𝒩_j · 2 J_I r^j × (chaîne de Pauli unitaire) .
```

Les chaînes de Pauli du développement anticommutent deux à deux, donc
`||Σ|| = (Σ r^{2(n−1)})^{1/2}` et

```
ε_edge(L, j) = 2 J_I r^j √( (1−r²) / (1−r^{2j}) ) + O(r^{L−j}) .
```

**Loi exigée, dérivée** : décroissance exponentielle en le support `j`
(pente `ln r`), **plate en `L`** dès `L > j`, exactement nulle à `h = 0`.
`Γ^{(j)}` reste neutre (`[Q, Γ] = 0`), impair sous `P_a` (un seul `X` par
terme), et son support est contigu depuis le bord  -  aucune non-localité
cachée. Sous perturbations symétriques génériques supplémentaires (termes
XYZ, §10), la même structure d'habillage existe dans la phase d'Ising ; son
audit est numérique (conjecture contrôlée, §14).

### 9.3 Lisibilité locale  -  statut corrigé

Il n'existe pas de doublet localisé au bord (les réduites des états `|±⟩`
coïncident sur toute sous-région stricte  -  Lemme 2), **mais** cette
indistinguabilité ne vaut que dans l'algèbre symétrique : `P X_j P = X̄` non
scalaire partout (§3.2 corrigé). L'information est portée par la parité
relative globale `P_a` (lecture `Z̄` non locale) et elle est **lisible /
basculable localement par tout opérateur P-impair**  -  c'est la définition
d'une mémoire SSB sous symétrie, pas d'un code topologique.

### 9.4 Ancre de bord  -  analyse sans complaisance

Trois ancres candidates, analysées par leurs commutateurs exacts.

**Ancre 1  -  potentiel de charge de bord** `A₁ = μ_e (n^a_0 + n^b_0)`.
Préserve `U(1)`, `P_a`, `P_b`. Commute avec `X_0`, avec tous les `Π^B`, tous
les `C_j`, tous les `n^d`. Sur le code : `A₁ P = μ_e P` (les états logiques
ont charge de barreau 1)  -  **scalaire sur le code** : ne fabrique aucune
dégénérescence, ne dégrade pas `ε_edge`, ne sélectionne rien logiquement.
Rôle légitime et unique : rigidifier la charge au bord le long de `H_micro`
(prévenir un état de bord de charge quand `U` est fini). **Admise.**

**Ancre 2  -  détuning de rail de bord** `A₂ = h₀ Z_0`.
Symétrique (`U(1)`, `P_a`, `P_b` ✓). Premier ordre nul sur le code
(`⟨Ω_±|Z_0|Ω_±⟩ = 0`  -  Lemme 2). Ne scinde le doublet qu'à l'ordre `L`
(retournement complet : splitting `∝ J_I r₀ r^{L−1}`), mais dégrade le mode
de bord au premier ordre : `[A₂, X_0] = 2i h₀ Y_0`, soit
`ε_edge(j=1) = 2h₀` immédiatement. Elle ne « sélectionne » aucun secteur de
parité (elle commute avec `P_a`). **Admise seulement transitoirement**
(pilotage), à `h₀ ≪ J`, et comptée dans le budget `ε_edge`.

**Ancre 3  -  ancre de conversion de bord** `A₃ = κ X_0` (transfert croisé de
barreau au bord). Conserve `U(1)` (`[Q, X_0] = 0`) mais **brise `P_a` et
`P_b`** (`X_0` est impair). Effet exact au point fixe : `X_0` agit comme
`diag(+1, −1)` dans la base `{Ω_+, Ω_−}` ⇒ splitting `2κ` **linéaire**  -  elle
détruit la protection en tant que terme statique. Elle sélectionne l'état
`|Ω_∓⟩` (base X̄), pas un secteur de parité. **Rejetée comme terme de
`H_bord`.** Reclassée explicitement comme **ressource externe active** :
impulsion d'initialisation / porte logique `X̄` (appliquée puis éteinte),
réalisable en détunant volontairement l'interféromètre à flux π du barreau de
bord. La sélection d'un secteur de `P_a` (base `Z̄`) se fait, elle, par
mesure de parité, pas par une ancre.

Conclusion dérivée (pas affirmée) : la quasi-conservation de `O_edge` ne
requiert **aucune** ancre ; l'ancre admissible (charge) est logiquement
inerte ; l'ancre qui « choisit » un état logique brise nécessairement la
symétrie protectrice  -  c'est un théorème de commutation, et c'est pourquoi
elle doit rester une ressource externe pulsée.

---

## 10. Réduction microscopique : de H_micro au point fixe, et explication du mur 6E

### 10.1 Élimination de Schrieffer–Wolff au second ordre

Secteur de basse énergie de `H_micro` : monomères (un par barreau),
`P_M` le projecteur associé. États virtuels : médiateur occupé (coût `Δ`),
barreau doublement occupé + trou (coût `E_v = 2U` via `C`, plus `2U_rung`
si l'on sépare le terme `n^a n^b`  -  noté `E_v` ci-dessous). Au second ordre :

- **canal médiateur** (`g`) : sur une liaison alignée,
  `|aa⟩ → |d⟩ → |aa⟩ ou |bb⟩` :

  ```
  H_g = −(g²/Δ) Σ_j ( |aa⟩ + |bb⟩ )( ⟨aa| + ⟨bb| )_j
      = −J_p Σ_j [ P^{aa} + P^{bb} + (σ⁺σ⁺ + σ⁻σ⁻) ]_j ,     J_p = g²/Δ ;
  ```

- **canal doublon** (`t`) : sur une liaison anti-alignée, deux chemins
  (`a` à droite puis `b` à gauche, ou l'inverse) :

  ```
  H_t = −J_ex Σ_j [ P^{ab} + P^{ba} + (σ⁺σ⁻ + σ⁻σ⁺) ]_j ,     J_ex = 2t²/E_v .
  ```

En notation Pauli de barreau (`σ⁺σ⁺ + h.c. = (XX − YY)/2`,
`σ⁺σ⁻ + h.c. = (XX + YY)/2`, `P^{aa}+P^{bb} = (1+ZZ)/2`,
`P^{ab}+P^{ba} = (1−ZZ)/2`), avec `δ ≡ J_p − J_ex` et le `V_leg` natif :

```
H_eff = Σ_j [ −((J_p+J_ex)/2) X_jX_{j+1} + (δ/2) Y_jY_{j+1} + ((V−δ)/2) Z_jZ_{j+1} ] + const .
```

### 10.2 Ligne d'Ising (condition de matching dérivée)

```
δ = 0  et  V = 0     ⇔     2 t² / E_v = g² / Δ ,   V_leg = 0
⇒   H_eff = −J̃ Σ_j X_j X_{j+1} + const = 2J̃ Σ_j Π^B_j + const ,    J̃ = g²/Δ .
```

Sur cette ligne, `H_micro` se projette **exactement** (au second ordre) sur
le parent dérivé, avec `J = 2g²/Δ`. Hors ligne, les termes `YY` et `ZZ` sont
des perturbations symétriques ; la phase d'Ising (donc le gap et la
dégénérescence) persiste pour `|δ|, |V| < O(J̃)` (stabilité perturbative,
statut conjecturé contrôlé  -  §14). Corrections d'ordre 4 : termes `XX` à
portée 2, symétriques, mêmes conclusions.

**Correction d'audit  -  validité restreinte à `θ = 0`.** Le calcul ci-dessus
vaut pour la représentation bosonique. À `θ = π`, le calcul de chemins avec
cordes rung-major (audit §4.2) montre que le superéchange **change de
signe** (`+2t²/E_v` : chaque chemin traverse une corde occupée une fois),
tandis que le canal médiateur garde le sien (cordes au carré). Sur la ligne
de matching, le parent fermionique est alors

```
H_eff(θ=π) = +J̃ Σ_j Y_j Y_{j+1} + const      (Ising en base Y, antiferro),
Π^B(θ=π)_j = ½( N_j N_{j+1} + Y_j Y_{j+1} ) ,   cibles = produits staggered d'états propres de Y.
```

Équivalent à `θ = 0` par unitaire on-site (`Π_{j impair} Z_j` puis
`e^{−iπ/4ΣZ_j}`), mais l'affirmation antérieure « même parent XX aux deux θ
via `d → −d` » était **fausse**. La famille est définie canoniquement à
`θ = 0` ; toute implémentation à `θ = π` doit utiliser la forme ci-dessus.

Enfin, la ligne de matching exige un superéchange uniforme : soit
`J1 = J2 = t` (écart déclaré à la dimérisation SSH du modèle gelé), soit un
médiateur dimérisé `g_i² = 2 J_i² Δ / E_v` liaison par liaison.

### 10.3 Explication analytique du mur 6E

Le scan 6E explorait le voisinage « médiateur pur » (`J_ex ≈ 0`, pas de
condition de matching). Alors `δ = J_p` et

```
H_eff = −(J_p/2) Σ ( XX − YY + ZZ )    ≅ (rotation de sous-réseau + Π_impairs σ^z)
        +(J_p/2) Σ ( XX + YY + ZZ ) :    chaîne de Heisenberg antiferromagnétique S=½.
```

Point **SU(2) caché, sans gap** (spectre de Bethe, échelles en `1/L`,
quasi-doublets non protégés). C'est exactement la phénoménologie observée en
6E : `split/gap = 1.96e−3` à `L=4` se dégradant à `3.53e−2` à `L=5`  -  un
comportement critique, pas une tendance protectrice. Le mur 6E n'était pas un
défaut de scan : le sous-espace scanné était **structurellement gapless**.
La sortie n'est pas « plus de points », c'est la condition de matching
`2t²/E_v = g²/Δ` qui brise le SU(2) émergent vers la phase d'Ising.

### 10.4 Fuite et rôle du médiateur

Au point fixe le médiateur est purement virtuel :
`⟨N_d⟩_micro ≈ (L−1) g²/Δ²` (amplitude `g/Δ` par liaison alignée, poids ½
par état logique, deux canaux cohérents  -  cf. la vérification de signe
§1.2). Le médiateur est l'**origine microscopique** du stabiliseur
(son canal croisé), pas un constituant de l'état cible : c'est la même
conclusion que 6D (« ingrédient local correct ») enfin munie de sa phase
many-body.

---

## 11. Finite-size predictions (L = 4, 6, 8, 10)

Point d'audit recommandé : `U = 4`, `Δ = 2`, `J = 1` (point fixe), et pour la
loi de bord `h = 0.1` (`J_I = 0.5`, `r = 0.2`). Secteur `Q = L`.

**Au point fixe exact (`h = 0`)**  -  tout écart au-delà de `1e−10` est un
rejet :

| observable | L=4 | L=6 | L=8 | L=10 |
|---|---|---|---|---|
| GSD OBC | 2 | 2 | 2 | 2 |
| GSD PBC | 2 | 2 | 2 | 2  -  **signature SSB, pas bord** |
| splitting logique | 0 | 0 | 0 | 0 |
| gap neutre | 1.0 | 1.0 | 1.0 | 1.0 (plat  -  contraste externe : 1.4621→0.2720) |
| gap de charge `E₀(L±1)−E₀(L)` | 4.0 | 4.0 | 4.0 | 4.0 (contraste externe : ≈ 0) |
| `max_j ‖P X_j P − ½tr(PX_jP)·I‖` | 1.0 | 1.0 | 1.0 | 1.0  -  **∀j : code NON localement indistinguable (algèbre complète)** |
| même test, opérateurs préservant `P_a,P_b` (≤ 2 barreaux) | 0 | 0 | 0 | 0 |
| `ε(X_0)` (bord) | 0 | 0 | 0 | 0 |
| `ε(X_{L/2})` (bulk) | 0 | 0 | 0 | 0  -  **aucune localisation au point fixe** |
| `⟨N_d⟩` | 0 | 0 | 0 | 0 |

**Sous `h = 0.1` (audit de la loi de protection)** :

```
ε(Γ^{(w)}_bord) = 2 J_I r^w √((1−r²)/(1−r^{2w})) :
   w=1 : 0.2000     w=2 : 0.03922     w=3 : 0.007839     w=4 : 0.001568
   (indépendant de L pour L > w, corrections O(r^{L−w}))
ε_min(bulk, w) : plateau Θ(2h) ≈ 0.2, indépendant de w (prédiction, audit §2.3)
P X_j P (bulk) ≈ ± m X̄,  m = (1−r²)^{1/8} ≈ 0.995  (action logique locale persistante)
```

- pente de `ln ε_edge` en `j` : `ln r = −1.609 ± 15 %` ;
- platitude en `L` : `|ε(10,j) − ε(8,j)| / ε(8,j) < 0.05` pour `j ≤ 3`  - 
  la signature exigée est le **couple** (décroissance en `j`, platitude en
  `L` à petite valeur), par opposition au générateur externe (plat en `L`
  mais à `6.31`) ;
- splitting logique `∝ r^L` : rapports `split(L+2)/split(L) = r² = 0.04
  ± 20 %` ;
- gap neutre `≈ 2(J_I − h) = 0.8` ;
- gap de charge `= 4.0 + O(h)` ;
- `d_local = O(r^{L−2})` (indétectable localement).

**Le long de `H_micro` sur la ligne d'Ising** (petit couplage,
`g, t ≪ Δ, U`) : mêmes signatures avec `J → 2g²/Δ`, fuite
`⟨N_d⟩ = (L−1) g²/Δ² ± 20 %`, gaps de charge `U − O(t²/U)` uniformes en `j`
(aucun mode de charge sous-gap au bord).

---

## 12. ANTLER implementation contract

1. **Modes** (aligné sur `antler/native_charge2_ladder.py`) :
   `mode_a(j) = 2j`, `mode_b(j) = 2j+1` pour `j = 0..L−1` ;
   `mode_d(bond) = 2L + bond` pour `bond = 0..L−2`. Total `3L − 1` modes
   hard-core ; base = masques binaires triés, secteur de charge pondérée
   `Q = Σ_{sites rails} n + 2 Σ_bonds n^d` (cf. `build_weighted_basis`).

2. **Règles de signe rung-major** : famille définie à `θ ∈ {0, π}`
   uniquement (Lemme 1). À `θ = 0` : aucune phase. À `θ = π` : appliquer aux
   opérateurs de paire les signes du §1.2
   (`a_ja_{j+1} → −(−1)^{n^b_j}`, `b_jb_{j+1} → −(−1)^{n^a_{j+1}}`, sur le
   masque intermédiaire, ordre normal « indice haut annihilé d'abord ») ;
   `X_j` sans corde ; `d` sans corde (défini JW-neutre). **Test de cohérence
   de signe obligatoire** : sur le secteur monomère, les canaux `aa` et `bb`
   doivent porter le même signe effectif (redéfinition `d → −d` admise, à
   consigner). **Attention (audit §4.2)** : à `θ = π` le superéchange change
   de signe ⇒ parent en base Y staggered (`Π^B(π) = ½(NN + YY)`) ; la table
   §11 et les états cibles du §3 valent tels quels **à `θ = 0` seulement**.

3. **Opérateurs à coder** (matrices creuses sur le secteur de charge) :
   `X_j, Y_j, Z_j, N_j` (§1.4) ; `q_cell(j)`, `C_j` (§1.5, §4.1) ;
   `Π^B_j` (§4.3) ; `S_j = X_jX_{j+1}` ; `P_a, P_b, Q` ; `H_fix` (§5.1) ;
   `H_micro` (§5.2) ; mode tronqué `Γ^{(j)}` (§9.2) ; ancres `A₁, A₂, A₃`
   (§9.4). Auditer chaque `Π` avec `local_projector_algebra`
   (hermiticité, idempotence, commutateurs disjoints/chevauchants) de
   `antler/phase7_parent_audit.py`.

4. **Observables obligatoires** :
   `ε_edge(L, j)` avec `P` = projecteur sur le doublet fondamental
   (vérifier d'abord `[Q, O_edge] = 0`  -  Théorème 1) ; gap neutre à
   `Q = L` ; gap de charge `E₀(L±1) − E₀(L)` ; `d_local` (max sur toutes les
   observables symétriques à support ≤ 2 barreaux) ; leakage `⟨N_d⟩` ;
   splitting logique ; dégénérescence OBC **et** PBC ; test d'injection de
   parité : ajouter `ε Σ_j X_j` et vérifier le splitting linéaire `2εL`
   (budget d'erreur logique, §13.2).

5. **Critères numériques de rejet** :
   - point fixe : toute ligne du tableau §11 violée au-delà de `1e−10` ;
   - commutateurs de symétrie ou `[Π^B_j, Π^B_{j+1}]` au-delà de `1e−12` ;
   - loi de bord : pente en `j` hors `ln r ± 15 %`, ou platitude en `L`
     violée (`> 5 %` pour `j ≤ 3`), ou `ε_edge(h=0) > 1e−10` ;
   - gap neutre non plat : `gap(L=10)/gap(L=4) < 0.9` ⇒ rejet (le parent
     externe donne `0.186`) ;
   - gap de charge `< U/2` à un `L` quelconque ⇒ rejet ;
   - `d_local > 1e−6` au point fixe, ou ne décroissant pas exponentiellement
     le long de la famille ⇒ rejet ;
   - dégénérescence PBC ≠ 2 ⇒ rejet de la présente analyse (et indice d'une
     physique différente à documenter) ;
   - leakage hors `(L−1)g²/Δ² ± 20 %` dans le régime perturbatif ⇒ rejet de
     la réduction §10.

---

## 13. No-gos et limites  -  liste exhaustive

1. **Pas d'ordre topologique  -  brisure spontanée Z₂ avec paramètre d'ordre
   local.** GSD = 2 en OBC et PBC (signature SSB) ; le paramètre d'ordre
   `X_j` est un opérateur physique local P-impair qui implémente `X̄`
   depuis n'importe quel site (`P X_j P` non scalaire ∀j  -  audit §1).
   La protection est celle d'une mémoire Ising sous symétrie exacte `P_a`
   (avec `P_b` esclave à `Q` fixé), dans le Hamiltonien **et** dans le
   couplage au bain. Aucun anyon non abélien, aucune tresse, aucun calcul
   topologique universel n'est revendiqué ni approché.
2. **Protection conditionnelle à l'exactitude des parités.** Un transfert
   croisé résiduel `ε_⊥ Σ_j X_j` (interféromètre imparfait) scinde le code
   **linéairement** : taux d'erreur `X̄` `≈ 2 ε_⊥ L` (élément de matrice
   exact au point fixe). Budget : avec la valeur auditée
   `ε_⊥ ≈ 3.5e−19`, négligeable ; mais la protection est de type
   « symétrie exacte », pas géométrique. À déclarer dans toute
   communication.
3. **Lecture logique `Z̄` non locale** (corde de parité `P_a` ou fusion
   bord-à-bord). Inhérent à la classe ; c'est le prix de
   l'indistinguabilité locale exacte (Lemme 2).
4. **`Δ → 0` interdit** : les pavages à médiateurs deviennent des modes
   zéro ⇒ dégénérescence extensive. La famille exige `Δ > 0` strictement ;
   le régime `g ≳ Δ` est hors du domaine dérivé (non exploré, non
   revendiqué).
5. **Secteurs dopés non protégés.** Hors `Q = L`, une lacune de monomère est
   mobile (canal de charge sans gap au sens du secteur) ; la mémoire n'est
   définie qu'au remplissage unité. Ce n'est pas une contrainte globale
   cachée : le coût local `U` (§8) rend le secteur `Q = L` auto-sélectionné
   pour `|μ| < U`.
6. **Pas d'auto-correction thermique** : murs de domaine 1D diffusifs à
   `T > 0` (taux logique `Z̄` ∼ Arrhenius `e^{−J/T}` puis marche
   aléatoire)  -  mémoire quantique 1D standard, pas une mémoire
   auto-corrigée.
7. **Restriction statistique** `θ ∈ {0, π}` (Lemme 1) : pour θ fractionnaire
   le vertex de médiation porte une corde  -  la famille Phase 7 ne définit
   pas d'extension anyonique.
8. **`H_micro` n'est pas commuting-projector** (§7.2) : toutes les
   propriétés exactes sont celles de `H_fix` ; celles de `H_micro` sont
   perturbatives d'ordre 2, avec stabilité conjecturée (§14).
9. **Opérateurs de bord chargés impossibles** (Théorème 1) : ferme
   définitivement la route « générateur Iemini tronqué » dans tout régime
   incompressible natif.

---

## 14. Claim boundary

**Prouvé (algébriquement, exact à tout L)** :
- Théorème 1 (mur U(1) : incompressibilité ⇒ `O_edge` nécessairement
  neutre) ;
- Lemme 1 (localité du vertex ⇔ `θ ∈ {0, π}`) et règles de signe ;
- symétries exactes de `H_fix` et `H_micro` (`Q`, `P_a`, `P_b`) ;
- `H_fix` : termes PSD mutuellement commutants, frustration-free, `ker = 𝒞`
  de dimension 2 (OBC/PBC), gap neutre `= min(J, Δ)`, gap de charge `= U`
  local et exact, compressibilité nulle ; forme commuting-projector stricte
  via la réécriture spectrale ou `K_j` (audit §3  -  `C_j` lui-même n'est PAS
  idempotent) ;
- indistinguabilité locale exacte des états logiques **sous l'algèbre
  symétrique seulement** (Lemme 2) ; dans l'algèbre complète,
  `P X_j P = X̄` non scalaire ∀j : code NON localement indistinguable
  (audit §1) ;
- `[H_fix, X_j] = 0` pour **tout** `j` au point fixe : `X_j` = paramètre
  d'ordre de bulk, « mode de bord » retiré ; sous champ `h`, mode fort de
  bord habillé avec loi `ε(w) = 2J_I r^w √((1−r²)/(1−r^{2w}))`, plate en
  `L`, contre plateau `Θ(h)` prédit en bulk ;
- ancres : `A₁` inerte, `A₂` symétrique mais coûteuse en `ε_edge`, `A₃`
  brise `P_a` et scinde linéairement (donc ressource pulsée, pas un terme) ;
- réduction SW d'ordre 2 de `H_micro` et condition de ligne d'Ising
  `2t²/E_v = g²/Δ`, `V_leg = 0` ;
- point « médiateur pur » = Heisenberg AFM caché, sans gap  -  explication
  structurelle du rejet 6E.

**Conjecturé (contrôlé, à auditer numériquement)** :
- stabilité de la phase d'Ising (gap, doublet, mode fort) hors ligne de
  matching pour `|δ|, |V|, |h| ≪ J̃`, et aux ordres SW ≥ 4 ;
- loi `ε_edge` sous perturbations symétriques génériques (XYZ) avec la même
  structure exponentielle en `j` ;
- valeurs de fuite et de gaps le long de `H_micro` (§11, régime
  perturbatif).

**Réfuté / fermé** :
- tout opérateur de bord chargé dans un bulk incompressible (Théorème 1) ;
- l'espoir qu'un scan supplémentaire de type 6E (voisinage médiateur pur)
  trouve une tendance protectrice : le sous-espace est gapless par symétrie
  émergente ;
- l'ancre de bord comme mécanisme statique de sélection de secteur : toute
  ancre sélective brise la symétrie protectrice (commutation, §9.4).

**Non revendiqué** : anyons, tresse, universalité, auto-correction
thermique, protection au-delà de l'exactitude microscopique des parités de
branche (Hamiltonien et bain), indistinguabilité locale dans l'algèbre
complète.

**Verdict final (après audit adversarial, `docs/PHASE7_ADVERSARIAL_AUDIT.md`)** :
**mémoire Ising/Z₂ sous symétrie exacte (SSB de la parité relative)**  - 
rejetée comme code localement indistinguable et comme candidat topologique ;
conservée comme benchmark auditable de protection sous symétrie, contrôle
positif/négatif du harnais `ε_edge` (bord vs bulk sous `h`), et support du
Théorème 1 (mur U(1)), qui contraint toute construction ultérieure : un futur
candidat devra rendre `max_j ‖P O_j P − scalaire‖ → 0` sur **tous** les
opérateurs locaux physiques, pas seulement les opérateurs symétriques.
