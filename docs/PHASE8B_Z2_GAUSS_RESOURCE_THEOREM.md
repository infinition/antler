# Phase 8B  -  théorème de ressource minimale pour une loi de Gauss Z2 de rail

**Révision 1.** Le gadget (T3) original est **réfuté par audit exact** :
avec un seul mode `μ` et un vertex unique `A = X_L + X_R + ηP_B`, tout
processus bas→bas est un mot pair en `A`, donc `H_eff = f(A²) ∈
span{1, X_LX_R, X_LP_B, P_BX_R}` ; le coefficient de `G_B` est nul
(`<2.8e-16` sur `λ/Δ_G = 0.025..0.20`). Le terme annoncé `K~λ³η/Δ_G²`
n'existe pas. Le §2 est remplacé par le §2R (gadget en chaîne Λ) ; la règle
de sélection et la primitive nouvelle requise sont au §2R.0.

Question : U(1) exacte + hopping `theta=pi` + médiateurs charge-2 + canaux de
paires + **au plus une nouvelle ressource locale par sommet** peuvent-ils
rendre le tunneling de barreau gauge-chargé, avec secteur bas local et gap
chargé fini ? Hors périmètre : tresse, jonction T, matrices logiques.

## 1. Verdict

Trois énoncés, du plus fort au conditionnel :

**(T1) Borne inférieure  -  zéro nouvelle ressource : IMPOSSIBLE.**
Toute charge locale conservée de l'algèbre de matière U(1) à `N` fixé est
fonction des `n_{j,l}` ; un générateur de Gauss `G_v` doit anticommuter avec
`c†_a c_b` (sinon il ne le charge pas) tout en commutant avec `H` et en étant
produit de facteurs locaux indépendants. Dans l'algèbre de matière seule, le
seul candidat est `(-1)^{n_{a,j}}` (ou produits), qui ne commute pas avec le
hopping intra-rail voulu : aucune assignation n'existe (c'est le no-go §7 du
verdict T reformulé). **Un degré de liberté de jauge par sommet est
indispensable.**

**(T2) No-go du gauging fin (par barreau) avec couplages ≤ 2-corps.**
Gauss fin : `G_j = τ^x_{l_j^-} (-1)^{n_{a,j}} τ^x_{l_j^+}` (qubit de jauge
`τ_l` charge 0 par lien). Le hopping intra-rail `a` voulu anticommute avec
`G_j, G_{j+1}` : sa version physique est le hopping habillé
`c†_{a,j} τ^z_l c_{a,j+1}`, **3-corps**. Or le hopping nu 2-corps est natif à
l'ordre 0 tandis que tout habillage dérivé par SW apparaît à un ordre
supérieur : le parasite violant Gauss **domine strictement** le terme voulu,
à tout point de la grammaire. C'est la hiérarchie 7C transposée au 3-corps,
ici fatale car le parasite est le terme cinétique lui-même. **NO-GO.**

**(T3) CONSTRUCTION CONDITIONNELLE RÉVISÉE  -  gauging à gros grain.**
Le no-go (T2) est contourné exactement par coarse-graining : charge de
matière par **bloc** de `b ≥ 2` barreaux, `P_B = (-1)^{N_{a,B}}`, qubit de
jauge sur les **frontières de blocs** seulement (≤ 1 ressource par sommet,
en fait 1 par `b` barreaux). Alors :
- le hopping intra-bloc commute avec tous les `G_B` (il conserve
  `N_{a,B} mod 2`) : **aucun habillage requis** ;
- le transport inter-blocs du rail `a` par saut simple est mis à zéro par le
  contrôle natif par lien, et remplacé par les **canaux de paires** déjà en
  grammaire (`ΔN_{a,B} = 0, ±2` : Gauss-pairs) ;
- le tunneling de barreau `c†_a c_b` **n'importe où** dans un bloc flippe
  `P_B` donc anticommute avec `G_B` : gauge-chargé (preuve §3).

## 2. Construction T3 originale  -  réfutée, conservée pour traçabilité

**Degrés de liberté par frontière de bloc** (la ressource nouvelle, unique,
déclarée) : un qubit `τ` charge 0 avec champs natifs `h_x τ^x`, `h_z τ^z`,
couplage densité natif `κ τ^α n_{a,j}` (2-corps), et un médiateur auxiliaire
hard-core charge 0 `μ_B` par bloc, détuning `Δ_G` (même famille que les
médiateurs existants, charge différente : composante de la même ressource
déclarée).

**Opérateur de Gauss exact :**
\[
G_B = \tau^x_{B-1,B}\; (-1)^{N_{a,B}}\; \tau^x_{B,B+1},\qquad [G_B,G_{B'}]=0 .
\]

**Hamiltonien microscopique local :**
\[
H = H_0^{\rm intra}
 + H_{\rm paires}^{\rm inter}
 + \Delta_G\sum_B \mu_B^\dagger\mu_B
 + \lambda\sum_B\left[\mu_B^\dagger\big(\tau^x_{B-1,B}+\tau^x_{B,B+1}
 + \eta\,(-1)^{N_{a,B}}\big)+h.c.\right]
 + H_{\rm contre}
\]
où `H_0^intra` est le ladder Phase 8 par bloc (hopping `t`, médiateurs
charge-2, canaux), `H_paires^inter` les canaux de paires de frontière
(natifs), et `(-1)^{N_{a,B}} = Π_j (1-2n_{a,j})` est couplé à `μ_B` via la
chaîne SW des couplages densité 2-corps déclarés (voir tableau ; pour `b=2`
c'est un terme d'ordre 2 dérivé, pas natif).

**SW/Floquet.** Tous les opérateurs `{τ^x, P_B, μ†μ}` engendrent une algèbre
**abélienne** ; l'élimination SW de `μ_B` reste dans cette algèbre :
- ordre 2, `λ²/Δ_G` : parasites `τ^x τ^x`, `τ^x P_B`, champs  -  tous **dans le
  commutant de tous les `G_B`** (ils ne violent jamais Gauss) ; diagonaux,
  compensés exactement par `H_contre` (champs `h_x`, potentiel `μ_a` natifs) ;
- ordre 3, `K = c₃ λ³η/Δ_G²` : le terme cible `-K G_B`. La hiérarchie
  d'ordres est ici **inoffensive** : les termes d'ordre inférieur commutent
  avec la cible (contrairement à 7C où ils la détruisaient) ;
- Floquet : `G_B` commute avec `H0` et `H1` intra-bloc et avec les canaux de
  frontière (tous Gauss-pairs) → Gauss exacte sur le cycle **si** les canaux
  gardent une parité `N_a` définie (même condition que §3.5 du verdict T :
  espèces de médiateurs séparées requises, ressource déjà déclarée là-bas).

**Tableau.**

| terme micro | ordre | terme effectif | parasite | statut Gauss |
|---|---|---|---|---|
| hopping intra-bloc `t` | 0 | lui-même |  -  | commute (pair) |
| canaux de paires (intra + frontière) | 2 | `U0`-termes, `NN+YY` | micromouvement `(g/Δ)²` | pairs |
| `λ μ†τ^x + h.c.` | 2 | `τ^xτ^x`, champs | diagonaux, commutant | compensés (`H_contre`) |
| `λη μ†P_B + h.c.` (chaîne densité) | 2 | `τ^xP_B`, `μ_a` effectif | diagonaux, commutant | compensés |
| croisés ordre 3 | 3 | `-K G_B`, `K~λ³η/Δ_G²` | corrections ordre 4 dans le commutant | **cible** |
| tunneling de barreau `ε c†_ac_b` | 0 |  -  |  -  | **anticommute avec `G_B` : chargé** |
| saut simple `a` inter-blocs résiduel `ε'` | 0 |  -  |  -  | chargé ; budget `ε'` = qualité du contrôle par lien |

**Preuve que le tunneling de rail viole `G_B`** : `c†_{a,j}c_{b,j}` change
`N_{a,B}` de `±1` donc `{c†_ac_b, P_B} = 0`, donc anticommute avec `G_B` et
envoie le secteur physique `G_B=+1` dans `G_B=-1`, orthogonal, à coût `2K`.
Il n'existe aucun habillage **local** le rendant invariant : toute corde
`τ^z` a deux extrémités et l'extrémité distante crée une autre violation.
Ce n'est pas une symétrie imposée : c'est une **contrainte énergétique
locale** ; en dessous de `2K`, la surselection est **émergente**.

**Conditions de gap et hiérarchie** (fenêtre pré-enregistrée) :
`Δ_G ≫ λ ≫` (échelles de compensation) ; `K = c₃λ³η/Δ_G²` doit dominer tout
résidu `Z2`-impair : `2K ≫ max(δφ·g²/Δ, ε'·t, flips de frontière Floquet)`.
`K` est petit (ordre 3) : c'est le prix de la dérivabilité ; la protection
est `exp(-2K/T_eff)`-thermique et `(ε/2K)²`-virtuelle (les processus
d'ordre 2 en `ε` retombent dans le commutant : Gauss-pairs, inoffensifs pour
la classe (c), splitting seulement).

## 2R. Réparation (T3) : gadget en chaîne Λ  -  VERDICT : CONSTRUCTION RÉPARÉE

### 2R.0 Règle de sélection et primitive requise

Règle qui interdisait le gadget original : si tous les vertex de conversion
partagent le même opérateur hermitien `A`, le secteur bas ne reçoit que des
mots pairs en `A` → jamais un produit impair de trois involutions distinctes.
Contournement minimal : un chemin bas→bas **connexe** traversant trois
vertex **distincts**. Cela exige une structure interne de médiateur
(≥ 2 espèces)  -  c'est la catégorie « couplage à trois facteurs natif »,
version minimale :

**Primitive nouvelle unique déclarée** : un « marcheur » neutre hard-core à
`b+1` espèces `μ_0..μ_b` par bloc (≤ 1 occupé), avec (i) création/annihilation
aux extrémités flippant le qubit de frontière : `λ_L μ_0†τ^x_{L} + h.c.`,
`λ_R μ_b†τ^x_{R} + h.c.` ; (ii) saut conditionné par parité locale :
`w_j μ_j†μ_{j-1}(1-2n_{a,j}) + h.c.`, `j=1..b`. Vertex 2- et 3-modes, même
classe de localité que `d†aa` existant, mais **non dérivée : nouvelle**.
Détunings `Δ_j > 0` sur chaque espèce. U(1) exacte (marcheur neutre, matière
touchée seulement en diagonal).

### 2R.1 Secteurs et vertex

Secteur bas : marcheur vide. Secteur haut : une espèce occupée, énergie
`Δ_j`. Chemin générateur (pour `b=2`) :
`vide →(λ_L, τ^x_L) μ_0 →(w_1, 1-2n_{a,1}) μ_1 →(w_2, 1-2n_{a,2}) μ_2
→(λ_R, τ^x_R) vide`, produisant l'opérateur
`τ^x_L (1-2n_{a,1})(1-2n_{a,2}) τ^x_R = X_L P_B X_R = G_B`.

### 2R.2 Non-annulation et coefficient

Le chemin est une chaîne de Raman à `b+2` vertex ; l'ordre SW exact est
`b+2` (sélection par le nombre d'occupation du marcheur : tout chemin plus
court doit rebrousser et donne un scalaire). Le terme de SW connexe d'ordre
`b+2` donne, pour `b=2`, `Δ_j = Δ_G` :
\[
K = -\,c\,\frac{\lambda_L w_1 w_2 \lambda_R}{\Delta_G^{3}} + O(\lambda^6/\Delta_G^5),
\]
`c > 0` combinatoire d'ordre 1. Les deux sens de parcours produisent le
**même** opérateur `G_B` (facteurs commutants) et s'ajoutent ; aucune
symétrie ne force `K=0` (contrairement au gadget réfuté, il n'existe pas de
réécriture `f(A²)` : les vertex sont distincts et ordonnés par l'occupation
du marcheur). Signe de `-K G_B` réglé par le signe de `w_1`.

### 2R.3 Tableau des termes d'ordre ≤ b+2 (bas secteur, `b=2`)

| ordre | chemins | opérateur | statut |
|---|---|---|---|
| 2 | créer/annihiler même espèce | `(τ^x)² = 1`, scalaires | inoffensif |
| 3 | aucun (parité de longueur : pas d'annihilation de `μ_1`) |  -  |  -  |
| 4 rebroussés | aller-retour partiels | `(1-2n_j)² = 1`, `(τ^x)² = 1` → scalaires | inoffensif |
| 4 connexe | chaîne complète | `-K G_B` | **cible** |
| blocs voisins | marcheurs indépendants | produits de scalaires | commutant |

**Démonstration de commutation** : tout terme non-cible est un mot rebroussé,
donc un produit de carrés d'involutions = scalaire ; les scalaires commutent
avec tous les `G_B`. Il n'y a **aucun parasite non trivial au même ordre ni
en dessous**  -  rien n'est masqué dans un contre-terme ; `H_contre` du §2
original devient inutile pour ce gadget.

**Pourquoi hors 7C** : 7C interdit de *dériver* un multi-corps sélectif
depuis la grammaire 2-corps figée (cible à ordre strictement supérieur aux
parasites agissant sur le code). Ici la primitive est *déclarée* (pas
dérivée), et la sélection par occupation du marcheur rend les termes
d'ordre inférieur scalaires : la compétition d'ordres de 7C est vide.

### 2R.4 Bloc ED minimal

`b=2` : 4 sites de matière (`N` fixé), 2 qubits `τ`, marcheur à 3 espèces
(≤ 1 occupé). Mesures pré-enregistrées : (1) coefficient de `G_B` en
`λ^{2}w^{2}/Δ_G^{3}` (exposants ajustés `2,2,3 ± 0.05`) ; (2) projection des
termes bas-secteur hors `span{1, G_B}` `< 1e-12` ; (3) gap `2K` entre
secteurs `G_B = ±1` ; (4) injection `ε c†_ac_b` : excitation détectable,
retour virtuel `(ε/2K)²` ; (5) deux blocs : indépendance des marcheurs et
commutation des deux `G_B` effectifs.

### 2R.5 Claim boundary du gadget

Affirmé : existence d'un chemin connexe non annulé produisant `G_B` à
l'ordre `b+2`, coefficient `K` non nul générique, absence de parasite non
scalaire aux ordres ≤ cible, au prix d'**une** primitive nouvelle déclarée
(marcheur à saut conditionné en parité). Non affirmé : la valeur de `c`,
la robustesse de `K` (ordre `b+2` : petit, décroît avec la taille de bloc),
toute phase, mémoire, fusion, indistinguabilité ou tresse ; le test ED
2R.4 est prérequis à toute suite.

## 3. Distinctions obligatoires

- **symétrie globale imposée** : `Z2` de rail actuelle  -  non protectrice
  (verdict T, classe (c)) ;
- **contrainte de jauge locale** : `G_B` ci-dessus, énergétique, dérivée avec
  Hamiltonien explicite  -  protectrice en dessous de `2K` ;
- **protection émergente** : à re-démontrer *après* : le test
  d'indistinguabilité doit être repassé dans le secteur `G=+1` (contrat §4) ;
  non acquise par la seule existence de `G_B` ;
- **post-sélection/mesure** : exclue ici ; la mesure de `G_B` serait une
  ressource distincte (mitigation, cf. verdict T §11.1).

## 4. Contrat ED original  -  obsolète, ne pas exécuter

Bloc : 2 blocs × `b=2` barreaux (8 sites fermioniques), 1 qubit `τ` de
frontière (+2 de bord gelés), 2 médiateurs `μ`, secteur de charge fixe.
1. coefficient de `G_B` : vérifier `K ~ λ^{3.0}/Δ_G^{2.0}` et l'annulation
   des parasites hors commutant (`<1e-12`) ;
2. compensation `H_contre` : splitting intra-secteur résiduel vs budget ;
3. spectre : gap `2K` mesuré entre `G=+1` et `G=-1` ; injecter
   `ε c†_ac_b` et vérifier l'excitation (pas d'action logique) pour
   `ε ≪ 2K`, et la loi `(ε/2K)²` du retour virtuel ;
4. transport inter-blocs par canaux de paires : Gauss-invariance exacte,
   survie de la phase appariée sans saut simple de frontière (question
   ouverte décisive : si la phase meurt, (T3) échoue et la ressource
   minimale remonte au hopping habillé 3-corps natif) ;
5. seulement ensuite : re-test d'indistinguabilité locale complète dans
   `G=+1`.

## 5. Claim boundary original  -  superseded by §2R.5

Affirmé : (T1) preuve qu'une nouvelle ressource par sommet est
indispensable ; (T2) no-go du gauging par barreau à couplages ≤2-corps
(le hopping nu domine son habillage) ; (T3) construction à gros grain
explicite dont la cible `G_B` émerge à l'ordre 3 dans un commutant abélien,
parasites compensables, tunneling de rail prouvé gauge-chargé.
Non affirmé : la survie de la phase appariée sans saut simple de frontière
(test 4), la valeur de `c₃`, une protection émergente démontrée, toute
indistinguabilité, mémoire, fusion, tresse ou jonction T. `K` d'ordre 3 est
petit : la fenêtre `2K ≫` résidus impairs est une hypothèse quantitative à
mesurer, pas un acquis.
