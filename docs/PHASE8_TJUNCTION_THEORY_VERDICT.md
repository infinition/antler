# Phase 8  -  jonction T Floquet à nombre conservé : verdict théorique

Document théorique pré-enregistré, antérieur à tout code de jonction T.
Il répond au mandat : établir ou réfuter une jonction en T Floquet à
conservation exacte de la charge, dérivable des ressources Phase 8, portant un
sous-espace de fusion de dimension 2 et deux échanges adjacents non
commutatifs.

**Révision 2 (2026-07-19).** Intègre l'audit indépendant de la réalisation
directe : les mêmes médiateurs physiques par lien sont réutilisés entre `H0`
et `H1`, et aucune assignation de charge ne fournit de `Z2` microscopique
commune aux deux demi-périodes (§3.5). La `Z2` exacte commune affirmée par la
première version de la table §3.4 est **fausse** pour cette réalisation ;
les §1, 3.4–3.5, 4, 5, 6, 9, 10, 11, 12 sont corrigés en conséquence. Le
no-go du §7 reste inchangé; la construction conditionnelle reste **non
dérivée**. La révision 2 ajoute le test de frontière Rabi-fermée : l'effet
projeté existe dans la réalisation partagée mais varie comme `(g/Delta)^5.93`,
non comme `(g/Delta)^2`; il ne peut donc pas être invoqué comme renforcement
perturbatif générique du no-go.

## 1. Verdict initial

**NO-GO CIBLÉ**, au sens strict des critères pré-enregistrés, assorti d'une
**construction conditionnelle** documentée ci-dessous au niveau effectif,
mais **non dérivée à ce jour** : l'audit des révisions 1--2 (§3.5) montre
qu'elle exige une ressource microscopique supplémentaire absente de la
réalisation directe actuelle.

Énoncé exact du no-go (périmètre délimité, cf. §7) :

> Avec les ressources autorisées 1–6 (hopping ANTLER à `theta=pi`, médiateurs
> charge-2 séparés positivement détunés, canaux de conversion de paires
> programmables, U(1) exacte à N fixé, graphe en T de largeur bornée,
> modulation Floquet des canaux déclarés), aucun sous-espace de fusion de
> dimension 2 ne peut satisfaire le test d'indistinguabilité locale **pour
> tout opérateur local conservant la charge** (critère pré-enregistré 4).
> Le nombre quantique qui distingue les deux états de fusion est la parité
> relative de rail `Z2 = (-1)^{N_a}` : c'est une symétrie d'ingénierie, non
> surselectionnée. Il existe des opérateurs locaux, hermitiens, U(1)-exacts et
> `Z2`-impairs  -  le tunneling de barreau `c†_a c_b + h.c.` en est le
> représentant minimal  -  dont l'action projetée sur le code est d'amplitude
> O(1) au voisinage d'un défaut, sans décroissance avec la séparation des
> défauts. Le blocage est **algébrique** (absence de surselection), non
> perturbatif, non spectral, non propre à la géométrie T.

Ce qui reste vrai et testable : la restriction de symétrie « `Z2` exacte »
échoue désormais à **deux niveaux distincts**  -  (i) contrôle : le budget
d'angle de canal mesuré `±0.003 rad` est un budget de violation de `Z2` ;
(ii) **structurel** : la réalisation directe actuelle réutilise les mêmes
médiateurs entre `H0` et `H1`, et aucune assignation de charge fixe ne rend
une `Z2` microscopique exacte sur le cycle complet (§3.5). La chaîne
microscopique → SW → Floquet → doublet de bord → protocole d'échange est
écrite ici **au niveau effectif**, avec lois d'erreur explicites ; elle ne
devient un candidat falsifiable de **mémoire protégée par symétrie** (niveau
SPT, pas mémoire topologique) qu'après acquisition de la ressource d'étage 1
du §11.0 (espèces de médiateurs à charges fixes, symétrie dynamique dérivée,
ou projection du secteur médiateur). Le mode neutre gapless n'est pas un
argument de rejet (cf. §8) : il dégrade la protection d'exponentielle à
polynomiale en temps de protocole, il ne l'interdit pas.

Le résultat 1 (« CONSTRUCTION VIABLE » au sens des neuf critères) est donc
réfuté pour les ressources actuelles ; le résultat 2 est établi avec une
identification **à deux étages** de ce qui manque (§11) : (é1) pour la
construction conditionnelle elle-même, des espèces de médiateurs distinctes
à charges `Z2` fixes, ou équivalent §11.0  -  ressource **nouvelle**, à
déclarer, jamais réputée présente ; (é2) pour lever le no-go proprement dit,
une **ancilla de jauge `Z2` avec loi de Gauss locale** (option 5, exposée à
la hiérarchie d'ordres du no-go Phase 7C), la **mesure de parité de rail**
(option 2) restant une **mitigation candidate** dont le protocole
measurement-based n'est pas dérivé  -  pas une ressource minimale suffisante.

## 2. Portée exacte, substrat et statut de la géométrie T

Trois mises au point de périmètre, obligatoires avant la dérivation.

**Substrat.** ANTLER v0.7 est un modèle abstrait de bosons hard-core en
espace de Fock avec corde de Jordan–Wigner fractionnaire définie sur un
ordre linéaire rung-major (`antler/model.py`). À `theta=pi`, il coïncide
avec des fermions spinless sur cet ordre (audit `2.45e-16`). Toute la
dérivation ci-dessous est menée dans le **langage fermionique abstrait**, où
un graphe trivalent est local sans obstruction. Le problème de corde JW sur
graphe trivalent n'est **pas** une obstruction native : il devient une
contrainte matérielle seulement si un substrat physique spin/qubit est choisi
plus tard, et devra alors être classé comme ressource à cette étape-là.
Dans la représentation native hard-core-boson, les liens de nœud des bras 2
et 3 portent une corde `(-1)^{N_{\rm intercalé}}` **par définition** de
l'extension de modèle (cf. orientation JW, §3).

**Géométrie T.** Le graphe en T, les couplages de nœud de degré 3 et leurs
médiateurs ne sont implémentés nulle part dans `antler/` (le `hop_list` est
strictement linéaire SSH). Ils sont classés ici **EXTENSION DE GRAPHE À
DÉRIVER** : le §3 en donne la définition complète, le contrat §10 en impose
l'audit avant toute promotion. Aucune capacité de nœud n'est réputée acquise.

**Chiffres de contrôle Phase 8.** Les tolérances mesurées (`±0.003 rad`,
crosstalk `0.3%`, rampes `1%`) sont des budgets de contrôle sur petits blocs
idéaux. Elles sont citées comme **corroboration numérique** de la sensibilité
au premier ordre prédite par le no-go, jamais comme sa preuve. La preuve est
l'argument algébrique et le calcul de `P O P` du §6, et le contrat §10 en
donne le test de réfutation direct.

## 3. Graphe T minimal et Hamiltonien microscopique

### 3.1 Modes

Trois bras `m ∈ {1,2,3}`, chacun un ladder à deux rails de `L` barreaux,
plus un barreau de nœud unique `0` partagé :

- sites fermioniques physiques : `(0,ℓ)` et `(m,j,ℓ)`, `j ∈ 1..L`,
  `ℓ ∈ {a,b}` ; total `M = 2(3L+1)` modes de charge 1 ;
- liens de l'arbre : intra-bras `e = (m,j)–(m,j+1)` et liens de nœud
  `e_m = 0–(m,1)`, `m = 1,2,3` ; total `3L` liens, arbre sans cycle ;
- médiateurs hard-core de charge 2 : deux par lien pour le segment `H0`
  (un par canal), idem pour le segment `H1` ; détuning positif `Δ` commun.
  Un médiateur est assigné à un lien et à un canal, jamais partagé
  (condition structurelle héritée du pont Phase 8 : elle annule le transfert
  de paires parasite au second ordre).

### 3.2 Orientation Jordan–Wigner

Ordre linéaire déclaré : barreau de nœud d'abord, puis bras 1 complet,
bras 2, bras 3, rung-major à l'intérieur de chaque bras :

```
k(0,ℓ) = ℓ ;  k(m,j,ℓ) = 2 + 2(m-1)L + 2(j-1) + ℓ .
```

Conséquences exactes à `theta=pi` :

- tous les liens intra-bras et le lien de nœud `e_1` ont la structure de
  phase du ladder 1D existant (saut de jambe = corde sur un site) ;
- les liens de nœud `e_2` et `e_3` portent, dans la représentation native,
  la corde `(-1)^{N_{\rm bras\,1}}` resp. `(-1)^{N_{\rm bras\,1}+N_{\rm bras\,2}}` ;
  dans l'espace de Fock fermionique (convention `_apply`), ces signes sont
  l'antisymétrie standard et aucun terme non local n'apparaît. C'est la
  **définition** de l'extension de graphe, à auditer par le contrat §10
  (invariance du spectre sous permutation de l'ordre des bras).

### 3.3 Hamiltonien microscopique

Pour chaque segment Floquet `s ∈ {0,1}` avec ensembles de canaux
orthonormés dans la base de paires `(aa, ab, ba, bb)` d'un lien :

```
s=0 : C_1 = (aa - bb)/sqrt(2),  C_2 = (aa + bb)/sqrt(2)
s=1 : C_1 = (aa - bb)/sqrt(2),  C_2 = (ab + ba)/sqrt(2)
```

\[
H^{(s)}_{\rm micro} = -t\sum_{e}\sum_{\ell}(c^\dagger_{e^-,\ell}c_{e^+,\ell}+h.c.)
 + \Delta\sum_{e,s'}d^\dagger_{e,s'}d_{e,s'}
 - g\sum_{e,s'}\left[d^\dagger_{e,s'}\,\hat C^{(s)}_{s'}(e)+h.c.\right]
\]

où `Ĉ^{(s)}_{s'}(e)` est l'opérateur d'annihilation de paire du canal `s'`
sur le lien `e`, y compris les trois liens de nœud. Les phases autorisées
sont les phases des coefficients de canaux (contrôle déclaré) et les phases
Peierls de jambe déjà auditées en Phase 7D.

### 3.4 Symétries exactes du modèle idéal

| symétrie | définition | statut |
| --- | --- | --- |
| U(1) totale | `N_tot = Σ n + 2 Σ n_d` | exacte par construction, tous segments |
| `Z2` globale de rail | `(-1)^{N_a,tot}`, extensions médiateur `Q_x` du §3.5 | **trois niveaux distincts, à ne jamais confondre (§3.5)** : exacte au niveau effectif projeté ; approximative au niveau stroboscopique ; **absente** au niveau microscopique commun dans la réalisation directe actuelle (révision 2) |
| `Z2_m` par bras | `(-1)^{N_a,m}` | brisée par les liens de nœud  -  voulu : le transport inter-bras l'exige ; seule la `Z2` globale survit |
| `S3` géométrique | permutation des bras à `L`, couplages égaux | exacte dans le modèle idéal ; diagnostic utile, non protectrice |
| conjugaison complexe | base de Fock réelle | vraie pour `H0` et pour `H1` aux points de canaux déclarés ; à vérifier numériquement, non porteuse |

Il n'existe **aucune autre** symétrie discrète locale disponible dans la
grammaire : c'est le point d'appui du no-go (§7).

### 3.5 Trois niveaux de `Z2`  -  correction d'audit (révision 2)

La famille de candidats est `Q_x = (-1)^{N_a + Σ_m x_m n_{d,m}}` avec
charges fixes `x_m ∈ {0,1}` par médiateur physique. Trois notions distinctes :

1. **`Z2` effective projetée.** Sur le secteur sans médiateur,
   `Q = (-1)^{N_a}` commute exactement avec `H_eff` : tous les termes
   effectifs des deux segments sont `Z2`-pairs (§4.1). Exacte, mais définie
   seulement **après** projection SW  -  elle n'existe pas comme symétrie de
   l'évolution microscopique.
2. **`Z2` stroboscopique approximative.** L'absence du niveau 3 permet en
   principe une violation projetée, mais son ordre ne se déduit pas de la
   seule occupation virtuelle. Le test explicite Rabi-fermé `L=3,N=2`
   (`U0=-1.5`, `g/Δ=0.20..0.0125`) mesure une pente `5.931` du flip projeté,
   et non `2`. Ainsi, la loi auparavant annoncée `O((g/Δ)^2)` est réfutée
   pour cette séquence idéale; l'accumulation pour une autre séquence reste
   à mesurer. Voir `PHASE8_MICROSCOPIC_PARITY_AUDIT.md`.
3. **`Z2` microscopique exacte commune.** Exigerait un `x` fixe commutant
   avec les deux demi-périodes. Audit exhaustif `L=3,N=2` de la réalisation
   directe : `H0` commute exactement pour `x=[0,0,0,0]` seulement, `H1`
   pour `x=[0,1,0,1]` seulement, **aucune assignation ne commute avec les
   deux** ; meilleur compromis, résidu normalisé `7.05e-2`. La cause est
   structurelle : le même médiateur physique porte le canal `(aa+bb)`
   (parité `N_a` paire) pendant `H0` et `(ab+ba)` (impaire) pendant `H1` ;
   une charge de symétrie ne peut pas changer avec le temps pour un même
   mode sans nouvelle structure physique. **Ce niveau est absent de la
   réalisation directe actuelle.**

Toute la suite du document repose au mieux sur le niveau 1. Le candidat
conditionnel (doublet SPT, protocole d'échange) exige le niveau 3, donc
l'une des ressources nouvelles de l'étage 1 du §11.0  -  jamais réputée
présente.

## 4. Réduction Schrieffer–Wolff et composition Floquet

### 4.1 SW ordre 2, lien par lien

Chaque médiateur étant séparé et positivement détuné, le générateur SW
standard donne dans le secteur sans médiateur, par lien et par canal :

\[
H^{(2)}_{\rm SW} = -\frac{g^2}{\Delta}\,\hat C^\dagger \hat C
 + O\!\left(\frac{g^4}{\Delta^3}\right).
\]

Sommes exactes des canaux (identités algébriques, vérifiées en Phase 8 à
`3.5e-16` sur le bloc local) :

- segment 0 : `Σ C†C = n_{a,j}n_{a,j+1} + n_{b,j}n_{b,j+1}` ; le terme
  croisé `(aa)†(bb)` s'annule exactement entre les deux canaux ;
- segment 1 : `Σ C†C = (1/2)[N_j N_{j+1} + Y_j Y_{j+1}]` avec
  `N = n_a+n_b`, `Y = -i c†_a c_b + i c†_b c_a` ; ce sont les termes de
  `P† H0 P` à `eta=pi/2`, sans impulsion `P` physique. Le développement de
  `YY` contient le pair-hopping inter-rail `a†a†bb + h.c.`
  (`ΔN_a = ±2`, `Z2`-pair) et l'échange `(ab)†(ba)` (`ΔN_a = 0`).

Sur les liens de nœud, la même réduction s'applique terme à terme : le
médiateur d'un lien `e_m` ne peut reconvertir sa paire que sur `e_m`
(médiateurs séparés), donc **aucun transfert de paire inter-bras n'apparaît
au second ordre**. Les premiers termes propres au nœud sont d'ordre
`g^4/Δ^3` : densités assistées à trois bras corrélées par l'occupation du
barreau 0, et boucles virtuelles à deux médiateurs de nœud. Ils sont
`U(1)`- et `Z2`-exacts mais brisent `S3` si les calibrations diffèrent ;
ils sont à borner numériquement (contrat §10, étape M1).

### 4.2 Composition Floquet

Alternance `H^{(0)}_{\rm eff}` (durée `αT`) / `H^{(1)}_{\rm eff}`
(durée `(1-α)T`) :

\[
H_F = \alpha H_0 + (1-\alpha)\,H_1
 - i\,\frac{T}{2}\,\alpha(1-\alpha)\,[H_1, H_0] + O(T^2),
\qquad H_1 = P^\dagger H_0 P .
\]

Le facteur `-i` est nécessaire : `[H_1,H_0]` est anti-hermitien, et c'est
`-i[H_1,H_0]` qui est hermitien (correction de la première version).

Points de contrôle explicites, hérités des audits enregistrés et repris
comme conditions du candidat :

- `U0 = -g^2/Δ` fixé (point de référence `-1.5`) ; profondeur SW
  `g/Δ ≤ 0.0125` (`Δ=9600`, `g=120`) ;
- fermeture Rabi virtuelle `T_m = 4πm/sqrt(Δ² + 4g²)` sur chaque
  demi-segment ; fenêtre de temps mesurée `±0.3%` ;
- rampes de canaux `sin²`, `≤ 1%` de période par transition ;
- erreur d'angle de canal `|δφ| ≤ 3e-3 rad`  -  noter que cette condition est
  un budget de **violation de `Z2`**, pas seulement de fuite ;
- crosstalk de lien voisin `ε ≤ 3e-3` ;
- le terme `-i[H_1,H_0]` est `Z2`-pair **au niveau effectif** (niveau 1 du
  §3.5, commutateur de deux termes pairs) : l'erreur Magnus du premier ordre
  ne brise pas la `Z2` projetée ; une composition Strang symétrique
  l'annule à cet ordre ;
- chauffage Floquet : `H_eff` locale bornée (hard-core), fréquence
  `ω = 2π/T ~ Δ/2m` très supérieure à la bande locale → borne de chauffage
  exponentiellement lente applicable ; compatible avec le point SW profond.

**Frontière de segment : correction numérique (révision 2).** L'absence
d'un `Q_x` microscopique commun est exacte (§3.5), mais elle ne fixe pas à
elle seule l'ordre du défaut stroboscopique. Le mécanisme de reconversion
croisée `(ab+ba) → (aa±bb)` est un chemin permis dans le modèle partagé; il
ne justifie cependant pas une amplitude projetée universelle
`O((g/Δ)^2)`. Pour le cycle Rabi-fermé enregistré, le calcul explicite donne
un flip bas-espace impair-vers-pair et un commutateur de parité polaire qui
varient tous deux comme `(g/Δ)^5.93`, pas comme `(g/Δ)^2`. Les répétitions
à `2` et `4` cycles, qui traversent aussi la frontière inverse, donnent les
pentes profondes `5.984` et `6.018`. Le contrôle à
quatre espèces, avec charges fixes `[0,0]` pour `H0` et `[0,1]` pour `H1`,
commute avec un `Q` microscopique à moins de `1.7e-17` et annule le flip à
l'arrondi. La séparation des espèces est donc un contrôle de symétrie valide;
la loi d'ordre deux et son interprétation comme action logique générique sont
retirées. Les séquences hors fermeture et la géométrie de nœud restent non
testées. Voir `PHASE8_MICROSCOPIC_PARITY_AUDIT.md`.

## 5. Tableau obligatoire

| terme microscopique | ordre SW | terme effectif | U(1) | parités | parasite associé | amplitude/échelle | condition de contrôle |
| --- | --- | --- | --- | --- | --- | --- | --- |
| saut de jambe intra-bras, `theta=pi` | 0 | `-t c†c + h.c.` intra-rail | exacte | `Z2` paire, `Z2_m` paire | aucun au même ordre | `t` | audit `2.45e-16` (acquis) |
| saut de jambe des liens de nœud `e_m` | 0 | idem, degré 3 au barreau 0 | exacte | `Z2` paire ; brise `Z2_m` (voulu) | corde `(-1)^N` en représentation native (définitionnel) | `t` | **extension de graphe à auditer** |
| détuning `Δ n_d` | 0 | projection secteur sans médiateur | exacte (charge 2) | selon canal | micromouvement virtuel, fuite `~n^2` par cycles | `Δ` | fermeture `T_m`, fenêtre `±0.3%` |
| conversion canaux `s=0` `(aa±bb)/√2`, tous liens | 2 | `U0 (n_a n_a + n_b n_b)` | exacte | `Z2` paire | `(aa)†(bb)` exactement nul ; résidu `~(g/Δ)^{1.97}` mesuré | `g²/Δ` | `g/Δ ≤ 0.0125` |
| conversion canaux `s=1` `{(aa-bb),(ab+ba)}/√2` | 2 | `(U0/2)[NN + YY]` | exacte | `Z2` **effective** paire (niveau 1) ; pas de `Z2` microscopique commune avec `s=0` sur médiateurs partagés (niveau 3 absent, §3.5) | erreur d'angle `δφ` → terme effectif `Z2`-impair `~ (g²/Δ)δφ` | `g²/Δ` | `|δφ| ≤ 3e-3 rad` |
| réutilisation des mêmes médiateurs entre segments | frontières de segment | chemins croisés permis; aucun `Q_x` microscopique commun | exacte | niveau 3 absent (`L=3,N=2`, meilleur résidu `7.05e-2`) | sous le cycle Rabi-fermé, flip projeté `~(g/Δ)^5.93`; aucune loi universelle ni action logique de classe (c) n'est déduite | séquence enregistrée seulement | **espèces séparées : contrôle exact; nœud/hors-fermeture encore à auditer (§11.0)** |
| conversion sur liens de nœud | 2 | idem sur `e_1,e_2,e_3` | exacte | `Z2` globale paire | ordre 4 : termes à trois bras assistés par `n_0`, boucles à deux médiateurs `~g⁴/Δ³` | `g²/Δ` | à borner (contrat M1) |
| crosstalk résiduel `εg` lien voisin | 2 | transfert de paires inter-liens `~ε g²/Δ` | exacte | `Z2`-impair possible si canaux de familles différentes voisinent (nœud !) | erreur cohérente logique, non détectable comme perte | `ε g²/Δ` | `ε ≤ 3e-3` ; au nœud : à re-borner |
| alternance Floquet | Magnus 1 | `-i(T/2)α(1-α)[H1,H0]` | exacte | `Z2` effective paire (niveau 1) | accumulation cohérente `n^{1.996}` mesurée | `O(T·U0·t)` par site | Strang + fermeture ; contre-termes ouverts |
| tunneling de barreau parasite `J⊥ c†_a c_b` | 0 | lui-même | exacte | **`Z2`-impair** | **opérateur logique `X` au défaut, O(1)** (§6) | `J⊥` résiduel | aucune suppression topologique disponible : c'est le no-go |

La dernière ligne est la ligne décisive : elle est au même ordre (zéro) que
les termes utiles, conserve exactement la charge, et aucune ressource
autorisée ne l'interdit ni ne la supprime avec la distance.

## 6. Sous-espace de fusion : dimension, gaps, localité

### 6.1 Secteur relatif et dimension

Bosonisation par bras (rails `a,b` → champs `φ_ρ,θ_ρ` totaux et
`φ_σ,θ_σ` relatifs ; validité au couplage faible, les conclusions étant
re-testées par ED indépendamment de cette validité) :

- le canal `YY` attractif génère `cos(2√2 θ_σ)` : il **épingle la phase
  relative** `θ_σ` sur deux minima images l'un de l'autre sous
  `Z2 : θ_σ → θ_σ + π/√2` ;
- le secteur relatif est donc gappé, de type Ising, avec longueur `ξ_σ` ;
- le secteur total `ρ` reste une ligne de Luttinger `c=1` : le **mode neutre
  gapless est un fait attendu de la classe**, conforme à la référence
  externe, et n'est pas un critère de rejet (garde-fou enregistré) ;
- il n'y a **pas** d'ordre à longue portée `Z2` : tout opérateur local
  contenant `e^{i√2θ_σ}` porte nécessairement un facteur dual
  `e^{±i√2φ_σ}` désordonné (le fermion est produit ordre × désordre), donc
  `⟨c†_a c_b⟩_{\rm bulk} = 0` et les corrélations de barreau décroissent
  exponentiellement. Le doublet est de type **bord/défaut**, pas un chat de
  bulk. C'est ce qui rend la construction conditionnelle non triviale.

Défauts : parois entre régions programmées topologiques (canaux `s=1`
actifs) et triviales (hopping dominant), positions contrôlées par les
amplitudes de canaux lien par lien  -  localité des défauts garantie par
construction, largeur de paroi fixée par le profil de rampe.

Comptage de dimension (correspondance secteur relatif ↔ chaîne de Kitaev
effective par bras, utilisée **comme cible de diagnostic**, conformément aux
interdictions) : trois bras topologiques joints au nœud → trois modes de
paroi externes plus un mode de nœud survivant génériquement à
l'hybridation de degré 3 → quatre modes, et à `(N, Z2_{\rm globale})` fixés :

\[
\dim \mathcal{H}_{\rm fusion} = 2 .
\]

Statut : **plausible, dérivé sous deux hypothèses explicites**  -  (i)
hybridation générique au nœud ne laissant qu'un mode (condition spectrale à
vérifier), (ii) `Z2` exacte **au niveau 3 du §3.5**, c'est-à-dire ressource
d'étage 1 du §11.0 acquise  -  ce qui n'est pas le cas de la réalisation
directe actuelle. La dimension est un prédicat du contrat §10 (étape T2),
pas un acquis.

### 6.2 Splitting interne et gap chargé

- splitting interne prédit : `~ exp(-L/ξ_σ)` avec préfacteur algébrique
  possible du secteur `ρ` (exposant fonction de `K_ρ`) ; la référence 1D
  externe donne `9.10e-3 → 1.36e-3 → 2.01e-4` pour `L=8,12,16`  -  compatible,
  loi exacte à mesurer sur le T ;
- gap chargé (excitations dangereuses = ajout/retrait de charge et
  quasi-particules chargées) : attendu fini, courbure `0.36–0.38` aux points
  publiés ; sur le T, à mesurer par `E0(N±1)` (contrat T2) ;
- gap neutre : `~ πv_ρ/L_tot`, ferme comme `1/L`  -  attendu, non disqualifiant,
  quantitativement contraignant pour l'adiabaticité (§8).

### 6.3 Localité et indistinguabilité : calcul de `P_code O P_code`

Norme déclarée : pour `O` local hermitien, `‖O‖=1`, support à distance `d`
du défaut le plus proche, séparation inter-défauts `D ≫ d` :

\[
\mathcal{D}(O) = \left\| P\,O\,P - c_O P \right\|_2,
\qquad c_O = \tfrac{1}{2}\mathrm{tr}(P O P).
\]

Trois classes, avec lois dérivées du §6.1 :

**(a) `O` `Z2`-pair** (densités, `n_a n_a`, pair-hopping, échange) :
l'élément hors-diagonal est nul par symétrie ; la différence diagonale ne
peut provenir que du contenu en secteur relatif, gappé :

\[
\mathcal{D}(O) \le C\, e^{-d/\xi_\sigma}\times\left(1 + O(L^{-\eta(K_\rho)})\right).
\]

Décroissance exponentielle en `d`, préfacteur algébrique neutre possible.
**Test satisfaisable.**

**(b) `O` `Z2`-impair, support en bulk** (`c†_a c_b` loin de toute paroi) :
tout terme de `O` porte le champ de désordre `e^{±i√2φ_σ}` ; connecter les
deux vacua épinglés exige le mode zéro, absent du support :

\[
\mathcal{D}(O) \le C'\, e^{-d/\xi_\sigma}.
\]

Décroissance également exponentielle. **Le bulk est sain.**

**(c) `O` `Z2`-impair, support sur une paroi** (`d ≲ ξ_σ`) : à la paroi, la
corde de désordre se termine ; l'opérateur de barreau local se réduit, avec
coefficient O(1), à l'opérateur de mode zéro habillé (l'analogue exact de
l'opérateur d'ordre de bord d'une chaîne d'Ising ouverte, dont l'élément de
matrice entre les deux quasi-fondamentaux est l'aimantation locale `m`,
d'ordre 1). Comme `c†_a c_b` est neutre en charge totale, le secteur `ρ`
gapless n'apporte **aucune suppression** :

\[
\mathcal{D}(O_{\rm paroi}) = O(1),
\quad \text{indépendant de } D \text{ et de } L.
\]

Ce n'est pas une hypothèse : c'est la phénoménologie même du mode zéro
(sans mode zéro de paroi à action O(1), pas de doublet du tout), corroborée
par le benchmark externe déjà archivé  -  transfert local de bord O(1)
(`0.200` à `L=4`) décroissant vers le bulk (`0.0769` à `L=8`),
`IEMINI_BRAID_AUDIT.md`. La prédiction quantitative complète (saturation en
classe (c), décroissance en classes (a),(b)) est pré-enregistrée comme test
T3 du contrat §10 : **si l'ED de la jonction T trouvait la classe (c)
décroissante avec `D`, le no-go serait réfuté** et la route promue.

Conclusion du critère 4 : l'indistinguabilité locale vaut pour la
sous-algèbre `Z2`-paire avec loi `exp(-d/ξ_σ)`, et **échoue irréductiblement
pour l'algèbre locale complète** à cause de la classe (c). Sous U(1) exacte
à `N` fixé, la parité totale fermionique est figée par `N` : il n'existe
aucune surselection résiduelle pour interdire la classe (c), et aucune autre
symétrie discrète locale n'est disponible dans la grammaire (§3.4) pour
recoder le doublet ailleurs que dans une parité non surselectionnée.

## 7. No-go ciblé

**Périmètre.** Ressources Phase 8 (1–6) + jonction T minimale de largeur
bornée + U(1) exacte à `N` fixé + exigence de protection contre **tous** les
opérateurs locaux conservant la charge. Ce no-go ne porte **pas** sur : les
systèmes U(1)-conservants en général, la classe 2D à ordre topologique
intrinsèque à charge fixée (le parent torique 2D de Phase 7B y échappe
précisément  -  c'est cohérent), ni la valeur de la route comme mémoire SPT
sous `Z2` exacte, ni la référence externe Defossez et al.

**Réponses aux questions imposées.**

- *Quel mécanisme manque ?* La surselection du nombre quantique logique. Le
  doublet est distingué par `(-1)^{N_a}`, symétrie d'ingénierie que tout
  opérateur local de mélange de rails viole légalement.
- *Impossibilité de quoi ?* Ni de la géométrie T (le nœud est dérivable,
  §3–4), ni de U(1) en soi, ni du protocole Floquet (l'alternance préserve
  `Z2` au premier ordre Magnus). C'est le **contenu microscopique**  -  deux
  rails discernables, aucun secteur de charge 2 condensé, aucune contrainte
  de jauge  -  qui ne fournit aucun mécanisme interdisant la classe (c).
- *Nature du blocage ?* **Algébrique** (structure de surselection de
  l'algèbre d'observables locale). Non perturbatif : aucune profondeur SW ne
  le réduit. Non spectral : il persiste à gap chargé fini et splitting
  exponentiellement petit. Non topologique de la géométrie : il est déjà
  présent sur le ladder 1D ; la jonction T n'ajoute que des parasites
  d'ordre 4 contrôlables.
- *Pourquoi ciblé et pas général ?* En 2D, un code à charge fixée peut
  stocker l'information dans des secteurs de boucles non locaux sans
  étiquette de symétrie locale (Phase 7B). En largeur bornée (arbre), les
  phases gappées n'offrent que des dégénérescences de type SSB (localement
  lisibles, rejetées) ou SPT (protégées par une symétrie non
  surselectionnée) : le T n'ouvre pas de troisième voie. La formalisation
  complète de cet argument de classification est le théorème Phase 7B
  encore en attente ; le présent no-go n'en dépend pas  -  il repose sur le
  calcul direct §6.3(c) pour l'encodage effectivement disponible.

## 8. Mode neutre gapless et adiabaticité : quantification, pas verdict

Conformément au garde-fou : le mode `c=1` n'est pas un no-go. Couplage et
échelles explicites, sous `Z2` exacte (niveau 3 du §3.5, c'est-à-dire
ressource §11.0 supposée acquise  -  hypothèse déclarée, non disponible) :

- une paroi en mouvement est un potentiel local dépendant du temps ; elle
  porte un profil de charge `δq = O(1)` et se couple linéairement à
  `∂_x φ_ρ` (densité). Les états logiques ne diffèrent que dans le secteur
  `σ` : le secteur `ρ` est **aveugle à la parité**  -  il ne lit pas le qubit,
  il ne fait que désadiabatiser et habiller ;
- taille finie : plus petit gap le long du chemin
  `Δ_n(L) = π v_ρ / L_tot` (critère 5 : le « gap minimal le long des
  chemins » est cette échelle, fermant en `1/L`  -  à déclarer comme tel, pas
  comme gap thermodynamique) ;
- rampe lisse `C^k` de durée `τ` : seuls les modes `ω ≲ 1/τ` absorbent ;
  aucun mode disponible dès `τ ≳ L_tot/v_ρ`, puis borne adiabatique standard

\[
\varepsilon_{\rm dia} \lesssim
\left[\frac{\|\partial_\lambda H\|\,\ell}{\tau\,\Delta_n(L)^2}\right]^2
\sim \left[\frac{\delta q\,\ell\, L_{\rm tot}^2}{v_\rho^2\,\tau}\right]^2 ,
\]

  soit `τ = poly(L)` : coût polynomial, pas exponentiel ;
- fenêtre de protocole : `L_tot/v_ρ ≪ τ ≪ ħ/ΔE_{\rm split} ~ e^{+L/ξ_σ}`  - 
  fenêtre large, le protocole existe ;
- phases abéliennes habillées par le secteur `ρ` : scalaires globaux, elles
  s'annulent exactement dans `‖[U1,U2]‖` et dans le résidu Yang–Baxter
  (les deux membres portent la même phase totale)  -  les métriques
  pré-enregistrées y sont insensibles ;
- budget Floquet composé : `n_cycles ~ τ/T` et fuite cohérente mesurée
  `~ n²(g/Δ)^{3.9}` à durée logique fixe ; la condition conjointe
  `ε_dia + P_leak < seuil` définit la profondeur SW requise en fonction de
  `L`  -  polynomiale également.

Conclusion : sous `Z2` exacte, l'échange adiabatique est **viable à coût
polynomial** malgré le mode neutre. La dégradation d'exponentiel à
polynomial doit être déclarée dans toute claim boundary, rien de plus.

## 9. Échanges adjacents : protocole défini, non promu

Le critère d'entrée de la section C (« seulement si le sous-espace existe et
est isolé ») **n'est pas satisfait** au sens du critère 4 (§6.3), et, depuis
les révisions 1--2, la restriction `Z2` déclarée n'est elle-même pas réalisable
sans la ressource d'étage 1 du §11.0 (§3.5). Aucun générateur logique,
aucune holonomie, aucune matrice d'échange n'est donc écrite ici comme
résultat. Ce qui suit est la **définition du protocole** qui serait mesuré
si le contrat §10 passait ses étapes T2–T3, sous la restriction `Z2`
déclarée **et** la ressource §11.0 déclarée et auditée :

- `U1` : échange des parois externes des bras 1 et 2 par la manœuvre en
  trois temps standard (parking sur le bras 3), réalisée exclusivement par
  interpolation lien-par-lien des amplitudes de canaux (ressource 6),
  rampes `sin²`, durée `τ` dans la fenêtre §8 ;
- `U2` : idem pour les bras 2 et 3 (parking bras 1) ;
- extraction : transport de Kato discrétisé du projecteur `P(λ)` le long du
  chemin (le harnais Phase 5 existe), matrices logiques **dérivées** de
  l'évolution, jamais imposées ;
- métriques pré-enregistrées : `‖[U1,U2]‖` brut ; résidu
  `‖U1U2U1 - U2U1U2‖` déclaré interprétable **uniquement si**
  `‖[U1,U2]‖ > 0.5` (la règle « commutateur quasi nul → YB non
  interprétable » reste en vigueur, cf. la falsification N=3 archivée) ;
  séparation erreur d'adiabaticité / fuite / erreur de projection / erreur
  de discrétisation par les quatre estimateurs déjà utilisés en Phase 5 ;
  gap minimal le long du chemin rapporté **par secteur** (chargé, relatif,
  neutre-taille-finie) et jamais agrégé.

## 10. Contrat ED/MPS pré-enregistré

Ordre impératif : géométrie T **effective d'abord, sans médiateurs
explicites** ; microscopique ensuite ; holonomies en dernier.

**T0  -  modèle effectif et symétries survivantes.** `H_eff` sur le T,
bras égaux `L ∈ {2,3,4}`, `M = 2(3L+1) ∈ {14,20,26}` sites. Secteurs :
`N` pair prioritaire (à `N` impair, la symétrie d'échange des rails force
l'égalité des secteurs de parité  -  caveat déjà enregistré) ; points
principaux `L=3, N ∈ {4,6}`, contrôles `L=2, N=4` et `L=4, N=8`
(`ν ≈ 0.29–0.31` apparié). Dimensions `C(20,6)=38 760` (ED dense),
`C(26,8)≈1.56e6` (Lanczos sparse). Vérifier : `[H, Z2] < 1e-12` ;
non-conservation effective des `Z2_m` ; multiplets `S3` ; réalité.
Audit de convention : spectre invariant sous permutation de l'ordre JW des
bras (`< 1e-12`), garde-fou de l'extension de graphe §3.2.

**T2  -  spectre.** Par secteur `(N, Z2)` : splitting du doublet, gap vers la
troisième valeur propre, **gap chargé** `E0(N+1)+E0(N-1)-2E0(N)`, gap
neutre vs `1/L`. Seuils de rejet du candidat conditionnel :
`split/gap ≥ 1e-2` à `L=3` ou non-décroissance en `L`. Le gap neutre ne
rejette rien (garde-fou §8) ; il calibre `v_ρ`.

**T3  -  base complète d'opérateurs locaux projetés.** Pour chaque barreau et
chaque lien : base hermitienne complète des opérateurs conservant la charge
(harnais exhaustif Phase 7 réutilisé), `D(O)` mesuré, classé par parité `Z2`
et distance `d` au défaut le plus proche. Prédictions pré-enregistrées :
classes (a),(b) décroissantes (ajustement exponentiel, `ξ_σ` extrait) ;
classe (c) plate et `≥ 0.3`. **Critère de réfutation du no-go : classe (c)
décroissante avec la séparation**  -  dans ce cas la route est promue et ce
document amendé.

**T4  -  holonomies (conditionnel : T2 et T3(a,b) passés, ressource §11.0
déclarée et auditée, restriction `Z2` déclarée).** Chemins de §9, transport
de Kato, `‖[U1,U2]‖`, résidu YB sous la règle d'interprétabilité, scaling en
`τ` de `ε_dia` (vérifier la loi §8), gaps par secteur le long du chemin.

**M1  -  microscopique de nœud (dernier).** Bloc explicite : barreau 0 + un
barreau par bras (8 sites) + 6 médiateurs de nœud, base pondérée à charge
fixée. Vérifier : absence de transfert de paire inter-bras à l'ordre 2
(`< 1e-12`), scaling `g^4` des termes de nœud d'ordre 4, budget de
crosstalk **inter-familles** au nœud (le crosstalk `Z2`-impair est une ligne
dangereuse du tableau §5  -  seuil à re-mesurer, le budget 1D `3e-3` n'est pas
transférable d'office). Correction révision 2 : le scan Rabi-fermé `L=3,N=2`
mesure le flip projeté des médiateurs partagés comme `(g/Δ)^5.931`, et non
comme `(g/Δ)^2`; la variante à espèces séparées commute avec `Q` à moins de
`1.7e-17`. M1 doit mesurer séparément les séquences hors-fermeture, les
rampes et le nœud; aucune extrapolation de cette pente locale n'est permise.

**MPS.** Ordre serpent bras-concaténés, couplages de nœud en termes MPO de
portée `~L` ; charges `U(1) × Z2` dans les tenseurs (infrastructure de
l'audit canonique réutilisée). `L` jusqu'à 12–16, mêmes observables T2–T3,
convergence en `χ` sérialisée (`χ ≥ 384`, audit `512`).

## 11. Ressources nouvelles nécessaires (classement explicite, révisé)

### 11.0 Étage 1  -  ressource requise par la construction conditionnelle
elle-même (révision 2)

Conséquence directe du §3.5 : sans l'une des trois ressources suivantes,
la construction conditionnelle n'est **pas dérivée**, même sous restriction
`Z2` déclarée. Aucune n'est présente dans la grammaire actuelle ; chacune
est une ressource **nouvelle** à déclarer comme telle.

1. **Espèces de médiateurs distinctes pour `H0` et `H1`, à charges `Z2`
   fixes** : quatre médiateurs par lien ; l'espèce `{(aa±bb)}` de charges
   `[0,0]` n'est couplée que pendant `H0`, l'espèce `{(aa-bb),(ab+ba)}` de
   charges `[0,1]` que pendant `H1` ; chaque mode garde une charge
   invariante dans le temps, et `Q = (-1)^{N_a + Σ n_{d,\rm impairs}}`
   devient exacte sur le cycle complet. Coût : doublement du nombre de
   médiateurs, de la calibration et du budget de crosstalk inter-espèces  - 
   à dériver et auditer (contrat M1).
2. **Symétrie dynamique spatio-temporelle** (combinant translation de
   demi-période et opération interne) rendant l'évolution stroboscopique
   `Z2`-covariante malgré le partage des médiateurs : simple possibilité
   logique, **aucune n'est dérivée** à ce jour.
3. **Projection ou mesure du secteur médiateur aux frontières de segment** :
   variante measurement-based du point précédent ; protocole non dérivé.

### 11.1 Étage 2  -  lever le no-go du §7 (protection contre l'algèbre
locale complète)

1. **Ancilla de jauge `Z2` avec loi de Gauss locale** (option 5 du menu)  - 
   seule voie identifiée vers une interdiction authentique de la classe
   (c) : gauger la parité de rail rend ces opérateurs non invariants de
   jauge, donc physiquement interdits. Coût : contraintes de Gauss
   multi-corps, exposées à la même hiérarchie d'ordres que le no-go
   Phase 7C ; une version Floquet de l'imposition de Gauss n'est pas
   dérivée. À ne poursuivre qu'avec un audit d'ordres préalable.
2. **Mesure de parité de rail** (option 2)  -  **mitigation candidate**, pas
   une ressource minimale suffisante démontrée. Une mesure projective
   répétée de `(-1)^{N_a}` de segment détecterait les flips de classe (c),
   sans briser U(1) ni exiger de terme Hamiltonien nouveau ; mais un
   protocole measurement-based complet reste à dériver : couplage dispersif
   ancilla–segment, fréquence de mesure contre le taux de flips (y compris
   ceux générés par l'horloge, §4.2), backaction sur le secteur `ρ`
   gapless, seuils de décision. Elle ne crée pas, à elle seule, une tresse.
3. **Réservoir/condensat de phase** (option 3)  -  restaure la structure de
   surselection de champ moyen ; route supraconductrice standard, **sort du
   mandat U(1) exacte** ; déclarée pour mémoire.
4. La restriction « `Z2` exacte du matériel » n'est pas une ressource du
   menu mais une **contrainte de symétrie non dérivée**, qui échoue
   désormais à deux niveaux (contrôle : `δφ`, crosstalk ; structurel :
   §3.5). Elle ne peut figurer dans aucune claim sans l'étage 1
   explicitement acquis.

## 12. Claim boundary

Est affirmé : (i) la jonction T minimale, ses liens de nœud et ses canaux
sont **dérivables** des ressources déclarées au sens des §3–4, comme
extension de graphe à auditer, avec parasites de nœud repoussés à l'ordre 4 ;
(ii) le doublet de fusion existerait au mieux comme doublet de bord protégé
par la symétrie `Z2` de rail, avec indistinguabilité exponentielle
restreinte à la sous-algèbre `Z2`-paire, **et seulement sous réserve de la
ressource d'étage 1 (§11.0)**, absente de la réalisation directe actuelle
(§3.5) ; (iii) le critère pré-enregistré 4 est inatteignable avec les
ressources 1–6 (no-go §7), sa levée exigeant l'étage 2 du §11.1  -  jauge
`Z2` locale, ou mesure de parité une fois dérivée en protocole
measurement-based complet ; (iv) le mode neutre gapless impose un coût
polynomial, non un blocage.

N'est pas affirmé : une phase thermodynamique, une mémoire topologique, une
tresse réalisée, une non-abélianité, un chiffre de `‖[U1,U2]‖`, une matrice
logique quelconque, l'universalité, la tolérance aux fautes, un no-go
au-delà du périmètre du §7, ni **aucune construction conditionnelle dérivée
en l'état** (blocage structurel §3.5 : pas de `Z2` microscopique commune
dans la réalisation directe à médiateurs partagés). Aucun résultat numérique
de jonction T n'existe à la date de ce document ; le contrat §10 est la
seule voie d'amendement, dans les deux directions.
