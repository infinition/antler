# Phase 8C  -  cartographie récente et décision d'architecture (2023--2026)

## Portée

Cette carte ne prétend pas couvrir toute la littérature. Elle couvre les
résultats primaires récents qui touchent directement les contraintes ANTLER :
U(1) exacte, lien de jauge local, code 2D, défauts/tresses et universalité.
Elle met à jour la décision de travail, sans requalifier les benchmarks
externes comme résultats natifs.

## Résultat de la recherche

### 1. La route `Z2` + twists reste la meilleure route courte vers une tresse

La construction de twists n'est pas une simple phase `pi` de hopping : elle
requiert une ligne de domaine qui échange réellement `e` et `m` dans un code
2D. C'est exactement le sens des défauts de Bombin. Des travaux récents
continuent de chercher des réalisations microscopiques de permutations
d'anyons dans des topological orders doubles, confirmant que le problème
physique est le design de la dislocation/check, pas un signe classique sur une
arête.

- H. Bombin, *Topological Order with a Twist: Ising Anyons from an Abelian
  Model* (2010), https://arxiv.org/abs/1004.1838
- G. M. Yoshitome, *Lattice Realization of Twist Defects in a
  Z2 x Z2 Topological Order* (2026, prépublication),
  https://arxiv.org/abs/2605.02039

**Décision ANTLER.** Conserver la route Phase 8C : lien neutre `Z2` -> code
de jauge pure -> mur `e<->m` -> deux paires de twists -> holonomies. Un twist
Ising est une étape non abélienne crédible, mais ses braids sont de Clifford
et ne suffisent pas à eux seuls à l'universalité.

### 2. La voie Floquet à nombre conservé est déjà un benchmark, pas une sortie

Defossez, Vanderstraeten, Peralta Gavensky et Goldman ont publié une échelle à
nombre conservé où une séquence Floquet vise
`H_eff=alpha H0+(1-alpha) P^dag H0 P` et un pair-hopping.

- A. Defossez *et al.*, *Dynamic Realization of Majorana Zero Modes in a
  Particle-Conserving Ladder* (2025), https://arxiv.org/abs/2412.14886

ANTLER a déjà reproduit le benchmark externe et sa décomposition directe
`H0/H1`. Le no-go spécifique de la réalisation enregistrée reste : les slots
de médiateur réutilisés n'ont pas de parité microscopique commune entre `H0`
et `H1`; les tests de pulse ont en outre séparé fermeture de population et
signal logique. Cette publication n'est donc pas une échappatoire à rouvrir
sans une nouvelle ressource qui réalise `P` ou des espèces séparées avec une
symétrie commune exacte. Voir `PHASE8_CANONICAL_MPS_AUDIT.md` et
`PHASE8_MICROSCOPIC_PARITY_AUDIT.md`.

### 3. La frontière universelle est un ordre non abélien et la fusion

Trois jalons matériels récents sont instructifs :

| Résultat | Ce qu'il démontre | Limite pertinente pour ANTLER |
|---|---|---|
| D4 sur processeur ions piégés (2024) | préparation adaptative et contrôle d'un ordre non abélien | circuit numérique, pas Hamiltonien ANTLER passif |
| Fibonacci sur processeur supraconducteur (2024) | fusion et braiding universel simulés | string-net/circuits imposés, pas lien microscopique U(1) dérivé |
| D(S3) sur processeur ions piégés (2026) | ensemble universel par **braiding + fusion**, état fondamental 54 qubits | réclame liens non abéliens à six états et Gauss non abélien, hors ressources actuelles |

Sources primaires :

- M. Iqbal *et al.*, *Non-Abelian Topological Order and Anyons on a
  Trapped-Ion Processor* (2024),
  https://www.nature.com/articles/s41586-023-06934-4
- S. Xu *et al.*, *Non-Abelian Braiding of Fibonacci Anyons with a
  Superconducting Processor* (2024),
  https://www.nature.com/articles/s41567-024-02529-6
- C. F. B. Lo *et al.*, *Universal Gates from Braiding and Fusing Anyons on
  Quantum Hardware* (2026),
  https://www.nature.com/articles/s41586-026-10709-y

**Décision ANTLER.** Séparer explicitement deux ambitions :

- **objectif proche :** tresse non commutative de twists `Z2`, sans promesse
  d'universalité par braiding seul ;
- **objectif universel à long terme :** soit fusion/ressource magique ajoutée
  à la route Ising, soit passage à un quantum double non abélien comme
  `D(S3)`. Ce dernier nécessite au minimum un lien neutre à six états (ou un
  encodage local équivalent), des contraintes de Gauss non abéliennes et une
  nouvelle dérivation microscopique. Il ne doit pas être présenté comme une
  simple extension du médiateur charge-2.

### 4. Les protocoles actifs par gauging sont une branche alternative

Des travaux de 2025 dérivent des rubans d'anyons et de défauts en utilisant des
circuits unitaires séquentiels ou, selon le cas, des mesures adaptatives.

- A. Lyons *et al.*, *Protocols for Creating Anyons and Defects via Gauging*
  (2025), https://arxiv.org/abs/2411.04181

Cette approche est très utile comme spécification d'opérateurs de ruban et de
tests de fusion. Elle est **hors mandat passif actuel** : elle autorise une
préparation/correction active par circuits ou mesures. Elle deviendra une
branche ANTLER distincte seulement si ces ressources sont déclarées et
autorisées, jamais comme une dérivation implicite du Hamiltonien gelé.

## Feuille de route falsifiable

### Track A  -  tresse non abélienne de twists (actif)

1. **T2 : code de jauge pure.** Sur un patch 2D, geler la matière ; établir
   stabilisateurs étoile/plaquette, gap de syndrome, distance et
   indistinguabilité locale exhaustive.
2. **T3 : mur de domaine.** Dériver un mur `e<->m` au niveau des checks. Une
   phase `pi`, un flux `B_p=-1` ou une ligne classique `eta` échoue cette
   porte.
3. **T4 : fusion.** Deux paires de twists séparées, dimension de fusion
   mesurée, opérateurs locaux projetés scalaires, gap et leakage enregistrés.
4. **T5 : braid.** Déformations adiabatiques microscopiquement définies ;
   mesurer `||[U1,U2]||`, puis Yang--Baxter seulement si ce commutateur est
   significativement non nul.
5. **T6 : pont matériel.** Dériver les liens neutres et les checks du modèle
   ANTLER étendu ; faire bruit, crosstalk, taille et stabilité de gap.

### Track B  -  universalité topologique (futur, séparé)

Après T5 seulement, évaluer quantitativement :

- une primitive de fusion/mesure ou une injection de magic state pour twists
  Ising ; ou
- une architecture `D(S3)` : lien six-états, opérateurs de multiplication de
  groupe, contrainte Gauss non abélienne, et fusion non abélienne.

Le choix B est une décision matérielle majeure. Aucun scan RL ne peut remplacer
le degré de liberté de lien requis ; l'optimisation ne devient utile qu'après
la déclaration d'une grammaire microscopique légale.

## Claim boundary

Affirmé : la littérature récente confirme la hiérarchie ANTLER `code ->
twist/fusion -> braid`, et distingue clairement la tresse non commutative de
l'universalité. Les résultats T0/T1 du dépôt qualifient seulement le premier
lien de cette chaîne.

Non affirmé : que ANTLER réalise un twist, un anyon non abélien, une tresse,
une fusion, `D(S3)`, une préparation active par gauging, ou une universalité.
