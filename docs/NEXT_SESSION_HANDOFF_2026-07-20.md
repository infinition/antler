# Handoff de reprise  -  2026-07-20

## Décision de reprise

**Les Phases 8C-T0 à T4 et T5a--T5g sont terminées à leur niveau de
référence.**
T0/T1 qualifient la grammaire de jauge et de flux ; T2 qualifie le code `Z2`
pur comme référence abélienne ; T3/T4 calibrent les twists de bord ; T5a
construit quatre twists intérieurs dans un formalisme de graphe externe.
T5f/T5g ferment deux boucles uniquement après l'ajout de mesures ponctuelles
externes et calibrent un commutateur logique non nul. La prochaine étape est
**la dérivation ou la réfutation d'un mécanisme physique pour ces mesures**,
et non une nouvelle recherche de pulses, un scan de paramètres ou une
simulation de tresse. Le bilan critique et les conditions de T3--T5 sont dans
`docs/EXTERNAL_Z2_PROPOSALS_CRITICAL_TRIAGE_2026-07-20.md`.

Chaîne de travail, dans cet ordre strict :

`lien Z2 neutre -> loi de Gauss exacte -> secteur physique -> espace de fusion -> jonction T -> holonomies`

Seuls les deux premiers maillons sont autorisés au départ. Une tresse non
abélienne, une protection topologique ou un espace de fusion ne doivent pas
être revendiqués avant leurs audits dédiés.

## État scientifique au moment du relais

- La primitive digitale de phase ANTLER est validée, mais sa famille logique
  actuelle est abélienne.
- Les réalisations dynamiques enregistrées du lien conditionnel Phase 8B sont
  fermées : fermeture Rabi, séquence `AB` rapide, écho lisse, groupe à quatre
  signes et baseline classique ne produisent pas de porte physique locale non
  triviale dans leurs boîtes de contrôle.
- L'identité contrôlée est explicite : moyenner les Hamiltoniens effectifs
  après Schrieffer-Wolff ne commute pas avec intégrer le Hamiltonien
  microscopique moyen. Le terme `XX` statique disparaît dans la limite rapide
  de la grammaire enregistrée.
- Le pont `XZ` à corde Jordan--Wigner ne fonctionne que si le médiateur
  charge-2 reçoit artificiellement un poids de corde impair. La convention
  rung-major gelée donne le poids `0` (le comptage de charge donne `2`) et
  donne exactement `c_XZ=0`. C'est un contrôle algébrique, pas une ressource
  native ANTLER.
- La ressource manquante identifiée est un **lien de jauge `Z2` neutre et
  impair**, distinct d'un médiateur moléculaire charge-2. C'est une nouvelle
  ressource déclarée, non dérivée du Hamiltonien gelé.

Lire dans cet ordre :

1. `docs/PHASE8B_JW_STRING_CONVENTION_AUDIT.md`
2. `docs/PHASE8B_FLOQUET_AVERAGING_OBSTRUCTION.md`
3. `docs/PHASE8B_Z2_GAUSS_RESOURCE_THEOREM.md`
4. `docs/PHASE8_TJUNCTION_THEORY_VERDICT.md`
5. `docs/RESULT_STATUS_MATRIX.md`

## Contrat de modèle Phase 8C-T0

Ne pas réutiliser la classe de médiateur charge-2. Créer un namespace séparé,
par exemple `experiments/phase8c/`, avec deux sous-espaces tensoriels
explicites :

- matière fermionique spinless sur les sommets `v` ;
- qubit de jauge neutre `tau_e` exclusivement sur chaque arête `e=(v,w)`.

Sur une arête, le saut autorisé est :

`h_vw = -t (c_v^dag tau^z_vw c_w + h.c.)`.

La correction essentielle au garde-fou initial est le générateur de Gauss :

`G_v = (-1)^n_v product_(e incident to v) tau^x_e`.

Le produit de liens seul `product tau^x_e` **ne commute pas** avec le saut
habillé : la parité de matière est indispensable. Cette égalité doit être
codée et testée comme identité matricielle, non inférée du spectre.

Le modèle de référence peut prendre :

`H_gauge = sum_(v,w) h_vw - K sum_v G_v - h sum_e tau^x_e`.

Toutes ces ressources et tous ces termes sont *insérés dans le modèle de
référence*. Il ne s'agit pas encore d'une dérivation microscopique ANTLER. Ne
pas ajouter de hopping nu `c_v^dag c_w+h.c.` : il viole Gauss. Toute future
dérivation doit expliquer pourquoi le terme habillé est présent et le terme
nu absent ou borné.

## Plus petit audit ED à lancer demain

Graphe : une étoile minimale, sommet central `v0` et trois feuilles
`v1,v2,v3`, avec trois arêtes. Il y a quatre sommets et trois liens.

- Hilbert complet : `2^4 * 2^3 = 128`.
- Secteur à charge fixée `N=2` : `C(4,2) * 2^3 = 48`.
- Ce bloc est assez petit pour diagonalisation exacte dense ; aucune charge
  GPU ni grande simulation n'est justifiée.
- Cette étoile **ne peut pas à elle seule certifier un espace de fusion ou une
  tresse**. Elle ne qualifie que la grammaire locale de jauge.

Script T0 exécuté :
`experiments/phase8c/run_phase8c_z2_star_gauss_audit.py`.

JSON associé : `results/phase8c/z2_star_gauss_audit.json`. Il passe tous les
commutateurs à précision machine (`128D` complet, `48D` à `N=2`, `6D`
physique) ; le hopping nu a une projection physique de `2.78e-17`, tandis que
la densité locale projetée reste non scalaire (`0.5`). T0 ne fait donc aucune
promotion de code ou de fusion.

Script T1 exécuté :
`experiments/phase8c/run_phase8c_z2_plaquette_preflight.py` est maintenant
exécuté. Le carré (`256D` complet, `96D` à `N=2`, `12D` physique) a
`max ||[H,G_v]||=0`, `max ||[B_p,G_v]||=0` et un gap statique de flux
`1.7491`, mais sa densité locale projetée vaut encore `0.5`. T1 n'est donc pas
un code ni une fusion.

Script T2 exécuté : `experiments/phase8c/run_phase8c_z2_code_patch_preflight.py`.
Son tore `3x3` de jauge pure (`18` liens) a rang stabilisateur `16`, GSD `4`,
distance `3`, gap de syndrome exact `4` et `1431` sondes de Pauli de poids 1--2
toutes scalaires ou nulles dans le code. Voir
`docs/PHASE8C_Z2_CODE_PATCH_PREFLIGHT.md` et
`results/phase8c/z2_code_patch_preflight.json`. La matière dynamique reste une
injection de défaut ultérieure, pas le support initial du code.

Le JSON doit sérialiser les conventions de base, la liste des sommets/arêtes,
la dimension complète et à charge fixée, toutes les normes ci-dessous, le
spectre par secteur de Gauss, et une `claim_boundary` explicite.

## Portes d'acceptation, dans l'ordre strict

### T0  -  algèbre avant spectre

**Terminé : PASS comme modèle de référence à liens neutres déclarés.** Les
tests ci-dessous restent les invariants obligatoires de tout prolongement.

Avant toute évolution temporelle, vérifier à précision machine :

- `||[H, G_v]|| / max(1, ||H||) < 1e-12` pour chaque sommet ;
- `||[G_v, G_w]|| < 1e-12` pour toute paire ;
- `||[H, N]|| / max(1, ||H||) < 1e-12` ;
- `G_v^2 = I` et hermiticité de `H` et de chaque `G_v` ;
- le saut nu a un commutateur Gauss non nul, tandis que le saut habillé le
  commute ;
- le projecteur physique `P_phys` sur le secteur choisi vérifie
  `P_phys G_v P_phys = P_phys`.

Un échec T0 clôt le candidat : aucune optimisation, dynamique ou tresse ne
doit être lancée pour cette définition.

### T1  -  secteur et spectre

Après T0 seulement : rapporter le splitting interne, le gap de charge, le gap
de syndrome/Gauss et le plus petit gap neutre. Les secteurs et la convention
de signe de `-K G_v` doivent être imprimés. Un gap isolé seul n'est pas une
preuve de code.

### T2  -  protection locale

**Terminé : PASS pour la référence de jauge pure.** Le tore `3x3` a un code
de dimension quatre, de distance trois, et aucune action logique locale de
poids inférieur à trois. C'est une porte de protection locale, pas une tresse
ni une dérivation ANTLER.

### T3  -  mur `e<->m`, conditionnel

**T3a terminé : calibrage seulement.** La dualité globale échange exactement
étoiles/plaquettes et une corde `Z`/`X`; flux, signe statique et Hadamard borné
échouent le test de support. Voir `docs/PHASE8C_EM_DUALITY_WALL_PREFLIGHT.md`.

**T3b-référence terminé :** le patch triangulaire non-CSS `[[7,1,3]]` possède
un check central `YZXIXII`, un gap `2J`, et des déformations explicites de
cordes `Z`/`X` vers des continuations mixtes. Il se termine au bord ; voir
`docs/PHASE8C_TRIANGLE_TWIST_REFERENCE.md`.

**T4 terminé comme calibration :** le patch de surface-code planaire à quatre
coins de type twist a un doublet de distance `3`, puis `5` lorsque la séparation
des coins croît de `2` à `4`; 1 089 525 sondes sous distance au point `d=5` ne
réalisent aucune action logique. Voir
`docs/PHASE8C_CORNER_TWIST_FUSION_CALIBRATION.md`.

**T5a terminé comme référence de graphe :** deux arêtes retirées d'un tore
forment quatre sommets impairs intérieurs. Les tailles `L=4,6` ont respectivement
`k=3`, GSD `8`, distance `3,4` et aucune sonde locale non scalaire sous cette
distance. Voir `docs/PHASE8C_INTERIOR_TWIST_GRAPH_PREFLIGHT.md`. Ce sont des
qubits de sommet externes, et non les liens de jauge Phase 8C.

**T5b terminé comme contrôle négatif :** l'interpolation hamiltonienne linéaire
reste gappée mais donne une lecture locale projetée de `0.2937` au milieu du
chemin. Elle est rejetée comme déplacement protégé ; voir
`docs/PHASE8C_INTERIOR_TWIST_LINEAR_DEFORMATION_AUDIT.md`.

**T5c terminé comme référence à mesures :** quatre mesures de checks finaux à
issue `+1` conservent rang `6`, GSD `8` et les 27 projections locales
scalaire/nulles à chaque étape. Voir
`docs/PHASE8C_INTERIOR_TWIST_MEASUREMENT_DEFORMATION.md`. L'appareil de mesure
et les frames issus des autres outcomes ne sont pas dérivés.

**T5d terminé comme référence :** les six Paulis logiques sont transportés et
les `16` frames d'outcome sont enregistrés. Voir
`docs/PHASE8C_INTERIOR_TWIST_LOGICAL_TRANSPORT.md`.

**T5e :** la grammaire des seuls checks cibles reste ouverte (la fermeture est
rejetée). **T5f/T5g :** une puis deux boucles ferment conditionnellement avec
des mesures locales externes et donnent un commutateur logique GF(2) non nul.
Voir `docs/PHASE8C_SECOND_AUXILIARY_HOLONOMY_COMMUTATOR.md`. Le prochain gate
est de dériver une mesure/commande symétrique et physique ; une matrice de
braid ne doit jamais être insérée à la main.

**T5h :** le contrôle local exhaustif jusqu'à six mesures à un sommet ne
ferme aucune boucle protégée. Il exclut un artefact local court, mais ne rend
pas les checks de graphe physiques ; voir
`docs/PHASE8C_LOCAL_MEASUREMENT_ONLY_CONTROL.md`.

**T5i :** les 19 checks requis ont une compilation exacte ancilla/CNOT, mais
le jeu de portes, la connectivité et la lecture sont des ressources externes.
La prochaine tâche est une réalisation locale et bruitée de ce contrat, ou un
no-go ciblé depuis le Hamiltonien ANTLER ; voir
`docs/PHASE8C_ANCILLA_CHECK_CIRCUIT_AUDIT.md`.

## Garde-fous anti-mirage

- Aucun terme de BdG/Majorana, pairing non conservant la charge, ou matrice de
  braid imposée.
- Aucun résultat de Hamiltonien effectif promu comme protocole sans séquence
  microscopique, convergence temporelle, signal relatif et leakage.
- Pas de RL/deep RL dans la grammaire Phase 8B fermée. Un optimiseur ne devient
  pertinent qu'après T0 pour une grammaire physique nouvelle et déclarée ; il
  ne peut pas créer le lien neutre manquant.
- Toute nouvelle primitive (qubit de lien neutre, terme à trois corps habillé,
  mesure de Gauss, ancilla) doit être déclarée dans le JSON, la matrice de
  statut et le résumé avant de l'utiliser.
- Mettre à jour `scripts/validate_saved_results.py`, le manifeste et les
  checksums à chaque résultat promu.

## Commandes de reprise et de clôture

Depuis la racine du dépôt :

```powershell
python scripts\validate_saved_results.py
python experiments\phase8c\run_phase8c_z2_star_gauss_audit.py
python experiments\phase8c\run_phase8c_z2_plaquette_preflight.py
python experiments\phase8c\run_phase8c_z2_code_patch_preflight.py
python experiments\phase8c\run_phase8c_auxiliary_closure_holonomy.py
python experiments\phase8c\run_phase8c_second_auxiliary_holonomy_commutator.py
python experiments\phase8c\run_phase8c_local_measurement_only_control.py
python experiments\phase8c\run_phase8c_ancilla_check_circuit_audit.py
python scripts\validate_saved_results.py
```

Après toute modification documentaire ou résultat promu, régénérer
`FILE_MANIFEST.txt`, `MANIFEST.sha256` et `SHA256SUMS.txt` avec le script de
scellage déjà utilisé dans cette session, puis vérifier chaque hash. Les
exclusions immuables sont `__pycache__`, `*.pyc`, les trois fichiers de
manifest eux-mêmes et `results/raw/phase4_7/`.

## Claim boundary de départ

Affirmé : une référence de jauge `Z2` avec liens neutres et saut habillé a été
auditée exactement sur une étoile 48D et une plaquette 96D à charge fixée ;
elle passe Gauss, U(1) et le test de flux. Son branchement de jauge pure passe
ensuite le patch de code torique `3x3` (GSD 4, distance 3 et localité sous la
distance). T3a calibre la dualité globale et rejette les faux murs ; T3b ajoute
une jonction mixte locale vers un bord. T4 calibre un doublet de code à quatre
coins de type twist. T5a fournit quatre twists intérieurs de graphe, mais pas
leur mouvement. T5b rejette l'interpolation linéaire, gappée mais localement
lisible. T5c fournit une déformation de référence à mesures, mais sans frames
ni transport logique. T5d a ensuite transporté la base logique et T5f/T5g
ont produit deux holonomies de référence conditionnelles avec un commutateur
non nul. Elles reposent sur des mesures locales externes et n'établissent pas
une tresse.

Non affirmé : que ce lien est dérivé d'ANTLER, qu'un espace de fusion existe,
qu'il y a protection topologique, ou qu'une tresse non abélienne est proche.
Le premier résultat utile peut parfaitement être un no-go ciblé ; il évitera
de consacrer du calcul à une jonction T ou à un faux twist sans permutation
`e<->m` cohérente.
