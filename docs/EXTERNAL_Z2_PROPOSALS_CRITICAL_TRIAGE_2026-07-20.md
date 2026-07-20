# Triage critique des propositions externes `Z2`  -  2026-07-20

Source examinée : `test.txt` (archivé à la racine). Le dossier contient une
bonne intuition architecturale, mais aussi plusieurs conclusions de fusion et
de tresse qui ne suivent pas de ses Hamiltoniens. Ce document sépare les
résultats exploitables des promotions à refuser.

## Verdict court

Le **lien de jauge neutre `Z2` sur les arêtes** est le seul pivot immédiatement
valide. Il remplace la corde impaire artificielle et contourne le no-go de
contrôle Phase 8B parce qu'il est une ressource locale nouvelle, pas un effet
Schrieffer--Wolff du médiateur charge-2.

Ce n'est pas, à lui seul, un breakthrough non abélien. Une étoile ou une
jonction en T sans plaquette ne possède ni flux `m`, ni ordre topologique, ni
espace de fusion protégé. La théorie `D(Z2)` ordinaire est abélienne : ses
canaux de fusion sont de multiplicité zéro ou un.

La cible non abélienne défendable est plus précise : des **défauts de twist
`e <-> m`** dans un code `Z2` 2D déjà qualifié. De tels défauts peuvent
reproduire les règles de fusion et de tresse d'Ising, mais demandent un mur de
domaine qui échange réellement les types `e` et `m`, donc une modification des
stabilisateurs étoile/plaquette ou une dislocation. Un simple signe classique
sur une ligne de hopping n'effectue pas cet échange.

Références primaires : [Bombin, 2010](https://arxiv.org/abs/1004.1838) montre
explicitement que les twists proviennent d'une ligne qui échange `e` et `m`,
et que les charges de `D(Z2)` seules ont une fusion abélienne. Le modèle
honeycomb de [Kitaev, 2006](https://arxiv.org/abs/cond-mat/0506438) atteint une
phase Ising non abélienne dans une construction de spins particulière, gappée
par champ magnétique ; ce n'est pas un simple hopping U(1)-conservant avec
liens `Z2`.

## Audit des propositions de `test.txt`

| Proposition | Statut | Raisonnement et action |
|---|---|---|
| Fermions sur sommets + qubit neutre `tau_e` sur chaque arête, `c_v^dag tau^z_e c_w+h.c.` | **Valide comme modèle de référence** | C'est exactement le couplage minimal d'une jauge `Z2`. Il donne une loi de Gauss exacte si `G_v=(-1)^n_v product tau^x_e`. Il introduit toutefois une ressource nouvelle, non dérivée du ladder ANTLER gelé. |
| Étoile/T à trois bras = espace de fusion 2D | **Rejeté** | Un arbre n'a aucune plaquette : aucun flux `m` et aucune holonomie de plaquette n'existent. Même sur un réseau 2D, `D(Z2)` a une fusion abélienne ; une doublet accidentel doit être rejeté s'il reste lisible localement. |
| Flux `pi` ou cocycle classique `eta_ij=-1` = twist Ising | **Rejeté sous cette forme** | Sur un arbre ou une région sans flux, le signe est une transformation de jauge/cobord. Sur une boucle, un produit `-1` est un vison de fond, pas un défaut qui permute `e` et `m`. Un twist exige une ligne de domaine qui transforme les opérateurs de corde/stabilisateurs `X` et `Z` l'un dans l'autre. |
| Réseau de Majorana « Kitaev-like » | **Benchmark externe possible, pas route ANTLER immédiate** | La construction honeycomb exacte repose sur un Hamiltonien de spins et un champ `Z2` statique émergent. Sa phase non abélienne est une ressource différente, et les Majorana/BdG ne respectent pas automatiquement le mandat U(1) à nombre fixé. Aucun Majorana, Hamiltonien BdG ou matrice de braid ne sera inséré dans ANTLER. |
| Higgs `U(1) -> Z2` en réutilisant le médiateur charge-2 | **Hypothèse d'architecture, non dérivée** | Un médiateur charge-2 avec U(1) globale n'est pas par lui-même un champ de Higgs d'une jauge compacte locale. Il faudrait définir un champ de jauge U(1), la loi de Gauss locale, un condensat Higgs et leur réduction contrôlée ; c'est une nouvelle grammaire matérielle. |
| Spins auxiliaires de sommet | **Indéterminé / secondaire** | Ils n'apportent rien sans Hamiltonien et Gauss exacts. Tout ajout doit d'abord passer les mêmes commutateurs et son comptage de ressources. |

## Résultat Phase 8C-T0 exécuté

Le premier test recommandé est maintenant terminé dans
`experiments/phase8c/run_phase8c_z2_star_gauss_audit.py`.

Géométrie : étoile `v0--v1`, `v0--v2`, `v0--v3`, quatre fermions de sommet et
trois liens neutres. Dimensions : `128` complète, `48` à `N=2`, `6` dans le
secteur `G_v=+1` à `N=2`.

| Test T0 | Valeur | Conclusion |
|---|---:|---|
| `max ||[H,G_v]||` | `0` | loi de Gauss exacte dans le modèle de référence |
| `max ||[G_v,G_w]||` | `0` | générateurs compatibles |
| `||[H,N]||` | `0` | U(1) totale exacte |
| commutateur maximal du hopping nu avec Gauss | `2` | le terme nu viole effectivement Gauss |
| norme de `P_phys H_nu P_phys` | `2.78e-17` | le hopping nu n'agit pas dans le secteur physique |
| non-scalarité du hopping habillé projeté | `1.303` | le transport habillé est une opération physique locale |
| non-scalarité de `P_phys n_v P_phys` | `0.5` | l'étoile n'est pas un code topologique |

Le test **passe T0 uniquement**. Son gap de secteur Gauss est un gap de
contrainte introduit par `-lambda sum_v G_v`; il ne doit pas être renommé gap
topologique. L'étoile a zéro plaquette et ne permet donc aucune conclusion sur
flux, déconfinement, fusion ou braid.

## Fork d'architecture résolu par T2

Le test T0/T1 contient de la matière mobile et impose
`G_v=(-1)^n_v A_v`. Il s'agit d'une théorie de jauge avec matière valide, mais
ce n'est pas automatiquement le même objet qu'un code torique qui stocke son
information sur les liens. Les deux chemins doivent désormais être distingués :

| Chemin | Secteur de départ | Objet que l'on peut tester | Ce qu'il ne faut pas conclure |
|---|---|---|---|
| **G  -  code de jauge pure** | matière vide/gelée, donc `G_v=A_v` | parent torique `-J_s sum A_s-J_p sum B_p`, GSD, distance, localité, puis murs de twists | ce code ne dérive pas encore le lien neutre ni le contrôle ANTLER |
| **M  -  matière dynamique** | `N>0`, `G_v=(-1)^n_v A_v` | transport gauge-invariant et création de charges/syndromes | une doublet ou un gap fini n'est pas un qubit topologique ; la matière peut lire ou hybrider le code |

T2 a donc commencé par le chemin **G**, avec matière strictement gelée, et a
passé le test d'indistinguabilité du code. Les charges ne pourront être
injectées que comme défauts contrôlés après l'audit du mur ; elles ne sont pas
le support initial de la mémoire topologique.

Le dépôt possède déjà deux calibrations utiles, qui ne doivent pas être
sur-vendues :

- `docs/PHASE7B_2D_REFERENCE_PREFLIGHT.md` : code torique de référence 3x3
  validé au niveau stabilisateur ; il n'est pas un Hamiltonien ANTLER natif.
- `docs/PHASE8B_2D_WALKER_STABILIZER_CLOSURE.md` et
  `docs/PHASE8B_TWIST_ENDPOINT_PREFLIGHT.md` : parents effectifs et briques de
  twist conditionnels, avec ressources de marcheur nouvelles et non dérivées.

T2 relie désormais proprement le nouveau lien neutre au chemin G. T3b doit
vérifier qu'un mur `e<->m` conserve ces portes de protection avant que la
matière contrôlée, la fusion ou les holonomies soient abordées.

## Feuille de route corrigée

1. **T0  -  fait :** qualifier l'algèbre locale du lien neutre sur l'étoile.
2. **T1  -  fait :** une plaquette carrée minimale (`256D` complet, `96D` à
   `N=2`, `12D` physique) passe Gauss/U(1), `B_p` commute avec Gauss et son
   gap statique de flux vaut `1.7491`. Une densité locale projetée de `0.5`
   ferme toute promotion de code ; une plaquette reste un préflight.
3. **T2  -  terminé, code abélien d'abord :** le tore de jauge pure `3x3` a
   rang stabilisateur `16`, GSD `4`, distance `3`, gap de syndrome `4` et
   aucune sonde locale non scalaire sous la distance. Voir
   `docs/PHASE8C_Z2_CODE_PATCH_PREFLIGHT.md`. Il établit l'ordre de référence
   avant de parler de défauts.
4. **T3a  -  calibrage terminé :** la dualité globale `e<->m` sur le tore
   échange exactement étoiles/plaquettes et cordes `Z`/`X`; le flux/signature
   `pi` et un changement de base borné échouent le test. Voir
   `docs/PHASE8C_EM_DUALITY_WALL_PREFLIGHT.md`.
5. **T3b  -  référence locale terminée :** un patch `[[7,1,3]]` possède une
   jonction non-CSS de twist-vers-bord et une déformation de corde X/Z-vers-
   mixte ; il n'a pas de défauts séparés. Voir
   `docs/PHASE8C_TRIANGLE_TWIST_REFERENCE.md`.
6. **T4  -  calibration de fusion terminée :** le code planaire CSS à quatre
   twists de bord passe aux distances `d=3,5` : doublet, gap de syndrome `2J`,
   distance `d` et aucune action logique de Pauli sous la distance. C'est une
   calibration externe de dimension d'espace de code et de protection locale,
   non une mesure de fusion physique ni des défauts mobiles.
7. **T5a  -  préflight intérieur terminé :** le formalisme de graphe donne quatre
   sommets impairs intérieurs, `k=3`, GSD `8` et un gap de syndrome sur les
   tores `L=4,6`. Ce repère externe n'est pas encore un lien-neutre ANTLER ni
   un mouvement de défaut.
8. **T5b  -  baseline fermée :** l'interpolation linéaire des checks reste gappée
   mais rend le sous-espace bas localement lisible (`0.2937`) ; elle est rejetée
   comme mouvement protégé, sans promotion de l'holonomie.
9. **T5c  -  référence de mesure terminée :** quatre mesures de checks conservent
   le rang, GSD et la localité à un qubit du petit code ; les outcomes et le
   transport logique ne sont pas encore calculés.
10. **T5d  -  véritable cible non abélienne :** dériver les frames et le transport
   de tous les Paulis logiques, puis composer deux mouvements adjacents avant
   tout commutateur ou résidu de Yang--Baxter.
11. **T6  -  pont matériel :** tenter de dériver le lien neutre et les checks
   depuis une ressource déclarée. Tout échec à cette étape est un no-go de
   réalisation ANTLER, pas une invalidation du benchmark de jauge.

## Décision sur les mots « breakthrough »

Le progrès concret est d'avoir remplacé une ambiguïté de contrôle par un
contrat algébrique testable, que T0 vient de passer. La direction twist est
une cible de recherche sérieuse et potentiellement compatible avec une tresse
non commutative, mais elle est connue dans la littérature et reste à établir
microscopiquement pour ANTLER. Le breakthrough du projet ne sera atteint que
lorsqu'un lien neutre déclaré, un code 2D protégé et un mur `e <-> m` auront
tous passé leurs audits indépendants.

## Claim boundary

Affirmé : le modèle de référence à liens neutres satisfait exactement Gauss et
U(1) sur l'étoile 48D, avec la sélection attendue entre hopping nu et habillé.

Non affirmé : toute dérivation depuis le Hamiltonien ANTLER gelé, une phase
déconfinée, ordre topologique, anyons, défauts de twist, fusion, tresse
non abélienne, universalité ou tolérance aux fautes.
