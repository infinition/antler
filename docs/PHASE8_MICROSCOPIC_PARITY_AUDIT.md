# Phase 8  -  audit de parité microscopique des canaux directs

## Verdict

La réalisation directe enregistrée, qui réutilise les deux mêmes médiateurs
de charge 2 par lien pendant `H0` et `H1`, ne possède pas de parité de rail
microscopique commune. C'est un obstacle d'architecture démontré sur le bloc
explicite `L=3, N=2`; ce n'est pas un no-go pour une réalisation utilisant des
espèces de médiateurs distinctes.

La correction importante est que, au cycle à fermeture Rabi enregistré, la
violation de parité effectivement projetée ne suit **pas** la loi
`O((g/Delta)^2)` auparavant avancée. Le fit mesuré est proche de l'ordre six.
L'absence de symétrie exacte subsiste, mais l'ordre de l'effet dynamique doit
être rapporté avec sa séquence et ne doit pas être transformé en no-go plus
fort sans un nouvel audit.

## 1. Symétrie commune : énumération exhaustive

On a énuméré les 16 opérateurs candidats

\[
Q_x=(-1)^{N_a + \sum_m x_m n_{d,m}}, \qquad x_m\in\{0,1\},
\]

sur le bloc direct `L=3,N=2`, avec les quatre slots
`[link0_slot0, link0_slot1, link1_slot0, link1_slot1]` :

- `H0`: canaux `(aa-bb)/sqrt(2), (aa+bb)/sqrt(2)`;
- `H1`: canaux `(aa-bb)/sqrt(2), (ab+ba)/sqrt(2)`.

Résultat exact :

| Condition | Assignation |
|---|---|
| commute avec `H0` | `[0, 0, 0, 0]` seulement |
| commute avec `H1` | `[0, 1, 0, 1]` seulement |
| commute avec les deux | aucune |
| meilleur compromis | résidu Frobenius normalisé `0.0705338` pour chaque segment |

Donc une charge `Z2` fixe ne peut pas être assignée aux slots physiques
réutilisés. Cette conclusion porte sur cette grammaire directe et ne décide ni
une jonction T, ni un espace de fusion, ni une tresse.

## 2. Test direct de la frontière temporelle

Le cycle testé est

\[
U=e^{-i(1-\alpha)T H_1}e^{-i\alpha T H_0},\quad
\alpha=1/2,\quad T=\frac{4\pi}{\sqrt{\Delta^2+4g^2}},
\]

à `U0=-g^2/Delta=-1.5`, pour `g/Delta` de `0.20` à `0.0125`. Il mesure les
blocs impair-vers-pair de l'évolution projetée et le commutateur de parité de
son unitaire polaire. Dans le modèle à médiateurs partagés, la norme
impair-vers-pair passe de `1.7851e-3` à `1.3374e-10`; la pente log-log est
`5.9312`. Le commutateur polaire donne la même pente `5.9314`. La frontière
inverse est également incluse en répétant le cycle : dans la fenêtre SW
profonde `g/Delta<=0.075`, les pentes du commutateur polaire pour `1,2,4,8`
cycles sont respectivement `5.983`, `5.984`, `6.018`, `6.319`. Il n'y a donc
aucun signal d'un terme d'ordre deux dans ce contrôle temporel.

La loi `O((g/Delta)^2)` est donc **réfutée pour ce protocole Rabi-fermé
idéal**, pas établie. Ce test ne borne pas les protocoles hors fermeture, les
rampes physiques ou les jonctions.

## 3. Contrôle à espèces séparées

Le contrôle utilise quatre médiateurs par lien : les deux canaux `H0` ont
charges `[0,0]`; les canaux `H1` ont charges `[0,1]`. Ainsi

\[
Q=(-1)^{N_a+\sum_{\text{slots }(ab+ba)}n_d}
\]

commute exactement avec les deux Hamiltoniens microscopiques. Son résidu de
commutateur est inférieur à `1.7e-17` dans tout le scan; les amplitudes de
flip projetées restent à l'arrondi (`<1.4e-17`). Le contrôle isole bien la
réutilisation des slots comme origine du défaut de symétrie, sans le confondre
avec la fuite.

## Claim boundary

Ce sont des calculs ED idéaux `L=3,N=2` à canaux imposés. Ils n'établissent
pas une symétrie de matériau, une jonction T, une phase thermodynamique, une
protection contre l'algèbre locale complète, une fusion, une tresse non
abélienne, l'universalité ou la tolérance aux fautes.

## Fichiers reproductibles

- `experiments/phase7/run_phase8_direct_channel_common_parity_audit.py`
- `results/phase7/direct_channel_common_parity_audit.json`
- `experiments/phase7/run_phase8_direct_channel_boundary_parity_scaling.py`
- `results/phase7/direct_channel_boundary_parity_scaling.json`
