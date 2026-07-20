# Phase 7C  -  environnement d'optimisation microscopique contraint

## Décision

Une boucle d'optimisation externe est faisable sur des blocs locaux, pas par ED microscopique globale. Un spectre seul n'est pas une ground truth suffisante : la récompense compare l'opérateur effectif, ses termes parasites et ses symétries au parent stabilisateur 2D.

`antler/phase7_microscopic_optimizer.py` instancie `MicroscopicCandidate2D` sur quatre barreaux. Sa grammaire contient seulement contrainte Mott, médiateurs hard-core de charge 2, conversions de paires, hoppings de même rail, densités `ZZ` et biais de rail. Aucun stabilisateur à quatre barreaux ne peut être inséré directement. La charge pondérée est exacte et chaque médiateur doit avoir une signature de parité de branche définie.

## Budget et ground truth

Le tore `3 x 3` a 36 modes de faible énergie. À charge 18, il contient déjà `C(36,18)=9 075 135 300` états sans médiateur, et 66 345 304 596 avec quatre médiateurs de charge 2. Une ED globale est donc exclue de la boucle.

Le bloc à quatre barreaux a 98, 127, 157 ou 188 états avec un à quatre médiateurs. L'environnement diagonalise exactement ce bloc, aligne les 16 états bas sur la base monomère par décomposition polaire, puis développe le Hamiltonien effectif dans les 256 Paulis. Les cibles locales sont `-J_s XXXX/2` et `-J_p ZZZZ/2`.

## Interface compacte

`compact_optimizer_observation` retourne le vecteur `S_t` suivant : recouvrement monomère minimal, gap haut/bas, gap normalisé, violation Mott, résidu de parité, coefficient cible, ratio à `-J/2`, alignement cible, deux résidus opérateurs et norme parasite. Il ajoute seulement l'histogramme des poids Pauli, les quatre parasites dominants, la connectivité et les échecs durs.

La récompense favorise le coefficient cible de bon signe et son alignement, mais pénalise parasites, perte de capture et gap insuffisant. Symétrie violée, gap nul ou capture inférieure à 0,90 imposent un rejet. L'interface d'optimisation retourne un seul payload JSON de candidat, jamais du code arbitraire ; `candidate_from_payload` le valide avant toute diagonalisation.

`OptimizationMonitor` journalise chaque audit en JSONL et CSV, affiche une ligne console et régénère un tableau Matplotlib avec loss/récompense, gap/capture, sélectivité opérateur et trajectoires de paramètres. Le script `experiments/phase7/run_phase7c_initial_observation.py` écrit le premier état complet dans `results/phase7/optimizer_s0_xxxx.json` et le tableau dans `results/phase7/optimizer_monitor/dashboard.png`.

## Contrôle négatif

`seeded_perturbative_candidate()` contient deux canaux de paires disjoints, avec `g/Delta=0.05`. Il a un manifold propre (recouvrement minimal `0.9998001`, gap `39.9751`, parité exacte) mais ses coefficients `XXXX` et `ZZZZ` sont nuls à environ `10^-15`. Il ne crée que `XX`, `YY` et `ZZ` de poids deux. La récompense ne confond donc pas un gap local avec un stabilisateur.

Le résultat est dans `results/phase7/microscopic_reward_baseline.json`.

Le premier anneau mixte soumis comme contrôle est archivé dans `results/phase7/mixed_species_ring_audit.json`. Son payload brut déclarait à tort `[0,0]` pour deux médiateurs `a-b`, alors que leurs termes imposent la signature `[1,1]`. Le parseur refuse désormais cette incohérence. Même après correction des seules métadonnées, le modèle conserve seulement des parités habillées et produit `XXXX` à l'ordre `10^-10`, non à l'ordre cible `-0.5`.

## Curriculum

1. Structure : un ou deux canaux pairs, sans hopping ; passer symétrie, capture et gap.
2. Connectivité : au plus quatre canaux et hoppings même rail ; connecter les quatre barreaux par graphe ou hypergraphe de médiateur et obtenir un coefficient cible de bon signe.
3. Sélectivité : optimiser `XXXX` et `ZZZZ` avec alignement `>=0.80` et norme parasite relative `<=0.25`.
4. Promotion : dérivation SW indépendante, contrôle de clusters chevauchants, puis seulement une tuile 2D scalable.

Un score local ne prouve ni Hamiltonien tuilable, ni phase topologique, ni tresse, ni non-abélianité, ni universalité, ni tolérance aux fautes.
