# Phase 5  -  statut vérifiable vers une opération non abélienne

Ce document sépare les résultats du ladder ANTLER gelé, les tests
d’implémentation d’extensions prescrites et les benchmarks d’architecture. Il
ne promeut aucune porte en calcul topologique universel sans sous-espace de
fusion protégé et tresses non commutatives démontrées.

## Résultats négatifs à conserver

### Ladder scalaire N=2

La famille digitale actuelle est diagonale de type `Z` et abélienne. Elle
fournit une excellente primitive de phase, mais pas deux générateurs de braid
non commutatifs.

### Médiateur N=3 statique : deux no-go indépendants

1. **No-go algébrique.** Pour les deux échanges locaux testés,
   `||[U_A,U_B]|| = 5.80e-8`, sous le seuil explicite `1e-5`. Les matrices
   sont quasi diagonales. Le résidu brut de Yang--Baxter (`5.79e-4`) est donc
   enregistré mais déclaré non interprétable : `braid_relation_residual=null`.
2. **No-go de contrôle adiabatique du protocole testé.** Le doublet présente
   un quasi-croisement interne (`split ≈ 1.02e-6`) alors que son gap vers le
   complément vaut seulement `≈1.25e-3`. Ce ne sont pas le même gap : le
   premier impose un transport de sous-espace, le second borne l’adiabaticité.
   Les liens de transport Kato ont une valeur singulière minimale `0.664`,
   sous le garde-fou `0.95`; le chemin n’est pas qualifié pour une porte à
   temps fini. Ce constat rejette le médiateur **épinglé** étudié, sans être
   un théorème d’impossibilité pour tout médiateur mobile ou tout graphe élargi.

Sources :

- `results/phase5/n3_local_exchange_holonomy.json` ;
- `results/phase5/n3_mediator_gap_optimization.json` ;
- `results/phase5/n3_pinned_mediator_no_go.json`.

## Ce que les extensions SU(2) ont réellement vérifié

Les liens `B1=exp(-i pi Z/4)` et `B2=exp(-i pi X/4)` ont été **écrits** sur
les liens synthétiques. Les audits à une particule puis à `N=2` avec corde
Jordan--Wigner montrent que le transport digital restitue fidèlement ces lignes
de Wilson SU(2) imposées, y compris la non-commutation prescrite. C’est un
test d’implémentation et de convention de corde utile ; ce n’est ni une
non-abélianité émergente, ni une preuve de protection topologique.

Le préflight complet L=14 renforce cette limite : malgré un gap spectral et
un doublet quasi dégénéré, une mesure locale de spin distingue le code avec
une norme `≈1.396`. Le spinor L=14 est donc explicitement **non protégé**.
Il ne faut pas lui ajouter des audits de leakage ou de bruit pour en déduire
une mémoire topologique ; il faut d’abord construire l’encodage défaut/fusion.

Sources :

- `results/phase5/synthetic_ed_link_transfer.json` ;
- `results/phase5/synthetic_n2_string_transfer.json` ;
- `results/phase5/synthetic_l14_code_preflight.json`.

## Premier pont local à conservation du nombre : préflight négatif du modèle minimal

Une famille additive à deux fils de fermions ordinaires a été ajoutée sans
modifier le Hamiltonien ANTLER gelé : hopping de jambe, transfert de paires
voisines entre fils et interaction densité-densité. Elle conserve le nombre
total et les parités de chaque fil. C'est un test **inspiré** de la voie
number-conserving de la littérature, pas une réimplémentation du modèle
exactement soluble d'Iemini et al.

Le scan honnête `L=4`, suivi des quatre meilleurs points à `L=6`, ne produit
aucun candidat de doublet protégé selon les critères prédéfinis. Le point le
plus proche (`w_pair=1`, `v_nn=-1`, `L=6`) a bien un split/gap `7.10e-5` et
une différence de densité locale maximale `8.35e-4`, mais son transfert local
reste massif dans le bulk : `max_bulk/mean_edge = 0.705`, très loin du seuil
exploratoire `0.1`. Le doublet n'est donc pas une paire de défauts localisés ;
on ne lui construit ni tresse ni porte.

Ce résultat **réfute seulement ce pont minimal et ce scan fini**. Il ne réfute
ni le modèle exactement soluble publié, ni une extension future avec
défauts/fusion explicitement construits.

Source : `results/phase5/pair_hopping_two_wire_preflight.json` ; code :
`antler/number_conserving_pairwire.py` et
`experiments/phase5/run_phase5_pair_hopping_two_wire_preflight.py`.

## Repère positif : parent exactement soluble d'Iemini reproduit indépendamment

Le Hamiltonien publié d'Iemini *et al.* (Eq. 3, ligne exactement soluble
`lambda=1`) a ensuite été implémenté séparément et soumis au même audit. C'est
le premier résultat positif de cette piste : les secteurs de parité relative
ont un doublet de sol exactement nul à l'erreur numérique (`<=1.0e-14`) pour
`L=4,6,8`, tandis que l'opérateur local `a_j^dag b_j` se localise aux bords.
Le ratio `max_bulk/mean_edge` décroît de `0.200` à `0.111` puis `0.0769` entre
`L=4,6,8` ; le dernier point est un calcul creux indépendant.

Ce résultat valide l'implémentation de l'audit et fournit un vrai **benchmark
à conservation du nombre** avec doublet de bord localisé. Il ne franchit pas
le mur ANTLER : ce Hamiltonien contient les interactions locales précises de
la référence, qui ne sont pas dérivées du hopping corrélé gelé. De plus, le
gap collectif à parité fixée diminue avec la taille, comme attendu pour cette
ligne ; il ne faut pas le présenter comme un bulk many-body gappé. Aucune
tresse dynamique ni universalité n'est encore établie ici.

Source : `results/phase5/iemini_exact_preflight.json` ; code :
`experiments/phase5/run_phase5_iemini_exact_preflight.py`.

## Audit de l'algèbre de tresse de la référence U(1)

Les générateurs de braid publiés ont été implémentés à support fini et projetés
sur le doublet exact. Le garde-fou de commutateur est franchi sans ambiguïté :
à `L=10`, support `j=4`, `||[R_aR,aL,R_aR,bR]||=1.281`. Le résidu
Yang--Baxter brut est donc interprétable et descend à `6.36e-2` (contre
`1.11e-1` à `L=8,j=3`), tandis que le défaut d'unitarité projeté descend à
`6.68e-2`.

Ce sont des tendances de troncature encourageantes, **pas encore une relation
de tresse établie à taille finie** : la fuite d'amplitude reste `0.223` au
meilleur point. L'audit démontre que la soufflerie détecte simultanément
localité de support, non-commutation et convergence vers Yang--Baxter sur une
référence U(1) externe. Il ne transforme pas ANTLER en ce modèle et ne prouve
pas une tresse dynamique.

Voir `docs/IEMINI_BRAID_AUDIT.md` et
`results/phase5/iemini_braid_scaling.json`.

## Benchmarks d’architecture, non dérivés d’ANTLER

Le réseau de Majorana effectif vérifie une tresse à temps fini et la jonction
Kitaev BdG à phase contrôlée possède un espace de parité fixe à deux
dimensions, quatre modes nuls, un gap et une indistinguabilité locale au point
doux. Ces chiffres servent de spécifications pour une extension future. Ils
nécessitent un appariement p-wave et un lien de jonction de phase contrôlée,
termes absents du Hamiltonien ANTLER gelé à nombre conservé.

Ils ne constituent donc pas une réalisation microscopique ANTLER, ni une
tresse BdG complète dérivée du ladder.

Sources :

- `results/phase5/tjunction_braid_dynamics.json` ;
- `results/phase5/tjunction_kitaev_preflight.json` ;
- `results/phase5/tjunction_kitaev_static_noise.json` ;
- `results/phase5/tjunction_kitaev_locality_audit.json`.

## Décision avant toute nouvelle dérivation

1. Fermer la Phase 4.7 sans promouvoir les points encore incomplets.
2. Utiliser le positionnement dans `docs/PHASE5_LITERATURE_POSITIONING.md`.
3. Choisir une extension locale précise à conservation du nombre ou une
   échelle trois-jambes Mott/chirale.
4. Établir ou réfuter un sous-espace de défaut/fusion localement
   indistinguable **avant** tout audit de porte, de leakage ou de bruit.
5. Ne tester la relation de braid que conjointement à une norme de
   commutateur non nulle et à un chemin spectralement continu.
