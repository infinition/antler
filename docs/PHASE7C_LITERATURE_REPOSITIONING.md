# Phase 7C  -  repositionnement bibliographique avant nouvelle exploration

## Pourquoi ce document

Les audits locaux ont rejeté la grammaire statique lien-par-lien et la famille
all-flip issue du médiateur croisé comme route directe vers un parent
stabilisateur commutant. Cette note évite de remplacer un scan aveugle par un
autre : elle distingue les mécanismes publiés qui changent réellement les
ressources physiques de ceux qui ne feraient que retuner les mêmes couplages.

## Résultats pertinents

1. **Pair-hopping par Floquet, conservant le nombre.** Defossez,
   Vanderstraeten, Peralta Gavensky et Goldman proposent une échelle à deux
   fils où une séquence de pulsations de hopping mono-particule engendre un
   Hamiltonien effectif avec pair-hopping, tout en préservant la parité aux
   temps stroboscopiques. Leur étude valide l'Hamiltonien effectif par
   dynamique exacte à petite taille, puis par bosonisation et MPS.
   Source primaire : https://arxiv.org/abs/2412.14886

2. **Modes de bord à nombre conservé ne suffisent pas au calcul topologique.**
   Thomas-Markarian, Agarwal et Martin étudient numériquement des modèles
   conservant le nombre avec interactions longues ; ils trouvent des signatures
   de modes de bord, mais pas une construction de tresse ou de code 2D.
   Source primaire : https://arxiv.org/abs/2509.00158

3. **Les parents 2D à deux corps demandent une nouvelle encodage/gadget, pas
   une simple compensation de couplages.** Brell, Flammia, Bartlett et Doherty
   construisent le toric code comme limite basse énergie de gadgets à deux
   corps avec encodage local supplémentaire. Ocko et Yoshida donnent une autre
   construction à gadget auxiliaire non perturbative. Ces travaux sont des
   précédents théoriques ; ils impliquent un coût en degrés de liberté et une
   dérivation complète, non un raccourci permettant d'appeler un terme imposé
   « émergent ».
   Sources primaires : https://arxiv.org/abs/1011.1942 et
   https://arxiv.org/abs/1107.2697

## Décision architecturale

La branche **médiateur charge-2 statique croisé** est gelée comme résultat
négatif ciblé. Le RL ne doit pas être lancé dans cette grammaire : les
contre-termes ne changent ni l'opérateur all-flip généré ni son algèbre de
recouvrement.

Deux hypothèses nouvelles, mutuellement distinctes, méritent une évaluation
préenregistrée :

| Hypothèse | Ressource nouvelle | Premier test obligatoire | Ce qu'elle ne prouve pas |
| --- | --- | --- | --- |
| **Phase 7D-F : compilateur Floquet de pair-hopping** | protocole périodique de hopping, canaux charge-2 même/opposé-rail sélectionnés dans le temps et observation stroboscopique | le bloc charge-2 positif compile les signes opposés de `ZZ` après pulse Rabi fermé, puis le compilateur quatre barreaux a une erreur Trotter quadratique ; il faut encore auditer crosstalk, fuites et tout le ladder | code 2D, tresse, non-abélianité |
| **Phase 7D-G : code-gadget 2D** | plusieurs degrés locaux/ancille par qubit logique, avec contraintes locales | algèbre parent et dérivation SW indépendante avant tout ED global | réalisation native par les médiateurs ANTLER actuels |

Dans les deux cas, le RL ne devient autorisé qu'après trois portes :

1. l'algèbre des termes effectifs voisins commute ou satisfait une relation
   explicitement visée ;
2. une dérivation indépendante fixe les paramètres et les termes parasites ;
3. un bloc exact valide la séparation bas/haut et la symétrie annoncée.

Le benchmark Floquet complet de petite taille ajoute une quatrième porte
pratique : il ne doit pas échouer au test de localisation de bord dans la
fenêtre de paramètres à optimiser. Le premier scan `L=4 -> 6` ne trouve aucun
candidat dans son domaine préenregistré. Il serait donc prématuré de lancer du
RL sur ce modèle avant une dérivation de protocole plus fidèle et une raison
physique de modifier la fenêtre explorée.

## Claim boundary

Cette note est un choix de stratégie inspiré de la littérature, pas une
implémentation de la proposition Floquet ni une dérivation d'un gadget ANTLER.
Les références 1D à nombre conservé ne démontrent pas une architecture de
calcul topologique non abélien.
