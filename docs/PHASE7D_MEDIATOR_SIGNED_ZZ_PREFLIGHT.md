# Phase 7D  -  préflight médiateur : interaction `ZZ` signée stroboscopique

## Objet

Ce test relie le compilateur Floquet abstrait de la Phase 7D à une extension
dynamique concrète de la brique ANTLER à médiateurs de charge 2. Il ne modifie
pas le Hamiltonien statique gelé : il ajoute une ressource de contrôle
temporelle, à savoir la sélection alternée de canaux de paire de même rail et
de rails opposés.

## Bloc exact à deux barreaux

Le bloc contient les modes `a0,b0,a1,b1` et deux médiateurs hard-core de
charge 2 par segment. Avec des détunings positifs identiques, les canaux

\[
(a_0a_1,b_0b_1) \quad\hbox{et}\quad (a_0b_1,b_0a_1)
\]

produisent, après une excursion Rabi complète, deux phases diagonales
logiques inverses l'une de l'autre. Les deux séquences sont donc une
réalisation pulsée de signes opposés de `ZZ`, sans introduire de détuning
négatif.

Paramètres enregistrés : `U=20`, `Delta=40`, `g=6`, donc `g/Delta=0.15`, et
temps de pulse complet `0.0776708795`.

## Résultats exacts du bloc

- fuite monomère, pour chacun des deux pulses : `1.26e-33` ;
- résidus de portes logiques même-rail et opposé-rail : `0` ;
- résidu de l'inversion des deux portes : `0` ;
- résidus de parité de rail nue dans le sous-espace logique, après pulse :
  `3.14e-16` ;
- angle de `ZZ` signé : `-0.03475747` ;
- résidu de compilation du pair-hopping après rotations de rails :
  `3.19e-16` ;
- après 45 portes compilées : probabilité de transfert `|bb>` vers `|aa>` =
  `0.99995498`.

Pendant le pulse de canaux opposés, la parité de rail nue du système complet
est habillée par le médiateur. Seule la porte logique après retour complet du
médiateur conserve la parité nue. Cette distinction stroboscopique est une
condition du protocole, pas une symétrie statique supplémentaire.

## Décision

La première porte matérielle est franchie : un bloc charge-2 à détuning
positif peut compiler exactement la primitive dynamique utilisée dans le
contrôle Floquet. C'est une **extension dynamique explicitement déclarée**,
pas une propriété déjà contenue dans la grammaire ANTLER statique gelée.

La suite obligatoire est l'audit de composition multi-liens, puis un audit de
fuite, crosstalk et corrections de pulse dans le système microscopique
étendu. Aucun gap de phase, qubit protégé, ordre 2D, tresse, non-abélianité,
universalité ou tolérance aux fautes ne découle de ce bloc.

Résultat machine : `results/phase7/mediator_signed_zz_preflight.json`.
