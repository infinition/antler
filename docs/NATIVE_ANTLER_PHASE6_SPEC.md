# Phase 6  -  découverte d’un Hamiltonien topologique natif ANTLER

## Principe

Iemini devient un benchmark externe gelé : il calibre les critères sans être
une cible à recopier. La Phase 6 cherche un Hamiltonien distinct, construit à
partir de la convention rung-major ANTLER, de termes locaux scalaires et de la
conservation stricte de la charge U(1).

Le Hamiltonien deux-jambes gelé reste la référence abélienne validée. Toute
extension est ajoutée dans un namespace distinct et ne modifie ni ses scripts
ni ses revendications.

## Première famille préenregistrée

La famille `native_threeleg` ajoute une troisième orbitale/jambe par rung et
conserve :

- le hopping corrélé scalaire induit par la corde Jordan--Wigner rung-major ;
- des hoppings de rung voisins ;
- des interactions de densité locales de rung et de jambe.

Elle exclut explicitement les liens matriciels SU(2), les potentiels de bord
qui fabriqueraient un qubit local, et les termes du parent Iemini.

## Portes de décision

Avant toute tresse, un candidat doit satisfaire simultanément :

1. un doublet isolé, avec séparation/gap décroissant avec la taille ;
2. indistinguabilité par un jeu croissant d’opérateurs locaux ;
3. réponse de bord localisée, sans support logique local dans le bulk ;
4. continuité spectrale sous les déformations prévues.

Le premier scan ne peut que réfuter vite un candidat localement lisible. Un
résultat positif ne donne droit qu’à un scan en taille et à un audit d’opérateurs
locaux plus complet. Les tests de commutateur, Yang--Baxter et dynamique ne
viennent qu’après.

## Mécanisme microscopique retenu pour le prochain test

La Phase 6C ne rajoute pas de terme effectif de paire à la main. Un bloc de
deux rails de médiateurs détunés à flux `pi` annule le cotunnelling à une
particule par interférence et laisse, grâce à l'interaction de médiateur, un
transfert de paire de quatrième ordre. Le contrôle sans interaction fait
disparaître ce dernier à la précision numérique.

Cette qualification est seulement locale. La prochaine famille doit tuiler
ces médiateurs explicites dans une échelle, conserver les parités de branche
et redériver la corde Jordan--Wigner globale. Le terme de paire effectif ne
doit pas être inséré dans le calcul de code : il devra ressortir de cette
extension microscopique.

## Retour du premier ladder à parité exacte

Le ladder tuilé à médiateurs de charge 2 a conservé exactement les deux
parités, mais aucun de ses 144 points préenregistrés ne satisfait le filtre de
protection. Le meilleur `L=4` manque le seuil de séparation/gap par un facteur
d'environ deux et se dégrade à `L=5`. Cette branche est donc conservée comme
un rejet de scan, sans nouvelle recherche de tresse.

La prochaine extension devra être dérivée depuis un invariant ou une limite
soluble qui fixe aussi les termes diagonaux induits, et non choisie par simple
balayage de couplages.

## Revendication interdite à ce stade

Ni la présence d’une troisième jambe, ni une quasi-dégénérescence de petit
système, ni une non-commutation écrite dans un lien ne sont une preuve de
calcul topologique. La découverte d’un Hamiltonien natif peut également se
conclure par un no-go propre ; ce résultat guiderait alors l’extension minimale
suivante au lieu de masquer l’obstacle.
