# ANTLER v0.7  -  Résumé français

## Résultat central

ANTLER produit désormais une véritable opération logique de phase dans le modèle simulé, pas seulement une phase mesurée sur un état isolé.

Pour le run de référence (`theta=0,3`, `R=4`, `T=20000`, Strang `dt=0,25`) :

- fuite pire cas : `9,38e-5` ;
- mélange logique parasite : `1,34e-6` ;
- phase logique impaire : `-0,2989278` ;
- pente : `-0,996426` ;
- fidélité moyenne nettoyée : `99,9999808 %`.

## Ce qui a été découvert ensuite

La navette gaussienne est très cohérente, mais sa phase dépend encore de la profondeur et de la largeur du puits. Elle constitue donc une porte géométrique calibrée, pas encore un invariant topologique strict.

La navette digitale séquentielle corrige ce point :

- sa phase ne dépend pratiquement plus de la distance parcourue ;
- elle converge vers `-theta` lorsque les particules deviennent fortement localisées ;
- l'erreur suit une loi empirique proche de `|D|^-2,5` ;
- elle reste robuste jusqu'à 20 % de désordre statique dans les tests terminés.

Dans la limite localisée idéale, le comptage exact des cordes donne :

- échange : une traversée non compensée, phase `e^{-i theta}` ;
- aller-retour : traversées opposées qui s'annulent, phase nulle.

Un lemme à deux niveaux vérifie aussi que chaque transfert adiabatique isolé porte exactement la phase du lien, indépendamment de la forme de la rampe, à environ `1e-14` près.

## Ce que cela implique

Le projet a établi une primitive de porte de phase anyonique/à hopping corrélé très haute fidélité. Il a aussi identifié la correction physique qui empêche la quantification parfaite à profondeur finie : les queues de délocalisation de la fonction d'onde.

Ce n'est pas encore un ordinateur quantique topologique complet. Le verrou suivant est la construction d'un espace de fusion non abélien et de deux générateurs de tressage qui ne commutent pas.

## Résultats non validés

Les logs de profondeur `D=8`, de composition multi-cycle et de désordre supérieur à 20 % n'ont pas produit de sortie finale complète. Ils sont conservés pour traçabilité, mais ne sont pas revendiqués comme résultats.

## Derniers compléments

- Le run digital `dt=0,125` confirme une fuite pire cas de seulement `3,02e-5` et un mélange parasite inférieur à `9e-8`.
- À `theta=0,6` et `0,9`, la famille reste diagonale et très haute fidélité, avec des pentes `-0,9772` et `-0,9826`.
- Le test d'algèbre ferme un mur conceptuel : les portes actuelles commutent et ne génèrent qu'une famille abélienne de rotations Z. Même parfaite, cette famille seule ne peut pas donner un calcul topologique universel.
- Une branche d'accélération Numba existe, mais reste exploratoire tant que sa convergence croisée avec le solveur principal n'est pas terminée.
