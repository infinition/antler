# Phase 7E  -  portes de validation pour une construction de code 2D native

## But réel

La Phase 7D fournit seulement une primitive locale de contrôle dynamique dans
un bloc à quatre barreaux. Elle ne définit ni un code topologique, ni une
phase, ni un braid. La Phase 7E ne pourra être promue que si elle relie cette
primitive à un Hamiltonien effectif 2D ayant un sous-espace de code vérifiable.

Le parent de référence est un modèle stabilisateur 2D à générateurs étoile et
plaquette. Ce parent sert de cible de diagnostic ; l'insérer directement dans
le Hamiltonien ANTLER ne compterait pas comme une réalisation microscopique.

## Portes à franchir, dans cet ordre

| porte | résultat exigé | rejet immédiat |
| --- | --- | --- |
| E0  -  grammaire de contrôle | les termes microscopiques autorisés, leurs symétries et les phases commandables sont écrits explicitement | un générateur de stabilisateur est ajouté à la main ou viole la conservation de charge déclarée |
| E1  -  dérivation indépendante | une réduction SW ou Magnus/Floquet contrôlée donne les coefficients connectés étoile/plaquette et les termes parasites | seul un fit spectral, une isospectralité, ou une projection non contrôlée est disponible |
| E2  -  sélectivité | le coefficient cible dépasse les termes non désirés selon un ratio préenregistré, avec une petite parameter `g/Delta` qui reste contrôlée | la cible apparaît à un ordre perturbatif plus élevé que des parasites dominants, comme dans le no-go statique Phase 7C |
| E3  -  code 2D | sur patchs croissants : gap neutre, séparation du code, indistinguabilité locale et boucles de Wilson sont mesurés dans le modèle microscopique ou son effectif contrôlé | un doublet est seulement global, lisible localement, ou son gap se ferme avec la taille |
| E4  -  opérations topologiques | des défauts ou une jonction rendent deux opérations réellement non commutatives ; le résidu Yang--Baxter est rapporté avec une norme de commutateur non nulle | une relation de braid est satisfaisfaite trivialement parce que les opérations commutent ou coïncident |

## Rôle précis du contrôle Peierls 7D

L'écho Peierls peut être une ressource pour compiler une séquence Floquet ou
supprimer un terme parasite. Il ne prouve pas à lui seul qu'un terme
stabilisateur à quatre corps est présent, et ne contourne pas le no-go Phase
7C pour sa génération statique perturbative. Une construction dynamique devra
donc exhiber son Hamiltonien effectif connecté, son régime de haute fréquence
et ses corrections de Magnus.

## Calculs à lancer seulement après E0

1. Construire le patch de référence et son suite d'audits indépendants.
2. Écrire le drive microscopique ANTLER admissible sans terme stabilisateur
   manuel.
3. Extraire les cumulants connectés de Floquet à l'ordre annoncé et comparer
   exactement opérateur par opérateur au parent de référence.
4. Refuser ou poursuivre avant toute optimisation de paramètres à grande
   échelle.

## Frontière de claim

Ce document est un contrat de recherche. Il ne fournit encore aucun
Hamiltonien 2D natif, phase protégée, anyon, braid non abélien, jeu universel
de portes ou tolérance aux fautes.
