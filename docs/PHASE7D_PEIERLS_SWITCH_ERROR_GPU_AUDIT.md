# Phase 7D  -  tolérance locale aux erreurs de commutation Peierls

## Question testée

Après la convergence 2/4/8 de la rampe finie, cet audit demande quelle erreur
déterministe de commande peut être tolérée par la primitive locale. Le même
bloc microscopique de 472 états est propagé exactement sur RTX 4070 Ti, en
`complex128`, avec 16 sous-cycles par impulsion médiateur et un crosstalk de
`0.01g`.

Deux défauts à un paramètre sont isolés :

1. une erreur systématique de la phase supposée égale à `pi` sur le plateau ;
2. un déséquilibre temporel signé, qui conserve la durée de chaque sous-cycle
   mais augmente les deux segments `phi=0` et diminue le segment `phi=pi`.

La cible locale préenregistrée est : fuite, résidu de parité et distance
logique au plus `1e-4`, avec valeur singulière minimale au moins `0.9999`.

## Résultats

| défaut isolé | dernier point testé passant | premier point testé échouant | interprétation correcte |
| --- | ---: | ---: | --- |
| décalage systématique du plateau `pi` | `5 deg` | `7.5 deg` | fenêtre locale encadrée, non interpolée |
| déséquilibre temporel signé | `5 %` | `10 %` | fenêtre locale encadrée, non interpolée |

Au dernier point de phase passant (`5 deg`), la fuite vaut `3.85e-5`, le plus
grand résidu de parité `3.74e-5`, la distance logique `4.72e-5` et la valeur
singulière minimale `0.9999807`. À `7.5 deg`, la distance logique devient
`1.42e-4`, donc dépasse la cible malgré une fuite encore sous `1e-4`.

Au dernier point temporel passant (`5 %`), la fuite vaut `4.83e-5`, le plus
grand résidu de parité `5.19e-5`, la distance logique `7.07e-5` et la valeur
singulière minimale `0.9999759`. À `10 %`, fuite, parité et distance logique
dépassent la cible.

## Décision

La primitive possède une marge locale non nulle contre ces deux erreurs
**systématiques et isolées**. Le résultat est volontairement rapporté comme
un intervalle de balayage (`<=5` testé passant, premier échec `7.5` ou `10`),
et non comme une spécification continue ou une tolérance expérimentale.

## Limites de promotion

Le test ne contient ni bande passante réelle, ni bruit gaussien ou corrélé,
ni jitter cycle-à-cycle, ni erreur différentielle entre liens/jambes, ni
calibration, ni dynamique du ladder complet. Il ne démontre donc pas une
implémentation matérielle, une phase protégée, un mode de bord, un code 2D,
un braid, une non-abélianité, l'universalité ou la tolérance aux fautes.

Résultat machine : `results/phase7/peierls_switch_error_gpu_audit.json`.
