# Phase 7D  -  audit microscopique à quatre barreaux des pulses charge-2

## Objet

Le compilateur multi-liens précédent agissait après fermeture de chaque pulse
dans l'espace logique. Cet audit traite explicitement les pulses : huit modes
fermioniques de rail, douze médiateurs hard-core de charge 2, contrainte de
charge totale quatre et terme de Mott sur chaque barreau. La dimension exacte
du bloc est 472.

La séquence est : couche paire `(0,1),(2,3)`, puis lien impair `(1,2)` ; sur
chaque couche, pulse même-rail, rotation de rail, pulse opposé-rail et retour.
Elle emploie le même point local enregistré : `U=20`, `Delta=40`, `g=6`,
`g/Delta=0.15`, temps Rabi complet `0.0776708795`.

## Modèle de crosstalk

Les canaux non sélectionnés ont reçu une amplitude cohérente résiduelle
`epsilon*g`. Ce n'est pas un modèle de bruit expérimental complet : il isole
l'erreur la plus directe de l'hypothèse « canal inactif ». Les valeurs testées
sont `epsilon=0, 1e-4, 1e-3, 3e-3, 1e-2`.

## Résultats

À `epsilon=0`, la fuite du schedule complet est `7.34e-29`, le défaut de
factorisation de la couche paire `9.93e-14`, l'erreur de composition de pulses
fermés `6.20e-14` et la plus petite valeur singulière logique vaut
`0.999999999999961`.

Sous le modèle de crosstalk déclaré, les erreurs suivent le comportement
quadratique attendu dans cette petite fenêtre. Au point le plus élevé,
`epsilon=0.01`, on obtient :

- fuite monomère finale : `5.67e-8` ;
- défaut de composition : `1.31e-8` ;
- défaut de factorisation de la couche paire : `1.62e-8` ;
- résidu de chaque parité logique : `1.60e-8` ;
- valeur singulière logique minimale : `0.9999999716`.

## Test séparé : hopping de jambe pendant les pulses

Le hopping intrajambe conserve la parité de rail en l'absence de médiateur,
mais les canaux opposé-rail ne conservent la parité nue qu'après une excursion
Rabi exactement fermée. Il fallait donc le simuler pendant les pulses plutôt
que présumer sa bénignité.

Le résultat est un rejet net de la séquence brute :

| `t_leg` pendant les pulses | fuite finale | défaut de parité logique | écart à la porte sans hopping |
| ---: | ---: | ---: | ---: |
| `0.1` | `2.02e-4` | `2.66e-4` | `4.02e-4` |
| `0.3` | `1.82e-3` | `2.39e-3` | `3.61e-3` |
| `1.0` | `1.98e-2` | `2.62e-2` | `3.98e-2` |

La plus petite valeur singulière tombe à `0.99003` au dernier point. Le
problème est cohérent : le hopping pendant l'excursion médiateur perturbe la
fermeture Rabi qui restaurait la parité de rail nue. Il ne s'agit pas d'un
simple défaut global de phase.

## Décision

Le contrat local du compilateur résiste numériquement au crosstalk cohérent
jusqu'à 1 % dans **ce modèle précis, avec hopping de jambe éteint**. En
revanche, le même audit falsifie la séquence brute dès que le hopping de jambe
est laissé actif. La prochaine dérivation doit donc fournir un mécanisme de
gel/refocalisation de ce hopping, ou des pulses recalibrés dans le Hamiltonien
complet, avant tout audit de phase.

Cette obligation a été satisfaite comme contrôle idéal dans
`PHASE7D_STAGGERED_ECHO_REFOCUSING_AUDIT.md` : une phase alternée de rail,
accompagnée d'une phase `pi` de médiateur, inverse uniquement le hopping de
jambe et restaure une fermeture convergente. Sa disponibilité physique reste
une hypothèse à dériver, pas une ressource ANTLER déjà démontrée.

Cela ne transforme pas la primitive en phase topologique et n'annule pas
l'échec de localisation observé dans le benchmark Floquet de petite taille.
Aucun résultat non abélien, universel ou tolérant aux fautes n'est déduit ici.

Résultat machine : `results/phase7/four_rung_microscopic_pulse_audit.json`.
