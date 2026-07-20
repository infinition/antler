# Phase 7D  -  écho par inversion directe de phase de Peierls

## Idée

L'écho de phase alternée précédent obtenait (H_{\rm leg}\to-H_{\rm leg})
par des kicks de potentiel. Une autre implémentation mathématiquement exacte
est une phase de Peierls synchrone `phi_leg=0` ou `pi` sur tous les hoppings de
jambe, sans modifier les conversions de paire charge-2.

Dans la matrice microscopique complète, cette réalisation directe est vérifiée
indépendamment :

\[
QH_+Q-H(\phi_{\rm leg}=\pi)=0
\]

en norme de Frobenius pour les canaux same-rail et opposite-rail. Elle produit
donc exactement le même Hamiltonien d'écho que le contrôle idéal, sans
potentiels de site de très grande amplitude.

## Résultat 472D

Le cas sévère reste `t_leg=1` et crosstalk `0.01g` :

| sous-cycles par pulse médiateur | fuite | résidu de parité | distance logique | bascules `pi` complètes |
| ---: | ---: | ---: | ---: | ---: |
| 16 | `6.97e-8` | `1.72e-5` | `3.12e-5` | 128 |
| 32 | `5.74e-8` | `4.29e-6` | `1.08e-5` | 256 |

La valeur singulière logique minimale au dernier point vaut `0.9999999713`.

## Décision

Cette voie qualifie un **pont de contrôle dynamique conditionnel** : si un
contrôle peut commuter de façon synchrone les phases de Peierls des jambes de
`0` à `pi` sans altérer la conversion charge-2, elle réalise la correction
locale nécessaire sans les kicks de potentiel qui ont échoué.

Ce contrat n'est pas encore une dérivation de matériel. Le précédent flux
`pi` local d'ANTLER montre que les phases de Peierls constituent une ressource
pertinente, mais ne démontre pas leur commutation globale rapide, leur bruit,
ni leur indépendance des canaux de paire. Les 256 changements de signe du
point le plus précis rendent cette exigence explicite.

Le premier audit CPU de rampes de phase finies est préservé comme incomplet,
mais sa complétion GPU est désormais disponible dans
`PHASE7D_PEIERLS_RAMP_GPU_AUDIT.md`. Les rampes courtes passent le bloc 472D
et leur représentation 2/4/8 segments est convergée. Cette validation est
strictement numérique : une dérivation du contrôle matériel, son bruit et son
synchronisme restent requis. Aucun résultat de phase protégée, code 2D, braid,
non-abélianité, universalité ou tolérance aux fautes n'est déduit ici.

Résultat machine : `results/phase7/direct_peierls_echo_audit.json`.
