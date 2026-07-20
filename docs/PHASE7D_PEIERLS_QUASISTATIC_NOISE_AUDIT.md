# Phase 7D  -  bruit quasi-statique global sur la primitive Peierls

50 réalisations déterministes par niveau ont été propagées exactement sur le
bloc 472 états (RTX 4070 Ti, `complex128`). Chaque réalisation tire une erreur
commune de plateau de phase et un déséquilibre de durée communs à la séquence
complète : c'est un modèle de dérive lente globale, ni du jitter cycle-à-cycle,
ni du bruit indépendant par lien.

| sigma | E[F] | Var(F) | E[fuite] | passages de la cible locale `1e-4` |
| ---: | ---: | ---: | ---: | ---: |
| `1%` | `0.99999839` | `2.45e-12` | `6.04e-6` | `50/50` |
| `3%` | `0.99998584` | `1.95e-10` | `5.36e-5` | `35/50` |
| `5%` | `0.99996086` | `1.49e-9` | `1.48e-4` | `20/50` |
| `10%` | `0.99984524` | `2.29e-8` | `5.82e-4` | `8/50` |
| `20%` | `0.99940648` | `3.15e-7` | `2.20e-3` | `2/50` |

La courbe donne un seuil de rupture local plutôt qu'une simple réussite à un
point : la cible stricte reste majoritairement satisfaite à `3%`, mais plus à
`5%`. Elle ne qualifie pas un matériel réel, un ladder long, une phase
protégée ou une tolérance aux fautes.

Résultat machine : `results/phase7/peierls_quasistatic_noise_gpu_audit.json`.
