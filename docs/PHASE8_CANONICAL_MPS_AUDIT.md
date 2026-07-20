# Phase 8  -  audit MPS canonique de la voie Floquet à nombre conservé

## Verdict actuel

La voie Floquet externe de Defossez *et al.* est maintenant reproduite avec
des MPS qui conservent **exactement** la charge totale `U(1)` et la parité de
branche `Z2`. Cela retire un faux positif important : avec `U(1)` seul, DMRG
glissait vers le secteur de parité fondamental et rendait le splitting
ininterprétable.

Le résultat est une validation sérieuse d'un **benchmark externe à modes de
Majorana à nombre conservé**, et non une réalisation native ANTLER, ni une
démonstration de braid ou d'ordinateur topologique.

## Contrôles numériques passés

À `U0=-2`, `alpha=0.5`, `N/(2L)=1/4`, les deux secteurs `P_A=+1,-1` sont
contraints dans les tenseurs.

| L | N | split de parité | plus petit gap neutre |
| ---: | ---: | ---: | ---: |
| 8 | 4 | `9.09936e-3` | `3.02086e-1` |
| 12 | 6 | `1.35878e-3` | `1.51653e-1` |
| 16 | 8 | `2.00566e-4` | `8.99732e-2` |

Le point `L=8` reproduit les deux énergies ED par secteur à mieux que
`2.7e-14`. Le splitting est stable lorsque `chi` passe de 256 à 384 :
changement `5.3e-15` à `L=12` et `9.95e-14` à `L=16`. Un fit descriptif aux
trois tailles donne une longueur de décroissance `xi_split ~= 2.10`; ce n'est
pas une extrapolation thermodynamique certifiée.

Le gap neutre diminue. Ce n'est **pas** un rejet de cette classe de modèle :
la référence prévoit un mode de charge total gapless (`c=1`) et un secteur
topologique relatif gappé. Employer le gap neutre comme critère de rejet
aurait donc été une erreur de diagnostic.

## Diagnostic adapté : gap chargé périodique

Sur un anneau, l'estimateur fini

\[
\Delta_{\rm charge}^{(L)}=
\frac{E(N+1)+E(N-1)-2E(N)}{2}
\]

évite les modes de bord OBC qui contaminent le changement de charge. Au point
mis en avant dans la référence, `U0=-1.5`, `alpha=0.5`, `nu=1/3`, il vaut :

| L | N | `Delta_charge^(L)` |
| ---: | ---: | ---: |
| 6 | 4 | `0.38114` |
| 9 | 6 | `0.35933` |

Deux tailles ne déterminent pas la limite infinie. Elles montrent seulement
qu'un diagnostic chargé, distinct du gap neutre gapless, reste non nul dans
le petit système périodique reproduit.

## Point publié et séquence de pulses réelle

Au même point publié, le splitting OBC à densité fixe est `0.21215`,
`0.07579`, `0.03575` pour `L=6,12,18`. Le point `L=6` est étalonné contre
l'ED à `2e-14`. Le calcul `L=18` atteint `chi=384`; une tentative de contrôle
`chi=512` a dépassé la limite d'exécution avant sérialisation et reste donc
**incomplète**.

Le cycle exact

\[
U(T)=P^\dagger e^{-i(1-\alpha)TH_0}P e^{-i\alpha TH_0}
\]

est aussi contrôlé, à `L=6,N=4`, contre `exp(-i T H_eff)` : son erreur par
période suit `T^1.976`. À l'angle idéal `eta=pi/2`, le commutateur de parité
stroboscopique est inférieur à `5e-16`; un décalage `eta-pi/2=0.1` rad donne
au contraire un résidu effectif `1.12e-1`. Le protocole dynamique est donc
algébriquement correct mais sa calibration est une exigence physique réelle.

## Portée et prochaines portes

Cette archive établit une reproduction contrôlée du mécanisme Floquet
externe : interaction attractive de jambe + rotation globale de rail pulsée
génèrent le pair-hopping au premier ordre effectif, sans l'impasse
perturbative statique `g^4` contre parasites `g^2` de Phase 7C.

Les conditions suivantes restent indépendantes et obligatoires avant toute
promotion ANTLER :

1. Convergence longue taille / grand `chi` et extraction du gap topologique
   chargé avec une méthode MPS infinie ou une alternative périodique contrôlée.
2. Dérivation microscopique ANTLER de l'interaction attractive de jambe et de
   la rotation de rail globale, avec les parasites et la fréquence finie.
3. Intégration des erreurs de pulse au ladder entier, pas seulement au bloc
   `L=6`.
4. Une jonction en T et un audit de braid non commutatif. Une échelle 1D ne
   suffit pas à elle seule.

Résultats machine :

- `results/phase7/tenpy_dmrg_parity_scaling.json`
- `results/phase7/tenpy_dmrg_bond_convergence.json`
- `results/phase7/tenpy_neutral_gap_audit.json`
- `results/phase7/pbc_charged_gap_ed.json`
- `results/phase7/tenpy_published_point_parity_scaling.json`
- `results/phase7/finite_pulse_stroboscopic_audit.json`

Référence primaire : Defossez, Vanderstraeten, Peralta Gavensky et Goldman,
“Dynamic Realization of Majorana Zero Modes in a Particle-Conserving Ladder”,
arXiv:2412.14886v2 (2025).
