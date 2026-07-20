# Phase 8B  -  audit de convention : pont conditionnel par corde JW

## Question

Un pont local dont la corde traverse le walker `d1` semble donner, a
`theta=pi`, le terme tres utile

```text
X_rail Z_walker.
```

Sur le petit bloc exact, ce terme est effectivement propre et sa propagation
peut etre excellente. Mais ceci n'est une prediction ANTLER que si la corde
compte `d1` avec une charge statistique impaire. Cette condition doit etre
verifiee contre les conventions gelees avant toute promotion.

## Controle de convention

Le modele gele fixe compte les occupations de rails unitaires dans une chaine
rung-major; un saut de barreau adjacent n'a **aucune** corde. Les mediateurs
charge-2 sont une extension separee et ne sont pas definis comme des sites
impairs de cette corde.

Trois choix explicites ont donc ete compares :

| poids de corde de `d1` | interpretation | `c_XZ` |
|---:|---|---:|
| 0 | convention gelee du saut de barreau | `0` |
| 2 | comptage de charge U(1) physique du mediateur moleculaire | `0` |
| 1 | mediateur hard-core artificiellement impair dans la corde | `-0.2` |

Le cas poids 1 donne bien trois portes locales passees a `Mott=120`, avec
leakage jusqu'a `3.41e-5`. Mais ce sont des **contre-factuels algebriques** :
ils ne sont ni une consequence de la corde rung-major gelee, ni une
realisation derivee d'un mediateur de charge 2.

## Resultat utile

Ce controle ne qualifie pas une porte ANTLER, mais identifie sans ambiguite
la ressource manquante : un lien de jauge `Z2` neutre (ou un degre de liberte
statistique impair explicitement derive) qui module le hopping de rail par
une phase `pi`.

Cette conclusion rejoint, sans le remplacer, le theoreme de ressource Gauss
Phase 8B : un degree de liberte de jauge supplementaire est necessaire. Ici,
le test precise pourquoi un mediateur moleculaire charge-2 ne peut pas jouer
ce role simplement par sa presence dans une corde JW.

## Consequence pour la suite

La prochaine construction defendable devra choisir explicitement l'une des
deux routes :

1. deriver un vrai lien de jauge `Z2` neutre et son terme de hopping habille,
   avec budget de parasites et loi de Gauss ; ou
2. produire un substrat ou une statistique de mediateur impair est physique,
   puis refaire l'audit depuis le Hamiltonien microscopique.

Ni RL ni une optimisation de pulses ne peuvent transformer le cas poids 0 ou
2 en cas poids 1 : la difference est une ressource, pas un parametre.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_jw_string_conditional_rung_bridge.py
```

## Claim boundary

Le `X Z` au poids 1 est garde comme controle algebrique d'une ressource
supplementaire. Il n'etablit pas une porte native du ladder gele, un espace de
code, une protection topologique, des defauts, une fusion, une tresse non
abelienne, l'universalite ou un ordinateur quantique topologique.
