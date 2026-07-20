# Phase 8B  -  obstruction de moyenne Floquet du lien shared-matter

## Resultat

Le lien conditionnel `XX` de l'audit shared-matter est valide **comme
moyenne arithmetique de deux Hamiltoniens deja downfoldes** :

```text
H_static = (H_eff(A) + H_eff(B)) / 2.
```

Cette egalite ne definit pas une porte Floquet. Un drive rapide realise au
premier ordre la moyenne microscopique

```text
H_micro = (H(A) + H(B)) / 2,
```

qu'il faut ensuite downfolder. Dans la grammaire enregistree
`A=(+J,+g,+g)`, `B=(-J,+g,-g)`, ces deux operations ne commutent pas :

```text
(H_eff(A)+H_eff(B))/2 != H_eff((H(A)+H(B))/2).
```

La premiere expression contient `XX`, tandis que la seconde a exactement
`c_XX = 0` et conserve `IZ`.

## Mesure exacte

A `g/Delta=0.05` :

| quantite | valeur |
|---|---:|
| `c_XX` de la moyenne des SW | `-1.25217e-4` |
| `c_XX` du SW de la moyenne microscopique | `0` |
| terme restant dominant | `c_IZ=-6.25511e-3` |
| rapport `|IZ|/|XX_static|` | `49.9542` |

Le `XX` statique suit bien la loi profonde `(g/Delta)^2`, mais il est absent
dans la limite Floquet rapide de ce mot de signes. La convergence est mesuree
sur le propagateur microscopique complet a duree totale 100 : la distance a
`exp(-i T H_micro)` tombe de `2.42e-3` (`tau=0.01`) a `6.03e-6`
(`tau=2.5e-4`). Le drive converge donc vers la mauvaise moyenne, et non vers
le compilateur statique.

## Extension equilibree a quatre signes

La meilleure correction interne du mot a deux segments utilise les quatre
signes

```text
(s_J,s_0,s_1) = (+,+,+), (-,+,-), (-,-,+), (+,-,-),
```

avec `s_J=s_0*s_1`. Chaque SW de segment porte alors le meme `XX`, et les
termes isoles s'annulent dans la moyenne statique. Surtout, **tous** les
controles microscopiques ont une moyenne nulle : contrairement au mot a deux
signes, le `IZ` rapide disparait aussi.

Le resultat exact est toutefois negatif pour la porte : la moyenne statique
reste `XX` propre, mais le SW de la moyenne microscopique est entierement
scalaire (`XX=0` et tout terme non scalaire nul a la precision machine). Les
30 protocoles lisses associes ne passent pas : meilleure erreur relative
`0.999239`, avec des points a faible fuite mais sans signal logique. Cette
extension ferme donc le raffinage par simple groupe de signes; elle ne
contredit pas une construction ou un terme conditionnel apparaitrait dans un
ordre Magnus derive et quantitativement controle.

## Tests de contournement honnetes

Une rampe `sin^2` zero-a-zero autour de chaque segment a ete testee sur 30
points exacts (`g/Delta=0.05,0.025`, cible `|theta_XX|=0.1`). Aucun ne passe
simultanement l'erreur relative au signal, la purete `XX` et la fuite. La
meilleure erreur relative reste `0.999278`; les points a faible fuite
retournent donc presque un scalaire.

Avant toute tentative RL, une baseline classique a egalement ete faite dans
la boite physique a deux maintiens longs (rampes fixees, durees libres
`[0,400]`, trois graines de differential evolution). Aucun run ne trouve la
porte `exp(-i 0.1 XX)` : les erreurs relatives sont `>= 1.000318`, quoique la
fuite puisse etre inferieure a `1e-5`. Le budget d'iterations est conserve
tel quel : c'est un baseline negatif reproductible, pas une preuve globale
d'absence de pulse.

## Decision de conception

Le probleme n'est pas encore une recherche de parametres : dans cette
algebre de controles, la limite rapide annule le processus virtuel que la
moyenne apres SW mettait en avant. Deep RL ne doit donc pas etre employe pour
« forcer » ce `XX` : il optimiserait une grammaire dont la limite physique
connue ne contient pas le signal cible.

La prochaine hypothese utile doit fournir, au niveau **microscopique**, une
sequence ou un cadre de controle dont la moyenne rapide conserve le terme
conditionnel et satisfait avant optimisation :

1. `c_XX != 0` dans `H_eff(H_micro,average)` ;
2. les termes hors cible sont bornes relativement a `c_XX` ;
3. U(1), localite et les ressources declarees sont preserves ;
4. le terme n'est pas introduit a la main.

Une fois cette grammaire derivee, une baseline classique reste prioritaire;
RL pourra alors raffiner une famille de pulses deja physiquement admissible,
notamment sous bruit et contraintes de bande passante.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_shared_matter_floquet_averaging_obstruction.py
python experiments/phase7/run_phase8b_shared_matter_adiabatic_echo_preflight.py
python experiments/phase7/run_phase8b_shared_matter_classical_control_baseline.py
python experiments/phase7/run_phase8b_shared_matter_four_sign_group_preflight.py
```

## Claim boundary

Ce document ferme uniquement la realisation Floquet rapide et les boites de
controles lisses explicitement testees pour le mot shared-matter `A/B`.
Il ne constitue pas un no-go pour une nouvelle ressource microscopique, un
kick/cadre derive, un autre controle Floquet, un code 2D, des defauts, la
fusion, une tresse non abelienne, l'universalite ou la tolerance aux fautes.
