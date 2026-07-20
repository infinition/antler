# Phase 8B  -  fermeture stabilisateur 2D par marcheurs neutres

## Verdict

La boucle virtuelle d'un marcheur neutre à quatre positions génère un mot de
quatre liens au quatrième ordre. Un marcheur conditionné par `X` donne une
étoile `A_s`; un marcheur conditionné par `Z` donne une plaquette `B_p`.
Leur chevauchement microscopique est calculé, puis leurs coefficients sont
pavés sur le tore de référence `3 x 3`.

Cela fournit une construction **effective et conditionnelle** d'un parent
torique abélien. Les deux types de marcheur sont de nouvelles primitives
déclarées, non dérivées du ladder ANTLER gelé.

## Chevauchement local exact

Le bloc contient six liens et deux marcheurs, donc 1 024 états. Le Schur
complément exact dans le sous-espace bas des marcheurs donne seulement
`I`, `A_s`, `B_p` et `A_s B_p` :

| quantité profonde (`g/Delta <= 0.075`) | valeur |
|---|---:|
| puissance de `c_A` et `c_B` | `4.0365` |
| puissance de `c_AB` | `8.0514` |
| `||[A_s,B_p]||_F` | `0` |
| coefficient hors algèbre stabilisateur | `< 5.6e-17` |

À `g/Delta=0.10`, `c_A=c_B=-2.17068e-3`; le produit commutant vaut
`-1.34665e-6` dans la convention stabilisateur. Il est donc plus petit d'un
facteur `6.2e-4`.

## Fermeture effective du tore `3 x 3`

En supposant des copies indépendantes des gadgets locaux, le parent
`sum c_A A_s + sum c_B B_p + sum c_AB A_s B_p` est énuméré exactement sur les
`2^16` syndromes indépendants. À tous les rapports testés :

- le syndrome `A_s=B_p=+1` est unique ;
- la dégénérescence logique vaut 4 ;
- le rang stabilisateur est 16, avec deux qubits encodés et distance 3 ;
- le gap de syndrome est non nul et suit la puissance `4.0372`.

Le gap effectif est `8.704e-3` à `g/Delta=0.10`, mais seulement `3.141e-5` à
`0.025`. La suppression des corrections hors ordre réduit donc fortement le
gap : ce n'est pas une protection gratuite.

## Limites

Ce résultat ne simule pas le pavage microscopique complet ni le crosstalk entre
de nombreux marcheurs. Il n'établit aucune réalisation ANTLER native, aucun
gap thermodynamique, anyon mobile, fusion, braid non abélien, universalité ou
tolérance aux fautes. Le code torique effectif est abélien.

## No-go de tresse pour la famille actuelle

La modulation temporelle des seuls coefficients existants ne contourne pas ce
point. Tous les Hamiltoniens accessibles ont la forme d'un polynôme dans les
mêmes générateurs `{A_s,B_p,A_sB_p}`; ils commutent donc à tous les temps et
se projettent comme `E_code(t) I_4` sur le code. Les 21 paires de points
auditées ont un commutateur effectif et logique exactement nul.

La ressource minimale suivante est une primitive de déformation de défaut qui
change le support des stabilisateurs et échange les caractères électrique et
magnétique sur une coupure de branche, ou un autre générateur logique non
commutant **dérivé**. Une porte Hadamard, BdG, Majorana ou matrice de braid
insérée à la main ne satisfait pas ce contrat.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_star_plaquette_walker_overlap_audit.py
python experiments/phase7/run_phase8b_effective_toric_tile_closure.py
```
