# Phase 8  -  pont microscopique natif vers le ladder Floquet

## Résultat

Le no-go Phase 7C concernait la création **statique** d'un stabilisateur à
quatre corps par médiateurs virtuels. Il ne bloque pas la ressource plus
modeste requise par le protocole Floquet : une attraction de densité de jambe
à deux corps. Un pont microscopique local vers le Hamiltonien de départ
Floquet est maintenant dérivé et audité.

Pour deux rails `a,b`, on vise

\[
H_0=-t\sum_{j,\ell}(c^\dagger_{j,\ell}c_{j+1,\ell}+h.c.)
 +U_0\sum_{j,\ell}n_{j,\ell}n_{j+1,\ell},\qquad U_0<0.
\]

Les paramètres ANTLER correspondants sont `theta=pi`, `J1=J2=t` et
`Jperp=0` pendant les segments libres. La matrice de hopping rung-major
ANTLER coïncide alors avec le hopping fermionique de référence avec une norme
de Frobenius `2.45e-16` sur le bloc `L=2,N=2`.

## Dérivation SW constructive

Chaque lien reçoit **deux médiateurs charge-2 indépendants**, un par rail :

\[
V=-g\sum_j\left[d^\dagger_{a,j}a_ja_{j+1}
                 +d^\dagger_{b,j}b_jb_{j+1}+h.c.\right],
\qquad H_d=\Delta\sum_{j,\ell}d^\dagger_{\ell,j}d_{\ell,j}.
\]

Dans le sous-espace sans médiateur,

\[
H_{\rm SW}^{(2)}=-\frac{g^2}{\Delta}
\sum_j\left(n_{a,j}n_{a,j+1}+n_{b,j}n_{b,j+1}\right).
\]

Le fait que les médiateurs soient séparés est décisif : un médiateur partagé
produirait aussi `a†a†bb+h.c.` au même ordre. Ici, le résidu de l'opérateur SW
visé et l'élément parasite `aa -> bb` sont exactement nuls dans l'audit.

À `U0=-1.5`, l'erreur spectrale exacte du bloc diminue comme
`(g/Delta)^1.966` : `3.23e-2`, `1.47e-2`, `3.73e-3` pour
`g/Delta=0.15,0.10,0.05`, tandis que la capture du sous-espace bas atteint
`0.979`, `0.990`, `0.998`.

Une impulsion de rung `Jperp=-1/2` de durée `pi/2` réalise exactement la
rotation globale `P=exp(-i pi Jx/2)` sur ce bloc. Le signe inverse ou une
impulsion Peierls équivalente est donc une ressource de contrôle nécessaire.

## Composition microscopique et fermeture Rabi

L'audit exact sur `L=3,N=2`, avec quatre médiateurs explicites, compare le
cycle microscopique à `exp(-i T H_eff)`. Une durée libre générique laisse une
micromotion de médiateur. À `g/Delta=0.05`, son enveloppe décroît comme
`(g/Delta)^1.961` à `U0` fixé : c'est une correction SW contrôlée, mais elle
ne doit pas être ignorée.

Le choix analytique

\[
T_m=\frac{4\pi m}{\sqrt{\Delta^2+4g^2}}
\]

ferme un nombre entier de cycles Rabi virtuels durant chaque demi-période.
Au plus court point (`m=1`, `T=0.0208400`) :

| métrique | valeur |
| --- | ---: |
| fuite microscopique pire cas | `2.16e-6` |
| distance polaire au cycle `H_eff` | `7.33e-5` |
| plus petite valeur singulière | `0.99999892` |

La fenêtre locale stricte (`fuite < 1e-4`, distance `< 1e-4`) est symétrique
sur les erreurs de période testées `[-0.3%, +0.3%]`. À `±1%`, la distance
cohérente échoue d'abord ; à `±3%`, fuite et distance échouent.

## Limite rencontrée, explicitement conservée

La primitive se compose, mais pas encore de façon innocente. Au point fermé
le plus court, la fuite évolue

\[
2.16\,10^{-6}\;\to\;8.64\,10^{-6}\;\to\;3.45\,10^{-5}
\;\to\;1.37\,10^{-4}
\]

pour `n=1,2,4,8`; l'ajustement donne `P_leak ~ n^1.996`, signe d'une
accumulation cohérente. L'erreur unitaire croît comme `n^0.997`.

Une composition Strang symétrique choisie pour fermer chaque sous-segment
réduit la distance cohérente à temps égal (`1.47e-4 -> 8.63e-5` au premier
cycle comparé), mais ne réduit pas la fuite virtuelle. Une séquence composite
ou un contre-terme/dressing stroboscopique est donc le prochain travail
microscopique, avant toute montée en taille.

## Second segment synthétisé directement et contrôle en taille

Une hypothèse testable était que l'accumulation provenait principalement de la
rotation physique abrupte `P` entre les deux demi-périodes. Elle est maintenant
écartée dans le bloc contrôlé. Le terme d'interaction de
`P^dag H0 P` est factorisé exactement en deux canaux de conversion de paires
cohérents : l'erreur de factorisation locale est `3.51e-16`. Les deux segments
Floquet sont donc construits directement dans le même espace nu, sans
impulsion `P` intermédiaire.

Cela ne donne pas une suppression qualitative de la fuite à ratio fixé. Sur
`L=3,N=2`, elle vaut `2.00e-6` par cycle et croît à `1.27e-4` après huit
cycles, très près du contrôle à impulsion. Au remplissage de référence
`nu=1/4`, le test exact direct donne `3.00e-6` sur `L=4,N=2` et
`1.57e-5` sur `L=6,N=3` par cycle. Après quatre cycles, ce dernier point
atteint `2.50e-4`. Le micromouvement virtuel des médiateurs, et non le seul
changement de frame, reste donc le mécanisme dominant dans cette réalisation.

Il existe néanmoins un paramètre de contrôle propre. À `L=6,N=3`, on fixe
`U0=-g^2/Delta=-1.5` et on compare des suites de fermetures Rabi dont les
durées totales ne diffèrent que de `0.47%` au plus :

| `g/Delta` | cycles fermés | `Delta`, `g` | fuite pire cas | distance polaire |
| ---: | ---: | ---: | ---: | ---: |
| `0.0500` | 1 | `600`, `30` | `1.57e-5` | `1.43e-4` |
| `0.0250` | 4 | `2400`, `60` | `4.06e-6` | `3.61e-5` |
| `0.0125` | 16 | `9600`, `120` | `1.02e-6` | `9.06e-6` |

Les ajustements donnent respectivement
`P_leak ~ (g/Delta)^1.970` et `distance ~ (g/Delta)^1.989` pour cette durée
logique approximativement fixée. C'est une convergence numérique utile : la
profondeur SW contrôle aussi l'erreur **accumulée**, pas seulement une fuite
par cycle. Ce n'est pas une solution matérielle gratuite : elle exige ici un
détuning multiplié par 16, un couplage multiplié par 4 et une cadence de
commutation multipliée par 16. Leur faisabilité, leur bande passante et leur
bruit ne sont pas dérivés.

## Premier budget de contrôle des canaux directs

Le point profond a ensuite été soumis à un changement **continu** des canaux
de paires, au lieu d'une commutation instantanée. La trajectoire conserve
exactement deux canaux orthonormés : elle tourne
`(aa+bb)/sqrt(2)` vers `(ab+ba)/sqrt(2)` et conserve
`(aa-bb)/sqrt(2)`. Les deux extrémités reproduisent les projecteurs voulus à
`3.16e-16` et `3.51e-16`.

Sur `L=3,N=2`, `g/Delta=0.0125` et 16 cycles (la durée logique du contrôle
profond), les rampes linéaires et `sin^2` jusqu'à `1%` de la période **par
transition** gardent simultanément `P_leak < 1e-4` et la distance logique
`<1e-4`. À `1%`, la rampe `sin^2` donne `P_leak=5.28e-5`, contre
`9.46e-5` pour la rampe linéaire. À `3%`, les deux échouent (`4.73e-4` pour
`sin^2`, `8.45e-4` pour la rampe linéaire). Le contrôle de discrétisation est
convergé : le passage de 16 à 32 segments de rampe modifie la fuite relative
de seulement `3.99e-6`.

Une erreur statique sur l'angle terminal du canal `H1` donne une contrainte
différente et plus sévère. Avec le même test, le seuil conjoint fuite,
distance et commutateur de parité ne passe que pour
`delta_phi in {-0.003, 0, +0.003}` rad, soit environ `+/-0.172 deg`.
À `+/-0.01` rad (`+/-0.573 deg`), la fuite est encore basse mais le
commutateur de parité vaut `1.62e-4` : ce n'est donc pas une porte protégée
au seuil enregistré. Ces nombres sont un **budget de calibration cohérente**
sur un petit bloc idéal, non une qualification de bruit matériel.

Le crosstalk est également injecté au niveau microscopique : un médiateur
assigné au lien `j` peut convertir une paire sur `j-1` ou `j+1` avec amplitude
résiduelle `epsilon*g`. Au point profond et à la durée logique enregistrée,
le scan déterministe passe jusqu'à `epsilon=0.003` sur `L=3,N=2` et
`L=4,N=2`, puis échoue sur la distance logique à `epsilon=0.01`. La fuite et
le commutateur de parité restent faibles dans ce modèle ; c'est donc une erreur
cohérente de Hamiltonien à calibrer, non une simple perte détectable.

Une campagne reproductible de 100 réalisations par niveau, sur `L=4,N=2`,
complète ce premier seuil. Chaque réalisation emploie une même erreur complexe
aléatoire sur toutes les conversions voisines, avec
`E|epsilon|^2=sigma^2` :

| `sigma` RMS de `epsilon/g` | passes stricts / 100 | erreur logique moyenne |
| ---: | ---: | ---: |
| `0.1%` | `100` | `1.70e-5` |
| `0.2%` | `100` | `3.58e-5` |
| `0.3%` | `89` | `5.30e-5` |
| `0.5%` | `72` | `7.98e-5` |
| `1.0%` | `39` | `1.48e-4` |

La distribution, les variances et les maximums sont sérialisés. Elle ne
représente qu'un crosstalk cohérent **totalement corrélé** entre liens ; les
erreurs indépendantes, les fluctuations temporelles et la combinaison avec
une rampe finie restent des tests séparés.

## Statut scientifique

Le résultat établit un **pont natif local contrôlé** : les ressources ANTLER
étendues par deux médiateurs charge-2 indépendants par lien, le hopping
rung-major à `theta=pi`, et une rotation de rung pulsée peuvent reproduire
les ingrédients du protocole Floquet à nombre conservé de Defossez *et al.*
au niveau microscopique fini.

Il ne démontre pas encore un Hamiltonien ANTLER étendu complet, un gap
thermodynamique, une mémoire topologique, une tresse, une non-abélianité,
l'universalité ou la tolérance aux fautes. Les prochaines portes sont :

1. dériver puis auditer une implémentation finie-bande passante des canaux de
   paires directs au point SW profond sur un bloc plus grand avec erreurs de
   lien indépendantes et rampes finies simultanées ;
2. tester les fluctuations temporelles, la bande passante et la calibration à
   cette cadence sur davantage de liens, sans confondre une fermeture Rabi
   idéale avec une porte robuste ;
3. injecter les corrections microscopiques dans le MPS à grande taille et
   vérifier le gap chargé ;
4. seulement alors, construire une jonction en T et tester une algèbre de
   braid non commutative.

Résultats machine :

- `results/phase7/native_h0_sw_bridge.json`
- `results/phase7/native_micro_floquet_l3.json`
- `results/phase7/native_micro_detuning_scaling.json`
- `results/phase7/native_micro_rabi_closure.json`
- `results/phase7/native_micro_rabi_timing_tolerance.json`
- `results/phase7/native_micro_rabi_composition.json`
- `results/phase7/native_micro_strang_closure.json`
- `results/phase7/native_direct_h1_closure.json`
- `results/phase7/native_direct_h1_size_audit.json`
- `results/phase7/native_direct_h1_fixed_time_scaling.json`
- `results/phase7/direct_channel_ramp_audit.json`
- `results/phase7/direct_channel_angle_audit.json`
- `results/phase7/direct_channel_crosstalk_audit.json`
- `results/phase7/direct_channel_crosstalk_ensemble.json`
