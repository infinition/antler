# Lemme de transfert digital  -  convention, composition et corrections finies

## Portée exacte

Ce lemme concerne la navette digitale séquentielle du Hamiltonien ANTLER :
deux bosons hard-core sur l’échelle SSH, avec saut corrélé dans l’ordre
Jordan–Wigner **rung-major**. Il explique la limite localisée et adiabatique
de cette convention de Hamiltonien. Il ne démontre ni une statistique anyonique
indépendante de convention, ni un espace de fusion non abélien, ni
l’universalité du calcul.

Les scripts de clôture Phase 4.7 associés sont :

```text
experiments/phase4_7/run_phase47_deep_limit.py
experiments/phase4_7/run_phase47_path_invariance.py
experiments/phase4_7/run_phase47_composition.py
```

Ils emploient sans le modifier le constructeur ED `L=14, N=2`, le code chat,
la reconstruction de la matrice logique et la décomposition polaire validés
dans les phases précédentes.

## Lemme local de transfert

Considérons d’abord un handoff isolé dans la base ordonnée `(source, cible)` :

\[
H_\varphi(t)=
\begin{pmatrix}
\epsilon_1(t) & -J\,e^{+i\varphi}\\
-J\,e^{-i\varphi} & \epsilon_2(t)
\end{pmatrix},
\qquad J>0.
\]

On suppose que `epsilon_1(0)<epsilon_2(0)`,
`epsilon_2(T)<epsilon_1(T)`, que le gap évité ne se ferme pas, et que le
transfert est adiabatique. Posons

\[
G_\varphi=\operatorname{diag}(e^{i\varphi},1).
\]

Alors

\[
H_\varphi(t)=G_\varphi H_0(t)G_\varphi^\dagger.
\]

Comme `G_phi` est indépendant du temps, le propagateur adiabatique obéit à la
même conjugaison. Si, pour `varphi=0`, le handoff envoie

\[
|\mathrm{source}\rangle\longmapsto
e^{i\gamma_{\rm dyn}}e^{i\gamma_0}|\mathrm{cible}\rangle,
\]

alors

\[
|\mathrm{source}\rangle\longmapsto
e^{-i\varphi}e^{i\gamma_{\rm dyn}}e^{i\gamma_0}
|\mathrm{cible}\rangle.
\]

`gamma_dyn` ne dépend que des valeurs propres, donc pas de `varphi`; `gamma_0`
est une phase de convention à `varphi=0`. Le facteur orienté est donc
exactement `exp(-i varphi)`, indépendamment de la forme de rampe, dès que les
hypothèses à deux niveaux et adiabatiques sont vérifiées.

### Attention indispensable sur le signe de la matrice écrite

La matrice mentionnée dans le cahier de route,

\[
\begin{pmatrix}
\epsilon_1 & J e^{-i\alpha}\\
J e^{+i\alpha} & \epsilon_2
\end{pmatrix},
\]

est, à un signe réel de `J` près, le cas `varphi=-alpha` de la matrice
précédente. Dans la base ordonnée littéralement `(1,2)`, elle donne donc
`|1> -> exp(+i alpha)|2>` (à phase dynamique commune près), et non
`exp(-i alpha)`. Cette différence est une convention de signe, mais elle doit
être écrite correctement dans une preuve.

Pour ANTLER, un saut longitudinal croissant reçoit la phase
`exp(-i theta n_mid)` dans la convention rung-major : le handoff pertinent
est donc le cas `varphi=theta`, d’où le facteur transporté `exp(-i theta)`.
Équivalemment, on peut conserver la matrice ci-dessus en définissant
`alpha=-theta`. Inverser l’orientation du lien inverse aussi le signe.

## Composition orientée des liens

Pour une succession de handoffs séparés, chaque lien `ell` porte une phase
matricielle `alpha_ell`. En notant `s_ell=+1` lorsque le transfert suit
l’orientation de référence et `s_ell=-1` lorsqu’il la remonte, les amplitudes
se multiplient :

\[
\phi_{\rm transport}=-\sum_\ell s_\ell\alpha_\ell
\pmod{2\pi}.
\]

Le signe de cette formule est celui du saut croissant ANTLER. Il change si
l’ordre de la corde ou l’orientation de référence change.

Dans la limite Fock strictement localisée de la navette séquentielle :

1. le premier saut longitudinal `0 -> 2` franchit l’occupation au site `1`
   et fournit `-theta` ;
2. les autres sauts longitudinaux ont une corde vide ;
3. les sauts de barreau ont `n_mid=0`.

Il reste ainsi

\[
\phi_{\rm ex}=-\theta.
\]

Pour l’aller-retour apparié, le saut retour est le conjugué du saut aller :

\[
\phi_{\rm rt}=-\theta+\theta=0.
\]

La différence idéale est donc `Delta phi=-theta`. Le comptage exécutable et
indépendant de toute propagation temporelle est conservé dans
`experiments/phase4_1/run_phase4_3_exact_path_count.py`.

## Pourquoi le résultat n’est pas encore exact à profondeur finie

À profondeur de piège finie, l’état transporté n’est pas un unique état de
Fock. Soit `P(t)` le sous-espace localisé visé par un handoff et `Q=1-P` les
configurations parasites. Une élimination de Schrieffer–Wolff donne, de façon
formelle,

\[
H_{\rm eff}=PHP+
PVQ\,\frac{1}{E-QHQ}\,QVP+\cdots .
\]

Le couplage du lien devient alors

\[
J_\ell^{\rm eff}e^{i\varphi_\ell^{\rm eff}}
=J_\ell e^{i\varphi_\ell}
+\sum_{c\in Q}
\frac{V_{1c}V_{c2}}{E-E_c}+\cdots .
\]

Les termes virtuels échantillonnent des cordes Jordan–Wigner additionnelles.
Ils déplacent simultanément le module du couplage, sa phase et la population
hors sous-espace logique. Si les dénominateurs parasites restent d’ordre `D`,
la norme d’amplitude parasite est typiquement `O(J/D)` et son poids
`w_Q=O((J/D)^2)`. Cela motive une petite correction de phase. Cela **ne fixe
pas** à lui seul l’exposant observé : la chaîne SSH, les handoffs et les
annulations entre chemins virtuels peuvent modifier l’exposant apparent.

La décomposition honnête est

\[
\Delta\phi(D,T,\Delta t)
=-\theta+\delta\phi_{\rm loc}(D)
+\delta\phi_{\rm ad}(D,T)
+\delta\phi_{\rm num}(\Delta t).
\]

Le protocole Phase 4.7 isole ces termes en prenant `T(D) propto D^2` et en
comparant, à `D=6` et `D=8`, `Delta t=0.5, 0.25, 0.125`. Après contrôle de
ces deux limites, l’ansatz empirique à tester est

\[
\frac{\Delta\phi}{\theta}
=-1+\frac{a}{D^p}+O(D^{-p-1}).
\]

Le script ajuste exactement la forme à limite fixée `-1`, publie les
résidus/covariances, et ne masque pas un point qui franchirait `-1`.

## Diagnostics exigés par les campagnes

Chaque point sauvegardé contient :

- pente `Delta phi/theta`, phase impaire et fidélité logique ;
- valeurs singulières du bloc logique, leakage pire cas et défaut d’unitarité ;
- norme du mélange hors diagonale ;
- gap d’isolation minimal pendant les handoffs ;
- matrices logiques brutes et unitaires, pour audit indépendant.

Le gap rapporté est la séparation spectrale entre la branche logique de
dimension deux, suivie par recouvrement au cours du temps, et son complément.
Ce n’est pas le splitting interne du doublet ; les deux quantités sont
enregistrées séparément. Cette définition évite d’appeler « gap protecteur »
un simple splitting logique.

## Limite de revendication

Une convergence réussie établira une quantification digitale asymptotique pour
le Hamiltonien de saut corrélé d’ANTLER, dans sa convention de corde
rung-major. La famille de portes résultante demeure diagonale et abélienne.
Il faudra toujours une opération avec commutateur non nul avant de conclure à
un calcul anyonique non abélien ou universel.
