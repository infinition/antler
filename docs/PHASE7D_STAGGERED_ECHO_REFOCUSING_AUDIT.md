# Phase 7D  -  écho de phase alternée contre le hopping de jambe

## Le problème résolu localement

Le pulse charge-2 brut ferme exactement lorsque le hopping de jambe est nul,
mais échouait nettement lorsqu'il restait actif. À `t_leg=1`, il donnait une
fuite `1.98e-2` et un résidu de parité logique `2.62e-2`. Cette note teste une
correction par contrôle dynamique dans exactement le même bloc à 472 états.

## Identité algébrique

Soit

\[
 Q=(-1)^{\sum_j j n_j+N_{\rm med}}.
\]

La phase alternée sur les rails change le signe de tout hopping de voisin
proche. Elle donnerait aussi un signe moins à toute paire sur un lien ; une
phase `pi` simultanée sur les médiateurs charge-2 compense précisément ce
signe. Pour chaque canal de paire de voisin proche,

\[
QH_{\rm med}Q=H_{\rm med},\qquad QH_{\rm leg}Q=-H_{\rm leg}.
\]

L'identité est vérifiée exactement dans la matrice 472D pour les canaux
same-rail et opposite-rail : résidu de Frobenius nul dans les deux cas.

Un sous-cycle symétrique applique

\[
e^{-iH_+\delta/4}e^{-iH_-\delta/2}e^{-iH_+\delta/4},
\qquad H_-=QH_+Q,
\]

ce qui annule le hopping de jambe au premier ordre sans annuler la conversion
de paire.

## Convergence exacte du contrôle

Pour `t_leg=1` sans crosstalk, le premier écho diminue la fuite de `1.98e-2`
à `3.29e-3`; avec 2, 4, 8 et 16 sous-cycles, elle devient respectivement
`1.27e-4`, `7.20e-6`, `4.39e-7`, `2.73e-8`.

Les fits à partir de deux sous-cycles donnent :

\[
P_{\rm leak}\propto n^{-4.0579},\qquad
\|U_n-U_0\|\propto n^{-2.0281},
\]

avec `R² > 0.99994`. C'est la loi attendue : l'erreur d'amplitude du schéma
symétrique est d'ordre deux et sa probabilité de fuite d'ordre quatre.

## Test combiné le plus sévère enregistré

Avec simultanément `t_leg=1`, couplage résiduel de canal inactif `0.01 g` et
32 sous-cycles par pulse :

- fuite hors sous-espace monomère : `5.74e-8` ;
- résidu de chacune des parités logiques : `4.29e-6` ;
- valeur singulière logique minimale : `0.9999999713` ;
- distance logique à la porte sans hopping : `1.08e-5`.

Le crosstalk crée ici le plancher principal de fuite ; augmenter les
sous-cycles ne le supprime donc pas indéfiniment.

## Décision

Le no-go « hopping de jambe laissé actif pendant le pulse brut » est réparé
dans le bloc microscopique par un protocole d'écho explicite et convergent.
C'est une avancée de **contrôle local conditionnelle** : elle qualifie la
primitive pulsée contre les deux erreurs cohérentes simulées, non une phase ou
un qubit topologique.

La nouvelle obligation physique est de dériver et réaliser les phases `pi`
alternées de rail et de médiateur, puis d'auditer leurs erreurs, les délais,
le crosstalk spatial, la dissipation et les tailles étendues. Le benchmark
Floquet de petite taille n'a toujours pas trouvé de mode de bord protégé. Il
n'y a donc aucun claim de phase topologique, code 2D, tresse, non-abélianité,
universalité ou tolérance aux fautes.

Le premier test de réalisation à durée finie est désormais archivé dans
`PHASE7D_FINITE_PHASE_KICK_AUDIT.md`. Il rejette les kicks de potentiel fort
avec Hamiltonien de pulse laissé actif sur la fenêtre `kappa=100..800` : le
protocole matériel demande donc encore une fenêtre de découplage, une séquence
composite plus économe, ou une justification de bande passante.

Résultat machine : `results/phase7/staggered_echo_refocusing_audit.json`.
