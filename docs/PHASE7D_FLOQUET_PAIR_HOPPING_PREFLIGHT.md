# Phase 7D-F  -  préflight Floquet de pair-hopping

## Motivation

La littérature récente montre qu'une échelle à nombre conservé peut obtenir du
pair-hopping par séquence périodique de pulsations, plutôt que par un terme
statique très faible d'ordre élevé. Le présent contrôle ne reproduit pas ce
protocole expérimental ; il vérifie la brique algébrique minimale que devrait
réaliser une extension dynamique d'ANTLER.

Référence d'inspiration : Defossez, Vanderstraeten, Peralta Gavensky et
Goldman, *Dynamic Realization of Majorana Zero Modes in a Particle-Conserving
Ladder*, arXiv:2412.14886 (version révisée 2025),
https://arxiv.org/abs/2412.14886.

## Ressource explicitement ajoutée

Sur deux barreaux, les segments périodiques sont

\[
H_x=UC+V X_0X_1,\qquad H_y=UC-V Y_0Y_1,
\]

avec un signe d'interaction modulable. Les rotations de rail compilent les
axes `X` et `Y`; la modulation de signe est une nouvelle ressource de contrôle,
non présente dans la grammaire statique charge-2 gelée.

Les deux segments commutent dans ce contrôle, donc sur une période `T=2 tau`,

\[
H_F=UC+\frac{V}{2}(X_0X_1-Y_0Y_1)
 = UC+V(S_0^+S_1^+ + S_0^-S_1^-).
\]

## Résultat exact sur le bloc à deux particules

Avec `U=20`, `V=1` et `tau=pi/40` :

- résidu de compilation d'une période : `2.40e-18` ;
- résidus de parité de rail stroboscopiques : `0` ;
- fuite hors sous-espace monomère par période : `0` ;
- coefficient de pair-hopping : `1.0` ;
- transfert `|bb>` vers `|aa>` après dix périodes : `1.0`.

Il s'agit d'un contrôle positif par construction : il démontre que le banc
d'essai reconnaît une primitive dynamique de pair-hopping propre. Il ne
démontre pas que cette primitive est déjà disponible dans le Hamiltonien natif
ANTLER.

## Suite et mise à jour du statut

Le premier pont charge-2 positif a depuis été vérifié sur un bloc exact : des
canaux même-rail et opposé-rail à détuning positif, sélectionnés dans le
temps, fournissent les deux signes de `ZZ` après excursion Rabi fermée. Le
compilateur à quatre barreaux montre également une erreur de Trotter
quadratique. Voir `PHASE7D_MEDIATOR_SIGNED_ZZ_PREFLIGHT.md` et
`PHASE7D_MULTILINK_COMPILER_AUDIT.md`.

Ces deux contrôles ne franchissent pas la porte de phase : le bridge minimal
à pair-hopping a été rejeté pour absence de localisation de bord, et un
contrôle séparé du Hamiltonien Floquet complet n'a pas trouvé de candidat sur
sa petite fenêtre `L=4 -> 6`. Voir
`PHASE7D_FULL_FLOQUET_LADDER_AUDIT.md`.

La route n'est donc ni fermée ni promue : elle demande une dérivation de
ladder à pulses complète, avec crosstalk et fuites, puis une étude de taille
thermodynamique. Une géométrie 2D/T-junction et un encodage non abélien restent
des problèmes distincts.
