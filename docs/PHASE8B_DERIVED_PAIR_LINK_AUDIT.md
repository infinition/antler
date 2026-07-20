# Phase 8B  -  lien de paire dérivé et intégration Gauss à deux blocs

## Verdict

Le terme de transport de paire qui était inséré dans l'audit initial à deux
blocs possède maintenant une dérivation locale explicite : un médiateur
charge-2 positivement détuné, partagé entre deux canaux de paires voisins,
génère le lien au deuxième ordre SW. Son intégration avec les deux marcheurs
Lambda conserve exactement les deux lois de Gauss grossières.

Il s'agit d'une avancée de mécanisme microscopique local. Ce n'est toujours ni
une phase à plusieurs paires, ni un code topologique, ni une tresse.

## Dérivation locale

Sur deux blocs `b=2`, avec les quatre modes de rail `a`, on introduit un
médiateur hard-core `m` de charge 2 :

\[
H=\Delta n_m + g m^\dagger(P_L+P_R)+h.c.,\qquad
P_L=a_0a_1,\quad P_R=a_2a_3.
\]

Dans la branche médiateur vide, à charge totale fixée, la réduction donne

\[
H_{\rm eff}^{(2)}=-\frac{g^2}{\Delta}
 (P_L^\dagger P_L+P_R^\dagger P_R+P_L^\dagger P_R+P_R^\dagger P_L).
\]

Les deux termes diagonaux ont le même coefficient que le transfert. Ils sont
scalaires pour une unique paire mobile sur ce lien, mais doivent être gardés
dans toute étude à plusieurs paires.

## Contrôle exact local (29 états)

| `g/Delta` | erreur relative du coefficient SW | transfert maximal | population virtuelle max de `m` |
|---:|---:|---:|---:|
| 0.10 | 1.92% | 0.99902 | 0.03704 |
| 0.05 | 0.495% | 0.99994 | 0.00980 |
| 0.025 | 0.125% | 0.999996 | 0.00249 |

L'erreur de coefficient se comporte comme `(g/Delta)^1.982`. Les deux parités
de bloc commutent exactement avec le Hamiltonien; le saut monoparticule de
contrôle a une projection nulle sur le secteur physique.

## Intégration exacte avec les marcheurs (3 712 états)

Le terme effectif imposé `P_L^dagger P_R+h.c.` de l'ancien contrôle a été
retiré et remplacé par le médiateur explicite. Les deux marcheurs neutres
Lambda et les trois qubits de frontière restent présents.

| `g/Delta` | population paire cible max | population lien `m` max | population marcheurs max | fuite Gauss max |
|---:|---:|---:|---:|---:|
| 0.10 | 0.96164 | 0.03936 | 0.14392 | `1.04e-12` |
| 0.05 | 0.99000 | 0.00998 | 0.03747 | `5.15e-13` |
| 0.025 | 0.99750 | 0.00250 | 0.00966 | `4.34e-11` |

Pour chaque point, `[H,G_0]=[H,G_1]=0` et le lien dérivé commute aussi avec
les deux générateurs. Le saut `a` monoparticule de frontière reste strictement
hors secteur physique.

## Ressources et limites

La brique emploie un médiateur de charge 2, mais avec une topologie de canal de
paires à quatre jambes (`m` couple aux deux paires locales). Cette topologie,
ainsi que les marcheurs Lambda conditionnés par densité, sont des ressources
nouvelles déclarées : elles ne sont pas dérivées de la grammaire gelée à deux
médiateurs par lien.

Le préflight à plusieurs paires reste négatif pour la continuation minimale
pair-hopping + répulsion : fluide compressible ou doublet CDW lisible. Il faut
donc encore dériver une interaction multi-blocs qui donne simultanément gap et
indistinguabilité locale avant toute jonction T.

## Claim boundary

Cette preuve couvre une brique locale et une navette à deux blocs, une paire,
à charge totale fixe. Elle n'établit aucun gap thermodynamique, aucune phase à
plusieurs paires, aucune indistinguabilité locale, aucun défaut/fusion, aucune
jonction T et aucune tresse non abélienne.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_shared_pair_link_sw_audit.py
python experiments/phase7/run_phase8b_two_block_derived_pair_link.py
```
