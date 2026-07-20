# Phase 5  -  positionnement bibliographique avant extension microscopique

Ce document fixe les comparaisons minimales avant de dériver de nouveaux
termes. Il ne transforme pas la similarité de forme entre Hamiltoniens en
équivalence topologique.

## 1. Hopping corrélé et anyon-Hubbard

Keilmann *et al.* ont montré qu’un modèle anyon-Hubbard 1D peut être réécrit
comme un Bose-Hubbard à hopping dépendant de l’occupation, avec angle
statistique contrôlable. C’est le précédent direct pour la phase de Peierls
conditionnelle du ladder ANTLER. Il établit la pertinence du mécanisme de
hopping corrélé ; il n’établit ni une dégénérescence de fusion, ni des tresses
non abéliennes dans une échelle scalaire 1D.

Référence primaire : [Keilmann *et al.*, *Nature Communications* 2, 361
(2011)](https://www.nature.com/articles/ncomms1353).

## 2. Majoranas avec conservation du nombre

La conservation du nombre n’interdit pas à elle seule des modes de bord de
type Majorana. Fidkowski, Lutchyn, Nayak et Fisher donnent une construction où
l’ordre supraconducteur de longue portée n’est pas requis ; Iemini *et al.*
construisent un modèle exact à deux fils, conservant le nombre, avec excitations
de bord Majorana-like. Le second résultat souligne aussi le piège essentiel :
un gap de particule unique peut coexister avec un Hamiltonien globalement
gapless. Il faut donc auditer séparément le gap du secteur logique, le mode de
charge et l’indistinguishabilité locale.

Références primaires :

- [Fidkowski *et al.*, *Physical Review B* 84, 195436
  (2011)](https://arxiv.org/abs/1106.2598) ;
- [Iemini *et al.*, *Physical Review Letters* 115, 156402
  (2015)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.156402) ;
- [Wang *et al.*, modèles à nombre conservé ayant des états fondamentaux BCS
  exacts (2017)](https://arxiv.org/abs/1703.01249).

**Conséquence pour ANTLER :** la piste crédible n’est pas d’ajouter une porte
SU(2) imposée au code L=14 actuel. C’est de chercher une extension locale à
conservation du nombre  -  par exemple un pair-hopping entre jambes  -  puis de
démontrer un sous-espace de défaut/fusion, son gap pertinent et sa localité.

## 3. Parafermions et géométries d’échelle

Les parafermions demandent des structures plus contraignantes que le simple
hopping scalaire à corde Jordan–Wigner. Calzona *et al.* donnent un Hamiltonien
fermionique 1D local et interactif exactement soluble avec modes de bord
Z4, protégé par la structure de symétrie correspondante. Santos et Béri relient
des isolants de Mott quasi-1D à trois jambes à des modèles d’horloge chiraux
Z3 et à des strong zero modes. Ces travaux rendent la troisième jambe
scientifiquement motivée, mais ils imposent un régime Mott/chiral et des
interactions spécifiques absents du Hamiltonien ANTLER gelé.

Références primaires :

- [Calzona *et al.*, *Physical Review B* 98, 201110(R)
  (2018)](https://arxiv.org/abs/1802.06061) ;
- [Santos & Béri, *Physical Review Letters* 125, 207201
  (2020)](https://arxiv.org/abs/2005.12288).

## Décision de conception

La prochaine étude n’est pas une nouvelle simulation de porte sur le spinor
L=14 : son audit montre déjà une distinguabilité locale d’environ `1.40`.
Elle doit être une **comparaison de Hamiltoniens** limitée et falsifiable :

1. choisir une extension locale précise (pair-hopping conservant le nombre,
   ou échelle à trois jambes en régime Mott chiral) ;
2. montrer par ED/DMRG qu’elle possède un sous-espace de défaut/fusion
   localement indistinguable ;
3. seulement alors construire deux contrôles d’échange et auditer
   commutateur, relation de braid, gap, leakage et bruit.

Le modèle BdG de jonction en T reste un benchmark architectural utile ; il ne
constitue pas une dérivation de l’extension ANTLER à nombre conservé.
