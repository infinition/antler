# Phase 8B  -  préflight à plusieurs paires

## Verdict

Le transport cohérent d'une paire ne se prolonge pas automatiquement en phase
protégée. Dans le modèle effectif le plus favorable, déjà projeté sur des
paires hard-core, les deux prolongements locaux simples échouent de façons
complémentaires :

- le fluide de paires mobile est compressible et sans doublet logique ;
- la répulsion forte forme un doublet CDW localement lisible.

Ce résultat ferme seulement la continuation *pair hopping + répulsion de
densité* de l'architecture Phase 8B. Il ne réfute pas une phase appariée plus
riche avec une interaction supplémentaire effectivement dérivée.

## Modèle de contrôle

On considère une chaîne périodique de `B` blocs `b=2`, projetée à une
occupation locale de paire hard-core `p_j=0,1`, à demi-remplissage :

\[
H=-J\sum_j(p_j^\dagger p_{j+1}+h.c.)+V\sum_j n_jn_{j+1}.
\]

Chaque saut change `N_a` de deux unités dans chacun des deux blocs concernés.
Il est donc compatible avec la loi de Gauss grossière déjà auditée. Cette
projection et le terme `V` sont toutefois insérés : ils ne sont pas une
dérivation microscopique du marcheur Lambda.

## Résultats

| Cas | Observation `B=4 → 12` | Lecture correcte |
|---|---|---|
| `V=0` | gap neutre `2.828 → 1.035`, courbure d'ajout de paire `1.657 → 0.527`; puissances `-0.916`, `-1.043` | fluide compressible, fondamental unique |
| `V=8J` | split du doublet `0.899 → 0.0217`, gap au-dessus `8.00 → 5.11` | doublet isolé en taille finie |
| `V=8J`, observateur local `n_j` | `||P n_j P - tr(Pn_jP)I/2|| = 0.4683` à `B=12` | le doublet CDW est localement lisible, donc pas un code topologique |

Le contrôle répulsif est important : un petit splitting et un grand gap vers le
bulk ne suffisent pas. La projection d'un opérateur local sur le doublet doit
être scalaire ; ici elle reste presque maximale.

## Conséquence architecturale

Le prochain filtre n'est ni une jonction T ni une matrice de tresse. Il faut
d'abord dériver, à partir d'une ressource microscopique déclarée, une
interaction à plusieurs blocs dont le sous-espace fondamental est à la fois
gappé et localement indistinguable. Sans cette étape, une T-junction ne ferait
que déplacer des paires ou des états CDW lisibles.

## Claim boundary

Ce préflight est un modèle effectif de paires, avec projection locale et
interaction `V` imposées. Il ne démontre aucune phase ANTLER native, aucun gap
thermodynamique protégé, aucune indistinguabilité locale, aucun défaut, aucune
fusion, aucune jonction T et aucune tresse non abélienne.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_pair_chain_manybody_preflight.py
```
