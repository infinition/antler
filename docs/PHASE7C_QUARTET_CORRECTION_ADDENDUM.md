# Phase 7C  -  addendum de correction : projection par monômes du quartet croisé

## Portée et reproductibilité

Addendum de correction aux audits Phase 7C du médiateur croisé
(`PHASE7C_CROSSED_MEDIATOR_AUDIT.md`, `PHASE7C_OVERLAP_AUDIT.md`,
`PHASE7C_ALL_FLIP_OVERLAP_ALGEBRA.md`). Calcul nouveau :
`experiments/phase7/run_phase7c_quartet_monomial_projection.py` →
`results/phase7/quartet_monomial_projection.json` (exécuté, sortie vérifiée
contre les coefficients archivés). La matrice de statut et les manifests ne
sont **pas** modifiés par cet addendum : leur mise à jour attend la
validation indépendante du nouveau JSON.

Toutes les affirmations d'identité d'opérateur ci-dessous sont limitées au
**Hamiltonien croisé de référence**  -  deux médiateurs charge-2
`d₁†(a₀a₁+b₂b₃)`, `d₂†(b₀b₁+a₂a₃)`, `U=20`, `Δ=10`, **sans** hopping de rail,
**sans** contre-termes  -  dans le **sous-espace monomère** downfoldé
exactement. Elles ne sont pas extrapolées à toute grammaire charge-2
imaginable.

## 1. Projection explicite : le quartet ne génère pas `XXXX`, il génère l'all-flip

`XXXX = Π_j(S_j⁺+S_j⁻)` contient seize monômes de flip ; sur les seize états
monomères, il connecte les **huit** paires antipodales (états différant sur
les quatre barreaux) avec des éléments égaux. Projection exacte du
Hamiltonien effectif traceless à `g=0.5` :

| observable | valeur |
| --- | ---: |
| élément `⟨bbbb\|H_eff\|aaaa⟩` (all-flip) | `−1.6244e-6` |
| plus grand des 7 autres éléments antipodaux | `2.31e-15` |
| norme hors-diagonale à distance de flip 1 / 2 / 3 | `1.4e-14` / `1.6e-14` / `1.7e-14` |
| norme hors-diagonale à distance de flip 4 | `2.2972e-6` (= √2·\|c_F\|) |
| loi de puissance de \|c_F\| sur `g=0.3..1.3` | `g^3.9984` |

**Dans ce bloc de référence, le contenu hors-diagonal monomère→monomère est
`c_F(|aaaa⟩⟨bbbb|+h.c.)` à la précision machine** : deux monômes sur seize,
pas `XXXX`. Raison mécanique, propre à cette instance de grammaire : chaque
canal lie des paires d'arêtes fixées à rail fixé ; depuis toute configuration
mixte (`|abab⟩`, `|aabb⟩`, …) aucune conversion n'est applicable ou aucun
chemin ne referme vers un monomère différent. Avec des hoppings de rail ou
d'autres termes, cette identité ne tient plus telle quelle  -  et ces termes
réintroduisent les parasites d'ordre deux du catalogue.

Identité de Pauli vérifiée numériquement (résidu de motif de signes
`8.0e-16`) :

```
|aaaa⟩⟨bbbb| + h.c. = (1/8) [ XXXX − (XXYY+XYXY+XYYX+YXXY+YXYX+YYXX) + YYYY ]
```

d'où `c_XXXX = c_F/8 = −1.6244e-6/8 = −2.030e-7`  -  **exactement** le
coefficient archivé dans `crossed_charge2_cancellation_audit.json`
(`−2.03e-7`) : rétro-validation croisée des deux calculs. La part des huit
mots de la famille all-flip dans tout le contenu quatre-corps hors-diagonal
est `1.0000` (aucun autre mot quatre-corps X/Y).

## 2. No-go reformulé : sélectivité, gaps et parasites, sans amplitude absolue

**Conventions d'abord.** Le JSON de l'optimiseur définit
`target_alignment = c_cible² / Σc²`  -  convention **au carré**, sur tous les
mots non-identité (diagonaux compris). En amplitude, le recouvrement
`F_p`↔`XXXX` vaut `1/√8 = 0.35355` ; au carré, `1/8 = 0.125`. Le `0.1135`
mesuré après contre-termes est en convention au carré.

Trois clauses structurelles, invariantes d'échelle :

1. **Plafond d'identité d'opérateur.** Pour le bloc croisé de référence, tout
   le contenu quatre-corps hors-diagonal est la famille de `F_p`
   (part `1.0000`), et les contre-termes admis (`rail_biases` → `Z`,
   `zz_couplings` → `ZZ`) sont diagonaux : ils ne peuvent ni créer les
   quatorze monômes manquants ni toucher les mots hors-diagonaux. Donc
   l'alignement au carré sur `XXXX` est borné par `1/8 = 0.125` **quels que
   soient** `g, Δ, U` et la calibration. Le `0.1135` mesuré en est 90.8 % ;
   le résidu vient des mots diagonaux d'ordre supérieur (`ZZZ`-type,
   `ZZZZ`) hors de portée de la grammaire de contre-termes. La porte de
   stage 2 exige `≥ 0.80` : **inatteignable structurellement**, sans aucun
   seuil d'amplitude.
2. **Hiérarchie d'ordres.** `|c_F| ∝ g^3.998` contre parasites `∝ g^1.999`
   (archivé) : la sélectivité varie en `(g/Δ)²`. La fenêtre SW enregistrée
   (`g/Δ ≤ 0.15`) est elle-même imposée par les portes dures de capture
   (`≥ 0.90`) et de gap : on ne peut pas y échapper en augmentant `g`. La
   sélectivité est donc bornée dans la fenêtre contrôlée et s'annule dans la
   limite contrôlée  -  énoncé indépendant de tout `|c| > 1e-7`.
3. **Portée des contre-termes.** La calibration ne peut que nettoyer une
   partie du diagonal (c'est ce que `0.1135` mesure) ; elle ne change ni
   l'identité (1) ni la hiérarchie (2).

Les critères d'amplitude absolue des audits antérieurs étaient suffisants
pour rejeter, mais pas nécessaires ; cette reformulation les remplace.

## 3. Plaquettes voisines : commutateurs, pavage, indistinguabilité

**Commutateur analytique.** Pour `F_p = Π_pS⁺ + Π_pS⁻`,
`F_q = Π_qS⁺ + Π_qS⁻`, recouvrement `O = p∩q ≠ ∅`, `p ≠ q` : sur `O`,
`S⁺S⁺ = 0`, `S⁺S⁻ = P^a`, `S⁻S⁺ = P^b`, donc avec
`W = (Π_{p∖O}S⁺)(Π_{q∖O}S⁻)` :

```
[F_p, F_q] = W ( Π_O P^a − Π_O P^b ) + W† ( Π_O P^b − Π_O P^a ) ,
```

de norme spectrale `1` (isométries partielles à supports orthogonaux),
indépendamment de `|O| = 1, 2, 3`. Cela **dérive** la table exacte archivée
`{0:0, 1:1, 2:1, 3:1, 4:0}` (`all_flip_overlap_algebra.json`) : tout pavage
à barreaux partagés est non commutant au niveau des termes, à hauteur de la
normalisation du terme lui-même.

**Pavage disjoint  -  fait nouveau.** Le spectre exact d'un terme idéal isolé
`−J F_p` sur les seize monomères : fondamental **unique**
`(|aaaa⟩+|bbbb⟩)/√2` à `−J`, quatorze états mixtes **sombres** à `0`,
partenaire cat à `+J` (`ground_degeneracy = 1`, `dark_state_count = 14`).
Un pavage disjoint du terme réellement engendré n'a donc **aucun espace de
code** (`k = 0`), avant même toute question d'indistinguabilité.

**Lisibilité  -  analyse conditionnelle du modèle idéal.** *Si* un doublet cat
dégénéré `span{|aaaa⟩, |bbbb⟩}` était fabriqué par des ressources
supplémentaires, une sonde à **un seul barreau** y serait déjà logique :
`‖P Z₀ P − scalaire‖ = 1.0` (idem `Z₁`), `X₀, Y₀ → 0`  -  la même classe de
lisibilité locale que le benchmark Ising archivé. Ceci est une analyse du
doublet idéalisé/fabriqué, **pas** du Hamiltonien effectif complet : son
isolement spectral d'un tel doublet n'est pas démontré ici (les énergies
diagonales effectives des monomères mixtes n'ont pas été analysées pour
l'isolement dans cet addendum).

## 4. Charge-4 : classement, sans conclusion anticipée

Un médiateur de charge 4 est un **nouvel opérateur de conversion à quatre
fermions**, hors de la chaîne d'audit charge-2 existante (6D). C'est une
**nouvelle ressource matérielle non auditée** : le contenu d'opérateur
qu'elle engendrerait dépend de sa structure de couplage microscopique
exacte, et n'est pas préjugé ici  -  ni en sa faveur, ni contre elle. Avant
tout usage, elle exigerait sa propre chaîne : signatures de parité,
annulation mono-particule, mesures d'ordre, algèbre de recouvrement, portes
7E. Elle ne doit pas être présentée comme une solution ANTLER de la branche
statique.

## 5. Verdict

**Branche statique charge-2 auditée : fermée structurellement.**

La grammaire auditée (catalogue lien-par-lien, quartet croisé, contre-termes
diagonaux) échoue par identité d'opérateur (all-flip, plafond d'alignement
au carré `1/8`), par hiérarchie d'ordres dans la fenêtre SW enregistrée, et
par algèbre de recouvrement/pavage (`k=0` en disjoint, non-commutation en
recouvrement). Cette fermeture ne dit **pas** « tous les Hamiltoniens
charge-2 possibles : impossibles » : une grammaire différente (hoppings,
autres topologies de canaux, ressources nouvelles) exigerait son propre
audit, à sa charge, via les portes 7E. Le front actif reste la Phase 8
dynamique.
