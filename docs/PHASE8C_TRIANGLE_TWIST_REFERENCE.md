# Phase 8C-T3b  -  Référence locale de twist vers un bord

## Verdict

**PASS comme référence finie de jonction mixte ; pas de fusion ni de tresse.** Un patch triangulaire à sept qubits neutres, trois plaquettes et trois boucles de bord fournit un contrôle local explicite de coupure `e<->m` qui se termine sur un bord. Le check central non-CSS est `YZXIXII`.

Le calcul est `experiments/phase8c/run_phase8c_triangle_twist_reference.py`; le JSON est `results/phase8c/triangle_twist_reference.json`.

## Contrat et géométrie

Les sept qubits sont étiquetés par les sommets du plus petit triangle
`r=s=t=2` :

```text
(0,0,0), (1,0,0), (0,1,0), (0,0,1),
(1,1,0), (1,0,1), (0,1,1).
```

Les checks inscrits sont deux plaquettes CSS (`XXIXIXI`, `ZIZZIIZ`), la plaquette mixte centrale `YZXIXII`, puis trois boucles de bord. Cette structure est une référence de stabilisateurs à liens neutres déclarés ; elle s'inspire du mécanisme de défaut non-CSS utilisé dans les codes à twists, mais n'est pas une dérivation du Hamiltonien ANTLER. Les twists des codes de surface exigent précisément une irrégularité non-CSS, et non un simple flux. [Bombín (2010)](https://arxiv.org/abs/1004.1838), [Yoder & Kim (2016)](https://arxiv.org/abs/1612.04795)

## Portes exactes

| Test | Résultat |
|---|---:|
| Checks commutants | oui (`0` anticommutation) |
| Rang stabilisateur | `6` |
| Qubits encodés / GSD | `1 / 2` |
| Distance | `3` |
| Gap de syndrome de `-J sum S_a` | `2J` |
| Paulis de poids 1--2 sondés | `210` |
| Actions logiques non scalaires sous distance | `0` |

## Témoin de conversion de corde

Les deux représentants logiques de bord sont

```text
Z-side = ZZZIIII
X-side = IIXXXII
```

et anticommuttent. Multiplier une corde par le check mixte ne change pas son
action dans le code mais la déforme à travers la jonction :

```text
ZZZIIII × YZXIXII = XIYIXII
IIXXXII × YZXIXII = YZIXIII
```

Le témoin est important : un représentant de type pur devient une continuation
mixte contenant le type Pauli opposé, avec le `Y` localisé à la jonction. C'est
la signature algébrique locale recherchée pour une coupure `e<->m`; elle va au-delà d'un signe `pi` ou de checks mixtes non reliés à une déformation de corde.

## Ce qui manque encore

Le patch contient une jonction qui se termine au bord, pas deux défauts
séparés dans une même géométrie 2D. Il ne possède donc pas d'espace de fusion
à mesurer. La prochaine porte T4 est : deux extrémités de twist explicites,
une distance géométrique variable, dimension du sous-espace de fusion,
indistinguabilité locale et gap. Les holonomies et le commutateur de deux
échanges ne seront définis qu'après ce passage.

## Claim boundary

Affirmé : un code de stabilisateurs local, non-CSS et de distance trois réalise
une déformation de corde X/Z-vers-mixte à travers une jonction vers un bord.

Non affirmé : un mur thermodynamique, une paire de twists, une fusion Ising,
une tresse non abélienne, une universalité, une dérivation ANTLER, ou une
réalisation matérielle.
