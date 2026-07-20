# Phase 8B - pont microscopique a matiere partagee vers un lien conditionnel

## Question

Le bloc minimal a reservoir disjoint est exactement separable : les phases de
conversion et le hopping de barreau ne peuvent y produire qu'un propagateur
`U_rail tensor U_walker`. Cette obstruction n'interdit pas un bloc dans lequel
le qubit de rail, une paire voisine et les mediateurs de charge deux partagent
des etats virtuels. Le present audit teste ce mecanisme minimal non separable.

## Bloc et controles enregistres

Le bloc conserve exactement la charge ponderee trois et contient

```text
q_a, q_b       : un fermion logique de rail,
r_a, r_b       : deux modes de reservoir de paire,
d0, d1         : deux mediateurs hard-core de charge deux.
```

Le Hamiltonien microscopique contient le hopping de barreau, les deux hoppings
correles de jambe a `theta=pi`, les conversions coherentes
`d_i^dag r_a r_b + h.c.`, un cout Mott qui fixe un fermion `q`, et le detuning
positif de la paire reservoir. Le cadre bas est les quatre etats
`|q_a/q_b> tensor |d0/d1>`.

Deux segments sont compares :

```text
A : J_perp=+J, g0=+g, g1=+g
B : J_perp=-J, g0=+g, g1=-g.
```

Le premier changement est un signe Peierls de barreau, le second une phase
relative `pi` entre canaux de paires. Ils sont des controles dynamiques
enregistres, et non des termes logiques ajoutes au Hamiltonien effectif.

## Resultat exact de downfolding

La moyenne des deux Hamiltoniens effectifs de Schrieffer--Wolff donne

```text
H_echo = c_I I + c_XX X_q tensor X_(d0,d1).
```

Les termes seuls `X_q` et `X_walker` changent de signe entre les segments et
s'annulent. Le terme conditionnel est pair sous les deux changements et reste.

| g/Delta | c_XX | capture minimale | gap bas-haut minimal |
|---:|---:|---:|---:|
| 0.20 | -2.0035e-3 | 0.98060 | 19.9590 |
| 0.10 | -5.0087e-4 | 0.99453 | 19.6677 |
| 0.05 | -1.2522e-4 | 0.99821 | 19.5935 |
| 0.025 | -3.1304e-5 | 0.99914 | 19.5748 |

Le fit profond donne `|c_XX| ~ (g/Delta)^2.0000`. Dans toutes les lignes,
le plus grand coefficient Pauli non scalaire hors cible est nul a la precision
serialisee et le residu d'hermiticite est inferieur a `9e-19`.

## Axes statiques et correction du test de pulse

La phase Peierls de barreau donne, dans le Hamiltonien effectif statique, la
famille continue :

```text
(-cos(phi) X_rail + sin(phi) Y_rail) tensor X_walker.
```

La rotation de cet axe est verifiee a la precision numerique et le coefficient
`YX` a `phi=pi/2` conserve le fit profond d'ordre `2.0000`.

Le troisieme axe ne requiert pas un terme logique insere : remplacer le
hopping de barreau par un biais de potentiel physique entre `q_a` et `q_b`,
dont le signe est echoe entre les segments, isole
`Z_rail tensor X_walker`. Son coefficient est aussi d'ordre `2.0000`.
Le registre local contient donc les trois liens derives :

```text
X_rail tensor X_walker : hopping de barreau reel
Y_rail tensor X_walker : hopping de barreau a phase Peierls pi/2
Z_rail tensor X_walker : biais de potentiel de rail signe
```

Pour les trois, aucun Pauli non scalaire hors cible ne depasse la precision
serialisee, et la capture a `g/Delta=0.025` est au moins `0.99913`.

La premiere lecture des pulses Rabi etait trop optimiste : une petite distance
absolue a la cible SW ne prouve rien lorsque la cible elle-meme est petite.
Au point `g/Delta=0.025`, pour `X` et `Y`, le signal SW non scalaire vaut
`3.93e-5`, le signal physique seulement `2.74e-7`, et l'erreur relative vaut
`1.0063`. Pour `Z`, les valeurs sont `7.86e-5`, `5.05e-7`, `1.0063`.
Le retour Rabi est donc presque scalaire : il ferme les occupations virtuelles
mais annule aussi la phase conditionnelle que le downfolding statique predit.

Verdict de pulse : **la famille de timings Rabi entiers est rejetee comme
porte**. Son petit leakage et son bracket de timing ne sont conserves que
comme controles de retour de population; ils ne qualifient aucune operation
logique. Une sequence off-resonante, habillee ou adiabatique doit etre derivee
et auditee avec l'erreur relative au signal, avant toute nouvelle promotion.

## Ce qui est franchi

Le verrou de factorisation du bloc a reservoir disjoint est donc leve dans une
extension explicite a matiere partagee. C'est un pont microscopique local vers
un lien conditionnel de marcheur, et non une matrice de braid inseree a la
main. L'echelle est favorablement d'ordre deux, contrairement aux stabilisateurs
statiques a quatre corps de la branche 7C.

## Ce qui reste obligatoire

Seule la compilation statique/downfolded du lien est etablie. Il faut ensuite :

1. compiler un marcheur neutre a quatre etats a partir de ces liens `X`, `Y`,
   `Z` et verifier ses crosstalks simultanes ;
2. integrer les liens dans le patch de controles mixtes et re-auditer les
   commutateurs, le gap et les termes parasites simultanes ;
3. construire des defauts mobiles et un espace de fusion ;
4. seulement alors calculer deux holonomies derivees, leur commutateur et la
   relation de Yang--Baxter.

## Claim boundary

Ce resultat est un compilateur local **statique/downfolded** dans une extension
declaree de la grammaire ANTLER (modes reservoir et controles coherents). Il
ne demontre pas une porte pulse, un marcheur complet, un patch natif, une phase
topologique, des twists, de la fusion, une tresse non abelienne, l'universalite
ou la tolerance aux fautes.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_shared_matter_conditional_link_sw_audit.py
python experiments/phase7/run_phase8b_shared_matter_pulse_closure_audit.py
python experiments/phase7/run_phase8b_shared_matter_pulse_timing_audit.py
python experiments/phase7/run_phase8b_shared_matter_phase_link_audit.py
python experiments/phase7/run_phase8b_shared_matter_xyz_link_audit.py
```
