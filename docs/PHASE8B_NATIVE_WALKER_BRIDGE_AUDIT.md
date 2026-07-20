# Phase 8B - audit du pont micro vers le marcheur neutre

## Question

Le patch de controles mixtes emploie un marcheur neutre dont les transitions
sont conditionnees par `X`, `Y` ou `Z` du qubit de rail. Cette ressource est
encore declaree. Peut-on l'obtenir simplement en combinant le hopping de
barreau natif et des conversions coherentes vers deux mediateurs charge-2 ?

## Bloc exact minimal

Le bloc a charge totale ponderee trois contient :

```text
q_a, q_b       qubit de rail, exactement une particule
r0, r1         paire-reservoir
d0, d1         deux mediateurs hard-core de charge 2
```

Les seuls termes testes sont ceux qui existent deja comme ingredients locaux :

1. `q_a^dag q_b + h.c.` ;
2. `d_i^dag r0 r1 + h.c.` avec phases coherentes independantes ;
3. detunings positifs des mediateurs.

Dans le sous-secteur Mott invariant, les six etats s'ecrivent exactement
`C^2_qubit tensor C^3_mediateur`. A chaque instant,

```text
H(t) = H_q(t) tensor I + I tensor H_w(t) - c(t) I.
```

## Resultat

| audit | residu de factorisation |
|---|---:|
| Hamiltonien statique, phases complexes | `2.72e-16` |
| pire segment Floquet | `2.72e-16` |
| propagateur de trois segments, phase globale retiree | `8.66e-16` |

Les phases des canaux mediateurs changent le Hamiltonien du marcheur, et le
hopping de barreau change le Hamiltonien du qubit, mais ils ne produisent
aucun terme connecte de forme `W tensor X`, `W tensor Y` ou `W tensor Z`.
Le propagateur reste exactement un produit `U_q tensor U_w` a une phase
globale pres.

## Consequence

La programmation Floquet de phases ne suffit pas a elle seule a transformer
des ingredients separes en lien conditionnel. Le premier pont microscopique
a essayer doit casser cette factorisation par un couplage non separable, par
exemple un hopping de mediateur conditionne par la densite/parite d'un mode de
rail, ou une boucle de hopping correle partageant explicitement les modes du
qubit et du reservoir.

Ce terme ne peut pas etre suppose : il devra etre derive du hopping correle
Jordan--Wigner et des canaux charge-2, puis passer les memes controles de
selectivite, fuite, symetrie et projection sur le code.

## Claim boundary

Il s'agit d'une obstruction de factorisation pour le bloc minimal a reservoir
disjoint. Elle ne refute pas toute construction partageant des modes de
matiere ou exploitant une boucle de hopping correle. Elle ne constitue ni une
derivation native positive, ni une refutation generale du marcheur, ni une
preuve de phase topologique, twist, fusion ou braid.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_native_walker_factorization_audit.py
```
