# Phase 8B  -  audit du bloc de Fock sous contrainte de Gauss

## Verdict

Le marcheur Lambda réparé s'embedde correctement dans un premier bloc de
Fock à charge fixe. C'est un résultat plus fort que l'audit de signes : le
Hamiltonien complet du bloc commute exactement avec le générateur de Gauss,
et le tunnel de barreau nu ne peut pas agir dans le secteur physique.

La primitive de saut conditionné du marcheur est cependant encore insérée à
la main. Ce n'est pas une dérivation native ANTLER ni une phase multi-blocs.

## Bloc audité

- deux barreaux (`b=2`), quatre modes fermioniques, secteur `N=2` ;
- deux qubits de jauge de frontière ;
- marcheur neutre `|vac>, |mu0>, |mu1>, |mu2>`, avec au plus une occupation ;
- générateur
  \[
  G_B=X_L(-1)^{N_{a,B}}X_R;
  \]
- saut de marcheur inséré :
  \[
  \mu_j^\dagger\mu_{j-1}(1-2n_{a,j})+h.c.
  \]
  pour les deux barreaux.

La dimension totale est `96`.

## Résultats exacts

| Test | Résultat |
|---|---:|
| `||[H,G_B]||_F` | 0 |
| `||{c_a†c_b+h.c.,G_B}||_F` | 0 |
| `||P_+(c_a†c_b+h.c.)P_+||` | 0 |
| norme de changement de secteur `P_- O P_+` | 1 |
| pente profonde du gap de secteurs | `3.973` |
| commutateur avec hopping intra-bloc `t_leg=0.3` | 0 |

À `coupling/Delta=0.025`, le gap entre `G_B=+1` et `G_B=-1` est
`1.5567e-5`; il est petit car il est d'ordre quatre. À ce stade, c'est une
échelle de contrainte démontrée, pas encore une fenêtre expérimentale viable.

## Ce que le résultat établit

Dans le bloc explicitement défini :

1. le mécanisme de boucle génère bien une pénalité de Gauss ;
2. le hopping intra-bloc utile respecte cette contrainte ;
3. le tunnel de rail parasite est une excitation de jauge, pas une opération
   logique interne au secteur `G_B=+1`.

## Prochaine porte obligatoire

Le prochain calcul doit mettre deux blocs côte à côte et ajouter le transport
inter-blocs par paires. Il doit vérifier simultanément :

1. `[H,G_0]=[H,G_1]=0` à précision machine ;
2. maintien de gaps de Gauss non nuls ;
3. absence de saut simple inter-blocs ;
4. survie d'une dynamique ou phase appariée.

En cas d'échec de cette quatrième porte, le marcheur ne donne pas une route
vers une jonction T utilisable malgré le succès du bloc local.

## Claim boundary

Pas de primitive native dérivée, pas de transport entre blocs, pas de phase
thermodynamique, pas d'indistinguabilité locale, pas de fusion ni de tresse.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_fock_gauss_block_audit.py
```
