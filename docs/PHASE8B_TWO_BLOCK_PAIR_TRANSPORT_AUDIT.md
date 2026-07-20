# Phase 8B  -  audit à deux blocs et transfert de paires

## Verdict

Le premier transport entre blocs passe la contrainte de Gauss : dans le
modèle à primitive de marcheur insérée, un transfert de deux particules du
rail `a` conserve exactement les deux générateurs de Gauss. Le saut simple
sur la même frontière les viole tous deux et ne possède aucune projection dans
le secteur physique.

Ce résultat valide l'algèbre de transport à deux blocs. Il ne démontre pas
encore une phase appariée ni une réalisation native.

## Modèle

- deux blocs de deux barreaux, huit modes de matière à `N=2` ;
- trois qubits de frontière, dans leur base propre de `X` ;
- un marcheur Lambda par bloc ;
- transfert admis :
  \[
  a_{0}^{\dagger}a_{1}^{\dagger}a_{3}a_{2}+h.c.,
  \]
  qui change chaque `N_{a,B}` de `0` ou `+/-2` ;
- contrôle négatif : saut simple `a_1^\dagger a_2+h.c.`.

La dimension exacte à charge fixe est `3584`.

## Résultats

| Test | Résultat |
|---|---:|
| `||[H,G_0]||_F`, `||[H,G_1]||_F` | 0 |
| commutateurs du transfert de paire avec `G_0,G_1` | 0 |
| norme du transfert de paire dans `G_0=G_1=+1` | `0.4` |
| anticommutateurs du saut simple avec `G_0,G_1` | 0 |
| projection physique du saut simple | 0 |
| pente profonde du plus petit gap de Gauss | `3.9733` |

Le gap est petit  -  `1.54e-5` au point le plus profond  -  car le gadget est
d'ordre quatre. Aucune fenêtre de bruit ou de matériel n'est donc revendiquée.

## Signification

Le projet possède maintenant un chemin cohérent au niveau des blocs :

1. contrainte de Gauss locale ;
2. bloc de Fock à charge fixe ;
3. transfert de paires inter-blocs qui conserve la contrainte ;
4. rejet exact du saut simple parasite.

## Portes suivantes

Il reste à établir :

1. la dérivation du marcheur conditionné et du transfert de paire à partir
   d'une ressource matérielle précise ;
2. une phase appariée sur plusieurs blocs, avec gap et stabilité ;
3. l'indistinguabilité locale de tout sous-espace de code ;
4. seulement après, une jonction T et deux holonomies non commutatives.

## Claim boundary

Le marcheur et le transfert de paire sont des termes effectifs insérés dans
ce contrôle. Pas de phase thermodynamique, mémoire, fusion, tresse ou
non-abélianité n'est démontrée.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_two_block_pair_transport_audit.py
```
