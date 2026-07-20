# Phase 8B - correction : fermeture Rabi n'est pas une porte

## Erreur de critere corrigee

Les premiers audits de pulse rapportaient une faible distance absolue entre
l'unitaire polaire microscopique et la cible Schrieffer--Wolff. Ce critere est
insuffisant si la cible non scalaire est elle-meme tres petite.

Le controle correct est

```text
erreur relative = d(U_phys, U_SW) / d(U_SW, I),
```

apres retrait de la phase globale. Une porte realisee exige une erreur relative
petite; une valeur proche de un signifie que le propagateur physique est reste
presque scalaire et n'a pas accumule le signal cible.

## Lien local

Au point profond `g/Delta=0.025`, la fermeture Rabi entiere de deux segments
donne :

| axe | signal SW | signal physique | erreur relative |
|---|---:|---:|---:|
| `X` | `3.9289e-5` | `2.7381e-7` | `1.0063` |
| `Y` | `3.9289e-5` | `2.7381e-7` | `1.0063` |
| `Z` | `7.8602e-5` | `5.0493e-7` | `1.0063` |

La fuite peut etre faible (`~5e-7`), mais le signal logique est egalement
supprime. La famille de timings Rabi entiers est donc un retour de population,
pas une implementation de l'Hamiltonien effectif statique.

## Echo C3

Le groupe de parite C3 est valide algebraiquement apres downfolding : il annule
les parasites et garde `XXX`. Mais son propagateur microscopique a la meme
pathologie. A `g/Delta=0.025`, le signal SW du groupe est `2.114e-7`, le
signal physique `3.30e-9`, et l'erreur relative `1.00012`. Aucun point du
scan de multiplicateurs Rabi `1,2,4` ne passe l'ecran pratique.

Un balayage exact independant de deux durees libres `t_A,t_B in [0,8]` pour
`g/Delta=0.10,0.05,0.025` confirme le diagnostic : aucune rotation `XX`
non triviale ne passe simultanement fuite `<1e-4` et fidelite de trace
`>0.999`. Les meilleurs points retombent sur l'angle minimal impose
`0.05` rad, c'est-a-dire presque l'identite.

## Preflight Floquet rapide abrupt

Une autre fenetre, distincte de la fermeture Rabi et de la boite a deux
durees, a ete verifiee exactement : repetitions du mot abrupt `A B`, avec
`g/Delta = 0.05, 0.025`, durees totales `20, 100` et pas de segment
`0.025 .. 0.4`. Les 20 lignes echouent. La meilleure erreur relative au
signal SW vaut `14.2956`, tres loin du seuil pre-enregistre `0.1`; les
fuites sont egalement superieures a `1e-4` dans toutes les lignes.

Ce resultat ferme **la fenetre de commutation rapide abrupte definie ici**.
Il ne constitue pas un no-go pour une commande derivee plus riche : rampe
compensee, preparation d'etat habille, pulse off-resonant/composite, ou une
nouvelle ressource microscopique doivent cependant etre definis avant toute
optimisation.

## Statut et suite

Les liens locaux `X/Y/Z` et l'echo de parite C3 restent des resultats
**statiques/downfolded**. Toutes les revendications de porte pulsee associees
a la fermeture Rabi sont retirees.

Le prochain probleme est bien defini : construire une sequence qui garde les
mediateurs virtuels mais accumule leur phase conditionnelle (pulse off-resonant
habille, rampe adiabatique, ou sequence composite derivee), puis mesurer
leakage, signal physique et erreur relative au signal.

## Reproduction

```powershell
python experiments/phase7/run_phase8b_shared_matter_pulse_closure_audit.py
python experiments/phase7/run_phase8b_shared_matter_xyz_link_audit.py
python experiments/phase7/run_phase8b_shared_matter_c3_parity_echo_pulse_audit.py
python experiments/phase7/run_phase8b_shared_matter_c3_parity_echo_duration_scan.py
python experiments/phase7/run_phase8b_shared_matter_two_segment_gate_search.py
python experiments/phase7/run_phase8b_shared_matter_fast_floquet_preflight.py
```
