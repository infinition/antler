# ANTLER 0.1 : Phase 0 (ED ladder SSH anyonique)

Implementation numerique de la specification ANTLER-0.1, avec correctifs
issus de la validation. Python 3.10+, numpy uniquement (scipy optionnel).

```
antler/
  basis.py     base hard-core, ordonnancement rung-by-rung, masques binaires
  model.py     H(lambda) exact, phases JW fractionnaires, hermiticite verifiee
  edge.py      extraction et rotation des 4 orbitales de bord (L1,L2,R1,R2)
  logical.py   etats 2 particules, base logique, projecteur P_L, leakage
run_phase0.py  criteres C1..C6, parametres lambda_ref
```

Usage : `python run_phase0.py`

## Parametres valides (lambda_ref)

L=14, J1=0.4, J2=1.0, Jperp=0.1, mu_edge=-0.35 sur les 4 sites extremes,
theta=0. Dimension exacte du secteur N=2 : 378.

## Verrous leves par la Phase 0 (a integrer dans la spec avant Phase 1)

### V1. Le doublet logique n'est PAS le fondamental

Les modes de bord SSH sont mid-gap. A 2 particules hard-core, le
fondamental est un etat de volume (E0 = -2.61 a L=6). Toute metrique
QGT calculee sur "l'etat fondamental obtenu par ED" (spec ANTLER-0.1,
section 5) mesure la geometrie du volume, pas celle du qubit. Correctif :
le QGT de la Phase 2 doit etre le QGT non abelien sur le sous-espace P_L
(formulation Wilczek-Zee), pas le QGT scalaire du fondamental.

### V2. Shell chiral E=0 : degenerescence massive avec le volume

La symetrie chirale rend le spectre 1p symetrique. Toute paire de modes
de volume (e, -e) donne un etat 2p a E=0, exactement la ou vivrait le
doublet logique naif. Le doublet serait noye dans un shell degenere de
taille O(L^2). Correctif adopte : biais de bord mu_edge = delta sur les
4 sites extremes. Les 4 orbitales de bord se decalent de delta, le
doublet logique migre a E ~ 2*delta, hors du shell. Verifie : capture
logique 1.99981/2 a delta=-0.35.

### V3. Definition figee du qubit : doublet chiral, pas doublet de jambe

Les etats du draft (|L1 R1>, combinaisons de jambe) ne sont pas quasi
degeneres : Jperp les etale sur une bande de largeur ~2*Jperp_eff
(verifie numeriquement : recouvrement disperse sur des etats a
E = -0.226, 0, +0.226 a L=6). Definition figee :

```
B_X = (X1 + X2)/sqrt(2)   bonding rung, energie delta - Jperp_eff
A_X = (X1 - X2)/sqrt(2)   antibonding,  energie delta + Jperp_eff

|0_L> = |B_L A_R>         |1_L> = |A_L B_R>
```

Les deux etats ont energie (delta - Jp) + (delta + Jp) = 2*delta :
degenerescence EXACTE par construction (verifie : split = 0 machine),
protegee tant que la symetrie de reflexion gauche-droite tient.
Gap interne de la bande de bord vers |B_L B_R> et |A_L A_R> : 2*Jperp.

### V4. Hybridation L-R : L=6 est trop court

t_LR(L=6) = 5.4e-2, du meme ordre que Jperp : les orbitales gauche/droite
ne sont pas separables proprement. Decroissance verifiee :

| L  | t_LR    |
|----|---------|
| 6  | 1.1e-1  |
| 8  | 4.3e-2  |
| 10 | 1.7e-2  |
| 14 | 2.7e-3  |

MVP fixe a L=14 (d=378, ED instantanee). L=6 ne doit plus etre utilise.

### V5. La corde JW brise la symetrie de jambe P a theta != 0

Avec l'ordonnancement rung-by-rung, la phase du saut leg 1 depend de
l'occupation (leg 2, rung i) mais celle du saut leg 2 depend de
(leg 1, rung i+1). P: leg1 <-> leg2 n'est donc PAS une symetrie de
H(theta != 0). Consequences :
- P ne peut servir qu'a definir la base logique a theta = 0 (fait)
- H(theta) est un modele de saut correle DEFINI par la convention
  d'ordonnancement ; sur une echelle, la statistique fractionnaire par
  corde 1D est dependante du chemin. Le papier doit le presenter ainsi,
  pas comme un anyon universel.

### V6. Robustesse statique en theta (P_L fige a theta = 0)

| theta  | split doublet | fuite statique |
|--------|---------------|----------------|
| 0.25   | 1.5e-5        | 0.21 %         |
| 0.50   | 5.8e-7        | 0.02 %         |
| 1.00   | 1.3e-4        | 0.41 %         |
| pi/2   | 1.8e-3        | 15.1 %         |

theta ~ pi/2 ouvre un canal de fuite statique majeur : le domaine de
theta exploitable par l'optimisation devra etre borne, ou la fuite
integree dans la fitness exactement comme prevu (terme alpha*P_leak).

### V7. Protection spectrale a reformuler

A L=14 le spectre 2p est dense (espacement moyen ~1.5e-2) : la distance
spectrale locale du doublet au reste est 2.4e-3. La protection reelle
est (i) la degenerescence interne exacte du doublet, (ii) le gap interne
2*Jperp de la bande de bord, (iii) la faiblesse des elements de matrice
vers le volume (localite), attestee par la capture 0.9999. Le terme
"Gap" de la fitness doit designer le gap interne de bande, pas un gap
spectral global qui n'existe pas.

## Sequence suivante

Phase 1 : ED systematique theta != 0 (spectres, scan delta x theta,
domaine de fuite). Phase 2 : QGT non abelien sur P_L + volume normalise
V_Q_bar. Phase 3 : boucles lambda(t) manuelles avec P_leak(t). Phase 4 :
optimisation bayesienne. L'IA n'entre qu'en Phase 4.

## Phase 1 : resultats (test d'hypothese statique)

Question testee : existe-t-il une region (theta, delta) avec courbure
non abelienne non nulle dans le sous-espace logique et P_leak ~ 0 ?

### V8. No-go de symetrie sur la variete symetrique

T = K o Rev (conjugaison complexe composee avec le renversement de
sites k -> M-1-k) est une symetrie anti-unitaire EXACTE de H(theta,
delta) pour tout theta tant que le biais de bord est symetrique
gauche-droite. Verifie : ||T H T^-1 - H|| = 0.0 machine. Consequence :
F identiquement nulle sur toute la variete (theta, delta_sym). Verifie
sur 50 plaquettes vertes (F = 0) et par plaquette fine (1.4e-12).
Un troisieme controle brisant T est obligatoire : delta_L != delta_R.

### V9. T brisee : la courbure apparente est un artefact de filaments

Scan (theta, delta_L), delta_R = -0.35 : ||F|| jusqu'a 1.76 dans la
zone P_leak < 1%. MAIS le test de convergence des plaquettes echoue :
|log W| ne scale pas en h^2 (champ lisse) mais de facon erratique
(ratios 0.8 a 448 selon le point). La courbure est concentree sur les
filaments d'anticroisements doublet-volume (gap local ~ 1e-3). La
courbure de fond lisse est < 1e-3 : inutilisable pour une holonomie
controlee (transiter pres d'un filament = machine a fuite).

### V10. Theoreme d'inertie : theta n'agit pas dans le code separe

Mesure directe : ||P_L (dH/dtheta) P_L|| = 8e-19 (zero machine), alors
que ||(1-P_L)(dH/dtheta) P_L|| ~ 2e-3 a 8e-3. Interpretation : les deux
particules logiques vivent aux bords opposes (P(meme moitie) = 3e-3),
les cordes JW ne sont jamais activees a l'interieur du code. theta est
donc un canal de FUITE PURE pour cet encodage, sans aucune action
logique. La tension architecturale est fondamentale : la separation
spatiale qui protege le qubit decouple exactement le qubit de la
statistique.

### VERDICT Phase 1

Hypothese statique FALSIFIEE. Aucune boucle adiabatique dans l'espace
(theta, delta_L, delta_R) ne peut produire d'holonomie logique avec
l'encodage a bords opposes : ce n'est pas un probleme d'optimisation,
c'est un zero structurel.

### Redefinition de la Phase 2

La statistique exige la proximite des particules. Trois voies :
  A. Protocole d'echange dynamique : navette d'une particule a travers
     le volume (potentiel mobile mu_i(t)), activation transitoire des
     cordes, holonomie de type tressage effectif. Voie principale.
  B. N = 3 : particule de controle mobile dont la position est un
     parametre de la variete.
  C. Encodage a recouvrement spatial (etats meme bord |B_X A_X>),
     au prix du cout d'interaction et de la perte de separation.

Le resultat V8-V10 est publiable en soi comme contrainte de design :
la protection par separation et le controle statistique statique sont
mutuellement exclusifs dans les echelles SSH anyoniques.

## Phase 2 (kickoff) : lemme, correction d'encodage, plan de controle

### V11. Lemme ANTLER-1 : certifie partiellement, mecanisme reel plus fort

Certifie exact (0.0 machine) : les sauts de rung sont Pi-pairs, les
sauts de jambe Pi-impairs, dH/dtheta est integralement Pi-impair.
Le lemme (P_L dans un secteur de Pi => P_L O_impair P_L = 0) est un
theoreme valide MAIS son hypothese n'est pas satisfaite par nos codes
(var(Pi) = 0.18 et 0.45) alors que le zero tient quand meme, jambe par
jambe (3e-19 et 6e-19, sans compensation croisee). Le mecanisme
operatoire est la PROXIMITE : un element intra-code de dH/dtheta exige
une amplitude sur des configurations a particules adjacentes dans le
bra ET le ket. Regle de design ANTLER-1' : l'action statistique au
premier ordre est bornee par le produit des poids de proximite des
etats logiques. Item ouvert : le mecanisme exact du zero profond
(1e-19 au lieu du ~1e-12 attendu par comptage de queues).

### V12. L'encodage |B_L A_L> etait invalide (interference bosonique)

Le produit symetrise de deux combinaisons orthogonales de jambes
s'annule sur la configuration dominante (rung 0 doublement occupe en
legs) : etat delocalise, capture 0.34, zeros factices dans le run
precedent. Base chat correcte : {|L1 L2>, |R1 R2>}.

### V13. Encodage chat corrige : le plan (theta, delta_asym) est utilisable

Resultats a lambda_ref (L=14) :
  - capture doublet = 1.691/2 (85%, a ameliorer par habillage du code)
  - ||P_L dH/dtheta P_L|| = 5.9e-2 : NON NUL par construction
    (particules adjacentes au meme bord, cordes actives)
  - action differentielle Z sous delta_asym = 0.10 :
    d(E_1 - E_0) entre theta=0 et 0.9 : -1.38e-2 (vs -3.4e-3 residuel
    a delta_asym = 0, artefact de tracking a elucider)
  - courbure F(theta, delta_asym) en (0.6, 0.05) : ||F|| ~ 17 a 21 avec
    scaling |log W| en h^2 (ratios 3.78, 3.42) : champ LISSE, ce que le
    plan (theta, delta_L) de l'encodage separe n'a jamais donne.

Verdict Phase 2 kickoff : l'hypothese renait avec l'encodage chat.
Compromis explicite : on echange la protection maximale (separation)
contre l'accessibilite du controle statistique, conformement a la
regle ANTLER-1'.

### Prochaines etapes

  1. Habillage du code (etats quasi-adiabatiques) pour porter la
     capture > 99% ; sinon reduire J1 ou augmenter le gap.
  2. Elucider le residu differentiel a delta_asym = 0 (symetrie Rev
     devrait l'annuler) et le mecanisme du zero profond V11.
  3. Premiere boucle holonomique manuelle (theta, delta_asym) fermee
     avec P_leak(t) et phase relative mesuree (Phase 3 originale).

## Phase 3 : boucle holonomique fermee (resultats)

### V14. Regime profond : isolation par liaison de bord (delta = -2.0)

Le doublet chat entre en collision avec le continuum a delta faible
(gap local 2e-4 a delta = -0.35). Correctif : biais profond delta = -2.0.
Les modes de bord deviennent des etats lies du puits, le doublet chat
{LL, RR} tombe sous les continuums 2p ((2,0) vs (1,1) : secteurs
quasi superselectionnes, tunneling e^{-L}). Resultat : split du doublet
1.1e-10, cadre localise a 99.5%, capture 0.94 a 0.9998 sur le domaine
de controle, P_leak dynamique max 3e-4 sur la boucle.

### V15. Chaine de mesure holonomique certifiee de bout en bout

Trois methodes independantes concordent :
  - courbure WZ par plaquettes : convergee (ratio scaling = 1.00),
    decomposition Pauli : composante X nulle (4e-9, tunneling LL-RR
    exponentiel), action purement Z + abelienne
  - transport de Wilson sur la boucle : phi = 6.2e-3 rad
  - evolution temporelle reelle (T = 800, 240 segments), aller/retour :
    phi_geo = (phi_a - phi_r)/2 = -6.2e-3 rad, phi_dyn = -1.52 rad,
    U_logical hors diagonale 6.6e-6, P_leak max 3e-4.
Accord Wilson / temps reel exact en amplitude (le signe releve de la
convention d'orientation du produit de Wilson). La phase geometrique
logique est MESUREE et reproductible : ~6 mrad par boucle
(theta_max = 0.8, da_max = 0.3), avec fz lineaire en da
(kappa ~ 0.085) et capture stable ~0.975 jusqu'a da = 0.30.

### V16. Frontiere courbure-fuite : compromis structurel

Scan (delta, J1) : aucun point avec capture > 0.99 ET fz > 0.01.
  delta=-2.0 : capture 0.94-0.96, fz ~ 5e-3 (propre mais faible)
  delta in [-1.2, -0.7] : fz ~ 1-2 mais capture 0.3-0.65 (detruit)
Raison physique : la courbure WZ et la fuite ont la MEME source
microscopique (elements de matrice hors code / gap). C'est ANTLER-1'
au niveau superieur : en controle statique, courbure forte et
sous-espace ferme sont antagonistes. Le regime exploitable est une
"holonomie milliradian" : exacte, protegee, faible.

### Bilan des trois nombres du reviewer

  ||P dH/dtheta P|| != 0 : ACQUIS (5.9e-2 encodage chat)
  capture logique : 0.975 a 0.9998 sur le domaine utile (cible > 0.99
    atteinte sur boucle restreinte, 0.975 min sur boucle elargie)
  holonomie utile : DEMONTREE et triple-certifiee, amplitude mrad ;
    insuffisante comme porte, valide comme signature mesurable

### Phase 4 recommandee

La phase dynamique differentielle controlee par theta (~1.5 rad/boucle
ici) est 250x plus grande que la phase geometrique : primitive de porte
Z realiste a court terme. Pour une holonomie geometrique O(1), la voie
statique est fermee (V16) : il faut l'echange REEL (navette mu_i(t)
d'une particule autour de l'autre, phase de tressage O(theta)
independante de l'aire de boucle). C'est l'experience decisive du
programme : braiding synthetique vs holonomie de courbure.
