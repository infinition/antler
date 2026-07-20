# Phase 8  -  route dynamique à nombre conservé : le pivot utile

## Le mur qui est réellement évité

La route statique ANTLER à médiateurs charge-2 cherchait un stabilisateur ou
un pair-hopping connecté dans une élimination perturbative. Elle est limitée
par une hiérarchie défavorable : le terme à quatre corps apparaît à un ordre
plus élevé que les parasites à deux corps. Le no-go Phase 7C reste valide pour
cette grammaire statique et ce régime perturbatif.

La route dynamique est différente. On part d'un Hamiltonien de deux jambes

\[
H_0=H_{\rm hop}+U_0\sum_{j,\ell}n_{j,\ell}n_{j+1,\ell}
\]

et d'une rotation globale de rail `P=exp(-i eta Jx)`. Un cycle Floquet donne,
à haute fréquence,

\[
H_{\rm eff}=\alpha H_0+(1-\alpha)P^\dagger H_0P+O(T).
\]

La rotation conjugue directement les interactions de densité en termes
d'échange et de transfert de paires. Le pair-hopping est donc d'ordre
principal dans `U0`, pas un terme faible d'ordre `g^4/Delta^3`. C'est le
mécanisme proposé pour une échelle à Majoranas à nombre conservé par Defossez
*et al.* (arXiv:2412.14886v2, 2025).

## Ce qu'ANTLER possède déjà, et ce qu'il doit encore dériver

| ressource | statut ANTLER | condition avant promotion |
| --- | --- | --- |
| rotations de rail à nombre conservé | réalisées localement dans les compilateurs Phase 7D | démontrer une rotation synchrone sur tout le ladder sans erreur différentielle de rail |
| modulation de phases Peierls | contrôle local convergé et tolérances locales mesurées | intégrer le contrôle dans la séquence Floquet complète |
| interaction attractive NN `U0` | pas encore reliée de façon native au Hamiltonien gelé | dériver son coefficient et ses parasites depuis les médiateurs ANTLER |
| pair-hopping dynamique | démontré comme mécanisme externe, pas ANTLER natif | comparer le cycle microscopique complet à `H_eff` à fréquence finie |

## Correction du benchmark précédent

L'ancien préflight Floquet ne clôt pas cette route, pour trois raisons
documentées :

1. son passage `L=4,N=2` à `L=6,N=4` changeait la densité de `1/4` à `1/3` ;
2. à charge totale impaire, la symétrie d'échange des rails force l'égalité des
   secteurs de parité de branche, donc leur split n'est pas un diagnostic ;
3. son ratio de transfert local bulk/bord n'est pas le couple de diagnostics
   (gap de bulk, spectre d'intrication) utilisé par la référence.

La réplication corrigée conserve la densité et sérialise les spectres de
Schmidt. Elle contient une séquence à charge paire, où le split de parité n'est
pas imposé par symétrie. Les candidats `U0=-2`, `alpha=0.25,0.5,0.75`, déjà
préenregistrés par le scan antérieur à `L=6,N=4`, sont réévalués avec ce
diagnostic ; aucune nouvelle fenêtre n'est choisie après coup.

## Mise à jour canonique

L'audit `PHASE8_CANONICAL_MPS_AUDIT.md` conserve désormais `U(1) x Z2`
directement dans les tenseurs. Il confirme le split de parité jusqu'à `L=16`
et le diagnostic de gap chargé périodique. La fermeture du gap neutre est
attendue pour le mode de charge total gapless ; ce n'est pas un critère de
rejet du secteur topologique relatif.

## Critères du prochain vrai jalon

La route devient un candidat ANTLER seulement si les quatre conditions
suivantes passent :

1. dérivation microscopique de `U0` et de la rotation de rail pulsée ;
2. accord contrôlé entre le cycle à fréquence finie et `H_eff`, avec termes de
   Magnus bornés ;
3. réplication indépendante à densité fixe et grande taille (MPS/DMRG) du gap
   **chargé** relatif et du spectre d'intrication, sans rejeter la branche sur
   le gap neutre volontairement gapless ;
4. géométrie de jonction qui permette un échange non homotope, après la phase
   1D établie.

Une réussite aux étapes 1--3 constituerait une avancée importante : une
plateforme ANTLER dynamique pour modes Majorana à nombre conservé. Elle ne
constituerait toujours pas, à elle seule, un ordinateur topologique ou une
preuve de braid non abélien.
