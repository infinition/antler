# Phase 7D  -  contrôle ED du ladder Floquet complet à nombre conservé

## Question testée

Après l'échec du pair-hopping minimal, le test suivant devait inclure les
termes supplémentaires du Hamiltonien effectif Floquet proposé par Defossez,
Vanderstraeten, Peralta Gavensky et Goldman : hopping intrajambe, attraction
intrajambe, densité inter-rail, échange et pair-hopping obtenus par rotation
globale des rails. Nous avons transcrit le contrôle haute fréquence

\[
H_{\rm eff}=\alpha H_0+(1-\alpha)P^\dagger H_0P,
\qquad P=e^{-i\pi J_x/2},
\]

sur notre convention fermionique rung-major. Cette étude est un **benchmark
externe** de petite taille ; elle n'est ni une reproduction complète du papier
ni une dérivation de pulses ANTLER à médiateur.

Référence : [Defossez et al., arXiv:2412.14886](https://arxiv.org/abs/2412.14886).

## Pré-enregistrement

Le scan `L=4,N=2` a couvert `u0=-0.5,-1,-1.5,-2` et
`alpha=0.25,0.5,0.75`, à hopping unitaire. Les quatre meilleurs points suivant
le ratio split/gap ont été prolongés à `L=6,N=4`. Un candidat devait satisfaire
simultanément :

1. commutateur de parité de branche inférieur à `1e-10` ;
2. `split/gap < 1e-2` ;
3. ratio du transfert local maximum dans le bulk sur le transfert de bord
   inférieur à `0.3`.

Le troisième test rejette explicitement un doublet faible mais localement
lisible.

## Résultats

La symétrie stroboscopique attendue est correctement reconstruite : les
commutateurs de parité valent au plus `2.33e-14` dans les quatre prolongements
`L=6`. Aucun ne passe le filtre.

Le meilleur split/gap prolongé est `0.05761`, pour `u0=-2, alpha=0.5`, mais son
ratio bulk/bord vaut `1.57065` : la réponse locale est plus grande dans le
bulk que sur les bords. Les trois autres ratios sont `1.58239`, `1.58239` et
`1.20589`. À `L=4`, aucun point n'avait déjà une réponse localisée (ratios
`1.066` à `1.295`).

## Décision et limite d'interprétation

Ce petit scan **ne trouve pas** de fenêtre candidate dans son domaine
préenregistré ; il ferme seulement cette fenêtre de benchmark. Il ne réfute
pas le résultat publié, lequel repose sur un régime et des analyses de taille
plus larges (dynamique exacte, bosonisation et MPS). Inversement, il interdit
de présenter notre compilateur de pair-hopping comme une phase protégée : sur
les tailles contrôlées ici, la localisation de bord échoue.

La voie Floquet reste donc une hypothèse matérielle mieux motivée que le
médiateur statique croisé, mais elle exige une dérivation fidèle de tout le
protocole, une exploration physique justifiée et un outil de taille
thermodynamique avant tout lancement RL. Aucun claim ANTLER natif,
topologique, 2D, de braid, non abélien, universel ou tolérant aux fautes n'est
établi.

Résultat machine : `results/phase7/full_floquet_ladder_preflight.json`.
