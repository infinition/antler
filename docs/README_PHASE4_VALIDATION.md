# ANTLER Phase 4  -  validation indépendante

## Verdict

La navette dynamique produit une signature impaire différentielle proche de la phase statistique attendue. Le run complet donne une pente moyenne de **-0.9774** sur θ = 0.3, 0.6, 0.9. Le contrôle round-trip reste proche de zéro.

Le contrôle de résolution montre toutefois que le résultat ne doit pas encore être présenté comme un tressage topologique universel :

- T=10 000 : pente -0.8655, fidélité minimale 0.802  -  régime non adiabatique, rejeté.
- T=20 000, N=5 000, R=6 : pente -0.9727, fidélité minimale 0.9945.
- T=30 000, N=7 500, R=6 : pente -0.9787, fidélité minimale 0.9904.
- T=20 000, N=5 000, R=4 : pente -0.9969, fidélité minimale 0.9951.

## Interprétation correcte

Le protocole démontre une **signature dynamique d'échange dans le modèle de saut corrélé défini par la convention Jordan–Wigner rung-major**. Il ne démontre pas encore :

1. une invariance topologique complète par rapport au chemin ;
2. une porte logique holonomique dans un sous-espace de qubit ;
3. une réalisation universelle d'anyons sur une échelle physique.

La dépendance résiduelle à R et T est faible mais mesurable. Elle doit être séparée en corrections non adiabatiques, phases géométriques ordinaires et effet de convention du modèle.

## Résultat scientifique actuel

Le résultat positif est plus précis que le claim initial : **l'encodage chat rend le canal statistique accessible et une navette fermée extrait une phase impaire presque linéaire en θ, alors que le trajet topologiquement trivial l'annule.**

## Go / No-Go

**GO conditionnel vers Phase 4.1.** La signature est robuste numériquement et mérite une suite. Le claim “braiding topologique démontré” reste NO-GO tant que les contrôles de déformation de chemin, de bruit et de projection logique ne sont pas passés.

## Phase 4.1 recommandée

- famille continue de chemins (R, largeur, forme du puits) à temps adiabatique comparable ;
- extrapolation T→∞ et Δt→0 ;
- Wilson loop / évolution projetée dans le sous-espace logique ;
- bruit statique et erreur de contrôle ;
- test d'une autre convention de représentation avec transformation unitaire cohérente.
