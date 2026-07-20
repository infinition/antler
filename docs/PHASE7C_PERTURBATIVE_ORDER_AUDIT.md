# Phase 7C  -  audit d'ordre perturbatif des canaux charge 2

## Verdict

Le diagnostic initial sur la hiérarchie des termes est qualitativement juste pour les canaux charge-2 **lien-par-lien**, mais son énoncé quantitatif doit être corrigé. Le calcul exact à quatre barreaux mesure `XXXX` à l'ordre six environ, pas huit, pour deux topologies connectées représentatives. Les parasites restent d'ordre deux et dominent dans toute la fenêtre enregistrée. Un médiateur croisé à support plaquette est hors du catalogue de ce document et fait l'objet d'un audit séparé.

## Catalogue discret

`experiments/phase7/run_phase7c_channel_topology_catalog.py` énumère 368 graphes connectés de trois ou quatre liens sur quatre barreaux, avec canaux lien-locaux `E=aa/bb` et `M=ab/ba`, à `U=20`, `Delta=10`, `g=0.5`. Les 31 topologies à parités de rail nues et les 337 topologies à parités habillées échouent toutes `|c_XXXX|>1e-7`. Le meilleur coefficient de bon signe vaut `1.372e-10`, face à la cible `-0.5`.

## Mesure d'ordre

`experiments/phase7/run_phase7c_fourbody_order_scaling.py` balaie `g=0.30..1.50`, donc `g/Delta<=0.15`, pour le meilleur motif nu du catalogue et l'anneau mixte à signatures corrigées de ses signatures. Les ajustements donnent `XXXX~g^5.9886` et parasite `~g^1.9976` pour le motif nu, puis `XXXX~g^5.9901` et parasite `~g^1.9980` pour l'anneau mixte ; tous les `R²` sont supérieurs à `0.999999`.

À `g=1.5`, les coefficients `XXXX` restent entre `-8.96e-8` et `-9.82e-8`, tandis que la norme parasite relative est déjà environ `0.155`. La sélectivité requise reste donc inaccessible dans cette fenêtre pour les motifs mesurés.

## Corrections du raisonnement externe

Compter quatre conversions comme ordre minimal pour huit opérateurs externes est une borne de comptage utile. Cela ne prouve pas l'ordre huit pour tout anneau indépendant : les contractions internes et les chemins connectés doivent être calculés, et ici le résultat exact est d'ordre six. La majoration `0.073` issue de `g^4/(Delta^2 U)` est une estimation dimensionnelle dans la boîte de recherche, non une borne rigoureuse sans préfacteurs ni dénominateurs. Enfin, avec `C2~g^2/Delta` et `C4~g^4/(Delta^2 U)`, le ratio d'échelle est `(U/Delta)(g/Delta)^-2`, pas seulement `(g/Delta)^-2`.

## Décision

Ne pas lancer une optimisation continue aveugle des amplitudes lien-par-lien : aucun motif discret testé n'offre un signal quatre-corps utile, et les représentants restent dominés par les parasites. La prochaine proposition doit modifier la primitive lien-locale  -  par exemple médiateur croisé à support plaquette, médiateur de charge plus élevée, contrainte auxiliaire ou cancellation démontrée  -  avant une nouvelle optimisation.

Il s'agit d'un no-go de la grammaire locale testée, pas d'un no-go universel de tout Hamiltonien ANTLER 2D, ni d'une conclusion sur une phase tuilée ou une tresse.
