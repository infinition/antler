# Phase 7D  -  audit du compilateur Floquet multi-liens

## Objet

Ce contrôle compose la primitive de pair-hopping à quatre barreaux logiques.
Il applique les liens pairs `(0,1)` et `(2,3)`, puis le lien impair `(1,2)`.
Chaque porte de lien est contractuellement le pulse charge-2 fermé audité dans
le préflight signé-`ZZ`.

La cible numérique est

\[
  \sum_j (S_j^+S_{j+1}^+ + S_j^-S_{j+1}^-).
\]

Le contrôle agit après fermeture de chaque pulse sur le sous-espace monomère :
il mesure donc l'erreur de Trotter du compilateur, pas encore le crosstalk
microscopique de plusieurs médiateurs actifs.

## Résultats

Le commutateur entre couches pairs et impairs est non nul (`||[H_e,H_o]||_F=4`),
comme attendu. Pour des angles de paire `0.008, 0.012, 0.018, 0.027, 0.040`, le
résidu de Trotter croît de `6.40e-5` à `1.60e-3`. L'ajustement donne

\[
  \epsilon_{\rm Trotter} \simeq 0.9994\,\phi^{1.99987},
\]

avec `R²=0.9999999992`. Les parités des deux rails restent exactement nulles
dans la précision numérique et le défaut d'unitarité maximal est
`3.35e-16`.

## Décision

La primitive se compose comme un compilateur numérique local, avec une erreur
quadratique mesurée et contrôlable par l'angle de pas. Cela autorise un audit
de Hamiltonien/dynamique de ladder complet, mais ne constitue pas un parent
commutant et ne prouve aucune phase topologique. Les audits encore requis sont
la fuite et le crosstalk pendant pulses, les erreurs de rotation de rail, les
termes de hopping de jambe, puis les diagnostics de gap et de protection.

Résultat machine : `results/phase7/multilink_floquet_compiler_audit.json`.
