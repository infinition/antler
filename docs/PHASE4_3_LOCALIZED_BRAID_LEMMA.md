# ANTLER Phase 4.3  -  Localized braid lemma

## Statement

Consider the rung-major correlated-hopping ladder with two hard-core particles,
initially localized on sites `{0,1}`. Under the sequential compact shuttle:

1. the mobile particle travels from site `0` to `2R` on leg 0;
2. the stationary particle moves `1 -> 0` across the left rung;
3. the mobile particle moves `2R -> 2R+1` across the remote rung;
4. the mobile particle returns to site `1` on leg 1.

In the strictly localized Fock-state limit, the product of fractional
Jordan–Wigner phases is

\[
\mathcal A_{\rm ex}/|\mathcal A_{\rm ex}|=e^{-i\theta}.
\]

For the matched round trip, where the mobile particle returns on leg 0 while
the stationary particle remains on site 1,

\[
\mathcal A_{\rm rt}/|\mathcal A_{\rm rt}|=1.
\]

Therefore the differential exchange phase is exactly

\[
\Delta\phi_{\rm odd}=-\theta\pmod{2\pi},
\]

independently of the turnaround distance `R`.

## Proof by string counting

In the implemented convention, an increasing-index hop `k -> l` contributes
`exp(-i theta n_mid)`, where `n_mid` is the occupation strictly between the two
sites. A decreasing-index hop contributes the Hermitian-conjugate phase.

During the exchange path, only the first hop `0 -> 2` has an occupied
intermediate site (`1`), hence exponent `-1`. Every later longitudinal hop has
an empty intermediate site, and rung hops have no intermediate site. The total
exponent is consequently `-1`.

During the round trip, the outbound `0 -> 2` hop has exponent `-1`, while the
return `2 -> 0` hop has exponent `+1`. They cancel exactly.

The executable enumeration is in
`experiments/phase4_1/run_phase4_3_exact_path_count.py`.

## Finite-localization correction

For finite trap depth, the transported eigenstate is not a single Fock
configuration. Virtual tails sample additional correlated-hopping links. The
measured quantization error for depths 3.5–6.0 is fitted by

\[
1-|\Delta\phi/\theta|
 = (0.84385\pm0.00426)|\Delta|^{-(2.4963\pm0.0037)},
\]

with `R^2 = 0.99999498`.

Allowing an asymptotic offset gives

\[
c=(-1.04\pm1.52)\times10^{-4},
\]

which is compatible with zero. This supports convergence to the exact
localized braid value rather than a persistent path-dependent offset.

This is currently an empirical scaling law, not yet an analytic derivation of
the exponent 2.496. A perturbative Schrieffer–Wolff treatment of the virtual
trap tails is the next theoretical task.
