# ANTLER Phase 4.2–4.3  -  Path invariance, digital shuttle and disorder

## Executive result

The full dressed logical doublet was propagated with a converged Strang
split-operator integrator. The Gaussian shuttle is an excellent calibrated
geometric phase gate, but its odd phase depends on the well profile even in the
adiabatic limit. It is therefore not topologically path-invariant.

Replacing the Gaussian by a compact site-to-site shuttle and separating the two
rung transfers removes the composite avoided crossing and yields a much cleaner
exchange primitive.

## Gaussian shuttle: decisive no-go for strict path invariance

At theta=0.3, dt=0.25, T=20000:

- R=3: slope -0.99175
- R=5: slope -0.97647
- A=2.3: slope -0.95921
- A=2.9: slope -1.16146
- w=0.8: slope -0.97187
- w=1.2: slope -1.04572

Repeating the extreme A and w cases at T=30000 leaves the slopes unchanged to
~1e-4. The variation is physical path dependence, not nonadiabatic error.

## Digital shuttle

A compact two-site cos²/sin² cross-fade was introduced. Executing the left and
remote rung swaps simultaneously creates a composite degeneracy and up to 70%
leakage. Sequentializing the two swaps reduces leakage to ~1e-4.

### Distance invariance

For depth D=-4:

- R=3: slope -0.973811
- R=4: slope -0.973511
- R=5: slope -0.973457

The phase is effectively independent of shuttle distance.

### Strong-localization convergence

For R=4:

- |D|=3.5: slope -0.963004
- |D|=4.0: slope -0.973511
- |D|=4.5: slope -0.980234
- |D|=5.0: slope -0.984783
- |D|=6.0: slope -0.990396

The deviation from the ideal exchange phase decreases approximately as a power
of inverse trap depth (empirical exponent about 2.5 over this finite range).
This supports the interpretation that the remaining error comes from finite
wavefunction overlap during the exchange, rather than from shuttle length.

At D=-6, T=30000, dt=0.125:

- sigma_min = 0.99996947
- worst-case leakage = 6.106e-5
- odd phase = -0.29670195 at theta=0.3
- odd slope = -0.98900648
- off-diagonal norm = 7.97e-9
- average gate fidelity to the ideal Z phase = 0.999998187

## Static disorder

The D=-4 digital sequential gate was tested with fixed Gaussian on-site disorder
using the same realization for exchange/round-trip and +/-theta.

For sigma/J2 from 0.01 through 0.20, across two seeds per level:

- worst leakage remained below 1.16e-4
- odd slope remained between -0.97297 and -0.97466
- off-diagonal norm remained below 1.7e-6
- average gate fidelity remained above 0.999989

The dominant systematic is therefore finite localization, not static diagonal
disorder.

## Scientific interpretation

The digital protocol produces a highly coherent Abelian exchange-phase gate.
It approaches the expected phase in the strong-localization limit and is
robust to path length and large static disorder. However, strict topological
invariance has not been proven at finite trap depth: the residual phase error is
controlled by localization.

This is sufficient for a publishable geometric/statistical gate primitive, but
not yet for a universal topological quantum computer. The next fundamental
requirements are:

1. direct composition of exchange cycles;
2. an inverse-cycle identity test;
3. a second non-commuting logical generator or a measurement-assisted route;
4. a two-qubit entangling operation;
5. a protected many-body/fusion-space encoding rather than a trapped two-state
   cat alone;
6. a hardware mapping and error-threshold study.
