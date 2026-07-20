# Project status and remaining walls

## Current scientific status

ANTLER has moved beyond a scalar exchange signature. The Phase 4.1 Strang runs propagate the complete two-dimensional logical frame, project the final evolution into the code, and separate loss from coherent rotation using a polar decomposition.

### Reference metrics

| Run | sigma_min | worst leakage | odd slope | off-diagonal norm | target F_avg |
|---|---:|---:|---:|---:|---:|
| theta=0.3, dt=0.25, R=4 | 0.9999531 | 9.38e-5 | -0.99643 | 1.34e-6 | 0.999999808 |
| theta=0.6, dt=0.25, R=4 | 0.9999523 | 9.53e-5 | -0.98394 | 9.90e-7 | 0.999984529 |
| theta=0.9, dt=0.25, R=4 | 0.9999526 | 9.49e-5 | -0.98109 | 7.19e-7 | 0.999951734 |
| theta=0.3, dt=0.25, R=6 | 0.9999579 | 8.42e-5 | -0.97386 | 5.72e-7 | 0.999989751 |

The coarse `dt=0.5` run fails the chosen leakage threshold (`1.91e-3`), which demonstrates that the tighter result depends on a converged integrator rather than optimistic projection.

## What has been demonstrated

1. A spatially separated code can be protected yet statistically inert.
2. Same-edge proximity restores first-order statistical control.
3. A shuttle exchange creates a phase odd in the hopping-statistics parameter.
4. The topologically trivial round trip provides a strong null control.
5. The complete logical subspace remains almost closed in the converged calculation.
6. The cleaned logical action is overwhelmingly diagonal and accurately approximates a phase gate.

## Remaining walls toward topological quantum computing

### Wall A  -  path-class invariance

Demonstrate that the odd logical phase is invariant under smooth deformations of the shuttle trajectory that preserve the exchange class. Vary well shape, timings, turn radius, transport distance, and local detours. Quantify the residual non-topological correction.

### Wall B  -  analytic effective theory

Derive an effective logical Hamiltonian or worldline/string-counting formula:

`Delta phi_odd = -theta + O(epsilon_adiabatic, finite-size, dressing)`.

This is needed to explain the numerical slope and identify correction terms.

### Wall C  -  realistic noise

Inject static disorder, temporal noise, calibration errors, timing jitter, hopping errors, and correlated noise. Report worst-case leakage and process fidelity distributions, not only means.

### Wall D  -  non-Abelian structure

The present phase is Abelian. A topological quantum computer requires a degenerate fusion space and non-commuting braid generators. Construct at least two independent exchanges `B1` and `B2` and test:

`B1 B2 != B2 B1`.

Without this, ANTLER is a geometric/anyon-inspired phase primitive, not universal topological computation.

### Wall E  -  universal logical gate set

Show arbitrary single-qubit control or a dense subgroup of SU(2), then a two-logical-qubit entangling operation. A protected Z-like phase alone is not universal.

### Wall F  -  scalability and fault tolerance

Move from one logical doublet to multiple encoded qubits, establish locality of controls, crosstalk bounds, error scaling, syndrome/readout strategy, and a threshold argument.

### Wall G  -  experimental mapping

Specify a realizable platform and derive how each Hamiltonian term is engineered. The rung-major fractional Jordan–Wigner convention defines a correlated-hopping model; claims must remain tied to that physical implementation unless representation independence is proven.

## Recommended next experiments

1. Automated homotopy/path-deformation scan at theta=0.3.
2. Noise Monte Carlo on the converged full logical doublet.
3. Effective phase derivation and finite-size scaling in L.
4. Two-exchange non-commutativity test on an expanded graph or synthetic dimension.
5. Candidate two-qubit encoding and entangling shuttle.
