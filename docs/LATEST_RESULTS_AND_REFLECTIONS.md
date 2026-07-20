# ANTLER v0.7  -  Latest Results and Scientific Reflections

**Status date:** 2026-07-18  
**Scope:** Phase 4.1 to Phase 4.5, including full logical-gate tomography, path-deformation tests, digital sequential exchange, localization scaling, static disorder, exact string counting, and the two-level handoff lemma.

## 1. Strongest fully validated logical-gate result

For the Gaussian shuttle at `theta=0.3`, `R=4`, `T=20000`, and converged Strang step `dt=0.25`:

- `sigma_min = 0.9999531172`
- worst-case leakage `= 9.3763e-5`
- coherent off-diagonal norm `= 1.3378e-6`
- odd logical phase `= -0.2989277507`
- slope `Delta phi_odd / theta = -0.9964258355`
- cleaned average gate fidelity `= 0.9999998084`

This validates a coherent logical phase gate in the model. It does not by itself establish topological path invariance.

## 2. Numerical wall removed: apparent recapture leakage

The early full-doublet test reported roughly 12% round-trip leakage. A dedicated recapture harness and a converged Strang split-operator showed that this was mainly a timestep/integrator artifact rather than a physical avoided-crossing catastrophe.

At converged resolution:

- pre-recapture adiabatic purity exceeds 99.99%;
- local recapture leakage is of order `1e-4` or below;
- the complete doublet remains almost closed.

The coarse `dt=0.5` gate run still reaches the correct phase but fails the chosen leakage threshold (`~1.9e-3`). The `dt=0.25` result is therefore the accepted reference.

## 3. Gaussian path deformation: high-fidelity but not topologically locked

At `theta=0.3`, smooth deformations of the Gaussian shuttle preserve the logical subspace but change the odd phase:

| deformation | odd phase | slope |
|---|---:|---:|
| `R=3` | `-0.29753` | `-0.99175` |
| `R=5` | `-0.29294` | `-0.97647` |
| `A=2.3` | `-0.28776` | `-0.95921` |
| `A=2.9` | `-0.34844` | `-1.16146` |
| `w=0.8` | `-0.29156` | `-0.97187` |
| `w=1.2` | `-0.31371` | `-1.04572` |

Repeating the extreme `A` and `w` cases at `T=30000` changes the phase only at the `1e-4` level. The path dependence is physical in the adiabatic limit. The Gaussian protocol is therefore a calibrated geometric/statistical gate, not a finite-depth topological invariant.

## 4. Digital sequential shuttle: isolation of the exchange class

A compact site-to-site shuttle was introduced to replace the broad Gaussian well. Performing the left and remote rung transfers simultaneously creates a composite near-degeneracy and can produce extreme leakage. Sequentializing these transfers removes that failure mode and restores leakage near `1e-4`.

### Distance invariance

At trap depth `D=-4`:

- `R=3`: slope `-0.973811`
- `R=4`: slope `-0.973511`
- `R=5`: slope `-0.973457`

The spread is only about `3.5e-4`. The digital sequential phase is essentially independent of shuttle distance, as predicted by string counting.

### Localization convergence

For `R=4`, the phase approaches the ideal exchange value as the trap is deepened:

| `|D|` | odd slope |
|---:|---:|
| 3.5 | -0.963004 |
| 4.0 | -0.973511 |
| 4.5 | -0.980234 |
| 5.0 | -0.984783 |
| 6.0 | -0.990396 |

A pure-power fit gives

`1 - |Delta phi / theta| = 0.84385 |D|^(-2.49627)`

with `R^2 = 0.99999498`. An offset fit yields an asymptotic offset compatible with zero. The empirical limit is approximately `-1.0001`.

The phase error also correlates almost perfectly with the delocalized orbital weight (`Pearson r = 0.9999988`). This identifies finite wavefunction tails as the dominant correction.

A slower and finer `D=-6`, `T=30000`, `dt=0.125` run gives:

- `sigma_min = 0.99996947`
- worst leakage `= 6.106e-5`
- odd phase `= -0.29670195` for `theta=0.3`
- slope `= -0.98900648`
- off-diagonal norm `= 7.97e-9`
- average target fidelity `= 0.999998187`

## 5. Exact localized path count

The Fock-path enumerator proves, in the strict localized limit:

- exchange exponent `= -1`, hence phase `exp(-i theta)`;
- matched round-trip exponent `= 0`, hence phase `1`.

Only the first longitudinal hop of the exchange sees the occupied intermediate site. The return path on the other leg has no compensating activated string. In the round trip, the outgoing and returning crossings cancel.

This establishes exact distance-independent string counting in the ideal localized path.

## 6. Exact two-level handoff lemma

For an isolated two-level transfer with link phase `phi`, adiabatic transport transfers the state with phase `-phi`, independently of the ramp shape.

The numerical lemma checks:

- `phi = 0.2, 0.6, 1.1`;
- `D/J = 2, 4, 8, 16`;
- linear, `sin^2`, and smoothstep ramps.

Maximum phase error is below `2.6e-14`.

This provides the local building block for the digital braid:

`phi_odd = -theta * sum_j eta_j n_j`.

At finite depth, additional virtual configurations invalidate the strictly isolated two-state reduction and generate the observed localization correction.

## 7. Dynamic odd contribution at theta = 0

At `theta=0`:

- the Hamiltonian can be chosen real;
- `dH/dtheta` is purely imaginary and antisymmetric;
- the expectation value on any instantaneous real eigenstate vanishes.

Thus the first-order odd dynamical phase is exactly zero. The measured odd exchange phase is a geometric/transport phase rather than an uncancelled first-order dynamical background.

A previously generated coarse response integration file (`results/phase4_4/response_baseline.json`) has large leakage and is retained for audit only. It is not used as evidence for the analytic statement above.

## 8. Static disorder

The sequential digital gate at `D=-4` was tested using fixed Gaussian on-site disorder, shared between exchange, round-trip, and `+/-theta` branches.

For disorder `sigma/J2 = 0.01, 0.03, 0.05, 0.10, 0.20`, with two seeds per level:

- worst leakage remains below `1.16e-4`;
- odd slope remains between `-0.97297` and `-0.97466`;
- coherent off-diagonal norm remains below `1.7e-6`;
- target fidelity remains above `0.999989`.

The dominant systematic in this range remains finite localization, not static diagonal disorder.

Logs for larger noise levels (`0.5` and above) are incomplete and are not treated as results.


## 9. Digital timestep and angle audit

The sequential digital gate was rerun at `D=-4`, `theta=0.3`, `dt=0.125`:

- `sigma_min = 0.9999849075`;
- worst leakage `= 3.0185e-5`;
- odd slope `= -0.97204098`;
- off-diagonal norm `= 8.99e-8`;
- target fidelity `= 0.9999882745`.

This confirms the gate closure and coherent diagonality at a finer timestep. The phase changes by roughly `1.5e-3` relative to the `dt=0.25` digital run, so the deepest localization-scaling numbers should still be interpreted with their stated finite-step uncertainty.

At `D=-4`, `dt=0.25`:

- `theta=0.6`: odd phase `-0.5863284`, slope `-0.9772140`, target fidelity `0.99996885`;
- `theta=0.9`: odd phase `-0.8843530`, slope `-0.9826144`, target fidelity `0.99995920`.

The digital family remains a high-fidelity diagonal phase family across the tested angle range, with finite-localization/higher-order corrections to exact linearity.

## 10. Abelian braid-algebra wall

The current logical gates are all effectively diagonal `Z` rotations. Their commutators vanish numerically, and the generated Lie algebra has rank one. A dedicated no-go comparison records:

- current Abelian commutator norm `= 0`;
- current generated Lie rank `= 1`;
- a target two-dimensional non-Abelian braid representation has nonzero commutator norm `sqrt(2)` while satisfying the braid relation.

This closes an important conceptual question: improving the fidelity or quantization of the present exchange phase cannot by itself produce universal all-topological computation. The architecture must acquire a genuinely degenerate internal/fusion space with matrix-valued, non-commuting exchanges.

## 11. Experimental acceleration branch

A Numba/local-Trotter Phase 4.6 prototype reduces the runtime of a full `T=20000`, `dt=0.25` digital run to roughly 24 seconds and reproduces the same qualitative phase gate with low leakage. Its slope differs from the primary sequential solver by several `1e-3`; it is retained as an acceleration prototype, not as the reference physics result, until cross-integrator convergence is completed.

## 12. Present claim boundary

### Demonstrated

- a structural protection/control trade-off;
- a static no-go for the spatially separated encoding;
- a high-fidelity Abelian logical phase gate;
- a strong round-trip null control;
- full logical-subspace closure and polar-decomposition audit;
- near-exact exchange phase in the localized digital limit;
- distance invariance of the digital path;
- strong static-disorder robustness up to `20%` in the tested model;
- an exact local handoff phase lemma;
- an exact ideal-path string-counting result.

### Not demonstrated

- finite-depth invariance under every smooth homotopy;
- a non-Abelian fusion space;
- non-commuting braid generators;
- a universal protected gate set;
- an entangling two-logical-qubit gate;
- a threshold theorem or scalable fault-tolerant architecture;
- an experimental realization of the exact rung-major correlated-hopping Hamiltonian.

## 13. Current deepest scientific interpretation

The Gaussian shuttle is a robust calibrated geometric gate. The sequential digital shuttle isolates a discrete exchange class. In the strictly localized limit, its phase is exactly the oriented Jordan-Wigner string count. Finite-depth deviations are not random calibration errors: they are controlled by wavefunction delocalization and vanish consistently with increasing localization.

The next fundamental wall is no longer the single-qubit phase. It is the transition from an Abelian exchange phase to a genuine non-Abelian protected computational space.
