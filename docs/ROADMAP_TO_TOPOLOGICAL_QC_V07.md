# ANTLER v0.7  -  Roadmap from the present gate to topological quantum computing

## Where the project stands

ANTLER currently has a high-fidelity **Abelian logical phase primitive**. In the sequential digital and strongly localized limit, the exchange phase is explained by an exact oriented string count. This is a strong result, but an Abelian phase gate alone is not a topological quantum computer.

## Wall 1  -  Many-body derivation of the finite-localization correction

Current evidence:

- the ideal localized path gives exactly `-theta`;
- finite-depth error correlates with delocalized weight;
- error fits approximately `|D|^-2.496` over the tested range.

Required:

1. perform a Schrieffer–Wolff reduction around the trapped transfer manifold;
2. derive the first virtual-tail correction to the link phase;
3. explain whether the observed exponent is asymptotic or pre-asymptotic;
4. add finite-size scaling in `L` and coupling scaling in `J1/J2`, `Jperp`.

Success criterion: an analytic bound or expansion for

`Delta phi_odd + theta` and `P_leak`.

## Wall 2  -  Composition and inverse-cycle algebra

Required tests:

- one exchange followed by its inverse must give identity;
- two identical exchanges must give twice the phase modulo `2pi`;
- error and leakage must scale predictably with cycle count;
- echo sequences must cancel even components without hiding odd coherent errors.

The current multi-cycle run log is incomplete and cannot be claimed.

## Wall 3  -  Non-Abelian code space

The current braid representation is one-dimensional/Abelian. A topological computer needs a degenerate fusion space of dimension greater than one, with matrix-valued exchanges.

Candidate routes:

1. extend the ladder with an internal synthetic dimension and at least three exchangeable defects;
2. construct parafermion-like or Fibonacci-inspired effective defects;
3. use measurement-assisted fusion channels;
4. introduce `N>=3` with a mobile mediator while preserving a protected logical doublet.

Required test:

`B1 B2 != B2 B1`

on the same protected subspace, with convergence and leakage audits.

## Wall 4  -  Universal logical operations

A protected `Z(theta)` family is not universal on its own.

Required:

- an independent protected non-commuting single-qubit operation, or a measurement-assisted equivalent;
- a two-logical-qubit entangling operation;
- process fidelities and worst-case leakage for arbitrary superpositions.

Target demonstration:

- dense single-qubit control in `SU(2)`;
- one entangling gate such as controlled phase;
- universality argument for the resulting set.

## Wall 5  -  Noise beyond static common-mode disorder

Completed: static diagonal disorder through 20% with two seeds per point.

Still required:

- time-dependent noise with controlled spectra;
- independent noise between `+theta`, `-theta`, exchange and null-control branches;
- timing jitter and ramp-shape errors;
- hopping calibration errors;
- correlated spatial noise;
- loss/dephasing channels using an open-system treatment.

Report distributions and worst-case tails, not only averages.

## Wall 6  -  Scalable encoding and error correction

Required:

- multiple logical qubits in one architecture;
- local addressability and crosstalk bounds;
- preparation and readout protocol;
- syndrome strategy or passive error-suppression argument;
- logical error scaling with system size;
- threshold estimate under a realistic noise model.

Without this, the work remains a gate primitive rather than a computer architecture.

## Wall 7  -  Experimental Hamiltonian mapping

The rung-major Jordan–Wigner construction defines a specific density-dependent correlated-hopping model.

A hardware proposal must specify:

- platform: cold atoms, Rydberg array, photonic synthetic dimension, superconducting analogue simulator, or another controlled lattice;
- implementation of density-dependent Peierls phases;
- SSH dimerization and rung coupling;
- moving/digital traps;
- preparation and readout of the cat logical states;
- accessible timescales relative to decoherence;
- parameter values corresponding to the validated dimensionless regime.

## Suggested paper sequence

### Paper A  -  Ready after finite-size and perturbative correction work

**Protection–control trade-offs and localized statistical phase gates in a correlated-hopping SSH ladder**

Core: no-go, dynamic rescue, full gate audit, digital exchange, localization scaling, disorder and local lemma.

### Paper B  -  Non-Abelian extension

Core: expanded protected fusion space, two non-commuting exchanges, matrix-valued holonomies.

### Paper C  -  Architecture

Core: multiple logical qubits, entangling gate, hardware mapping and fault-tolerance analysis.
