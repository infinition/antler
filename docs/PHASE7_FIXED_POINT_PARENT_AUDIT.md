# Phase 7 -- independent audit of the fixed-point parent

## Scope

This audit implements only the open-boundary fixed-point parent `H_fix` from
`PHASE7_STATE_TO_HAMILTONIAN_DERIVATION.md`. It does not test the proposed
Schrieffer--Wolff microscopic Hamiltonian, a periodic-boundary extension,
braiding, or a non-Abelian phase.

The implementation uses the established weighted-charge basis of the
charge-two mediator ladder and verifies it exactly at `L=4`.

## What the matrix audit confirms

With `U=4`, `Delta=2`, and `J=1`, the proposed OBC parent has the expected
two-dimensional zero-energy space and the expected finite-size spectral gaps:

| L | neutral gap above doublet | addition energy | removal energy |
|---:|---:|---:|---:|
| 4 | 1.0 | 4.0 | 4.0 |
| 6 | 1.0 | 4.0 | 4.0 |
| 8 | 1.0 | 4.0 | 4.0 |

All supplied local terms commute numerically in the `L=4,Q=4` audit. The
code frame is annihilated by the parent to machine precision, and weighted
charge plus both branch parities commute with the OBC parent.

The U(1) observation in the derivation is also retained as an analytic filter:
under its stated spectral hypotheses, a charge-changing edge operator cannot
have a vanishing projected commutator in an incompressible system.

## Corrections required before any promotion

### The cell term is not a projector

`C_j=(q_cell(j)-1)^2` is Hermitian, positive semidefinite, local, and
commuting, but it is not idempotent: its eigenvalues include `4` and `9`.
The measured `L=4` maximum term idempotency residual is `16.97`. Therefore
the fixed-point model is a **commuting positive-semidefinite,
frustration-free parent**, not a strict commuting-projector Hamiltonian unless
`C_j` is replaced by an actual local occupancy projector and the construction
is re-audited.

### The alleged edge operator is a bulk local logical operator

For the supplied code frame,

\[
P X_0 P=P X_2 P=\mathrm{diag}(1,-1).
\]

Both `X_0` and a bulk `X_2` have zero leakage and zero projected commutator.
Consequently `X_0` is not localized at an edge: every rung carries the same
local logical action. The result is the standard Ising ordered/cat doublet,
not a protected edge mode.

The full local-probe audit gives a non-scalar norm `1.4142` for `X_0`. Only
after excluding branch-parity-odd probes does the restricted set `{Z_0,n^a_0}`
give zero. Thus the valid statement is symmetry-restricted protection against
an idealized parity-preserving noise algebra, not full local
indistinguishability or intrinsic topological protection.

### PBC has not yet been defined

The delivered mode layout is OBC, with `L-1` bond mediators. A PBC claim needs
an extra wraparound mediator, a periodic cell-charge definition, and a new
algebra audit. It cannot be inferred from the OBC implementation.

## Status

The fixed point is retained as a useful, exactly soluble **symmetry-restricted
Ising/cat benchmark** with a real local charge gap. It is not promoted to a
native topological memory, a locally indistinguishable code, or a physical
edge-mode platform. The microscopic mapping and the requested adversarial
a revised microscopic construction remains a prerequisite for a new candidate.

## Reproducibility

- `antler/phase7_ising_parent.py`;
- `experiments/phase7/run_phase7_fixed_point_parent_audit.py`;
- `results/phase7/fixed_point_parent_independent_audit.json`.
