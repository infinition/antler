# Phase 7 -- contract for a theory-supplied parent Hamiltonian

Phase 7 separates analytic construction from numerical auditing. A candidate
arrives from the theory side as explicit local terms and is not promoted by a
small spectral split alone.

## Required theory payload

For each finite length and charge sector, the implementation needs:

1. a rung-major basis and weighted charge assigned to every mode;
2. a complete matrix or matrix-free action for every local term `Pi_j`;
3. the physical support of each `Pi_j`;
4. diagonal labels for total charge and every claimed parity symmetry;
5. a normalized proposed code frame `G` with columns spanning the claimed
   low-energy code;
6. a growing family of local probes, including bulk and both boundaries;
7. each proposed edge operator with its finite support;
8. the boundary-anchor term separated from the bulk parent.

## Mandatory gates

### Local algebra

If the candidate is called a commuting-projector parent, every local term must
be Hermitian, idempotent, and commute with every overlapping term. A failure of
any of these tests is not silently tolerated: the theory must then supply a
different gap proof for a noncommuting parent.

### Exact symmetries

The finite matrix must commute with weighted U(1) charge and every claimed
branch parity. A fixed-charge diagonalization is not evidence of an exact
microscopic U(1) symmetry by itself.

### Protection, not merely edge conservation

For every edge proposal, the audit reports

\[
\epsilon_{\rm edge}=\frac{\|(1-P)[H,O_{\rm edge}]P\|_F}{\|O_{\rm edge}P\|_F}.
\]

It also reports logical leakage and the non-scalar projection of every local
probe. A code fails if an individual boundary probe distinguishes its logical
states, even when `epsilon_edge` is exactly zero. The Phase 7 self-test is a
commuting-projector classical doublet with this exact false-positive pattern.

### Complete physical local algebra

For each proposed mode support, the audit must not use a curated list of
densities alone. It enumerates every local weighted-charge-conserving
occupation matrix unit and both Hermitian quadratures of every allowed local
transition. Rail modes carry their Fock signs; declared charge-two mediators
are hard-core commuting modes. This necessarily includes the local rail
transfer `a^dagger b + b^dagger a` whenever it is physically allowed.

The maximum non-scalar code projection over this complete basis must vanish
with the claimed protection parameter. A restriction to a symmetry-preserving
subalgebra may be reported as a separate SPT/SSB diagnostic, but cannot pass
the full-local-code gate.

For a charge-frozen two-dimensional edge-qubit parent, an exact binary
stabilizer certificate may complement the dense audit: the implementation must
prove that the projected full physical charge-conserving algebra is generated
by the local Pauli algebra, and enumerate every Pauli below the claimed code
distance. This validates only the specified parent. It is not evidence that a
low-body ANTLER Hamiltonian generates that parent.

### Solver-proposed microscopic candidates

An external optimizer may propose numeric parameters only through a registered
local-term grammar. The environment must reject a proposal before spectrum
calculation if it violates weighted U(1), a declared mediator-parity label,
finite support, or the registered perturbative window. Its reward must include
the target Pauli coefficient, all leading unwanted Pauli coefficients, monomer
capture, symmetry residuals and the eliminated-state gap. A spectral distance
alone is insufficient and no local reward promotes a candidate to a tiled
two-dimensional Hamiltonian without an independent controlled reduction.

### Independent spectral questions

The following must be reported separately:

- gap above the requested code multiplicity within a fixed charge/parity
  sector;
- addition and removal energies in neighbouring weighted-charge sectors;
- dependence of the neutral gap on length;
- dependence of `epsilon_edge` on *both* length and support, with one axis
  fixed while the other changes.

A global charging term is an external resource, not a substitute for a local
bulk charge gap. It must be reported as such.

### Boundary anchor

The anchor is tested once absent and once present. The audit reports which
symmetries it changes, its effect on the charge-sector spectrum, its effect on
local distinguishability, and its effect on `epsilon_edge`. A boundary term
that merely writes a locally readable qubit is rejected.

## Provided implementation

`antler/phase7_parent_audit.py` provides dense small-system checks for:

- local projector Hermiticity, idempotency and pairwise commutators;
- diagonal-symmetry commutators;
- code-projected edge metrics;
- local indistinguishability;
- sector-resolved finite-size gaps.

`antler/phase7_local_algebra.py` adds exhaustive charge-conserving local probe
generation. Its self-test rejects the archived Ising/cat benchmark already on
one rung, identifying the physical rail-transfer quadrature as the local
logical operator.

`experiments/phase7/run_phase7_harness_selftest.py` is the negative software
control. It must remain passing whenever the harness evolves.

`experiments/phase7/run_phase7_full_local_algebra_selftest.py` is the
complete-local-algebra negative control and must also remain passing.

## Claim boundary

Passing this contract could establish a finite-size candidate worthy of
scaling. It does not establish intrinsic topological order, non-Abelian braid
statistics, universal computation, or experimental feasibility.
