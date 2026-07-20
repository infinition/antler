# ANTLER  -  Walls between the current result and a topological computer

## Passed

1. **State → logical subspace.** The full dressed doublet is propagated.
   Worst-case leakage is measured through singular values, not a single state.
2. **Coherent gate extraction.** Polar decomposition separates loss from the
   closest logical unitary. Off-diagonal logical mixing is ~1e-6 or lower.
3. **Exchange-vs-trivial control.** The phase is extracted by exchange minus
   round trip and odd-in-theta subtraction.
4. **Localized braid count.** The sequential compact path has exactly one JW
   string crossing; the matched round trip has zero net crossing.
5. **Distance deformation.** R=3,4,5 gives the same digital-gate slope to
   roughly 3.5e-4 at finite trap depth.
6. **Static calibrated disorder.** The gate remains closed and coherent through
   sigma_mu=0.5; a breakdown is observed by sigma_mu=1.0.

## In progress

7. **Integrator convergence at moderate trap depth.** Fine-dt D=4 run.
8. **Gate composition.** Verify two exchange cycles give twice the phase and
   bounded leakage.
9. **Linearity.** Verify theta=0.6 and 0.9 in the compact protocol.
10. **Disorder threshold.** Locate breakdown between sigma_mu=0.5 and 1.0.

## Fundamental walls not yet passed

11. **Passive topological memory.** The cat states |LL> and |RR> are locally
    distinguishable. A local potential can measure or split the logical state.
    The current object is a robust exchange-controlled phase gate, not yet a
    topologically protected memory.
12. **Non-Abelian action.** The present braid is Abelian: it produces a scalar
    phase. All such braids commute. Universal all-topological computation needs
    a degenerate fusion space with noncommuting braid matrices, or a hybrid
    architecture with non-topological resources and error correction.
13. **Two-qubit entanglement.** A controlled phase between two encoded qubits
    has not been implemented.
14. **Scaling.** Four-particle sparse dynamics, initialization/readout, thermal
    errors, time-dependent noise and fault-tolerance thresholds remain open.

## Recommended architecture after Phase 4.3

Store information in a spatially separated, locally indistinguishable sector;
use a third mobile ancilla to activate the correlated statistical channel only
during a gate. In parallel, explore a two-species/synthetic-dimension extension
with matrix-valued hopping, seeking noncommuting Wilson loops rather than a
single Abelian phase.
