# Phase 7 -- Ising edge-versus-bulk calibration

## Scope

The revised Phase 7 Ising/cat construction is not a topological candidate. It
is nevertheless a useful controlled calibration for the projected-commutator
diagnostic: an edge strong-zero-mode recurrence can be compared with the same
one-sided recurrence started in the bulk.

The calculation uses the open transverse-field Ising reduction

\[
H=\frac{J}{2}\sum_{j=0}^{L-2}(1-X_jX_{j+1})+h\sum_{j=0}^{L-1}Z_j,
\qquad J=1,\quad h=0.1,\quad r=2h/J=0.2,
\]

at `L=8`.

## Sign convention correction

For this displayed Hamiltonian, the left strong-zero-mode recurrence is

\[
\Gamma_w=\mathcal N_w\sum_{n=0}^{w-1}(-r)^n
 \left(\prod_{m=0}^{n-1}Z_m\right)X_n,
\qquad
\mathcal N_w=\sqrt{\frac{1-r^2}{1-r^{2w}}}.
\]

The alternating sign is required to cancel the `2ihY_0` commutator at the
first step. It does not change the previously stated norm formula:

\[
\epsilon(\Gamma_w)=2(J/2)r^w
\sqrt{\frac{1-r^2}{1-r^{2w}}},
\]

up to the finite-length correction.

## Matrix result

| support w | left epsilon | formula | same recurrence in bulk |
|---:|---:|---:|---:|
| 1 | 0.2000000 | 0.2000000 | 0.2000000 |
| 2 | 0.0392232 | 0.0392232 | 0.2037721 |
| 3 | 0.0078386 | 0.0078386 | 0.2001907 |
| 4 | 0.0015677 | 0.0015677 | 0.2000084 |

The formula agrees to better than `1.3e-6` relative error through support 4.
The one-sided bulk recurrence retains nontrivial logical action but does not
become quasi-conserved as its support grows.

## Boundary

This validates a useful edge-versus-bulk diagnostic on the symmetry-restricted
Ising benchmark only. It does not remove the bulk local logical operator at
the fixed point, restore full local indistinguishability, establish a native
microscopic realization, or support a braid claim.

## Sources

- `experiments/phase7/run_phase7_ising_edge_bulk_calibration.py`;
- `results/phase7/ising_edge_bulk_calibration.json`.
