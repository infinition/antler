# Reproducibility notes

## Model convention

Sites are ordered rung-by-rung as `k = 2*i + sigma`. The fractional Jordan–Wigner string is therefore a definition of the correlated-hopping ladder model. Results should not be described as representation-independent anyon physics unless an explicit unitary equivalence is shown.

## Numerical hierarchy

- Exact diagonalization is used for static spectra and logical-frame construction.
- Sparse propagation was used for early shuttle exploration.
- The final Phase 4.1 gate audit uses a Strang split propagator:
  - exact hopping exponential, precomputed per theta;
  - exact midpoint diagonal-potential half steps;
  - timestep convergence explicitly checked.

## Reference environment

- Python 3.10+
- NumPy
- SciPy
- Matplotlib

## Saved-result validation

Run:

```bash
python scripts/validate_saved_results.py
```

This checks the archived reference run against the declared leakage, coherent-mixing, and fidelity thresholds.
