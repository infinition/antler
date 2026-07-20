# Incomplete or rejected runs

These files are preserved for audit but must not be cited as completed evidence.

## Incomplete executions

- Phase 4.7 deep-limit extrapolation: all current path and composition runs,
  plus the `D=12, dt=0.125` refinement, are serialized. What remains
  incomplete is a promoted deep physical fit whose retained points have a
  separately bounded timestep error; see `docs/PHASE4_7_CAMPAIGN_STATUS.md`.

- `results/raw/phase4_3/D80_seq_T40.log`
- `results/raw/phase4_3/composition/cycles2.log` and any related clean log: no final serialized two-cycle metrics
- noise logs above the curated `sigma=0.20` dataset, including partial `0.50`, `0.70`, `0.75`, `0.85`, and `1.00` runs
- partial adiabatic-response tracking logs
- `experiments/phase7/run_phase7d_peierls_phase_ramp_audit.py`: smooth finite
  Peierls-ramp stress test; the full 472D execution exceeded the local run
  limit before serializing final metrics. It is not evidence for or against
  the direct Peierls echo. Its CUDA replacement completed the corresponding
  2/4/8-segment numerical audit; see `PHASE7D_PEIERLS_RAMP_GPU_AUDIT.md`.

### Phase 8 published-point MPS convergence

- `experiments/phase7/run_phase8_tenpy_published_l18_convergence.py` was
  launched for `L=18,N=12` with `chi=256,384,512`, but the complete series
  reached its execution limit before its JSON was serialized. No partial
  result is promoted. The retained `chi=384` size-scaling point reaches its
  bond-dimension ceiling, so it is not by itself a long-size convergence
  certificate. A rerun should checkpoint every `chi` or resume from a saved
  MPS.

## Rejected diagnostics

- `results/phase4_4/response_baseline.json`: coarse propagation with large final leakage; retained only to document the failed direct-response integration strategy.
- early recapture files with coarse timesteps: useful for diagnosing integrator error, not for physical claims.
- simultaneous digital rung-swap experiments: rejected because the composite transfer creates severe leakage.

## Policy

A result is promoted to the curated results only if:

1. the run completes and serializes its final metrics;
2. the logical frame is orthonormal;
3. convergence or an appropriate null control is available;
4. leakage is explicitly reported;
5. the result is not contradicted by a finer integration step.
