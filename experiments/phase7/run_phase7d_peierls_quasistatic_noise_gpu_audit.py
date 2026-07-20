"""Statistical local-control audit for quasi-static Peierls switching errors.

Each realization uses one common phase-plateau error and one signed timing
imbalance for the whole four-pulse schedule.  This represents slow global
control drift, not independent link noise or cycle-to-cycle jitter.  The full
472-state block is propagated exactly on CUDA for every realization.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import numpy as np
import torch
import torch


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import FRAME, pair_gate, projected
from run_phase7d_peierls_switch_error_gpu_audit import (
    DEVICE, LOCAL_TARGET, direct_pair_gate, metrics, pulse_eigens, to_gpu,
)


SAMPLES_PER_LEVEL = 50
SIGMAS = (0.01, 0.03, 0.05, 0.10, 0.20)
SEED = 20260719
OUTPUT = ROOT / "results" / "phase7" / "peierls_quasistatic_noise_gpu_audit.json"


def subspace_average_overlap_fidelity(logical: np.ndarray, reference: np.ndarray) -> float:
    """Haar average target-state overlap on the 16D monomer frame, including loss."""
    dimension = logical.shape[0]
    survival = float(np.trace(logical.conj().T @ logical).real)
    coherent_overlap = abs(np.trace(reference.conj().T @ logical)) ** 2
    return float((survival + coherent_overlap) / (dimension * (dimension + 1)))


def scalar_summary(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(array.mean()),
        "variance": float(array.var(ddof=0)),
        "minimum": float(array.min()),
        "maximum": float(array.max()),
    }


def write_output(rows: list[dict], state: str, started: float) -> None:
    out = {
        "schema": "antler.phase7d.peierls-quasistatic-noise-gpu-audit.v1",
        "method": "exact 472-state CUDA spectral propagation, complex128, seeded common quasi-static control offsets",
        "device": torch.cuda.get_device_name(0),
        "samples_per_noise_level": SAMPLES_PER_LEVEL,
        "seed": SEED,
        "noise_model": {
            "phase_plateau_error": "common normal offset with std = sigma*pi radians on every pi plateau in the four-pulse schedule",
            "signed_time_imbalance": "common normal offset with std = sigma, preserving each subcycle duration",
            "correlation": "the two offsets are independent; each offset is common to all legs, links and subcycles within one realization",
        },
        "registered_local_target": LOCAL_TARGET,
        "rows": rows,
        "run_state": state,
        "wall_clock_seconds": time.perf_counter() - started,
        "claim_boundary": (
            "This is a local quasi-static global-control noise test only. It does not establish bandwidth, independent link noise, "
            "cycle-to-cycle noise, calibration, full-ladder dynamics, a protected phase, edge modes, a 2D code, braid, "
            "non-Abelian statistics, universality or fault tolerance."
        ),
    }
    OUTPUT.write_text(json.dumps(out, indent=2), encoding="utf-8")


def main() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for this registered statistical audit")
    torch.cuda.synchronize()
    started = time.perf_counter()
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    frame = to_gpu(FRAME)
    standard_draws = np.random.default_rng(SEED).normal(size=(SAMPLES_PER_LEVEL, 2))
    rows: list[dict] = []
    write_output(rows, "in_progress", started)

    for sigma in SIGMAS:
        fidelities: list[float] = []
        leakages: list[float] = []
        deviations: list[float] = []
        parities: list[float] = []
        singular_mins: list[float] = []
        passed: list[bool] = []
        phase_errors_deg: list[float] = []
        timing_errors_percent: list[float] = []
        for phase_z, time_z in standard_draws:
            phase_error = float(phase_z * sigma * np.pi)
            timing_imbalance = float(time_z * sigma)
            eigens = pulse_eigens(phase_error)
            even = direct_pair_gate(frame.clone(), (0, 2), eigens, timing_imbalance)
            complete = direct_pair_gate(even, (1,), eigens, timing_imbalance)
            output = complete.detach().cpu().numpy()
            logical = projected(output)
            record = metrics(output, reference)
            fidelities.append(subspace_average_overlap_fidelity(logical, reference))
            leakages.append(float(record["monomer_leakage"]))
            deviations.append(float(record["logical_deviation_from_zero_leg_schedule"]))
            parities.append(max(float(record["logical_parity_a_residual"]), float(record["logical_parity_b_residual"])))
            singular_mins.append(float(record["logical_singular_value_min"]))
            passed.append(bool(record["passes_registered_local_target"]))
            phase_errors_deg.append(float(np.rad2deg(phase_error)))
            timing_errors_percent.append(100.0 * timing_imbalance)
        row = {
            "sigma_fraction": sigma,
            "sigma_percent": 100.0 * sigma,
            "phase_error_degrees": scalar_summary(phase_errors_deg),
            "time_imbalance_percent": scalar_summary(timing_errors_percent),
            "average_subspace_fidelity": scalar_summary(fidelities),
            "monomer_leakage": scalar_summary(leakages),
            "logical_deviation": scalar_summary(deviations),
            "worst_logical_parity_residual": scalar_summary(parities),
            "logical_singular_value_min": scalar_summary(singular_mins),
            "registered_target_pass_fraction": float(np.mean(passed)),
            "registered_target_pass_count": int(sum(passed)),
        }
        rows.append(row)
        write_output(rows, "in_progress", started)
        print(f"completed sigma={100.0 * sigma:.1f}% ({len(rows)}/{len(SIGMAS)})", flush=True)

    write_output(rows, "complete", started)
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)), "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
