"""Seeded ensemble of correlated complex nearest-link mediator crosstalk.

One complex residual epsilon is shared by all adjacent-link conversions in a
sample.  This is a fully correlated coherent selectivity error, deliberately
separate from the deterministic scan and from unimplemented independent-link
noise.  Sigma is defined by E|epsilon|^2=sigma^2.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import site_index
from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation
from run_phase8_native_micro_floquet_l3 import ALPHA, TARGET_U0, phase_aligned_distance, polar_unitary
from run_phase8_native_direct_h1_closure import build_micro
from run_phase8_direct_channel_ramp_audit import channel_path


LENGTH, PARTICLE_NUMBER = 4, 2
RATIO, LOGICAL_CYCLES = 0.0125, 16
SIGMAS = (0.001, 0.002, 0.003, 0.005, 0.01)
SAMPLES, SEED = 100, 20260719


def metrics_for_epsilon(epsilon: complex, g: float, detuning: float, period: float,
                        target: np.ndarray, target_states: np.ndarray) -> dict:
    h0, states, frame = build_micro(
        LENGTH, PARTICLE_NUMBER, channel_path(0.0), g, detuning, nearest_link_crosstalk=epsilon,
    )
    h1, states_h1, frame_h1 = build_micro(
        LENGTH, PARTICLE_NUMBER, channel_path(np.pi / 2.0), g, detuning, nearest_link_crosstalk=epsilon,
    )
    if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
        raise RuntimeError("crosstalk segments have incompatible bases")
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("low frame and target basis disagree")
    cycle = expm(-1j * (1.0 - ALPHA) * period * h1) @ expm(-1j * ALPHA * period * h0)
    unitary = np.linalg.matrix_power(cycle, LOGICAL_CYCLES)
    raw = frame.conj().T @ unitary @ frame
    logical = polar_unitary(raw)
    projector = frame @ frame.conj().T
    pa = np.diag([
        -1.0 if sum((int(state) >> site_index(rung, 0)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in low_states
    ])
    return {
        "epsilon_real": float(epsilon.real), "epsilon_imag": float(epsilon.imag),
        "epsilon_abs": float(abs(epsilon)),
        "low_frame_leakage_worst": float(np.linalg.norm((np.eye(h0.shape[0]) - projector) @ unitary @ frame, ord=2) ** 2),
        "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
        "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
    }


def summarize(samples: list[dict], sigma: float) -> dict:
    keys = ("low_frame_leakage_worst", "polar_logical_vs_target_distance", "logical_branch_parity_commutator_normalized")
    passing = [
        row for row in samples
        if row["low_frame_leakage_worst"] < 1e-4
        and row["polar_logical_vs_target_distance"] < 1e-4
        and row["logical_branch_parity_commutator_normalized"] < 1e-4
    ]
    return {
        "sigma_rms_complex_epsilon_over_g": sigma,
        "samples": len(samples),
        "pass_count_strict": len(passing),
        "metrics": {
            key: {
                "mean": float(np.mean([row[key] for row in samples])),
                "variance": float(np.var([row[key] for row in samples])),
                "maximum": float(np.max([row[key] for row in samples])),
            }
            for key in keys
        },
        "raw_samples": samples,
    }


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    period = 4.0 * np.pi / np.sqrt(detuning**2 + 4.0 * g**2)
    target_h0, target_p, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    target = expm(-1j * LOGICAL_CYCLES * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    rng = np.random.default_rng(SEED)
    rows = []
    for sigma in SIGMAS:
        errors = sigma / np.sqrt(2.0) * (rng.normal(size=SAMPLES) + 1j * rng.normal(size=SAMPLES))
        samples = [metrics_for_epsilon(complex(error), g, detuning, period, target, target_states) for error in errors]
        rows.append(summarize(samples, sigma))
    out = {
        "schema": "antler.phase8.direct-channel-crosstalk-ensemble.v1",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "filling": PARTICLE_NUMBER / (2.0 * LENGTH),
            "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g,
            "closure_period": period, "logical_cycles": LOGICAL_CYCLES,
            "total_duration": LOGICAL_CYCLES * period, "sigmas_rms_complex_epsilon_over_g": list(SIGMAS),
            "samples_per_sigma": SAMPLES, "seed": SEED,
            "correlation_model": "one complex epsilon shared by every nearest-link crosstalk conversion in each realization",
        },
        "rows": rows,
        "decision": "Seeded statistical correlated-crosstalk control on the direct deep-SW finite block.",
        "claim_boundary": "This is only one fully correlated, coherent nearest-link crosstalk distribution on L=4,N=2. It excludes independent-link errors, temporal noise, ramps, pulse-angle error, hardware calibration, thermodynamic protection, braiding, non-Abelian statistics, universality and fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "direct_channel_crosstalk_ensemble.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"parameters": out["parameters"], "summary": [
        {key: value for key, value in row.items() if key != "raw_samples"} for row in rows
    ]}, indent=2))


if __name__ == "__main__":
    main()
