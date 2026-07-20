"""Detuning scaling of the explicit-mediator L=3 Floquet bridge.

The L=3 composition audit finds physical micromotion leakage.  Here we keep
the desired U0=-g^2/Delta fixed and vary g/Delta, so a fitted suppression law
can distinguish a controlled SW correction from an uncontrolled parasite.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation
from run_phase8_native_micro_floquet_l3 import (
    ALPHA, ETA, LENGTH, PARTICLE_NUMBER, TARGET_U0, build_micro_h0,
    phase_aligned_distance, polar_unitary,
)


REFERENCE_PERIOD = 0.10
PERIOD_SAMPLES = (0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20)
RATIOS = (0.0125, 0.025, 0.05, 0.075)


def main() -> None:
    h_target, p_target, _, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    rows = []
    for ratio in RATIOS:
        detuning = abs(TARGET_U0) / ratio**2
        g = ratio * detuning
        h_micro, jx_micro, _, frame = build_micro_h0(g, detuning)
        pulse = expm(-1j * ETA * jx_micro)
        projector = frame @ frame.conj().T
        samples = []
        for period in PERIOD_SAMPLES:
            cycle = (
                pulse.conj().T
                @ expm(-1j * (1.0 - ALPHA) * period * h_micro)
                @ pulse
                @ expm(-1j * ALPHA * period * h_micro)
            )
            raw = frame.conj().T @ cycle @ frame
            target = expm(-1j * period * (ALPHA * h_target + (1.0 - ALPHA) * (p_target.conj().T @ h_target @ p_target)))
            samples.append({
                "period": period,
                "low_frame_leakage_worst": float(np.linalg.norm((np.eye(h_micro.shape[0]) - projector) @ cycle @ frame, ord=2) ** 2),
                "polar_logical_vs_target_distance": phase_aligned_distance(polar_unitary(raw), target),
                "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
            })
        reference = next(sample for sample in samples if sample["period"] == REFERENCE_PERIOD)
        rows.append({
            "g_over_detuning": ratio,
            "detuning": detuning,
            "g": g,
            "at_reference_period": reference,
            "max_low_frame_leakage_over_period_samples": float(max(sample["low_frame_leakage_worst"] for sample in samples)),
            "max_polar_logical_distance_over_period_samples": float(max(sample["polar_logical_vs_target_distance"] for sample in samples)),
            "samples": samples,
        })
    leakage_slope = float(np.polyfit(
        np.log([row["g_over_detuning"] for row in rows]),
        np.log([row["max_low_frame_leakage_over_period_samples"] for row in rows]),
        1,
    )[0])
    distance_slope = float(np.polyfit(
        np.log([row["g_over_detuning"] for row in rows]),
        np.log([row["max_polar_logical_distance_over_period_samples"] for row in rows]),
        1,
    )[0])
    out = {
        "schema": "antler.phase8.native-micro-detuning-scaling.v2",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "reference_period": REFERENCE_PERIOD, "period_samples": list(PERIOD_SAMPLES), "target_u0": TARGET_U0, "alpha": ALPHA, "eta": "pi/2"},
        "rows": rows,
        "max_leakage_loglog_slope": leakage_slope,
        "max_logical_distance_loglog_slope": distance_slope,
        "decision": "Fixed-U0 detuning scaling of the explicit-mediator Floquet composition-error envelope.",
        "claim_boundary": "This is an ideal-pulse L=3 scaling control, not a many-link or hardware robustness result.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_detuning_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
