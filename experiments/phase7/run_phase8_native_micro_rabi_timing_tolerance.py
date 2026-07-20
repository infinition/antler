"""Timing tolerance around the shortest native virtual-Rabi closure."""
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
    ALPHA, ETA, LENGTH, PARTICLE_NUMBER, RATIO, TARGET_U0, build_micro_h0,
    phase_aligned_distance, polar_unitary,
)


RELATIVE_TIMING_ERRORS = (-0.05, -0.03, -0.01, -0.003, -0.001, 0.0, 0.001, 0.003, 0.01, 0.03, 0.05)


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    closure_period = 4.0 * np.pi / omega
    h_micro, jx_micro, states, frame = build_micro_h0(g, detuning)
    pulse = expm(-1j * ETA * jx_micro)
    h_target, p_target, _, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    projector = frame @ frame.conj().T
    identity = np.eye(h_micro.shape[0], dtype=complex)
    pa = np.diag([
        -1.0 if sum((int(state) >> (2 * rung)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in [states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])]
    ])
    rows = []
    for relative_error in RELATIVE_TIMING_ERRORS:
        period = closure_period * (1.0 + relative_error)
        cycle = (
            pulse.conj().T
            @ expm(-1j * (1.0 - ALPHA) * period * h_micro)
            @ pulse
            @ expm(-1j * ALPHA * period * h_micro)
        )
        target = expm(-1j * period * (ALPHA * h_target + (1.0 - ALPHA) * (p_target.conj().T @ h_target @ p_target)))
        raw = frame.conj().T @ cycle @ frame
        logical = polar_unitary(raw)
        rows.append({
            "relative_period_error": relative_error,
            "period": float(period),
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ cycle @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
            "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
        })
    strict = [row for row in rows if row["low_frame_leakage_worst"] < 1e-4 and row["polar_logical_vs_target_distance"] < 1e-4]
    out = {
        "schema": "antler.phase8.native-micro-rabi-timing-tolerance.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g, "closure_period": closure_period, "alpha": ALPHA, "eta": "pi/2"},
        "rows": rows,
        "strict_local_target": {"leakage_below": 1e-4, "logical_distance_below": 1e-4, "passing_relative_timing_errors": [row["relative_period_error"] for row in strict]},
        "decision": "Local timing-tolerance bracket around the analytically registered shortest virtual-Rabi closure.",
        "claim_boundary": "This timing scan is an ideal-pulse L=3 coherent control, not a material-noise, many-cycle, many-link, thermodynamic, or braiding qualification.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_rabi_timing_tolerance.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
