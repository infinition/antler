"""Does a deeper SW regime improve *accumulated* direct-bridge leakage?

At fixed U0=-g^2/Delta, decreasing r=g/Delta shortens the Rabi-closed
Floquet period as r^2.  A one-cycle leakage improvement is therefore not a
useful control claim: a fixed logical evolution needs more cycles.  This
audit compares the closest integer-closure sequence durations on L=6,N=3.
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

from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation
from run_phase8_native_micro_floquet_l3 import ALPHA, TARGET_U0, phase_aligned_distance, polar_unitary
from run_phase8_native_direct_h1_closure import build_micro, rotated_channels


LENGTH, PARTICLE_NUMBER = 6, 3
RATIOS_AND_CYCLES = ((0.05, 1), (0.025, 4), (0.0125, 16))


def main() -> None:
    channels_h1, factorization = rotated_channels()
    channels_h0 = np.asarray(((1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0)), dtype=complex)
    target_h0, target_p, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    rows = []
    for ratio, cycles in RATIOS_AND_CYCLES:
        detuning = abs(TARGET_U0) / ratio**2
        g = ratio * detuning
        omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
        period = 4.0 * np.pi / omega
        h0_micro, states, frame = build_micro(LENGTH, PARTICLE_NUMBER, channels_h0, g, detuning)
        h1_micro, states_h1, frame_h1 = build_micro(LENGTH, PARTICLE_NUMBER, channels_h1, g, detuning)
        if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
            raise RuntimeError("direct segments have incompatible bases")
        low_states = np.asarray([
            states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])
        ], dtype=np.int64)
        if not np.array_equal(low_states, target_states):
            raise RuntimeError("low frame does not match target physical basis")
        cycle = expm(-1j * (1.0 - ALPHA) * period * h1_micro) @ expm(-1j * ALPHA * period * h0_micro)
        unitary = np.linalg.matrix_power(cycle, cycles)
        target = expm(-1j * (cycles * period) * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
        raw = frame.conj().T @ unitary @ frame
        projector = frame @ frame.conj().T
        rows.append({
            "g_over_detuning": ratio,
            "cycles_at_closure": cycles,
            "detuning": detuning,
            "g": g,
            "closure_period": period,
            "total_duration": cycles * period,
            "low_frame_leakage_worst": float(np.linalg.norm((np.eye(h0_micro.shape[0]) - projector) @ unitary @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(polar_unitary(raw), target),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        })
    reference_duration = rows[0]["total_duration"]
    for row in rows:
        row["duration_over_r005_reference"] = row["total_duration"] / reference_duration
    log_ratio = np.log([row["g_over_detuning"] for row in rows])
    leakage_power = float(np.polyfit(log_ratio, np.log([row["low_frame_leakage_worst"] for row in rows]), 1)[0])
    distance_power = float(np.polyfit(log_ratio, np.log([row["polar_logical_vs_target_distance"] for row in rows]), 1)[0])
    out = {
        "schema": "antler.phase8.native-direct-h1-fixed-time-scaling.v1",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "filling": PARTICLE_NUMBER / (2.0 * LENGTH),
            "target_u0": TARGET_U0, "alpha": ALPHA,
            "protocol": "integer virtual-Rabi closure with the nearest number of cycles to retain the r=0.05 total duration",
        },
        "rotated_channel_factorization_frobenius_residual": factorization["factorization_frobenius_residual"],
        "closest_fixed_duration_loglog_powers": {
            "accumulated_leakage_vs_g_over_detuning": leakage_power,
            "polar_logical_distance_vs_g_over_detuning": distance_power,
        },
        "rows": rows,
        "decision": "Fixed-logical-duration control: lower r is credited only if the accumulated, not one-cycle, leakage improves.",
        "claim_boundary": "This is exact L=6 finite-block evidence only. Integer Rabi closures change the total duration by the reported small mismatch. Neither a positive nor negative result proves a thermodynamic topological phase or hardware implementability.",
    }
    path = ROOT / "results" / "phase7" / "native_direct_h1_fixed_time_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
