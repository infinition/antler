"""Multi-cycle composition of the shortest native virtual-Rabi closure."""
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


CYCLES = (1, 2, 4, 8)


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 4.0 * np.pi / omega
    h_micro, jx_micro, _, frame = build_micro_h0(g, detuning)
    pulse = expm(-1j * ETA * jx_micro)
    one_cycle = (
        pulse.conj().T
        @ expm(-1j * (1.0 - ALPHA) * period * h_micro)
        @ pulse
        @ expm(-1j * ALPHA * period * h_micro)
    )
    h_target, p_target, _, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    target_cycle = expm(-1j * period * (ALPHA * h_target + (1.0 - ALPHA) * (p_target.conj().T @ h_target @ p_target)))
    projector = frame @ frame.conj().T
    identity = np.eye(h_micro.shape[0], dtype=complex)
    rows = []
    for count in CYCLES:
        full = np.linalg.matrix_power(one_cycle, count)
        target = np.linalg.matrix_power(target_cycle, count)
        raw = frame.conj().T @ full @ frame
        logical = polar_unitary(raw)
        rows.append({
            "cycles": count,
            "total_time": float(count * period),
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ full @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        })
    out = {
        "schema": "antler.phase8.native-micro-rabi-composition.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g, "closure_period": period, "alpha": ALPHA, "eta": "pi/2"},
        "rows": rows,
        "leakage_loglog_cycle_slope": float(np.polyfit(
            np.log([row["cycles"] for row in rows]),
            np.log([row["low_frame_leakage_worst"] for row in rows]),
            1,
        )[0]),
        "logical_distance_loglog_cycle_slope": float(np.polyfit(
            np.log([row["cycles"] for row in rows]),
            np.log([row["polar_logical_vs_target_distance"] for row in rows]),
            1,
        )[0]),
        "decision": "Multi-cycle coherent composition of the shortest explicit-mediator Rabi-closure primitive.",
        "claim_boundary": "This is not a logical topological-gate composition: it is a finite L=3 microscopic Floquet/H_eff bridge under ideal pulses.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_rabi_composition.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
