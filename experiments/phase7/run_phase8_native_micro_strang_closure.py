"""Symmetric (Strang) native Floquet cycle at simultaneous Rabi closure.

The first-order cycle accumulates a coherent O(T^2) Floquet discrepancy.  At
alpha=1/2, T=8*pi/Omega makes the outer H0/4 pieces and central H0/2 piece
integer virtual Rabi cycles.  This compares that symmetric sequence directly
with the first-order cycle at the same total period and checks n=1,2,4,8
composition.
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
    ALPHA, ETA, LENGTH, PARTICLE_NUMBER, RATIO, TARGET_U0, build_micro_h0,
    phase_aligned_distance, polar_unitary,
)


CYCLES = (1, 2, 4, 8)


def metrics(unitary: np.ndarray, target: np.ndarray, frame: np.ndarray, projector: np.ndarray, identity: np.ndarray) -> dict:
    raw = frame.conj().T @ unitary @ frame
    return {
        "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
        "polar_logical_vs_target_distance": phase_aligned_distance(polar_unitary(raw), target),
        "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
    }


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 8.0 * np.pi / omega
    h_micro, jx_micro, _, frame = build_micro_h0(g, detuning)
    pulse = expm(-1j * ETA * jx_micro)
    outer = expm(-1j * period * h_micro / 4.0)
    central = expm(-1j * period * h_micro / 2.0)
    strang = outer @ pulse.conj().T @ central @ pulse @ outer
    first_order = pulse.conj().T @ central @ pulse @ expm(-1j * period * h_micro / 2.0)
    h_target, p_target, _, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    target_cycle = expm(-1j * period * (ALPHA * h_target + (1.0 - ALPHA) * (p_target.conj().T @ h_target @ p_target)))
    projector = frame @ frame.conj().T
    identity = np.eye(h_micro.shape[0], dtype=complex)
    rows = []
    for count in CYCLES:
        target = np.linalg.matrix_power(target_cycle, count)
        rows.append({
            "cycles": count,
            "first_order": metrics(np.linalg.matrix_power(first_order, count), target, frame, projector, identity),
            "strang": metrics(np.linalg.matrix_power(strang, count), target, frame, projector, identity),
        })
    out = {
        "schema": "antler.phase8.native-micro-strang-closure.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g, "omega_pair_mediator": omega, "period": period, "alpha": ALPHA, "eta": "pi/2", "closure": "outer free segment: one virtual Rabi cycle; central free segment: two"},
        "rows": rows,
        "decision": "Symmetric-cycle finite-block comparison against the same H_eff target.",
        "claim_boundary": "A successful ideal-pulse Strang control would still require finite-bandwidth, timing-noise, many-link and thermodynamic audits before any ANTLER phase or computing claim.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_strang_closure.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
