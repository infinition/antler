"""Stroboscopic virtual-mediator closure scan for the native L=3 bridge.

For alpha=1/2, choosing T = 4*pi*m/sqrt(Delta^2+4g^2) makes each free
half-period an integer pair--mediator Rabi cycle in the isolated-link limit.
This scan tests whether that analytically motivated closure suppresses the
explicit-mediator leakage in the composed three-rung protocol.
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


M_VALUES = tuple(range(1, 9))


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    h_micro, jx_micro, _, frame = build_micro_h0(g, detuning)
    pulse = expm(-1j * ETA * jx_micro)
    h_target, p_target, _, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    projector = frame @ frame.conj().T
    identity = np.eye(h_micro.shape[0], dtype=complex)
    rows = []
    for m in M_VALUES:
        period = 4.0 * np.pi * m / omega
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
            "rabi_cycles_per_half_period": m,
            "period": float(period),
            "half_period_omega_over_2pi": float((period * omega / 2.0) / (2.0 * np.pi)),
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ cycle @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        })
    best = min(rows, key=lambda row: row["low_frame_leakage_worst"])
    out = {
        "schema": "antler.phase8.native-micro-rabi-closure.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g, "omega_pair_mediator": omega, "alpha": ALPHA, "eta": "pi/2"},
        "rows": rows,
        "best_registered_closure": best,
        "decision": "Analytically registered virtual-Rabi closure scan in the composed explicit-mediator Floquet block.",
        "claim_boundary": "A local coherent closure point does not establish robustness to timing errors, pulse bandwidth, many-link crosstalk, a thermodynamic phase, braid, non-Abelian statistics, universality, or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_rabi_closure.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
