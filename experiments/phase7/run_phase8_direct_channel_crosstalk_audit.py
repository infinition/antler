"""Nearest-link mediator crosstalk audit for the direct Phase-8 bridge.

Each mediator nominally addresses one pair bond.  This control adds a residual
amplitude epsilon*g to the immediately adjacent pair bonds of that same
mediator.  It is a microscopic selectivity error, not an effective logical
perturbation.  The audit holds U0, the deep SW ratio and logical duration
fixed, then checks leakage, target distance and branch-parity restoration.
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


SYSTEMS = ((3, 2), (4, 2))
RATIO, LOGICAL_CYCLES = 0.0125, 16
EPSILONS = (0.0, 0.0003, 0.001, 0.003, 0.01, 0.03)


def audit_one(length: int, particle_number: int, epsilon: float, g: float, detuning: float, period: float) -> dict:
    h0, states, frame = build_micro(
        length, particle_number, channel_path(0.0), g, detuning, nearest_link_crosstalk=epsilon,
    )
    h1, states_h1, frame_h1 = build_micro(
        length, particle_number, channel_path(np.pi / 2.0), g, detuning, nearest_link_crosstalk=epsilon,
    )
    if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
        raise RuntimeError("crosstalk segments have incompatible bases")
    target_h0, target_p, target_states, _ = build_h0_and_rotation(length, particle_number, TARGET_U0)
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("bare low frame and target basis disagree")
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    target = expm(-1j * LOGICAL_CYCLES * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    cycle = expm(-1j * (1.0 - ALPHA) * period * h1) @ expm(-1j * ALPHA * period * h0)
    unitary = np.linalg.matrix_power(cycle, LOGICAL_CYCLES)
    raw = frame.conj().T @ unitary @ frame
    logical = polar_unitary(raw)
    projector = frame @ frame.conj().T
    pa = np.diag([
        -1.0 if sum((int(state) >> site_index(rung, 0)) & 1 for rung in range(length)) % 2 else 1.0
        for state in low_states
    ])
    return {
        "L": length, "N": particle_number, "filling": particle_number / (2.0 * length),
        "microscopic_dimension": int(h0.shape[0]), "bare_low_frame_dimension": int(frame.shape[1]),
        "nearest_link_crosstalk_over_g": epsilon,
        "low_frame_leakage_worst": float(np.linalg.norm((np.eye(h0.shape[0]) - projector) @ unitary @ frame, ord=2) ** 2),
        "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
        "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
        "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
    }


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    period = 4.0 * np.pi / np.sqrt(detuning**2 + 4.0 * g**2)
    rows = [
        audit_one(length, number, epsilon, g, detuning, period)
        for length, number in SYSTEMS for epsilon in EPSILONS
    ]
    passing = [
        {"L": row["L"], "N": row["N"], "epsilon": row["nearest_link_crosstalk_over_g"]}
        for row in rows
        if row["low_frame_leakage_worst"] < 1e-4
        and row["polar_logical_vs_target_distance"] < 1e-4
        and row["logical_branch_parity_commutator_normalized"] < 1e-4
    ]
    out = {
        "schema": "antler.phase8.direct-channel-crosstalk-audit.v1",
        "parameters": {
            "systems": [{"L": length, "N": number} for length, number in SYSTEMS],
            "target_u0": TARGET_U0, "g_over_detuning": RATIO,
            "detuning": detuning, "g": g, "closure_period": period,
            "logical_cycles": LOGICAL_CYCLES, "total_duration": LOGICAL_CYCLES * period,
            "nearest_link_crosstalk_scan_over_g": list(EPSILONS),
        },
        "rows": rows,
        "strict_local_target": {
            "leakage_below": 1e-4, "distance_below": 1e-4,
            "logical_parity_commutator_below": 1e-4, "passing_rows": passing,
        },
        "decision": "Microscopic nearest-link pair-conversion selectivity audit at deep SW.",
        "claim_boundary": "This is a deterministic nearest-link coherent crosstalk model on finite L=3,4 blocks. It does not establish stochastic noise tolerance, material selectivity, a thermodynamic phase, an edge qubit, a braid, non-Abelian statistics, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "direct_channel_crosstalk_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
