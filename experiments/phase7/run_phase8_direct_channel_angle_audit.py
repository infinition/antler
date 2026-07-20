"""Static calibration-offset audit for the direct rotated pair channels.

The second direct Floquet segment needs a pi/2 rotation in the two-particle
channel space.  This records how a coherent endpoint-angle offset affects the
same deep-SW, 16-cycle logical duration used by the finite-ramp control.
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


LENGTH, PARTICLE_NUMBER = 3, 2
RATIO, LOGICAL_CYCLES = 0.0125, 16
ANGLE_OFFSETS = (-0.05, -0.03, -0.01, -0.003, 0.0, 0.003, 0.01, 0.03, 0.05)


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 4.0 * np.pi / omega
    h0_micro, states, frame = build_micro(LENGTH, PARTICLE_NUMBER, channel_path(0.0), g, detuning)
    target_h0, target_p, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("low frame and target basis disagree")
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    target = expm(-1j * LOGICAL_CYCLES * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    projector = frame @ frame.conj().T
    identity = np.eye(h0_micro.shape[0], dtype=complex)
    pa = np.diag([
        -1.0 if sum((int(state) >> site_index(rung, 0)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in low_states
    ])
    rows = []
    for offset in ANGLE_OFFSETS:
        h1_micro, states_h1, frame_h1 = build_micro(LENGTH, PARTICLE_NUMBER, channel_path(np.pi / 2.0 + offset), g, detuning)
        if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
            raise RuntimeError("miscalibrated channel changed the basis")
        cycle = expm(-1j * (1.0 - ALPHA) * period * h1_micro) @ expm(-1j * ALPHA * period * h0_micro)
        unitary = np.linalg.matrix_power(cycle, LOGICAL_CYCLES)
        raw = frame.conj().T @ unitary @ frame
        logical = polar_unitary(raw)
        rows.append({
            "channel_angle_offset_rad": offset,
            "channel_angle_offset_deg": float(np.rad2deg(offset)),
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
            "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        })
    passing = [
        row["channel_angle_offset_rad"] for row in rows
        if row["low_frame_leakage_worst"] < 1e-4
        and row["polar_logical_vs_target_distance"] < 1e-4
        and row["logical_branch_parity_commutator_normalized"] < 1e-4
    ]
    out = {
        "schema": "antler.phase8.direct-channel-angle-audit.v1",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0,
            "g_over_detuning": RATIO, "detuning": detuning, "g": g,
            "closure_period": period, "logical_cycles": LOGICAL_CYCLES,
            "total_duration": LOGICAL_CYCLES * period,
            "ideal_direct_h1_channel_angle_rad": float(np.pi / 2.0),
        },
        "rows": rows,
        "strict_local_target": {
            "leakage_below": 1e-4,
            "distance_below": 1e-4,
            "logical_parity_commutator_below": 1e-4,
            "passing_channel_angle_offsets_rad": passing,
        },
        "decision": "Coherent direct-channel endpoint calibration control at fixed deep-SW logical duration.",
        "claim_boundary": "This models a static coherent channel-angle offset only. It does not model stochastic amplitude/phase noise, bandwidth, crosstalk, thermal effects, a thermodynamic phase, a braid, non-Abelian statistics, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "direct_channel_angle_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
