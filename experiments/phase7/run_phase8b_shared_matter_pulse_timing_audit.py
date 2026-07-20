"""Timing-error bracket for the registered shared-matter pulse closure."""
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

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    DELTA_PAIR,
    RUNG_HOPPING,
    build_segment,
    code_indices,
    schur_effective,
)
from run_phase8b_shared_matter_pulse_closure_audit import (
    phase_aligned_distance,
    polar_unitary,
)


RATIO = 0.025
TIMING_OFFSETS = (-0.004, -0.003, -0.002, -0.001, 0.0, 0.001, 0.002, 0.003, 0.004)
STRICT_LEAKAGE = 1e-5
STRICT_DISTANCE = 1e-4


def main() -> None:
    coupling = RATIO * DELTA_PAIR
    h_a, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
    h_b, states_b, positions_b = build_segment(-RUNG_HOPPING, coupling, -coupling)
    if not np.array_equal(states, states_b) or positions != positions_b:
        raise RuntimeError("the echoed segments do not share a basis")
    low = code_indices(positions)
    h_eff_a, _, _ = schur_effective(h_a, low)
    h_eff_b, _, _ = schur_effective(h_b, low)
    frame = np.zeros((len(states), len(low)), dtype=complex)
    frame[low, np.arange(len(low))] = 1.0
    projector = frame @ frame.conj().T
    identity = np.eye(len(states), dtype=complex)
    nominal_duration = float(2.0 * np.pi / np.sqrt(DELTA_PAIR**2 + 4.0 * coupling**2))
    rows = []
    for offset in TIMING_OFFSETS:
        duration = nominal_duration * (1.0 + offset)
        cycle = expm(-1j * duration * h_b) @ expm(-1j * duration * h_a)
        raw = frame.conj().T @ cycle @ frame
        logical = polar_unitary(raw)
        target = expm(-1j * duration * h_eff_b) @ expm(-1j * duration * h_eff_a)
        leakage = float(np.linalg.norm((identity - projector) @ cycle @ frame, ord=2) ** 2)
        distance = phase_aligned_distance(logical, target)
        rows.append({
            "relative_segment_timing_offset": offset,
            "segment_duration": duration,
            "low_frame_leakage_worst": leakage,
            "polar_distance_to_piecewise_sw": distance,
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
            "passes_strict_local_target": bool(leakage < STRICT_LEAKAGE and distance < STRICT_DISTANCE),
        })
    output = {
        "schema": "antler.phase8b.shared-matter-pulse-timing-audit.v1",
        "parameters": {
            "coupling_over_detuning": RATIO,
            "nominal_rabi_multiplier_per_segment": 1,
            "nominal_segment_duration": nominal_duration,
            "strict_local_target": {
                "low_frame_leakage_worst": STRICT_LEAKAGE,
                "polar_distance_to_piecewise_sw": STRICT_DISTANCE,
            },
        },
        "rows": rows,
        "passing_offsets": [row["relative_segment_timing_offset"] for row in rows if row["passes_strict_local_target"]],
        "decision": "A registered local timing bracket for the one-cycle shared-matter Rabi closure.",
        "claim_boundary": "This one-block deterministic timing bracket is not a pulse-bandwidth, noise, many-link, code, defect, fusion or braid qualification.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_pulse_timing_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
