"""Test the intermediate-time Floquet window of the shared-matter echo.

Integer Rabi returns erase the conditional phase.  The other natural attempt
is rapid A/B switching, with a segment time between the mediator scale and the
target effective scale.  This exact finite-block preflight compares the
physical repeated product with the repeated Schrieffer--Wolff product using a
relative-to-signal metric, so a large unwanted physical rotation cannot be
mistaken for successful conditional transport.
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

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    DELTA_PAIR,
    RUNG_HOPPING,
    build_segment,
    code_indices,
    schur_effective,
)
from run_phase8b_shared_matter_pulse_closure_audit import polar_unitary, phase_aligned_distance


RATIOS = (0.05, 0.025)
TOTAL_DURATIONS = (20.0, 100.0)
SEGMENT_DURATIONS = (0.025, 0.05, 0.10, 0.20, 0.40)
RELATIVE_ERROR_TARGET = 0.1
LEAKAGE_TARGET = 1e-4


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        h_a, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        h_b, states_b, positions_b = build_segment(-RUNG_HOPPING, coupling, -coupling)
        if not np.array_equal(states, states_b) or positions != positions_b:
            raise RuntimeError("Floquet segment bases disagree")
        low = code_indices(positions)
        h_eff_a, _, _ = schur_effective(h_a, low)
        h_eff_b, _, _ = schur_effective(h_b, low)
        frame = np.zeros((len(states), len(low)), dtype=complex)
        frame[low, np.arange(len(low))] = 1.0
        for desired_total in TOTAL_DURATIONS:
            for duration in SEGMENT_DURATIONS:
                repetitions = int(round(desired_total / (2.0 * duration)))
                total = float(2.0 * duration * repetitions)
                physical_cycle = expm(-1j * duration * h_b) @ expm(-1j * duration * h_a)
                target_cycle = expm(-1j * duration * h_eff_b) @ expm(-1j * duration * h_eff_a)
                physical = np.linalg.matrix_power(physical_cycle, repetitions)
                target = np.linalg.matrix_power(target_cycle, repetitions)
                raw = frame.conj().T @ physical @ frame
                logical = polar_unitary(raw)
                target_signal = phase_aligned_distance(target, np.eye(target.shape[0], dtype=complex))
                physical_signal = phase_aligned_distance(logical, np.eye(logical.shape[0], dtype=complex))
                error = phase_aligned_distance(logical, target)
                leakage = float(np.linalg.norm(physical @ frame - frame @ raw, ord=2) ** 2)
                rows.append({
                    "coupling_over_detuning": ratio,
                    "requested_total_duration": desired_total,
                    "segment_duration": duration,
                    "repetitions": repetitions,
                    "actual_total_duration": total,
                    "piecewise_sw_signal_distance_from_identity": target_signal,
                    "physical_polar_signal_distance_from_identity": physical_signal,
                    "relative_piecewise_sw_error": float(error / target_signal) if target_signal > 1e-15 else None,
                    "low_frame_leakage_worst": leakage,
                    "passes_signal_screen": bool(error / target_signal < RELATIVE_ERROR_TARGET and leakage < LEAKAGE_TARGET),
                })
    output = {
        "schema": "antler.phase8b.shared-matter-fast-floquet-preflight.v1",
        "parameters": {
            "coupling_ratios": list(RATIOS),
            "total_durations": list(TOTAL_DURATIONS),
            "segment_durations": list(SEGMENT_DURATIONS),
            "screen": {"relative_piecewise_sw_error": RELATIVE_ERROR_TARGET, "low_frame_leakage_worst": LEAKAGE_TARGET},
        },
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_signal_screen"]],
        "decision": "No registered fast A/B switching row realizes the SW conditional signal: every row must be judged by relative signal error and leakage, not by bare population return.",
        "claim_boundary": "This rejects the stated abrupt two-segment fast-switching window only. It does not refute a pulse sequence with additional derived controls, smooth switch compensation, dressed-state preparation, or a different microscopic resource grammar.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_fast_floquet_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
