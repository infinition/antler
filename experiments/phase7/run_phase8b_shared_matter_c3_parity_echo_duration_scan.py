"""Duration scan for the finite-pulse C3 parity-echo compiler.

The first equal-duration group pulse shows a speed/leakage tension.  This scan
keeps the same microscopic eight-segment group, changes only the integer
virtual-Rabi multiplier per segment, and reports the accumulated phase and a
linear leakage-budget diagnostic for a pi/4 logical rotation.  It is a search
for a concrete timing sweet spot, not an extrapolation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_c3_parity_echo_audit import PARITY_EVEN_SIGNS
from run_phase8b_shared_matter_c3_walker_preflight import (
    PAIR_DETUNING,
    build_segment,
    code_indices,
    pauli_coefficients,
    schur_zero,
)
from run_phase8b_shared_matter_pulse_closure_audit import phase_aligned_distance, polar_unitary


RATIOS = (0.05, 0.0375, 0.025)
MULTIPLIERS = (1, 2, 4)
STRICT_LEAKAGE = 1e-4
STRICT_DISTANCE = 1e-4
MAX_LINEAR_BUDGET = 1.0


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * PAIR_DETUNING
        physical_segments: list[np.ndarray] = []
        effective_segments: list[np.ndarray] = []
        frame: np.ndarray | None = None
        grouped_effective = np.zeros((8, 8), dtype=complex)
        for signs in PARITY_EVEN_SIGNS:
            h_a, _, positions = build_segment(coupling, echoed=False, link_signs=signs)
            h_b, _, _ = build_segment(coupling, echoed=True, link_signs=signs)
            low = code_indices(positions)
            if frame is None:
                frame = np.zeros((h_a.shape[0], len(low)), dtype=complex)
                frame[low, np.arange(len(low))] = 1.0
            h_eff_a, h_eff_b = schur_zero(h_a, low), schur_zero(h_b, low)
            physical_segments.extend((h_a, h_b))
            effective_segments.extend((h_eff_a, h_eff_b))
            grouped_effective += h_eff_a + h_eff_b
        if frame is None:
            raise RuntimeError("C3 frame was not initialized")
        grouped_effective /= len(effective_segments)
        xxx = pauli_coefficients(grouped_effective)["XXX"]
        rabi_duration = float(2.0 * np.pi / np.sqrt(PAIR_DETUNING**2 + 4.0 * coupling**2))
        for multiplier in MULTIPLIERS:
            duration = multiplier * rabi_duration
            evolved = frame.copy()
            target = np.eye(frame.shape[1], dtype=complex)
            for physical, effective in zip(physical_segments, effective_segments):
                evolved = expm_multiply((-1j * duration) * csr_matrix(physical), evolved)
                target = expm(-1j * duration * effective) @ target
            raw = frame.conj().T @ evolved
            logical = polar_unitary(raw)
            leakage = float(np.linalg.norm(evolved - frame @ raw, ord=2) ** 2)
            distance = phase_aligned_distance(logical, target)
            phase = float(duration * len(physical_segments) * abs(xxx))
            groups = float(np.pi / 4.0 / phase)
            budget = float(groups * leakage)
            rows.append({
                "coupling_over_detuning": ratio,
                "rabi_multiplier_per_segment": multiplier,
                "segment_duration": duration,
                "xxx_phase_magnitude_per_group": phase,
                "ideal_groups_to_pi_over_4": groups,
                "low_frame_leakage_worst": leakage,
                "polar_distance_to_piecewise_sw": distance,
                "linear_leakage_budget_estimate_to_pi_over_4": budget,
                "passes_registered_practical_screen": bool(
                    leakage < STRICT_LEAKAGE and distance < STRICT_DISTANCE and budget < MAX_LINEAR_BUDGET
                ),
            })
    output = {
        "schema": "antler.phase8b.shared-matter-c3-parity-echo-duration-scan.v1",
        "parameters": {
            "coupling_ratios": list(RATIOS),
            "rabi_multipliers": list(MULTIPLIERS),
            "screen": {
                "leakage_per_group": STRICT_LEAKAGE,
                "polar_distance_to_piecewise_sw": STRICT_DISTANCE,
                "linear_leakage_budget_to_pi_over_4": MAX_LINEAR_BUDGET,
            },
        },
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_registered_practical_screen"]],
        "decision": "No registered integer-Rabi duration supplies a practical C3 gate: the screen has no passing row, and the companion signal audit shows that the physical operation is nearly scalar at these closure timings.",
        "claim_boundary": "The linear leakage budget is a screening diagnostic, not a composition proof. Passing it would still require explicit multi-group composition, timing/noise/crosstalk audits and a C4 code construction.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_c3_parity_echo_duration_scan.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
