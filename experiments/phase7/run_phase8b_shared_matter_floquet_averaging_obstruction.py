"""Audit the noncommuting SW and fast-Floquet averages of the shared-matter echo.

The static compiler result is the arithmetic average
``(H_eff(A) + H_eff(B))/2``.  A realizable rapid A/B drive instead tends to
the Schrieffer--Wolff reduction of the microscopic average ``(H(A)+H(B))/2``.
Those operations need not commute.  This exact finite-block audit records the
difference, then verifies that the physical rapid product converges to the
microscopic average rather than to the desired average-of-SW Hamiltonian.
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
    RATIOS,
    RUNG_HOPPING,
    build_segment,
    code_indices,
    pauli_coefficients,
    schur_effective,
)
from run_phase8b_shared_matter_pulse_closure_audit import phase_aligned_distance


FAST_TAUS = (0.01, 0.002, 0.0005, 0.00025)
FAST_TOTAL_DURATION = 100.0


def largest_non_scalar(coefficients: dict[str, float], allowed: set[str]) -> tuple[str, float]:
    candidates = {key: value for key, value in coefficients.items() if key not in allowed}
    return max(candidates.items(), key=lambda item: abs(item[1]))


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        h_a, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        h_b, states_b, positions_b = build_segment(-RUNG_HOPPING, coupling, -coupling)
        if not np.array_equal(states, states_b) or positions != positions_b:
            raise RuntimeError("segment bases differ")
        low = code_indices(positions)
        effective_a, _, _ = schur_effective(h_a, low)
        effective_b, _, _ = schur_effective(h_b, low)
        average_of_sw = 0.5 * (effective_a + effective_b)
        microscopic_average = 0.5 * (h_a + h_b)
        sw_of_average, _, _ = schur_effective(microscopic_average, low)
        avg_sw_coefficients = pauli_coefficients(average_of_sw)
        sw_avg_coefficients = pauli_coefficients(sw_of_average)
        unwanted_label, unwanted_value = largest_non_scalar(sw_avg_coefficients, {"II", "XX"})
        rows.append({
            "coupling_over_detuning": ratio,
            "average_of_sw_xx": avg_sw_coefficients["XX"],
            "sw_of_microscopic_average_xx": sw_avg_coefficients["XX"],
            "sw_of_microscopic_average_largest_unwanted_pauli": unwanted_label,
            "sw_of_microscopic_average_largest_unwanted": unwanted_value,
            "averaging_noncommutation_frobenius": float(np.linalg.norm(average_of_sw - sw_of_average, ord="fro")),
            "fast_limit_unwanted_over_static_xx": float(abs(unwanted_value / avg_sw_coefficients["XX"])),
        })

    # The first-order Trotter product must approach the full microscopic
    # average.  Measure it directly rather than assuming the asymptotic law.
    ratio = 0.05
    coupling = ratio * DELTA_PAIR
    h_a, _, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
    h_b, _, _ = build_segment(-RUNG_HOPPING, coupling, -coupling)
    h_bar = 0.5 * (h_a + h_b)
    fast_rows = []
    for tau in FAST_TAUS:
        repetitions = int(round(FAST_TOTAL_DURATION / (2.0 * tau)))
        total = 2.0 * tau * repetitions
        exact_fast = np.linalg.matrix_power(expm(-1j * tau * h_b) @ expm(-1j * tau * h_a), repetitions)
        microscopic_average_propagator = expm(-1j * total * h_bar)
        fast_rows.append({
            "coupling_over_detuning": ratio,
            "segment_duration": tau,
            "repetitions": repetitions,
            "total_duration": total,
            "full_space_distance_to_microscopic_average": phase_aligned_distance(exact_fast, microscopic_average_propagator),
        })

    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    static_xx_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["average_of_sw_xx"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-floquet-averaging-obstruction.v1",
        "parameters": {
            "segments": "A=(+J,+g,+g), B=(-J,+g,-g)",
            "fast_taus": list(FAST_TAUS),
            "fast_total_duration": FAST_TOTAL_DURATION,
        },
        "rows": rows,
        "deep_static_xx_power": static_xx_power,
        "fast_limit_rows": fast_rows,
        "decision": (
            "For the registered two-segment shared-matter grammar, average-of-SW and SW-of-microscopic-average differ: "
            "the desired static XX is absent from the rapid-drive microscopic average, which instead retains a lower-order "
            "single-walker term. The exact fast product converges to that microscopic average."
        ),
        "claim_boundary": (
            "This establishes a fast-Floquet obstruction for this A/B sign grammar and its high-frequency limit. It does not "
            "prove a no-go for a different control algebra, a derived kick/frame resource, a higher-order Floquet construction, "
            "a distinct mediator species, a complete walker, a code, defects, fusion or a non-Abelian braid."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_floquet_averaging_obstruction.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
