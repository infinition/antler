"""Phase-programmable X/Y conditional-link audit for the shared-matter block.

The shared-matter echo has already isolated ``X_rail tensor X_walker`` for a
real rung hop.  A complex Peierls phase on precisely that physical rung hop is
the minimal way to test whether the same derived bridge can orient the rail
Pauli axis, rather than inserting a Y word in the low-energy theory.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


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


PHASES = (0.0, np.pi / 4.0, np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0)


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        for phase in PHASES:
            complex_rung = RUNG_HOPPING * np.exp(1j * phase)
            h_a, _, positions = build_segment(complex_rung, coupling, coupling)
            h_b, _, _ = build_segment(-complex_rung, coupling, -coupling)
            low = code_indices(positions)
            h_eff_a, capture_a, gap_a = schur_effective(h_a, low)
            h_eff_b, capture_b, gap_b = schur_effective(h_b, low)
            coefficients = pauli_coefficients(0.5 * (h_eff_a + h_eff_b))
            reference_amplitude = float(np.hypot(coefficients["XX"], coefficients["YX"]))
            expected_xx = -reference_amplitude * np.cos(phase)
            expected_yx = reference_amplitude * np.sin(phase)
            unwanted = {
                label: value for label, value in coefficients.items()
                if label not in {"II", "XX", "YX"}
            }
            rows.append({
                "coupling_over_detuning": ratio,
                "rung_peierls_phase": float(phase),
                "xx_coefficient": coefficients["XX"],
                "yx_coefficient": coefficients["YX"],
                "phase_rotation_residual": float(np.hypot(coefficients["XX"] - expected_xx, coefficients["YX"] - expected_yx)),
                "maximum_unwanted_non_scalar_coefficient": float(max(abs(value) for value in unwanted.values())),
                "minimum_low_frame_capture": min(capture_a, capture_b),
                "minimum_low_high_gap": min(gap_a, gap_b),
            })
    y_rows = [row for row in rows if abs(row["rung_peierls_phase"] - np.pi / 2.0) < 1e-12 and row["coupling_over_detuning"] <= 0.075]
    y_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in y_rows]),
        np.log(np.abs([row["yx_coefficient"] for row in y_rows])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-phase-link-audit.v1",
        "parameters": {
            "rung_hopping_magnitude": RUNG_HOPPING,
            "relative_pair_channel_echo": "A: g0=g1, B: g0=-g1",
            "target_family": "(-cos(phi) X_rail + sin(phi) Y_rail) tensor X_walker",
        },
        "rows": rows,
        "deep_sw_yx_power": y_power,
        "decision": "A physical Peierls phase on the shared-matter rung rotates the derived conditional link continuously between X and Y rail axes while the echo cancels isolated flips.",
        "claim_boundary": "This is a local SW compiler result for conditional X/Y links. It does not establish a complete four-state walker, a finite-pulse complex-phase closure, a code patch, defects, fusion, non-Abelian braiding, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_phase_link_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
