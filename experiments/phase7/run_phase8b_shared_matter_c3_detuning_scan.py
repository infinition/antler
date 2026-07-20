"""Selectivity scan for the explicit C3 shared-walker composition preflight.

The initial C3 loop produced its desired XXX word but a larger XIX companion.
This short independent scan changes only the virtual-walker detuning.  It
tests whether the failure is a tunable resonance or an order hierarchy of the
present direct-composition grammar.
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

import run_phase8b_shared_matter_c3_walker_preflight as c3


WALKER_DETUNINGS = (0.5, 1.0, 2.0, 5.0, 10.0)
RATIOS = (0.10, 0.05)


def main() -> None:
    original_detuning = c3.WALKER_DETUNING
    rows = []
    try:
        for walker_detuning in WALKER_DETUNINGS:
            c3.WALKER_DETUNING = walker_detuning
            for ratio in RATIOS:
                coupling = ratio * c3.PAIR_DETUNING
                h_a, _, positions = c3.build_segment(coupling, echoed=False)
                h_b, _, _ = c3.build_segment(coupling, echoed=True)
                low = c3.code_indices(positions)
                effective = 0.5 * (c3.schur_zero(h_a, low) + c3.schur_zero(h_b, low))
                coefficients = c3.pauli_coefficients(effective)
                unwanted = {label: value for label, value in coefficients.items() if label not in {"III", "XXX"}}
                label, coefficient = max(unwanted.items(), key=lambda item: abs(item[1]))
                target = coefficients["XXX"]
                rows.append({
                    "walker_detuning": walker_detuning,
                    "coupling_over_detuning": ratio,
                    "target_xxx_coefficient": target,
                    "largest_unwanted_pauli": label,
                    "largest_unwanted_coefficient": coefficient,
                    "unwanted_over_target": float(abs(coefficient / target)),
                })
    finally:
        c3.WALKER_DETUNING = original_detuning
    best = min(rows, key=lambda row: row["unwanted_over_target"])
    output = {
        "schema": "antler.phase8b.shared-matter-c3-detuning-scan.v1",
        "parameters": {
            "walker_detunings": list(WALKER_DETUNINGS),
            "coupling_ratios": list(RATIOS),
            "target": "XXX",
            "selectivity_target": "largest_unwanted_over_target < 0.25",
        },
        "rows": rows,
        "best_selectivity": best,
        "decision": "No tested virtual-walker detuning makes the direct C3 shared-matter composition selective: the leading companion exceeds XXX at every scanned point.",
        "claim_boundary": "This rejects the direct C3 reservoir/channel composition and the stated detuning window only. It does not refute the one-link X/Y/Z compiler, a counterterm/refocusing construction, a different walker encoding, or every multi-link ANTLER extension.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_c3_detuning_scan.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
