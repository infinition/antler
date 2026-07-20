"""Parity-group echo for the rejected direct C3 shared-walker composition.

The direct C3 loop has a desired closed XXX return and a lower-order XIX
companion.  The four sign words with even parity

    (+++), (+--), (-+-), (--+)

all retain the triple product but average every one- and two-link sign product
to zero.  This audit applies that control group to the *microscopic* C3
Hamiltonian, downfolds every segment independently, and only then averages.
It therefore tests whether the group argument survives the real shared-matter
virtual paths rather than being assumed at a low-energy link level.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_c3_walker_preflight import (
    PAIR_DETUNING,
    RATIOS,
    build_segment,
    code_indices,
    pauli_coefficients,
    schur_zero,
)


PARITY_EVEN_SIGNS = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))


def effective_for_signs(coupling: float, signs: tuple[int, int, int]) -> np.ndarray:
    h_a, _, positions = build_segment(coupling, echoed=False, link_signs=signs)
    h_b, _, _ = build_segment(coupling, echoed=True, link_signs=signs)
    low = code_indices(positions)
    return 0.5 * (schur_zero(h_a, low) + schur_zero(h_b, low))


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * PAIR_DETUNING
        grouped = sum((effective_for_signs(coupling, signs) for signs in PARITY_EVEN_SIGNS), start=np.zeros((8, 8), dtype=complex)) / len(PARITY_EVEN_SIGNS)
        coefficients = pauli_coefficients(grouped)
        unwanted = {label: value for label, value in coefficients.items() if label not in {"III", "XXX"}}
        label, value = max(unwanted.items(), key=lambda item: abs(item[1]))
        baseline = effective_for_signs(coupling, (1, 1, 1))
        baseline_coefficients = pauli_coefficients(baseline)
        baseline_unwanted = {
            name: coefficient for name, coefficient in baseline_coefficients.items()
            if name not in {"III", "XXX"}
        }
        baseline_label, baseline_value = max(baseline_unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "echoed_xxx_coefficient": coefficients["XXX"],
            "maximum_unwanted_non_scalar_coefficient": abs(value),
            "largest_unwanted_non_scalar_pauli": label,
            "unwanted_over_target": float(abs(value / coefficients["XXX"])),
            "baseline_largest_unwanted_pauli": baseline_label,
            "baseline_largest_unwanted_coefficient": baseline_value,
            "baseline_unwanted_over_target": float(abs(baseline_value / baseline_coefficients["XXX"])),
            "schur_hermiticity_residual": float(np.linalg.norm(grouped - grouped.conj().T, ord="fro")),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.10]
    target_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["echoed_xxx_coefficient"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-c3-parity-echo-audit.v1",
        "parameters": {
            "sign_group": [list(signs) for signs in PARITY_EVEN_SIGNS],
            "group_property": "each sign product s0*s1*s2=+1; all single and pair averages vanish",
            "target": "XXX",
            "segment_count": 8,
        },
        "rows": rows,
        "deep_sw_xxx_power": target_power,
        "decision": "Microscopic parity-group echo audit; selectivity must be read from the serialized unwanted/target ratios rather than inferred from the sign algebra.",
        "claim_boundary": "This is an average of independently downfolded static segments. A positive result would still require a finite-pulse group sequence, pulse leakage/timing/crosstalk audits, a C4 geometry, code protection and defect/braid construction. A negative result would reject this parity group only.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_c3_parity_echo_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
