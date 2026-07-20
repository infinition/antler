"""Exhaustive real crossed-hyperedge baseline before any Phase 7C RL run.

Each hard-core charge-two mediator coherently converts two disjoint pairs on a
four-rung plaquette.  The catalog covers all three rung pairings, bare and
mediator-dressed branch parities, and a relative sign for the two conversions.
It then evaluates every two-mediator candidate at a fixed perturbative point.
No direct four-body term, arbitrary counterterm fit, or continuous optimizer
is used: this is a discrete grammar baseline for deciding whether RL has a
credible local direction to refine.
"""
from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    MicroscopicCandidate2D,
    PairConversionChannel,
    evaluate_local_candidate,
    rail_mode,
)


U, DELTA, G = 20.0, 10.0, 0.5
PARTITIONS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def pair_modes(pair: tuple[int, int], kind: str) -> tuple[int, int]:
    first, second = pair
    rails = {
        "aa": (0, 0), "bb": (1, 1), "ab": (0, 1), "ba": (1, 0),
    }
    left, right = rails[kind]
    return rail_mode(first, left), rail_mode(second, right)


def channel_patterns() -> list[dict]:
    patterns = []
    for partition_index, partition in enumerate(PARTITIONS):
        for parity_family, kinds in (("bare", ("aa", "bb")), ("dressed", ("ab", "ba"))):
            for first_kind, second_kind in product(kinds, repeat=2):
                for sign in (1.0, -1.0):
                    patterns.append({
                        "name": f"{parity_family}_p{partition_index}_{first_kind}_{second_kind}_{'plus' if sign > 0 else 'minus'}",
                        "family": parity_family,
                        "partition": partition,
                        "kinds": (first_kind, second_kind),
                        "relative_sign": sign,
                    })
    return patterns


def materialize(pattern: dict, name_suffix: str) -> PairConversionChannel:
    pair_a, pair_b = pattern["partition"]
    first, second = pair_modes(pair_a, pattern["kinds"][0]), pair_modes(pair_b, pattern["kinds"][1])
    return PairConversionChannel(
        name=f"{pattern['name']}_{name_suffix}", detuning=DELTA, coupling=G, phase=0.0,
        pair_terms=((first[0], first[1], 1.0), (second[0], second[1], pattern["relative_sign"])),
    )


def row_for(first: dict, second: dict) -> dict:
    candidate = MicroscopicCandidate2D(
        mott_u=U, maximum_sw_ratio=0.15,
        channels=(materialize(first, "a"), materialize(second, "b")),
    )
    audit = evaluate_local_candidate(candidate, "XXXX", include_pauli_coefficients=True)
    coefficients = audit["full_traceless_pauli_coefficients"]
    fourbody_square = sum(
        value * value for label, value in coefficients.items() if label != "IIII" and label.count("I") == 0
    )
    xxxx = coefficients["XXXX"]
    fourbody_alignment = float(xxxx * xxxx / fourbody_square) if fourbody_square > 1e-28 else 0.0
    state = audit["state_vector"]
    return {
        "channels": [first["name"], second["name"]],
        "parity_type": candidate.parity_type(),
        "xxxx_coefficient": xxxx,
        "fourbody_xxxx_alignment": fourbody_alignment,
        "fourbody_norm": float(np.sqrt(fourbody_square)),
        "full_operator_alignment": state["target_alignment"],
        "fixed_scale_residual": state["fixed_scale_spectral_algebraic_residual"],
        "minimum_monomer_overlap": state["minimum_monomer_overlap_singular_value"],
        "low_high_gap": state["low_to_high_gap"],
        "hard_failures": audit["hard_failures"],
        "top_unwanted_paulis": audit["top_unwanted_paulis"],
    }


def main() -> None:
    patterns = channel_patterns()
    rows = [row_for(first, second) for first, second in combinations(patterns, 2)]
    signal_rows = [row for row in rows if abs(row["xxxx_coefficient"]) >= 1e-9]
    by_signal = sorted(signal_rows, key=lambda row: (-row["xxxx_coefficient"], -row["fourbody_xxxx_alignment"]))
    by_alignment = sorted(signal_rows, key=lambda row: (-row["fourbody_xxxx_alignment"], -row["xxxx_coefficient"]))
    known = next(
        row for row in rows
        if set(row["channels"]) == {"bare_p0_aa_bb_plus", "bare_p0_bb_aa_plus"}
    )
    out = {
        "schema": "antler.phase7c.crossed-hypergraph-catalog.v1",
        "parameters": {"U": U, "Delta": DELTA, "g": G, "g_over_delta": G / DELTA},
        "grammar": {
            "channel_patterns": len(patterns), "two_mediator_candidates": len(rows),
            "discrete_choices": "three disjoint rung pairings, bare/dressed parity family, rail orientation, relative sign",
            "excluded": "continuous phases, rail hopping, static counterterm fitting, direct four-body interactions",
        },
        "signal_threshold": 1e-9,
        "signal_candidate_count": len(signal_rows),
        "best_sign_correct_signal": by_signal[:12],
        "best_fourbody_xxxx_alignment": by_alignment[:12],
        "known_crossed_reference": known,
        "decision": (
            "This is a pre-RL discrete baseline. A candidate would need both a sign-correct, non-noise XXXX coefficient "
            "and high four-body XXXX alignment before continuous optimization becomes informative."
        ),
        "claim_boundary": (
            "The finite catalog covers real two-mediator crossed hyperedges only. It does not rule out phase-engineered, "
            "higher-mediator, different-degree-of-freedom, or nonperturbative primitives."
        ),
    }
    path = ROOT / "results" / "phase7" / "crossed_hypergraph_catalog.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "grammar": out["grammar"], "signal_candidate_count": len(signal_rows),
        "best_sign_correct_signal": by_signal[:3], "best_fourbody_xxxx_alignment": by_alignment[:3],
    }, indent=2))


if __name__ == "__main__":
    main()
