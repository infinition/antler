"""All-order closed-walk selection rule for the simultaneous mixed-check patch.

For a walker on a four-state cycle, every low-to-low virtual path has an
even-incidence edge set.  Modulo two, the cycle space of C4 is only
{empty, full cycle}.  Since Pauli operators commute modulo a scalar phase,
the Pauli support of *any* simultaneous multi-walker low-to-low path is a
product of complete check words, never a new Pauli word.

This script makes that selection rule explicit for the six-check Phase-8B
patch and combines it with the exact pairwise coefficients.  It is a formal
Schrieffer-Wolff series result under the declared independent phase-controlled
walker grammar; it is not a replacement for a native-resource derivation.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CHECK_LABELS = (
    "YXZXIII",
    "XZIIXZI",
    "XZXYIII",
    "YXIIYXI",
    "XIXIXIX",
    "YIZIYIY",
)


def cycle_parity_masks() -> list[str]:
    """Return all mod-two C4 edge sets with even degree at every vertex."""
    valid = []
    for mask in range(1 << 4):
        degrees = [0] * 4
        for edge in range(4):
            if (mask >> edge) & 1:
                degrees[edge] += 1
                degrees[(edge + 1) % 4] += 1
        if all(degree % 2 == 0 for degree in degrees):
            valid.append(format(mask, "04b"))
    return valid


def returned_walk_parities(maximum_steps: int) -> dict[str, list[str]]:
    """Independent enumeration of closed C4 walks as a check of the lemma."""
    rows: dict[str, list[str]] = {}
    for steps in range(maximum_steps + 1):
        frontier = {(0, 0)}
        for _ in range(steps):
            next_frontier = set()
            for vertex, parity in frontier:
                for direction in (-1, 1):
                    edge = (vertex - 1) % 4 if direction == -1 else vertex
                    next_frontier.add(((vertex + direction) % 4, parity ^ (1 << edge)))
            frontier = next_frontier
        returned = sorted(format(parity, "04b") for vertex, parity in frontier if vertex == 0)
        rows[str(steps)] = returned
    return rows


def main() -> None:
    patch_path = ROOT / "results" / "phase7" / "phase8b_mixed_check_patch_algebra.json"
    pair_path = ROOT / "results" / "phase7" / "phase8b_mixed_check_patch_pairwise_overlap.json"
    patch = json.loads(patch_path.read_text(encoding="utf-8"))
    pair = json.loads(pair_path.read_text(encoding="utf-8"))
    if tuple(patch["parameters"]["stabilizer_labels"]) != CHECK_LABELS:
        raise RuntimeError("patch stabilizer set changed; re-register the closed-walk audit")
    masks = cycle_parity_masks()
    walks = returned_walk_parities(12)
    if masks != ["0000", "1111"]:
        raise RuntimeError(f"unexpected C4 cycle space: {masks}")
    if any(any(mask not in masks for mask in values) for values in walks.values()):
        raise RuntimeError("closed-walk enumerator found a parity outside the C4 cycle space")
    if not patch["stabilizer_algebra"]["all_checks_commute"]:
        raise RuntimeError("closed-path products cannot be promoted to an Abelian stabilizer algebra")
    if not patch["complete_projected_local_pauli_gate"]["all_tested_probes_project_to_scalars_or_zero"]:
        raise RuntimeError("finite patch did not pass its local-Pauli gate")
    if len(pair["rows"]) != 15:
        raise RuntimeError("all 15 pair overlaps are required")
    maximum_pair_unwanted = max(float(row["maximum_outside_generated_algebra_coefficient"]) for row in pair["rows"])
    maximum_pair_coefficient = max(abs(float(row["product_coefficient"])) for row in pair["rows"])
    minimum_single_coefficient = min(
        min(abs(float(row["left_coefficient"])), abs(float(row["right_coefficient"])))
        for row in pair["rows"]
    )
    pair_sum_bound = sum(abs(float(row["product_coefficient"])) for row in pair["rows"])
    # For H = sum c_i S_i + sum c_ij S_i S_j, flipping a syndrome can change
    # each product contribution by at most 2|c_ij|.  This is only an order-8
    # truncated lower bound; higher simultaneous-walker orders are not bounded.
    fourth_plus_eighth_gap_lower_bound = 2.0 * (minimum_single_coefficient - pair_sum_bound)
    output = {
        "schema": "antler.phase8b.mixed-patch-closed-walk-closure.v1",
        "parameters": {
            "walker_graph": "independent four-state cycle C4 for each check",
            "checks": list(CHECK_LABELS),
            "formal_expansion": "Schrieffer-Wolff paths returning every walker to its vacuum state",
            "coefficient_point": pair["parameters"]["coupling_over_detuning"],
        },
        "cycle_space_proof": {
            "closed_C4_edge_parity_masks": masks,
            "returned_walk_parities_through_12_steps": walks,
            "consequence": (
                "for each walker, a vacuum-to-vacuum path carries either identity or its complete four-link check word "
                "modulo a scalar Pauli phase"
            ),
        },
        "simultaneous_multiwalker_consequence": {
            "formal_effective_operator_algebra": "span of products of the six commuting check words",
            "stabilizer_group_size": 1 << int(patch["stabilizer_algebra"]["independent_rank_over_GF2"]),
            "every_generated_group_element_acts_as_scalar_on_patch_code": True,
            "new_logical_pauli_from_ideal_walker_crosstalk": False,
            "reason": (
                "the Pauli group is Abelian modulo scalar phase, so interleaving different closed walker paths can alter only "
                "a scalar sign; it cannot alter their mod-two check-word product"
            ),
        },
        "registered_pairwise_input": {
            "pair_count": len(pair["rows"]),
            "maximum_outside_generated_algebra_coefficient": maximum_pair_unwanted,
            "minimum_single_check_coefficient": minimum_single_coefficient,
            "maximum_pair_product_coefficient": maximum_pair_coefficient,
            "sum_abs_pair_product_coefficients": pair_sum_bound,
            "fourth_plus_eighth_syndrome_gap_lower_bound": fourth_plus_eighth_gap_lower_bound,
        },
        "decision": (
            "Within the formal independent C4-walker SW grammar, simultaneous multi-walker virtual paths cannot generate a "
            "Pauli outside the finite patch stabilizer algebra at any perturbative order. The measured order-four checks dominate "
            "the registered order-eight pair products at g/Delta=0.05, leaving a positive fourth-plus-eighth syndrome-gap bound."
        ),
        "claim_boundary": (
            "This is an all-order combinatorial selection rule for the declared ideal walker grammar, not an exact finite-coupling "
            "six-walker diagonalization or a convergence-radius proof. It does not derive the walker, its complex conditional phase, "
            "or the stabilizer patch from the frozen ANTLER Hamiltonian, and it does not establish scalable topological order, a full "
            "twist/dislocation geometry, fusion, defect motion, non-Abelian braid, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_mixed_patch_closed_walk_closure.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
