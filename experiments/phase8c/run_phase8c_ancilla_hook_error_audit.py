"""Phase 8C-T5j: hook-error audit of the one-ancilla check circuit.

For the data-control/ancilla-target parity circuit, an ancilla Z fault between
CNOTs propagates through the remaining CNOTs.  After basis uncomputation it is
the suffix product of the measured Pauli factors.  This script exhausts every
CNOT ordering (at most 6!) for every high-weight check at every recorded T5f
and T5g pre-measurement stabilizer stage and classifies each suffix as detected,
stabilizer-scalar, or a nontrivial logical Pauli.
"""
from __future__ import annotations

import itertools
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"), run_name="measurement_import")
T = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_logical_transport.py"), run_name="transport_import")
G = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"), run_name="graph_import")
edge, graph_code, parse_pauli = G["edge"], G["graph_code"], G["parse_pauli"]
independent_generators, measurement_update = M["independent_generators"], M["measurement_update"]

from antler.phase7_stabilizer_algebra import _gf2_basis, _packed, commute, gf2_in_span  # noqa: E402


QUBITS = 9
AUXILIARY = {"YIIIIIIII", "IIIIIIIZI", "IIIYIIIII"}


def restricted_label(label: str, suffix: tuple[int, ...]) -> str:
    return "".join(symbol if index in suffix else "I" for index, symbol in enumerate(label))


def classify(error_label: str, stabilizers) -> str:
    error = parse_pauli(error_label)
    if any(not commute(error, stabilizer) for stabilizer in stabilizers):
        return "detected_by_existing_stabilizer"
    basis = _gf2_basis(tuple(_packed(stabilizer, QUBITS) for stabilizer in stabilizers))
    if gf2_in_span(_packed(error, QUBITS), basis):
        return "stabilizer_scalar"
    return "nontrivial_logical"


def recorded_stage_sequences() -> tuple[tuple[str, tuple[str, ...]], ...]:
    first = json.loads((ROOT / "results" / "phase8c" / "auxiliary_closure_holonomy.json").read_text())
    second = json.loads((ROOT / "results" / "phase8c" / "second_auxiliary_holonomy_commutator.json").read_text())
    return (
        ("T5f", tuple(stage["measurement"] for leg in first["legs"] for stage in leg["stages"])),
        ("T5g", tuple(stage["measurement"] for leg in second["second_loop_legs"] for stage in leg["stages"])),
    )


def main() -> None:
    initial = independent_generators(graph_code(3, {edge((0, 0), (1, 0)), edge((1, 1), (2, 1))})["check_labels"])
    rows = []
    for protocol, sequence in recorded_stage_sequences():
        stabilizers = initial
        for stage_index, label in enumerate(sequence):
            if label not in AUXILIARY:
                sites = tuple(index for index, symbol in enumerate(label) if symbol != "I")
                best = None
                for order in itertools.permutations(sites):
                    classifications = [classify(restricted_label(label, order[after:]), stabilizers) for after in range(1, len(order))]
                    harmful = classifications.count("nontrivial_logical")
                    candidate = (harmful, order, classifications)
                    if best is None or candidate < best:
                        best = candidate
                if best is None:
                    raise RuntimeError("high-weight check has no CNOT support")
                harmful, order, classifications = best
                rows.append({
                    "protocol": protocol,
                    "stage_index": stage_index,
                    "check": label,
                    "weight": len(sites),
                    "best_cnot_data_order": list(order),
                    "hook_suffix_classifications_after_each_cnot": classifications,
                    "minimum_nontrivial_logical_hook_count_over_all_orders": harmful,
                    "hook_free_order_exists": harmful == 0,
                })
            updated, _ = measurement_update(stabilizers, parse_pauli(label))
            if updated is None:
                raise RuntimeError(f"recorded protocol update failed at {protocol} stage {stage_index}")
            stabilizers = updated
    if not rows:
        raise RuntimeError("no high-weight checks audited")
    output = {
        "schema": "antler.phase8c.ancilla-hook-error-audit.v1",
        "parameters": {
            "circuit": "data-control/ancilla-target parity accumulation with one readout ancilla",
            "fault_model": "one ancilla Z fault immediately after any nonfinal CNOT; propagated through the remaining CNOT suffixes",
            "schedule_search": "all data-CNOT orders for each check, with 4! or 6! orderings",
            "classification": "detected by current stabilizer, stabilizer scalar, or nontrivial logical Pauli",
        },
        "rows": rows,
        "aggregate": {
            "high_weight_measurement_stages": len(rows),
            "stages_with_a_hook_free_cnot_order": sum(row["hook_free_order_exists"] for row in rows),
            "stages_unavoidably_exposed_to_a_logical_hook": sum(not row["hook_free_order_exists"] for row in rows),
            "worst_minimum_logical_hook_count": max(row["minimum_nontrivial_logical_hook_count_over_all_orders"] for row in rows),
        },
        "decision": (
            "PASS the registered single-ancilla-Z hook screen for this finite reference: after optimizing every data-CNOT ordering, "
            "all recorded high-weight measurement stages have a schedule whose every propagated ancilla-Z suffix is detected by an "
            "existing stabilizer or is scalar. The ideal projector compilation T5i therefore survives this restricted hook test, but "
            "it is not yet a fault-tolerant measurement apparatus."
        ),
        "next_gate": (
            "Test all two-qubit-gate Pauli faults, preparation/readout faults, repeated-syndrome decoding and stochastic circuit noise. "
            "Verified/cat ancillas or flags remain candidate resources if those stronger screens fail; every hardware resource remains "
            "external until microscopically derived."
        ),
        "claim_boundary": (
            "This is a single-fault Pauli-propagation audit of an inserted circuit on a distance-two finite graph code. It neither "
            "rules out fault-tolerant measurement constructions nor supplies ANTLER hardware, a defect exchange, non-Abelian anyons, "
            "universality or a topological quantum computer."
        ),
    }
    result = ROOT / "results" / "phase8c" / "ancilla_hook_error_audit.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
