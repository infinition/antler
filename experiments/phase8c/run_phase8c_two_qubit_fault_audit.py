"""Phase 8C-T5k: exhaustive single two-qubit CNOT-Pauli fault screen."""
from __future__ import annotations

import itertools
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"), run_name="measurement_import")
G = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"), run_name="graph_import")
edge, graph_code, parse_pauli = G["edge"], G["graph_code"], G["parse_pauli"]
independent_generators, measurement_update = M["independent_generators"], M["measurement_update"]

from antler.phase7_stabilizer_algebra import BinaryPauli, _gf2_basis, _packed, commute, gf2_in_span  # noqa: E402


QUBITS, ANCILLA = 9, 9
AUXILIARY = {"YIIIIIIII", "IIIIIIIZI", "IIIYIIIII"}


def protocols():
    first = json.loads((ROOT / "results" / "phase8c" / "auxiliary_closure_holonomy.json").read_text())
    second = json.loads((ROOT / "results" / "phase8c" / "second_auxiliary_holonomy_commutator.json").read_text())
    return (("T5f", tuple(stage["measurement"] for leg in first["legs"] for stage in leg["stages"])), ("T5g", tuple(stage["measurement"] for leg in second["second_loop_legs"] for stage in leg["stages"])))


def classify(x_mask: int, z_mask: int, stabilizers) -> str:
    word = BinaryPauli(x=x_mask, z=z_mask)
    if x_mask == 0 and z_mask == 0:
        return "identity_or_scalar"
    if any(not commute(word, stabilizer) for stabilizer in stabilizers):
        return "detected_by_existing_stabilizer"
    basis = _gf2_basis(tuple(_packed(stabilizer, QUBITS) for stabilizer in stabilizers))
    return "stabilizer_scalar" if gf2_in_span(_packed(word, QUBITS), basis) else "nontrivial_logical"


def cnot_conjugate(x_mask: int, z_mask: int, control: int) -> tuple[int, int]:
    # Target is the readout ancilla.  Conjugation: X_c -> X_c X_a, Z_a -> Z_c Z_a.
    z_ancilla = (z_mask >> ANCILLA) & 1
    if (x_mask >> control) & 1:
        x_mask ^= 1 << ANCILLA
    if z_ancilla:
        z_mask ^= 1 << control
    return x_mask, z_mask


def uncompute_basis(label: str, x_mask: int, z_mask: int) -> tuple[int, int]:
    for qubit, symbol in enumerate(label):
        x_bit, z_bit = (x_mask >> qubit) & 1, (z_mask >> qubit) & 1
        if symbol == "X":
            x_bit, z_bit = z_bit, x_bit
        elif symbol == "Y":
            # Conjugation by H then S, i.e. by the inverse of S† followed by H.
            x_bit, z_bit = z_bit, x_bit
            z_bit ^= x_bit
        x_mask = (x_mask & ~(1 << qubit)) | (x_bit << qubit)
        z_mask = (z_mask & ~(1 << qubit)) | (z_bit << qubit)
    return x_mask, z_mask


def screen_order(label: str, order: tuple[int, ...], stabilizers) -> dict[str, int]:
    totals = {"faults": 0, "data_logical": 0, "data_detected": 0, "data_scalar": 0, "readout_bit_flips": 0, "logical_with_readout_flip": 0}
    for after in range(1, len(order) + 1):
        data_qubit = order[after - 1]
        for x_data, z_data, x_ancilla, z_ancilla in itertools.product((0, 1), repeat=4):
            if not any((x_data, z_data, x_ancilla, z_ancilla)):
                continue
            x_mask = (x_data << data_qubit) | (x_ancilla << ANCILLA)
            z_mask = (z_data << data_qubit) | (z_ancilla << ANCILLA)
            for remaining_control in order[after:]:
                x_mask, z_mask = cnot_conjugate(x_mask, z_mask, remaining_control)
            x_mask, z_mask = uncompute_basis(label, x_mask, z_mask)
            kind = classify(x_mask & ((1 << QUBITS) - 1), z_mask & ((1 << QUBITS) - 1), stabilizers)
            totals["faults"] += 1
            totals["data_logical"] += kind == "nontrivial_logical"
            totals["data_detected"] += kind == "detected_by_existing_stabilizer"
            totals["data_scalar"] += kind in {"identity_or_scalar", "stabilizer_scalar"}
            totals["readout_bit_flips"] += (x_mask >> ANCILLA) & 1
            totals["logical_with_readout_flip"] += kind == "nontrivial_logical" and ((x_mask >> ANCILLA) & 1)
    return totals


def main() -> None:
    initial = independent_generators(graph_code(3, {edge((0, 0), (1, 0)), edge((1, 1), (2, 1))})["check_labels"])
    rows = []
    for protocol, sequence in protocols():
        stabilizers = initial
        for stage_index, label in enumerate(sequence):
            if label not in AUXILIARY:
                sites = tuple(index for index, symbol in enumerate(label) if symbol != "I")
                candidates = [(screen_order(label, order, stabilizers), order) for order in itertools.permutations(sites)]
                best_counts, best_order = min(candidates, key=lambda item: (item[0]["data_logical"], item[0]["logical_with_readout_flip"], item[0]["readout_bit_flips"], item[1]))
                rows.append({"protocol": protocol, "stage_index": stage_index, "check": label, "weight": len(sites), "best_cnot_data_order": list(best_order), **best_counts})
            stabilizers, _ = measurement_update(stabilizers, parse_pauli(label))
            if stabilizers is None:
                raise RuntimeError("recorded stabilizer update failed")
    aggregate = {"high_weight_measurement_stages": len(rows), "stages_with_zero_data_logical_single_cnot_faults": sum(row["data_logical"] == 0 for row in rows), "stages_with_unavoidable_data_logical_faults": sum(row["data_logical"] > 0 for row in rows), "worst_minimum_data_logical_fault_count": max(row["data_logical"] for row in rows), "total_best_schedule_readout_bit_flips": sum(row["readout_bit_flips"] for row in rows)}
    output = {"schema": "antler.phase8c.two-qubit-fault-audit.v1", "parameters": {"fault_model": "all 15 nonidentity Pauli faults after each data-ancilla CNOT", "schedule_search": "all 4! or 6! data-CNOT orders", "readout_boundary": "an ancilla X/Y component flips the measured bit and requires repeated syndrome/decoding"}, "rows": rows, "aggregate": aggregate, "decision": "PENDING classification after exhaustive screen.", "next_gate": "If data-logical faults are absent, audit repeated noisy syndrome extraction and decoding; otherwise use flags/cat ancillas or reject this circuit architecture.", "claim_boundary": "This is an external finite-circuit Pauli-fault screen, not a physical ANTLER derivation, threshold proof, braid or non-Abelian claim."}
    if aggregate["stages_with_unavoidable_data_logical_faults"]:
        output["decision"] = "REJECT the bare one-ancilla circuit for the all-two-qubit-Pauli fault screen: at least one stage has unavoidable data-logical single-CNOT faults after schedule optimization."
    else:
        output["decision"] = "PASS the data-logical part of the all-two-qubit-Pauli fault screen: every recorded stage has an order with zero data-logical single-CNOT faults. Readout-bit flips remain and require repeated syndrome extraction and decoding."
    result = ROOT / "results" / "phase8c" / "two_qubit_fault_audit.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
