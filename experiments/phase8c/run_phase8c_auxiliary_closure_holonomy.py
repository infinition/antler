"""Phase 8C-T5f: conditional closure with one declared local auxiliary check."""
from __future__ import annotations

import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"), run_name="measurement_import")
T = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_logical_transport.py"), run_name="transport_import")
G = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"), run_name="graph_import")
edge, graph_code, parse_pauli = G["edge"], G["graph_code"], G["parse_pauli"]
independent_generators, measurement_update, pauli_label, local_gate, symplectic_parity = M["independent_generators"], M["measurement_update"], M["pauli_label"], M["local_gate"], M["symplectic_parity"]
xor, logical_basis, coefficient_solver, span_key = T["xor"], T["logical_basis"], T["coefficient_solver"], T["span_key"]


LOOP_SEQUENCES = (
    ("IIIYYIXZI", "IIIZIYZIX", "IZXIXZIXZ", "XZIIIIZXI"),
    ("IIIXIYZIX", "IIIZYIXZI", "IZYIXZIXZ", "YXIXZIIII"),
    ("IZXIXZIXZ", "XZXIIIYXZ"),
    ("IYXIXZIXZ", "YIZZIXIII", "YIIIIIIII", "XZIXZIZXI"),
)
AUXILIARY_LABEL = "YIIIIIIII"


def main() -> None:
    fixed = edge((1, 1), (2, 1))
    rotating = (edge((0, 0), (1, 0)), edge((0, 0), (0, 1)), edge((0, 0), (2, 0)), edge((0, 0), (0, 2)), edge((0, 0), (1, 0)))
    graphs = [graph_code(3, {item, fixed}) for item in rotating]
    stabilizers = independent_generators(graphs[0]["check_labels"])
    pairs = logical_basis(stabilizers)
    initial_logicals = tuple(word for pair in pairs for word in pair)
    logicals = initial_logicals
    legs = []
    for index, sequence_labels in enumerate(LOOP_SEQUENCES):
        target = independent_generators(graphs[index + 1]["check_labels"])
        stages = []
        for label in sequence_labels:
            measured = parse_pauli(label)
            anticommuting = [position for position, check in enumerate(stabilizers) if symplectic_parity(check, measured)]
            if not anticommuting:
                raise RuntimeError("registered measurement is not rank-preserving")
            pivot = stabilizers[anticommuting[0]]
            logicals = tuple(xor(word, pivot) if symplectic_parity(word, measured) else word for word in logicals)
            stabilizers, count = measurement_update(stabilizers, measured)
            if stabilizers is None or not local_gate(stabilizers)["all_single_qubit_probes_scalar_or_zero"]:
                raise RuntimeError("loop step fails the registered local-protection gate")
            stages.append({"measurement": label, "anticommuting_generators": count, "local_gate": local_gate(stabilizers)})
        if span_key(stabilizers) != span_key(target):
            raise RuntimeError("leg does not reach its target graph code")
        legs.append({"measurements": list(sequence_labels), "stages": stages})
    if span_key(stabilizers) != span_key(independent_generators(graphs[0]["check_labels"])):
        raise RuntimeError("loop does not return to the initial graph code")
    solve = coefficient_solver(tuple(stabilizers) + initial_logicals)
    matrix = [solve(word)[len(stabilizers):] for word in logicals]
    if not (all(stage["local_gate"]["all_single_qubit_probes_scalar_or_zero"] for leg in legs for stage in leg["stages"]) and any(row != [int(i == j) for j in range(6)] for i, row in enumerate(matrix))):
        raise RuntimeError("auxiliary loop is not both protected and logically nontrivial")
    output = {
        "schema": "antler.phase8c.auxiliary-closure-holonomy.v1",
        "parameters": {"new_declared_resource": "one single-vertex Pauli-Y measurement YIIIIIIII used only on the closing leg", "logical_order": ["X1", "Z1", "X2", "Z2", "X3", "Z3"]},
        "legs": legs,
        "returns_to_initial_stabilizer_span": True,
        "closed_loop_logical_symplectic_map": matrix,
        "decision": "PASS as a conditional reference closed-loop holonomy: the auxiliary Y measurement closes the four-leg code-deformation cycle, every stage passes the one-qubit local gate, and the resulting logical symplectic map is nonidentity. It is a single conditional holonomy, not a second adjacent loop or a non-Abelian braid.",
        "next_gate": "A second independent closed loop using declared resources must be built on this same code and its commutator with this map must be nonzero before Yang-Baxter or non-Abelian claims.",
        "claim_boundary": "The single-site Y measurement is a newly declared external resource, not derived from ANTLER or a physical measurement apparatus. No exchange protocol, non-Abelian braid, universality, noise robustness or fault tolerance is established.",
    }
    path = ROOT / "results" / "phase8c" / "auxiliary_closure_holonomy.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
