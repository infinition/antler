"""Phase 8C-T5e: audit closure of the first interior-twist reference cycle."""
from __future__ import annotations

import itertools
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"), run_name="measurement_import")
G = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"), run_name="graph_import")
T = runpy.run_path(str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_logical_transport.py"), run_name="transport_import")
edge, graph_code, parse_pauli = G["edge"], G["graph_code"], G["parse_pauli"]
independent_generators, measurement_update, pauli_label, local_gate = M["independent_generators"], M["measurement_update"], M["pauli_label"], M["local_gate"]
span_key = T["span_key"]


def find_sequence(source: dict, target: dict, pool: tuple, target_key: frozenset[int]):
    start = independent_generators(source["check_labels"])
    tested = 0
    for length in range(1, len(pool) + 1):
        for sequence in itertools.permutations(pool, length):
            tested += 1
            current = start
            valid = True
            for measured in sequence:
                current, _ = measurement_update(current, measured)
                if current is None or not local_gate(current)["all_single_qubit_probes_scalar_or_zero"]:
                    valid = False
                    break
            if valid and span_key(current) == target_key:
                return sequence, tested
    return None, tested


def main() -> None:
    fixed = edge((1, 1), (2, 1))
    rotating_edges = (edge((0, 0), (1, 0)), edge((0, 0), (0, 1)), edge((0, 0), (2, 0)), edge((0, 0), (0, 2)))
    graphs = [graph_code(3, {rotating, fixed}) for rotating in rotating_edges]
    legs = []
    for index in range(4):
        source, target = graphs[index], graphs[(index + 1) % 4]
        all_target = tuple(parse_pauli(label) for label in sorted(target["check_labels"]))
        sequence, tested = find_sequence(source, target, all_target, span_key(independent_generators(target["check_labels"])))
        legs.append({
            "from_edge": sorted([list(vertex) for vertex in rotating_edges[index]]),
            "to_edge": sorted([list(vertex) for vertex in rotating_edges[(index + 1) % 4]]),
            "target_check_candidates": len(all_target),
            "permutations_tested": tested,
            "protected_measurement_sequence": None if sequence is None else [pauli_label(word) for word in sequence],
            "passes_registered_single_qubit_gate": sequence is not None,
        })
    if not (all(leg["passes_registered_single_qubit_gate"] for leg in legs[:3]) and not legs[3]["passes_registered_single_qubit_gate"]):
        raise RuntimeError("the registered reference cycle did not reproduce its three-leg pass / closing-leg fail pattern")
    output = {
        "schema": "antler.phase8c.interior-twist-cycle-closure-audit.v1",
        "parameters": {
            "reference_geometry": "3x3 periodic vertex-qubit graph-code; one degree-three endpoint is moved through four incident-edge orientations",
            "search": "all nonrepeating orders of all seven target face checks, accepting only rank-preserving stabilizer updates with the registered one-qubit local gate",
        },
        "legs": legs,
        "decision": "REJECT closure of this first four-leg reference cycle in the registered target-check measurement grammar. Three legs admit protected sequences, but the final S-to-E leg has no accepted ordering even with all seven target checks. Thus no closed holonomy, exchange, commutator or braid is reported.",
        "next_gate": "A new declared auxiliary measurement/check or a different interior-twist cellulation must be derived and pass the same local gate before a closed loop can be attempted.",
        "claim_boundary": "This is a finite grammar closure result only. It does not rule out all twist-deformation protocols and does not demonstrate a braid, non-Abelian statistics, ANTLER realization, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase8c" / "interior_twist_cycle_closure_audit.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
