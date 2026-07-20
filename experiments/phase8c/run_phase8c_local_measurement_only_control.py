"""Phase 8C-T5h: bounded local-measurement false-holonomy control.

The conditional T5f/T5g loops use high-weight graph checks plus declared local
auxiliaries.  This control asks whether arbitrary single-vertex Pauli
measurements alone already close a protected stabilizer loop on the same
initial code.  It is a finite-depth rejection test, not a no-go theorem for
arbitrarily long measurement protocols.
"""
from __future__ import annotations

from collections import deque
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
M = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"),
    run_name="measurement_import",
)
T = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_logical_transport.py"),
    run_name="transport_import",
)
G = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"),
    run_name="graph_import",
)
edge, graph_code = G["edge"], G["graph_code"]
independent_generators = M["independent_generators"]
measurement_update, local_gate, pauli_label = M["measurement_update"], M["local_gate"], M["pauli_label"]
span_key = T["span_key"]

from antler.phase7_stabilizer_algebra import local_paulis  # noqa: E402


QUBITS = 9
MAXIMUM_DEPTH = 6


def main() -> None:
    initial_graph = graph_code(3, {edge((0, 0), (1, 0)), edge((1, 1), (2, 1))})
    initial = independent_generators(initial_graph["check_labels"])
    initial_key = span_key(initial)
    candidates = tuple(local_paulis(QUBITS, 1))

    # A future stabilizer update depends only on the current group and on the
    # candidate measurement.  The last label is retained solely to avoid an
    # immediate remeasurement, which would have zero deformation content.
    frontier = deque([(initial, None)])
    visited = {(initial_key, None)}
    rows = []
    for depth in range(1, MAXIMUM_DEPTH + 1):
        following = deque()
        closure_count = 0
        accepted_transitions = 0
        for stabilizers, last_label in frontier:
            for measured in candidates:
                label = pauli_label(measured)
                if label == last_label:
                    continue
                updated, anticommutes = measurement_update(stabilizers, measured)
                if updated is None or anticommutes == 0:
                    continue
                gate = local_gate(updated)
                if not gate["all_single_qubit_probes_scalar_or_zero"]:
                    continue
                accepted_transitions += 1
                if span_key(updated) == initial_key:
                    closure_count += 1
                    continue
                state_key = (span_key(updated), label)
                if state_key not in visited:
                    visited.add(state_key)
                    following.append((updated, label))
        rows.append({
            "depth": depth,
            "accepted_rank_preserving_local_transitions_from_frontier": accepted_transitions,
            "new_unique_protected_states": len(following),
            "closed_nonempty_loops_to_initial_span": closure_count,
        })
        frontier = following

    if any(row["closed_nonempty_loops_to_initial_span"] for row in rows):
        raise RuntimeError("a bounded single-vertex-only closed loop was found; this control must not be promoted")
    output = {
        "schema": "antler.phase8c.local-measurement-only-control.v1",
        "parameters": {
            "reference_geometry": "same 3x3 periodic vertex-qubit graph code used by T5f/T5g",
            "measurement_vocabulary": "all 27 nonidentity one-vertex Pauli measurements only",
            "outcome_branch": "+1 stabilizer branch",
            "maximum_nontrivial_measurement_depth": MAXIMUM_DEPTH,
            "retained_gate": "rank-preserving update and all 27 one-qubit probes scalar or zero after every step",
        },
        "initial_code": {
            "rank": len(initial),
            "encoded_qubits": QUBITS - len(initial),
            "ground_space_degeneracy": 1 << (QUBITS - len(initial)),
        },
        "rows": rows,
        "decision": (
            "PASS as a bounded false-holonomy control: no nonempty protected loop returning to the initial stabilizer span was found "
            "through six arbitrary single-vertex Pauli measurements. Thus the T5f/T5g closed loops cannot be reproduced by this "
            "short local-measurement-only grammar; their high-weight graph-check deformations remain operationally essential."
        ),
        "next_gate": (
            "This finite-depth control does not establish topological defect motion. A physical, symmetry-preserving realization of "
            "the required high-weight checks and of the declared auxiliary measurements must still be derived before any exchange claim."
        ),
        "claim_boundary": (
            "The search is exhaustive only through depth six in an imposed stabilizer measurement vocabulary. It does not rule out "
            "longer local circuits, does not validate the external checks physically, and does not demonstrate non-Abelian anyons, "
            "braiding, an ANTLER Hamiltonian, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "local_measurement_only_control.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
