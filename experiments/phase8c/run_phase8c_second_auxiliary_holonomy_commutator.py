"""Phase 8C-T5g: a second conditional loop and logical-map commutator.

This is deliberately a reference stabilizer calculation.  It searches neither
for an ANTLER Hamiltonian nor for a braid: it asks the narrower question whether
two explicitly declared, outcome-conditioned measurement loops on the same
3x3 graph code induce noncommuting logical symplectic maps.
"""
from __future__ import annotations

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
edge, graph_code, parse_pauli = G["edge"], G["graph_code"], G["parse_pauli"]
independent_generators = M["independent_generators"]
measurement_update, local_gate, symplectic_parity = M["measurement_update"], M["local_gate"], M["symplectic_parity"]
xor, logical_basis, coefficient_solver, span_key = T["xor"], T["logical_basis"], T["coefficient_solver"], T["span_key"]


# Existing T5f loop about (0, 0), imported as a fixed independently recorded map.
FIRST_MAP = (
    (1, 1, 0, 1, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (0, 1, 1, 1, 0, 0),
    (0, 0, 0, 1, 0, 0),
    (0, 0, 0, 0, 1, 0),
    (0, 0, 0, 0, 0, 1),
)

# Second loop rotates the other missing edge around (1, 1): E -> N -> W -> E.
# The two marked single-vertex checks are deliberately external resources.
SECOND_LOOP_SEQUENCES = (
    ("IIIXIZZIX", "IIIZYXXXZ", "IXZIIIIYX", "IYXIXZIII"),
    ("IIIYIZZIX", "IIIIIIIZI", "IIIIYXIXZ", "IXZIIIIZX"),
    ("XZIXZIZXI", "IIIXIYZIX", "IYXIXZIXZ", "IIIYIIIII", "IIIZYIXZI"),
)
EXTERNAL_AUXILIARIES = ("IIIIIIIZI", "IIIYIIIII")


def matmul(left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(sum(left[row][mid] * right[mid][col] for mid in range(6)) % 2 for col in range(6)) for row in range(6))


def identity() -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(row == col) for col in range(6)) for row in range(6))


def inverse(matrix: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    augmented = [list(row) + [int(row_index == col) for col in range(6)] for row_index, row in enumerate(matrix)]
    for col in range(6):
        pivot = next((row for row in range(col, 6) if augmented[row][col]), None)
        if pivot is None:
            raise RuntimeError("logical map is singular over GF(2)")
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        for row in range(6):
            if row != col and augmented[row][col]:
                augmented[row] = [left ^ right for left, right in zip(augmented[row], augmented[col])]
    return tuple(tuple(row[6:]) for row in augmented)


def symplectic(matrix: tuple[tuple[int, ...], ...]) -> bool:
    form = tuple(tuple(int((row // 2 == col // 2) and ((row + col) % 2 == 1)) for col in range(6)) for row in range(6))
    return matmul(matmul(matrix, form), tuple(tuple(matrix[col][row] for col in range(6)) for row in range(6))) == form


def run_loop() -> tuple[tuple[tuple[int, ...], ...], list[dict[str, object]]]:
    fixed = edge((0, 0), (1, 0))
    rotating = (
        edge((1, 1), (2, 1)),
        edge((1, 1), (1, 2)),
        edge((1, 1), (0, 1)),
        edge((1, 1), (2, 1)),
    )
    graphs = [graph_code(3, {fixed, item}) for item in rotating]
    stabilizers = independent_generators(graphs[0]["check_labels"])
    initial_logicals = tuple(word for pair in logical_basis(stabilizers) for word in pair)
    logicals = initial_logicals
    legs: list[dict[str, object]] = []
    for index, labels in enumerate(SECOND_LOOP_SEQUENCES):
        target = independent_generators(graphs[index + 1]["check_labels"])
        stages = []
        for label in labels:
            measured = parse_pauli(label)
            anticommuting = [position for position, check in enumerate(stabilizers) if symplectic_parity(check, measured)]
            if not anticommuting:
                raise RuntimeError(f"measurement {label} is rank-changing")
            pivot = stabilizers[anticommuting[0]]
            logicals = tuple(xor(word, pivot) if symplectic_parity(word, measured) else word for word in logicals)
            stabilizers, count = measurement_update(stabilizers, measured)
            gate = local_gate(stabilizers) if stabilizers is not None else None
            if stabilizers is None or gate is None or not gate["all_single_qubit_probes_scalar_or_zero"]:
                raise RuntimeError(f"measurement {label} fails the local-protection gate")
            stages.append({
                "measurement": label,
                "is_declared_external_auxiliary": label in EXTERNAL_AUXILIARIES,
                "anticommuting_generators": count,
                "local_gate": gate,
            })
        if span_key(stabilizers) != span_key(target):
            raise RuntimeError(f"leg {index} does not reach the intended graph-code span")
        legs.append({"measurements": list(labels), "stages": stages})
    if span_key(stabilizers) != span_key(independent_generators(graphs[0]["check_labels"])):
        raise RuntimeError("second loop does not return to its initial stabilizer span")
    solve = coefficient_solver(tuple(stabilizers) + initial_logicals)
    matrix = tuple(tuple(solve(word)[len(stabilizers):]) for word in logicals)
    if not symplectic(matrix) or matrix == identity():
        raise RuntimeError("second loop is not a nonidentity symplectic logical map")
    return matrix, legs


def main() -> None:
    second_map, legs = run_loop()
    first_then_second = matmul(FIRST_MAP, second_map)
    second_then_first = matmul(second_map, FIRST_MAP)
    commutator = matmul(matmul(matmul(FIRST_MAP, second_map), inverse(FIRST_MAP)), inverse(second_map))
    noncommuting = first_then_second != second_then_first
    if not (symplectic(FIRST_MAP) and noncommuting and commutator != identity()):
        raise RuntimeError("the two conditional loops do not furnish a noncommuting reference-map pair")
    output = {
        "schema": "antler.phase8c.second-auxiliary-holonomy-commutator.v1",
        "parameters": {
            "reference_geometry": "3x3 periodic vertex-qubit graph-code reference",
            "logical_order": ["X1", "Z1", "X2", "Z2", "X3", "Z3"],
            "external_resources": [
                "T5f first-loop single-vertex Pauli-Y measurement YIIIIIIII",
                "T5g second-loop single-vertex Pauli-Z measurement IIIIIIIZI",
                "T5g second-loop single-vertex Pauli-Y measurement IIIYIIIII",
            ],
        },
        "second_loop_legs": legs,
        "returns_to_initial_stabilizer_span": True,
        "first_closed_loop_map_T5f": [list(row) for row in FIRST_MAP],
        "second_closed_loop_map_T5g": [list(row) for row in second_map],
        "both_maps_symplectic_over_GF2": True,
        "first_then_second_map": [list(row) for row in first_then_second],
        "second_then_first_map": [list(row) for row in second_then_first],
        "maps_noncommute": noncommuting,
        "group_commutator_map": [list(row) for row in commutator],
        "decision": (
            "PASS as a conditional stabilizer-reference noncommutativity benchmark: two distinct closed, outcome-conditioned "
            "measurement loops return to the same code, pass the registered one-qubit local gate at every stage, and induce "
            "noncommuting logical symplectic maps. This is not a non-Abelian exchange or braid: the maps require declared "
            "single-site external measurements and no defect worldlines, microscopic ANTLER derivation, physical apparatus, "
            "adiabatic holonomy, fusion rule, or Yang-Baxter relation has been established."
        ),
        "next_gate": (
            "Before any braid language, replace every external measurement by an explicitly specified, symmetry-preserving physical "
            "measurement/control mechanism; derive it from a microscopic Hamiltonian; then certify defect locality, a protected "
            "two-dimensional fusion sector, gap/leakage and an outcome-resolved exchange protocol."
        ),
        "claim_boundary": (
            "This exact finite-dimensional calculation calibrates the audit stack only. It does not demonstrate non-Abelian anyons, "
            "topological braiding, an ANTLER-native Hamiltonian, fault tolerance, scalability, or a topological quantum computer."
        ),
    }
    path = ROOT / "results" / "phase8c" / "second_auxiliary_holonomy_commutator.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
