"""Phase 8C-T5b: reject or certify the naive graph-check interpolation.

The two graph codes registered by T5a differ in five local checks.  This audit
does not assume that a nonzero spectral gap is enough to call that replacement a
protected defect motion.  On the exact 3x3 reference it diagonalizes

    H(s) = -(1-s) sum_f S_f^(A) - s sum_f S_f^(B)

and separately tests whether every one-qubit Pauli is scalar on the full
eight-dimensional low-energy band.  The latter is the protection gate: a
gapped but locally readable band must not be promoted to an adiabatic braid or
even a protected code deformation.
"""
from __future__ import annotations

import json
import runpy
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parents[2]
GRAPH_MODULE = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"),
    run_name="antler_phase8c_graph_import",
)
edge = GRAPH_MODULE["edge"]
graph_code = GRAPH_MODULE["graph_code"]
parse_pauli = GRAPH_MODULE["parse_pauli"]

SIZE = 3
QUBITS = SIZE * SIZE
DIMENSION = 1 << QUBITS
GROUND_DIMENSION = 8
PATH_POINTS = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
LOCAL_SCALAR_TOLERANCE = 1e-8


def pauli_matrix(label: str) -> np.ndarray:
    pauli = parse_pauli(label)
    y_count = (pauli.x & pauli.z).bit_count()
    matrix = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    for column in range(DIMENSION):
        row = column ^ pauli.x
        phase = (1j ** y_count) * (-1 if (pauli.z & column).bit_count() % 2 else 1)
        matrix[row, column] = phase
    return matrix


def parent(check_labels: list[str]) -> np.ndarray:
    return sum((-pauli_matrix(label) for label in check_labels), np.zeros((DIMENSION, DIMENSION), dtype=complex))


def one_qubit_readability(frame: np.ndarray) -> tuple[float, str]:
    worst = -1.0
    worst_label = ""
    for qubit in range(QUBITS):
        for symbol in "XYZ":
            label = "I" * qubit + symbol + "I" * (QUBITS - qubit - 1)
            projected = frame.conj().T @ pauli_matrix(label) @ frame
            nonscalar = projected - np.eye(GROUND_DIMENSION) * np.trace(projected) / GROUND_DIMENSION
            value = float(np.linalg.norm(nonscalar, ord=2))
            if value > worst:
                worst, worst_label = value, label
    return worst, worst_label


def main() -> None:
    fixed_pair = edge((1, 1), (2, 1))
    initial_edge = edge((0, 0), (1, 0))
    final_edge = edge((0, 0), (0, 1))
    initial = graph_code(SIZE, {initial_edge, fixed_pair})
    final = graph_code(SIZE, {final_edge, fixed_pair})
    h_initial = parent(initial["check_labels"])
    h_final = parent(final["check_labels"])
    rows = []
    for parameter in PATH_POINTS:
        eigenvalues, eigenvectors = eigh((1.0 - parameter) * h_initial + parameter * h_final)
        frame = eigenvectors[:, :GROUND_DIMENSION]
        readability, witness = one_qubit_readability(frame)
        rows.append({
            "s": parameter,
            "ground_energy": float(eigenvalues[0]),
            "ground_band_width": float(eigenvalues[GROUND_DIMENSION - 1] - eigenvalues[0]),
            "gap_above_eight_state_band": float(eigenvalues[GROUND_DIMENSION] - eigenvalues[GROUND_DIMENSION - 1]),
            "maximum_single_qubit_projected_nonscalar_norm": readability,
            "worst_single_qubit_witness": witness,
        })

    midpoint = next(row for row in rows if row["s"] == 0.5)
    endpoint_rows = (rows[0], rows[-1])
    if not (
        all(row["ground_band_width"] < 1e-10 for row in rows)
        and min(row["gap_above_eight_state_band"] for row in rows) > 1.0
        and all(row["maximum_single_qubit_projected_nonscalar_norm"] < LOCAL_SCALAR_TOLERANCE for row in endpoint_rows)
        and midpoint["maximum_single_qubit_projected_nonscalar_norm"] > 0.1
    ):
        raise RuntimeError("linear-deformation control did not realize the registered gapped-but-readable negative control")

    output = {
        "schema": "antler.phase8c.interior-twist-linear-deformation-audit.v1",
        "parameters": {
            "reference_geometry": "3x3 periodic graph-code negative control; exact dense Hilbert space",
            "qubits": QUBITS,
            "full_hilbert_dimension": DIMENSION,
            "tracked_low_band_dimension": GROUND_DIMENSION,
            "hamiltonian": "H(s)=-(1-s) sum S_A - s sum S_B",
            "path_points": list(PATH_POINTS),
            "local_scalar_tolerance": LOCAL_SCALAR_TOLERANCE,
        },
        "initial_static_code": {
            "encoded_qubits": initial["encoded_qubits"],
            "ground_space_degeneracy": initial["ground_space_degeneracy"],
            "minimum_logical_weight": initial["minimum_logical_weight"],
        },
        "final_static_code": {
            "encoded_qubits": final["encoded_qubits"],
            "ground_space_degeneracy": final["ground_space_degeneracy"],
            "minimum_logical_weight": final["minimum_logical_weight"],
        },
        "rows": rows,
        "minimum_tracked_band_gap": min(row["gap_above_eight_state_band"] for row in rows),
        "maximum_path_single_qubit_readability": max(row["maximum_single_qubit_projected_nonscalar_norm"] for row in rows),
        "decision": (
            "REJECT the naive linear interpolation as a protection-preserving interior-twist deformation. The eight-state band "
            "stays exactly degenerate and remains separated on this small reference, but a one-qubit Pauli becomes non-scalar "
            "inside the band (registered midpoint witness). A spectral gap alone is therefore insufficient; no protected motion, "
            "holonomy or braid is claimed."
        ),
        "next_gate": (
            "A future T5c schedule must keep projected local operators scalar throughout the deformation, not only at its endpoints. "
            "It must be derived as a measurement-resolved or microscopic gadget protocol before any logical transport is calculated."
        ),
        "claim_boundary": (
            "This is an exact negative control in an imposed nine-vertex graph code. It neither rules out every code-deformation "
            "protocol nor derives a measurement/gadget schedule, and it does not demonstrate any physical e<->m transport, fusion, "
            "non-Abelian braid, ANTLER realization, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "interior_twist_linear_deformation_audit.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
