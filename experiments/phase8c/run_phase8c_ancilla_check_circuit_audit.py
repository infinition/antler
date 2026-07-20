"""Phase 8C-T5i: explicit ancilla circuit audit for the measured graph checks.

Each high-weight Pauli check used by the conditional T5f/T5g reference loops
is compiled into basis changes, data-to-ancilla CNOTs, and one ancilla Z
readout.  The code verifies the exact +/− projective instrument on reproducible
statevectors.  The CNOT/readout gate set is a newly declared external resource.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA_QUBITS = 9
ANCILLA = DATA_QUBITS
TOTAL_QUBITS = DATA_QUBITS + 1
DIM_DATA = 1 << DATA_QUBITS
DIM_TOTAL = 1 << TOTAL_QUBITS
H = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / np.sqrt(2.0)
S = np.array([[1.0, 0.0], [0.0, 1.0j]], dtype=complex)
S_DAGGER = np.conjugate(S.T)


def support(label: str) -> list[int]:
    return [index for index, symbol in enumerate(label) if symbol != "I"]


def torus_distance(left: int, right: int) -> int:
    x0, y0 = left % 3, left // 3
    x1, y1 = right % 3, right // 3
    dx, dy = abs(x0 - x1), abs(y0 - y1)
    return min(dx, 3 - dx) + min(dy, 3 - dy)


def apply_one_qubit(state: np.ndarray, qubit: int, matrix: np.ndarray) -> None:
    bit = 1 << qubit
    for index in range(state.size):
        if index & bit:
            continue
        partner = index | bit
        low, high = state[index], state[partner]
        state[index] = matrix[0, 0] * low + matrix[0, 1] * high
        state[partner] = matrix[1, 0] * low + matrix[1, 1] * high


def apply_cnot_to_ancilla(state: np.ndarray, control: int) -> None:
    control_bit, ancilla_bit = 1 << control, 1 << ANCILLA
    for index in range(state.size):
        if (index & control_bit) and not (index & ancilla_bit):
            partner = index | ancilla_bit
            state[index], state[partner] = state[partner], state[index]


def apply_pauli(label: str, vector: np.ndarray) -> np.ndarray:
    x_mask = sum(1 << index for index, symbol in enumerate(label) if symbol in {"X", "Y"})
    z_mask = sum(1 << index for index, symbol in enumerate(label) if symbol in {"Z", "Y"})
    y_count = label.count("Y")
    output = np.zeros_like(vector)
    for index, amplitude in enumerate(vector):
        phase = (1.0j ** y_count) * (-1 if ((z_mask & index).bit_count() % 2) else 1)
        output[index ^ x_mask] += phase * amplitude
    return output


def measure_pauli_with_ancilla(label: str, vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    state = np.zeros(DIM_TOTAL, dtype=complex)
    state[:DIM_DATA] = vector
    for qubit, symbol in enumerate(label):
        if symbol == "X":
            apply_one_qubit(state, qubit, H)
        elif symbol == "Y":
            apply_one_qubit(state, qubit, S_DAGGER)
            apply_one_qubit(state, qubit, H)
    for qubit in support(label):
        apply_cnot_to_ancilla(state, qubit)
    for qubit, symbol in enumerate(label):
        if symbol == "X":
            apply_one_qubit(state, qubit, H)
        elif symbol == "Y":
            apply_one_qubit(state, qubit, H)
            apply_one_qubit(state, qubit, S)
    return state[:DIM_DATA].copy(), state[DIM_DATA:].copy()


def test_vectors() -> tuple[np.ndarray, ...]:
    generator = np.random.default_rng(20260720)
    vectors = []
    for _ in range(3):
        vector = generator.normal(size=DIM_DATA) + 1.0j * generator.normal(size=DIM_DATA)
        vectors.append(vector / np.linalg.norm(vector))
    return tuple(vectors)


def main() -> None:
    first = json.loads((ROOT / "results" / "phase8c" / "auxiliary_closure_holonomy.json").read_text())
    second = json.loads((ROOT / "results" / "phase8c" / "second_auxiliary_holonomy_commutator.json").read_text())
    auxiliary = {"YIIIIIIII", "IIIIIIIZI", "IIIYIIIII"}
    labels = {
        stage["measurement"]
        for leg in first["legs"]
        for stage in leg["stages"]
    } | {
        stage["measurement"]
        for leg in second["second_loop_legs"]
        for stage in leg["stages"]
    }
    checks = sorted(label for label in labels if label not in auxiliary)
    vectors = test_vectors()
    rows = []
    for label in checks:
        plus_residual = minus_residual = completeness_residual = 0.0
        for vector in vectors:
            ancilla_zero, ancilla_one = measure_pauli_with_ancilla(label, vector)
            pauli_vector = apply_pauli(label, vector)
            plus_target = 0.5 * (vector + pauli_vector)
            minus_target = 0.5 * (vector - pauli_vector)
            plus_residual = max(plus_residual, float(np.linalg.norm(ancilla_zero - plus_target)))
            minus_residual = max(minus_residual, float(np.linalg.norm(ancilla_one - minus_target)))
            completeness_residual = max(completeness_residual, abs(float(np.vdot(ancilla_zero, ancilla_zero).real + np.vdot(ancilla_one, ancilla_one).real) - 1.0))
        sites = support(label)
        diameter = max((torus_distance(left, right) for left in sites for right in sites), default=0)
        rows.append({
            "check": label,
            "weight": len(sites),
            "support_vertices": [[site % 3, site // 3] for site in sites],
            "torus_support_diameter": diameter,
            "data_to_ancilla_cnot_count": len(sites),
            "single_qubit_basis_operations": 2 * label.count("X") + 4 * label.count("Y"),
            "max_plus_instrument_residual": plus_residual,
            "max_minus_instrument_residual": minus_residual,
            "max_probability_completeness_residual": completeness_residual,
        })
    if not rows or max(max(row["max_plus_instrument_residual"], row["max_minus_instrument_residual"], row["max_probability_completeness_residual"]) for row in rows) >= 1e-12:
        raise RuntimeError("ancilla circuit does not reproduce the requested Pauli projectors")
    output = {
        "schema": "antler.phase8c.ancilla-check-circuit-audit.v1",
        "parameters": {
            "data_qubits": DATA_QUBITS,
            "readout_ancillas_per_check": 1,
            "declared_external_gate_set": "arbitrary data-ancilla CNOT, single-qubit H/S basis rotations, ancilla Z preparation/reset/readout",
            "checked_outcome_instrument": "M_+ = (I+P)/2 and M_- = (I-P)/2 on the data after uncomputing the basis rotation",
            "test_vectors": 3,
        },
        "declared_single_vertex_auxiliary_measurements": sorted(auxiliary),
        "high_weight_graph_checks": rows,
        "aggregate": {
            "unique_high_weight_checks": len(rows),
            "maximum_check_weight": max(row["weight"] for row in rows),
            "maximum_data_to_ancilla_cnot_count": max(row["data_to_ancilla_cnot_count"] for row in rows),
            "maximum_torus_support_diameter": max(row["torus_support_diameter"] for row in rows),
            "worst_instrument_residual": max(max(row["max_plus_instrument_residual"], row["max_minus_instrument_residual"], row["max_probability_completeness_residual"]) for row in rows),
        },
        "decision": (
            "PASS as an explicit external measurement-apparatus compilation: every high-weight Pauli check used by T5f/T5g is "
            "realized by one readout ancilla, basis rotations and data-to-ancilla CNOT parity accumulation, reproducing both "
            "projective branches below 1e-12 on reproducible statevectors."
        ),
        "next_gate": (
            "A physical architecture must provide geometrically local, fault-tolerant versions of these CNOT/readout operations and "
            "derive them from a microscopic Hamiltonian. In particular, this audit does not solve hook errors, connectivity, timing, "
            "noise, fault tolerance, or the ANTLER resource gap."
        ),
        "claim_boundary": (
            "The CNOT/readout gate set is inserted as a new external circuit resource. This is not an ANTLER derivation, a locality- or "
            "noise-qualified measurement apparatus, a physical defect exchange, non-Abelian anyon evidence, universality, or a topological quantum computer."
        ),
    }
    result = ROOT / "results" / "phase8c" / "ancilla_check_circuit_audit.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
