"""Phase 8C-T5d: logical transport and outcome-frame audit for T5c.

The measured checks of T5c are replayed algebraically.  At each +1-outcome
measurement, logical Paulis that anticommute with the measured check are
multiplied by the pre-measurement pivot stabilizer.  This is the stabilizer
transport rule; no braid matrix is inserted.  The same pivots give the Pauli
frame that would correct a -1 outcome before the next measurement.
"""
from __future__ import annotations

import itertools
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MEASUREMENT_MODULE = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_measurement_deformation.py"),
    run_name="antler_phase8c_measurement_import",
)
GRAPH_MODULE = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"),
    run_name="antler_phase8c_graph_import",
)
edge = GRAPH_MODULE["edge"]
graph_code = GRAPH_MODULE["graph_code"]
parse_pauli = GRAPH_MODULE["parse_pauli"]
independent_generators = MEASUREMENT_MODULE["independent_generators"]
measurement_update = MEASUREMENT_MODULE["measurement_update"]
pauli_label = MEASUREMENT_MODULE["pauli_label"]
local_gate = MEASUREMENT_MODULE["local_gate"]

from antler.phase7_stabilizer_algebra import (  # noqa: E402
    BinaryPauli,
    _gf2_basis,
    _packed,
    commute,
    gf2_in_span,
    local_paulis,
    symplectic_parity,
)


SIZE = 3
QUBITS = 9
MEASUREMENT_LABELS = ("IIIYYIXZI", "IIIZIYZIX", "IZXIXZIXZ", "XZIIIIZXI")


def xor(left: BinaryPauli, right: BinaryPauli) -> BinaryPauli:
    return BinaryPauli(x=left.x ^ right.x, z=left.z ^ right.z)


def span_key(generators: tuple[BinaryPauli, ...]) -> tuple[int, ...]:
    """Return the canonical GF(2) row-space representative.

    ``_gf2_basis`` is a row-echelon *membership* representation.  Its rows
    depend on insertion order, so comparing its raw values can reject two
    identical stabilizer spans.  Eliminate every lower pivot to obtain a
    reduced, ordered basis before using it as an equality key.
    """
    basis = _gf2_basis(tuple(_packed(word, QUBITS) for word in generators))
    for pivot in sorted(basis):
        row = basis[pivot]
        for lower_pivot in range(pivot):
            lower_row = basis.get(lower_pivot)
            if lower_row is not None and ((row >> lower_pivot) & 1):
                row ^= lower_row
        basis[pivot] = row
    return tuple(basis[pivot] for pivot in sorted(basis))


def logical_basis(stabilizers: tuple[BinaryPauli, ...]) -> tuple[tuple[BinaryPauli, BinaryPauli], ...]:
    """Return a deterministic symplectic quotient basis, ordered by Pauli weight."""
    span = _gf2_basis(tuple(_packed(word, QUBITS) for word in stabilizers))
    candidates: list[BinaryPauli] = []
    for weight in range(1, QUBITS + 1):
        for candidate in local_paulis(QUBITS, weight):
            if all(commute(candidate, stabilizer) for stabilizer in stabilizers) and not gf2_in_span(_packed(candidate, QUBITS), span):
                candidates.append(candidate)
                span = _gf2_basis(tuple(_packed(word, QUBITS) for word in (*stabilizers, *candidates)))
                if len(candidates) == 6:
                    break
        if len(candidates) == 6:
            break
    if len(candidates) != 6:
        raise RuntimeError("could not construct the three logical Pauli pairs")
    remaining = candidates[:]
    pairs: list[tuple[BinaryPauli, BinaryPauli]] = []
    while remaining:
        x_word = remaining.pop(0)
        z_index = next(index for index, word in enumerate(remaining) if symplectic_parity(x_word, word))
        z_word = remaining.pop(z_index)
        updated = []
        for word in remaining:
            updated.append(
                BinaryPauli(
                    x=word.x ^ (x_word.x if symplectic_parity(word, z_word) else 0) ^ (z_word.x if symplectic_parity(word, x_word) else 0),
                    z=word.z ^ (x_word.z if symplectic_parity(word, z_word) else 0) ^ (z_word.z if symplectic_parity(word, x_word) else 0),
                )
            )
        remaining = updated
        pairs.append((x_word, z_word))
    return tuple(pairs)


def logical_contract(logicals: tuple[BinaryPauli, ...], stabilizers: tuple[BinaryPauli, ...]) -> dict[str, object]:
    stabilizer_span = _gf2_basis(tuple(_packed(word, QUBITS) for word in stabilizers))
    matrix = [[symplectic_parity(left, right) for right in logicals] for left in logicals]
    expected = [[0] * 6 for _ in range(6)]
    for pair in range(3):
        expected[2 * pair][2 * pair + 1] = 1
        expected[2 * pair + 1][2 * pair] = 1
    return {
        "all_commute_with_current_stabilizers": all(all(commute(logical, stabilizer) for stabilizer in stabilizers) for logical in logicals),
        "all_outside_current_stabilizer_span": all(not gf2_in_span(_packed(logical, QUBITS), stabilizer_span) for logical in logicals),
        "symplectic_matrix": matrix,
        "canonical_symplectic_form": expected,
        "is_canonical_symplectic_basis": matrix == expected,
    }


def coefficient_solver(generators: tuple[BinaryPauli, ...]):
    rows: dict[int, tuple[int, int]] = {}
    for index, generator in enumerate(generators):
        vector = _packed(generator, QUBITS)
        coefficients = 1 << index
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = (vector, coefficients)
                break
            prior_vector, prior_coefficients = rows[pivot]
            vector ^= prior_vector
            coefficients ^= prior_coefficients
    if len(rows) != len(generators):
        raise RuntimeError("coordinate generators are not independent")

    def solve(word: BinaryPauli) -> list[int]:
        vector = _packed(word, QUBITS)
        coefficients = 0
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in rows:
                raise RuntimeError("word is outside the declared centralizer basis")
            row_vector, row_coefficients = rows[pivot]
            vector ^= row_vector
            coefficients ^= row_coefficients
        return [(coefficients >> index) & 1 for index in range(len(generators))]

    return solve


def outcome_frame_table(pivots: tuple[BinaryPauli, ...]) -> list[dict[str, object]]:
    rows = []
    for bits in itertools.product((0, 1), repeat=len(pivots)):
        frame = BinaryPauli(x=0, z=0)
        for bit, pivot in zip(bits, pivots):
            if bit:
                frame = xor(frame, pivot)
        rows.append({
            "minus_outcome_bits": list(bits),
            "chronological_pauli_frame_mod_global_phase": pauli_label(frame),
        })
    return rows


def main() -> None:
    fixed_pair = edge((1, 1), (2, 1))
    initial_graph = graph_code(SIZE, {edge((0, 0), (1, 0)), fixed_pair})
    final_graph = graph_code(SIZE, {edge((0, 0), (0, 1)), fixed_pair})
    stabilizers = independent_generators(initial_graph["check_labels"])
    target_stabilizers = independent_generators(final_graph["check_labels"])
    logical_pairs = logical_basis(stabilizers)
    logicals = tuple(word for pair in logical_pairs for word in pair)
    initial_logicals = logicals
    measured = tuple(parse_pauli(label) for label in MEASUREMENT_LABELS)
    stages = [{
        "measurement_plus_outcome": None,
        "logical_labels": [pauli_label(word) for word in logicals],
        "logical_contract": logical_contract(logicals, stabilizers),
    }]
    pivots = []
    for word in measured:
        anticommuting = [index for index, stabilizer in enumerate(stabilizers) if symplectic_parity(stabilizer, word)]
        if len(anticommuting) != 2:
            raise RuntimeError("registered measurement lost its two-generator replacement structure")
        pivot = stabilizers[anticommuting[0]]
        pivots.append(pivot)
        logicals = tuple(xor(logical, pivot) if symplectic_parity(logical, word) else logical for logical in logicals)
        stabilizers, measured_anticommutations = measurement_update(stabilizers, word)
        if stabilizers is None or measured_anticommutations != 2:
            raise RuntimeError("invalid stabilizer update")
        stages.append({
            "measurement_plus_outcome": pauli_label(word),
            "minus_outcome_frame_generator_before_measurement": pauli_label(pivot),
            "logical_labels": [pauli_label(logical) for logical in logicals],
            "logical_contract": logical_contract(logicals, stabilizers),
            "single_qubit_local_gate": local_gate(stabilizers),
        })

    if span_key(stabilizers) != span_key(target_stabilizers):
        raise RuntimeError("measurement sequence did not reach the T5c final stabilizer span")
    final_pairs = logical_basis(target_stabilizers)
    final_logicals = tuple(word for pair in final_pairs for word in pair)
    solve = coefficient_solver(tuple(stabilizers) + final_logicals)
    coordinates = [solve(logical) for logical in logicals]
    logical_coordinates = [row[len(stabilizers):] for row in coordinates]
    if not (
        all(stage["logical_contract"]["is_canonical_symplectic_basis"] for stage in stages)
        and all(stage["logical_contract"]["all_commute_with_current_stabilizers"] for stage in stages)
        and all(stage["logical_contract"]["all_outside_current_stabilizer_span"] for stage in stages)
        and all(stage["single_qubit_local_gate"]["all_single_qubit_probes_scalar_or_zero"] for stage in stages[1:])
    ):
        raise RuntimeError("logical transport failed the commutation, quotient or local-protection gate")

    output = {
        "schema": "antler.phase8c.interior-twist-logical-transport.v1",
        "parameters": {
            "reference_geometry": "3x3 periodic vertex-qubit graph-code reference",
            "outcome_convention": "+1 branch transports logical representatives; each -1 outcome is corrected before the next measurement by its recorded pivot Pauli",
            "logical_order": ["X1", "Z1", "X2", "Z2", "X3", "Z3"],
            "final_coordinate_order": ["X1_final", "Z1_final", "X2_final", "Z2_final", "X3_final", "Z3_final"],
        },
        "measurement_sequence": list(MEASUREMENT_LABELS),
        "initial_logical_pairs": [[pauli_label(x_word), pauli_label(z_word)] for x_word, z_word in logical_pairs],
        "final_canonical_logical_pairs": [[pauli_label(x_word), pauli_label(z_word)] for x_word, z_word in final_pairs],
        "stages": stages,
        "outcome_frame_table": outcome_frame_table(tuple(pivots)),
        "transported_initial_logicals": [pauli_label(word) for word in logicals],
        "transport_coordinates_in_final_logical_basis": logical_coordinates,
        "decision": (
            "PASS T5d as a +1-outcome logical-transport and Pauli-frame reference. The six transported representatives remain a "
            "canonical symplectic logical basis at every stabilizer stage, terminate in the final-code centralizer, and have explicit "
            "coordinates in a declared final basis. All sixteen classical outcome patterns have a recorded chronological Pauli-frame "
            "word modulo global phase. This is one code deformation, not two exchanges or a non-Abelian braid."
        ),
        "next_gate": (
            "T5e must construct a second adjacent interior-twist deformation on the same declared logical basis, concatenate the two "
            "outcome-resolved transports, and only then evaluate a commutator. Yang-Baxter is gated on a nonzero commutator."
        ),
        "claim_boundary": (
            "The result is an imposed stabilizer-measurement reference. It does not derive a physical measurement apparatus, ANTLER "
            "microscopic Hamiltonian, e<->m transport, fusion readout, an exchange path, non-Abelian braid statistics, universality, "
            "noise robustness or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "interior_twist_logical_transport.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
