"""Phase 8C-T5c: measurement-resolved interior-twist code deformation.

T5b rejected the direct linear Hamiltonian interpolation because it makes the
low band locally readable.  This script searches only the registered final
face checks as projective Pauli measurements.  A measurement that anticommutes
with current stabilizers replaces one generator and multiplies the other
anticommuting generators by that pivot, preserving the stabilizer rank.

The output is a code-space reference protocol with recorded +1 outcomes (other
outcomes require a Pauli frame).  It is not a physical measurement apparatus,
a microscopic ANTLER control, or a braid calculation.
"""
from __future__ import annotations

import itertools
import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GRAPH_MODULE = runpy.run_path(
    str(ROOT / "experiments" / "phase8c" / "run_phase8c_interior_twist_graph_preflight.py"),
    run_name="antler_phase8c_graph_import",
)
edge = GRAPH_MODULE["edge"]
graph_code = GRAPH_MODULE["graph_code"]
parse_pauli = GRAPH_MODULE["parse_pauli"]

from antler.phase7_stabilizer_algebra import (  # noqa: E402
    BinaryPauli,
    _gf2_basis,
    _packed,
    commute,
    gf2_in_span,
    gf2_rank,
    local_paulis,
    symplectic_parity,
)


SIZE = 3
QUBITS = SIZE * SIZE
EXPECTED_RANK = 6
EXPECTED_CODE_DIMENSION = 8


def pauli_label(word: BinaryPauli) -> str:
    symbols = []
    for qubit in range(QUBITS):
        x_bit = (word.x >> qubit) & 1
        z_bit = (word.z >> qubit) & 1
        symbols.append("Y" if x_bit and z_bit else "X" if x_bit else "Z" if z_bit else "I")
    return "".join(symbols)


def independent_generators(labels: list[str]) -> tuple[BinaryPauli, ...]:
    generators: list[BinaryPauli] = []
    basis: dict[int, int] = {}
    for label in sorted(labels):
        word = parse_pauli(label)
        packed = _packed(word, QUBITS)
        if not gf2_in_span(packed, basis):
            generators.append(word)
            basis = _gf2_basis(tuple(_packed(item, QUBITS) for item in generators))
    return tuple(generators)


def basis_key(generators: tuple[BinaryPauli, ...]) -> tuple[int, ...]:
    """Return an insertion-order-independent GF(2) row-space key.

    Raw row-echelon vectors are sufficient for membership tests but are not a
    canonical description of a stabilizer group.  The final back-elimination
    makes equality tests sound for deformation endpoints.
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


def measurement_update(
    generators: tuple[BinaryPauli, ...],
    measured: BinaryPauli,
) -> tuple[tuple[BinaryPauli, ...] | None, int]:
    """Return the +1-outcome stabilizer update and its anticommuting count."""
    anticommuting = [index for index, generator in enumerate(generators) if symplectic_parity(generator, measured)]
    if not anticommuting:
        in_span = gf2_in_span(
            _packed(measured, QUBITS),
            _gf2_basis(tuple(_packed(word, QUBITS) for word in generators)),
        )
        return (generators if in_span else None), 0
    pivot_index = anticommuting[0]
    pivot = generators[pivot_index]
    updated = []
    for index, generator in enumerate(generators):
        if index == pivot_index:
            updated.append(measured)
        elif symplectic_parity(generator, measured):
            updated.append(BinaryPauli(x=generator.x ^ pivot.x, z=generator.z ^ pivot.z))
        else:
            updated.append(generator)
    return tuple(updated), len(anticommuting)


def local_gate(generators: tuple[BinaryPauli, ...]) -> dict[str, object]:
    basis = _gf2_basis(tuple(_packed(word, QUBITS) for word in generators))
    counts = {"projects_to_zero": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
    for candidate in local_paulis(QUBITS, 1):
        if any(not commute(candidate, generator) for generator in generators):
            counts["projects_to_zero"] += 1
        elif gf2_in_span(_packed(candidate, QUBITS), basis):
            counts["stabilizer_scalar"] += 1
        else:
            counts["nontrivial_logical"] += 1
    counts["tested_nonidentity_paulis"] = sum(counts.values())
    counts["all_single_qubit_probes_scalar_or_zero"] = counts["nontrivial_logical"] == 0
    return counts


def stage_record(
    generators: tuple[BinaryPauli, ...],
    *,
    measurement: str | None,
    anticommuting_generators: int | None,
) -> dict[str, object]:
    rank = gf2_rank(tuple(_packed(word, QUBITS) for word in generators))
    return {
        "measurement_plus_outcome": measurement,
        "anticommuting_generator_count_before_measurement": anticommuting_generators,
        "generator_labels": [pauli_label(word) for word in generators],
        "pair_anticommutations": sum(
            not commute(left, right)
            for index, left in enumerate(generators)
            for right in generators[index + 1:]
        ),
        "independent_rank_over_GF2": rank,
        "encoded_qubits": QUBITS - rank,
        "ground_space_degeneracy": 1 << (QUBITS - rank),
        "local_protection": local_gate(generators),
    }


def main() -> None:
    fixed_pair = edge((1, 1), (2, 1))
    initial = graph_code(SIZE, {edge((0, 0), (1, 0)), fixed_pair})
    final = graph_code(SIZE, {edge((0, 0), (0, 1)), fixed_pair})
    current = independent_generators(initial["check_labels"])
    target = independent_generators(final["check_labels"])
    target_key = basis_key(target)
    if not (len(current) == len(target) == EXPECTED_RANK):
        raise RuntimeError("unexpected static-code rank")

    solution: tuple[BinaryPauli, ...] | None = None
    for length in range(1, len(target) + 1):
        for candidate_sequence in itertools.permutations(target, length):
            trial = current
            accepted = True
            for measured in candidate_sequence:
                trial, _ = measurement_update(trial, measured)
                if trial is None or not local_gate(trial)["all_single_qubit_probes_scalar_or_zero"]:
                    accepted = False
                    break
            if accepted and basis_key(trial) == target_key:
                solution = candidate_sequence
                break
        if solution is not None:
            break
    if solution is None:
        raise RuntimeError("no protected measurement sequence found in the registered final-check vocabulary")

    stages = [stage_record(current, measurement=None, anticommuting_generators=None)]
    for measured in solution:
        current, anticommutes = measurement_update(current, measured)
        if current is None:
            raise RuntimeError("search returned an invalid measurement")
        stages.append(stage_record(current, measurement=pauli_label(measured), anticommuting_generators=anticommutes))

    if not (
        basis_key(current) == target_key
        and len(solution) == 4
        and all(stage["pair_anticommutations"] == 0 for stage in stages)
        and all(stage["independent_rank_over_GF2"] == EXPECTED_RANK for stage in stages)
        and all(stage["ground_space_degeneracy"] == EXPECTED_CODE_DIMENSION for stage in stages)
        and all(stage["local_protection"]["all_single_qubit_probes_scalar_or_zero"] for stage in stages)
        and all(stage["anticommuting_generator_count_before_measurement"] == 2 for stage in stages[1:])
    ):
        raise RuntimeError("measurement deformation failed a registered protection gate")

    output = {
        "schema": "antler.phase8c.interior-twist-measurement-deformation.v1",
        "parameters": {
            "reference_geometry": "3x3 periodic vertex-qubit graph-code reference",
            "measurement_outcome_convention": "+1 outcomes; other outcomes require tracked Pauli-frame corrections",
            "measurement_vocabulary": "independent final face checks only; no hand-inserted logical or braid operator",
            "local_protection_gate": "all one-qubit Paulis must project to zero or a stabilizer scalar at every stabilizer stage",
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
        "measurement_sequence": [pauli_label(measured) for measured in solution],
        "stages": stages,
        "final_stabilizer_span_matches_target": basis_key(current) == target_key,
        "decision": (
            "PASS T5c as an external measurement-resolved code-deformation reference. Four recorded final-check measurements map "
            "the initial stabilizer span to the final one while retaining rank six, GSD eight, commutation and a scalar-or-zero "
            "projection for all 27 one-qubit Paulis at every stabilizer stage. This closes the direct-linear-interpolation loophole, "
            "but does not yet calculate logical transport, outcomes/Pauli frames, physical measurement errors or a braid."
        ),
        "next_gate": (
            "T5d must propagate a complete logical Pauli basis through the four measurements, resolve all outcome-dependent Pauli "
            "frames, and compare two adjacent such moves before any commutator or Yang-Baxter calculation."
        ),
        "claim_boundary": (
            "This is a stabilizer-measurement reference, not a Hamiltonian-time evolution and not a measurement apparatus derived "
            "from ANTLER. It does not establish physical e<->m transport, fusion readout, non-Abelian braid statistics, universality, "
            "experimental noise robustness or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "interior_twist_measurement_deformation.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
