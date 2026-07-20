"""Exhaustive pairwise walker-overlap audit for the seven-link mixed patch.

The finite patch in ``run_phase8b_mixed_check_patch_algebra.py`` has six
commuting four-link checks.  Algebraic commutation alone does not exclude
microscopic two-walker crosstalk.  This script therefore downfolds every one
of its 15 check pairs in the explicit phase-controlled walker Hamiltonian at
one deep Schrieffer-Wolff point, and performs complete seven-qubit Pauli
tomography of each effective block.

It is intentionally a pairwise closure gate, not a simultaneous six-walker
tiling calculation.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


ROOT = Path(__file__).resolve().parents[2]

DETUNING = 10.0
COUPLING_OVER_DETUNING = 0.05
N_QUBITS = 7
GAUGE_DIMENSION = 1 << N_QUBITS
WALKER_STATES = 4
DIMENSION = GAUGE_DIMENSION * WALKER_STATES * WALKER_STATES
PAULI_NAMES = "IXYZ"
CHECK_LABELS = (
    "YXZXIII",
    "XZIIXZI",
    "XZXYIII",
    "YXIIYXI",
    "XIXIXIX",
    "YIZIYIY",
)
PAULI_CODE = {name: code for code, name in enumerate(PAULI_NAMES)}


def index(gauge: int, left_walker: int, right_walker: int) -> int:
    return ((left_walker * WALKER_STATES + right_walker) * GAUGE_DIMENSION) + gauge


def z_eigenvalue(gauge: int, qubit: int) -> float:
    return -1.0 if (gauge >> qubit) & 1 else 1.0


def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
    matrix[row, column] += value
    matrix[column, row] += np.conj(value)


def link_transition(gauge: int, qubit: int, pauli: int, coupling: float) -> tuple[int, complex]:
    if pauli == 1:  # X
        return gauge ^ (1 << qubit), coupling
    if pauli == 2:  # i X Z = Y
        return gauge ^ (1 << qubit), 1j * coupling * z_eigenvalue(gauge, qubit)
    if pauli == 3:  # Z
        return gauge, coupling * z_eigenvalue(gauge, qubit)
    raise ValueError(f"identity cannot be a four-step walker link ({pauli=})")


def active_pattern(label: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    entries = tuple((qubit, PAULI_CODE[letter]) for qubit, letter in enumerate(label) if letter != "I")
    if len(entries) != WALKER_STATES:
        raise ValueError(f"expected a four-link check, got {label}")
    return tuple(entry[0] for entry in entries), tuple(entry[1] for entry in entries)


def add_walker_transitions(
    matrix: sparse.lil_matrix,
    gauge: int,
    left_walker: int,
    right_walker: int,
    coupling: float,
    edges: tuple[int, ...],
    pattern: tuple[int, ...],
    *,
    which: str,
) -> None:
    walker = left_walker if which == "left" else right_walker
    next_walker = (walker + 1) % WALKER_STATES
    next_gauge, amplitude = link_transition(gauge, edges[walker], pattern[walker], coupling)
    if which == "left":
        row = index(next_gauge, next_walker, right_walker)
    else:
        row = index(next_gauge, left_walker, next_walker)
    add_symmetric(matrix, row, index(gauge, left_walker, right_walker), amplitude)


def build_hamiltonian(left: str, right: str, coupling: float) -> sparse.csr_matrix:
    left_edges, left_pattern = active_pattern(left)
    right_edges, right_pattern = active_pattern(right)
    hamiltonian = sparse.lil_matrix((DIMENSION, DIMENSION), dtype=complex)
    for gauge in range(GAUGE_DIMENSION):
        for left_walker in range(WALKER_STATES):
            for right_walker in range(WALKER_STATES):
                column = index(gauge, left_walker, right_walker)
                if left_walker:
                    hamiltonian[column, column] += DETUNING
                if right_walker:
                    hamiltonian[column, column] += DETUNING
                add_walker_transitions(
                    hamiltonian, gauge, left_walker, right_walker, coupling,
                    left_edges, left_pattern, which="left",
                )
                add_walker_transitions(
                    hamiltonian, gauge, left_walker, right_walker, coupling,
                    right_edges, right_pattern, which="right",
                )
    return hamiltonian.tocsr()


def pauli_action(codes: tuple[int, ...], state: int) -> tuple[int, complex]:
    target = state
    phase = 1.0 + 0.0j
    for qubit, code in enumerate(codes):
        bit = (state >> qubit) & 1
        if code == 1:
            target ^= 1 << qubit
        elif code == 2:
            target ^= 1 << qubit
            phase *= 1j if bit == 0 else -1j
        elif code == 3:
            phase *= -1.0 if bit else 1.0
    return target, phase


def pauli_coefficients(hamiltonian: np.ndarray) -> dict[str, complex]:
    coefficients: dict[str, complex] = {}
    for codes in itertools.product(range(4), repeat=N_QUBITS):
        value = 0.0j
        for column in range(GAUGE_DIMENSION):
            row, phase = pauli_action(codes, column)
            value += np.conj(phase) * hamiltonian[row, column]
        coefficients["".join(PAULI_NAMES[code] for code in codes)] = value / GAUGE_DIMENSION
    return coefficients


def multiply_labels(left: str, right: str) -> tuple[str, complex]:
    table = {
        ("I", "I"): ("I", 1), ("I", "X"): ("X", 1), ("I", "Y"): ("Y", 1), ("I", "Z"): ("Z", 1),
        ("X", "I"): ("X", 1), ("Y", "I"): ("Y", 1), ("Z", "I"): ("Z", 1),
        ("X", "X"): ("I", 1), ("Y", "Y"): ("I", 1), ("Z", "Z"): ("I", 1),
        ("X", "Y"): ("Z", 1j), ("Y", "X"): ("Z", -1j),
        ("X", "Z"): ("Y", -1j), ("Z", "X"): ("Y", 1j),
        ("Y", "Z"): ("X", 1j), ("Z", "Y"): ("X", -1j),
    }
    letters, phase = [], 1.0 + 0.0j
    for first, second in zip(left, right):
        letter, local_phase = table[first, second]
        letters.append(letter)
        phase *= local_phase
    return "".join(letters), phase


def downfold_at_zero(hamiltonian: sparse.csr_matrix) -> tuple[np.ndarray, float]:
    low = np.asarray([index(gauge, 0, 0) for gauge in range(GAUGE_DIMENSION)], dtype=int)
    mask = np.ones(DIMENSION, dtype=bool)
    mask[low] = False
    high = np.flatnonzero(mask)
    h_lh = hamiltonian[low][:, high]
    h_hh = hamiltonian[high][:, high].tocsc()
    h_hl = hamiltonian[high][:, low].toarray()
    effective = -np.asarray(h_lh @ sparse_linalg.spsolve(h_hh, h_hl))
    return effective, float(np.linalg.norm(effective - effective.conj().T, ord="fro"))


def main() -> None:
    rows = []
    coupling = COUPLING_OVER_DETUNING * DETUNING
    for left_index, right_index in itertools.combinations(range(len(CHECK_LABELS)), 2):
        left, right = CHECK_LABELS[left_index], CHECK_LABELS[right_index]
        product_label, product_phase = multiply_labels(left, right)
        if abs(product_phase.imag) > 1e-12:
            raise RuntimeError(f"noncommuting target product for {left}, {right}: {product_phase}")
        effective, hermiticity = downfold_at_zero(build_hamiltonian(left, right, coupling))
        coefficients = pauli_coefficients(effective)
        allowed = {"I" * N_QUBITS, left, right, product_label}
        unwanted = {label: value for label, value in coefficients.items() if label not in allowed}
        largest_label, largest_value = max(unwanted.items(), key=lambda item: abs(item[1]))
        shared = [qubit for qubit, (a, b) in enumerate(zip(left, right)) if a != "I" and b != "I"]
        anticommutes = sum(a != b for a, b in zip(left, right) if a != "I" and b != "I")
        rows.append({
            "check_indices": [left_index, right_index],
            "left_check": left,
            "right_check": right,
            "shared_links": shared,
            "single_link_anticommutations": anticommutes,
            "product_label": product_label,
            "product_phase_in_left_times_right": float(np.real_if_close(product_phase)),
            "left_coefficient": float(np.real_if_close(coefficients[left])),
            "right_coefficient": float(np.real_if_close(coefficients[right])),
            "product_coefficient": float(np.real_if_close(coefficients[product_label])),
            "maximum_outside_generated_algebra_coefficient": float(abs(largest_value)),
            "largest_outside_generated_algebra_pauli": largest_label,
            "zero_energy_schur_hermiticity_residual": hermiticity,
        })
    output = {
        "schema": "antler.phase8b.mixed-check-patch-pairwise-overlap.v1",
        "parameters": {
            "gauge_links": N_QUBITS,
            "detuning": DETUNING,
            "coupling_over_detuning": COUPLING_OVER_DETUNING,
            "checks": list(CHECK_LABELS),
            "pair_count": len(rows),
            "new_primitives": "two neutral four-state phase-controlled X/Y/Z-conditioned walkers",
            "downfolding": "exact zero-energy Schur complement onto both walkers low",
        },
        "dimensions": {"total": DIMENSION, "low_gauge_subspace": GAUGE_DIMENSION, "high_walker_subspace": DIMENSION - GAUGE_DIMENSION},
        "rows": rows,
        "worst_case": {
            "maximum_outside_generated_algebra_coefficient": max(row["maximum_outside_generated_algebra_coefficient"] for row in rows),
            "maximum_hermiticity_residual": max(row["zero_energy_schur_hermiticity_residual"] for row in rows),
            "minimum_target_check_coefficient": min(min(abs(row["left_coefficient"]), abs(row["right_coefficient"])) for row in rows),
        },
        "decision": (
            "At the registered deep-SW point, every pair of checks in the finite mixed patch jointly closes on its generated "
            "commuting Pauli algebra; no pairwise walker crosstalk word is resolved."
        ),
        "claim_boundary": (
            "This closes only all pairwise overlaps at one deep point. It does not evaluate the simultaneous six-walker "
            "Hamiltonian, higher-body crosstalk, scaling, physical U(1) embedding, a complete twist/dislocation geometry, "
            "fusion, defect motion, non-Abelian braid, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_mixed_check_patch_pairwise_overlap.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
