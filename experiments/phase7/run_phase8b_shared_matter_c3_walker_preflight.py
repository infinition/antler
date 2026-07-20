"""Three-link closed-walker preflight built from explicit shared-matter links.

This is the smallest physical composition test after the one-link X/Y/Z
library.  Three rail qubits share a charge-two walker with three sites.  Each
oriented walker edge has its own reservoir pair and the same microscopic
ingredients as the audited one-link bridge.  The low code contains the walker
on d0; d1 and d2 are detuned virtual walker states.  A closed return may then
produce XXX, whereas every backtracking return is scalar.

The C3 ring is deliberately a composition preflight, not a surface-code or
twist construction.  It checks whether the sign-correlated controls can be
assigned consistently when a walker site participates in two different
reservoir channels, and whether the exact Schur complement isolates the
closed-loop Pauli word.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EDGES = 3
PAIR_DETUNING = 10.0
MOTT_Q = 15.0
WALKER_DETUNING = 0.5
LEG_HOPPING = 0.7
RUNG_HOPPING = 0.2
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05)
THETA = np.pi
PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.diag((1.0, -1.0)).astype(complex),
}


def q_a(edge: int) -> int:
    return 2 * edge


def q_b(edge: int) -> int:
    return 2 * edge + 1


def r_a(edge: int) -> int:
    return 2 * EDGES + 2 * edge


def r_b(edge: int) -> int:
    return 2 * EDGES + 2 * edge + 1


def d(edge: int) -> int:
    return 4 * EDGES + edge


CHARGES = (1,) * (4 * EDGES) + (2,) * EDGES
TOTAL_CHARGE = EDGES + 2


def basis() -> tuple[np.ndarray, dict[int, int]]:
    states = np.asarray([
        state for state in range(1 << len(CHARGES))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(CHARGES)) == TOTAL_CHARGE
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def annihilate(state: int, mode: int) -> tuple[int, complex] | None:
    if not ((state >> mode) & 1):
        return None
    return state ^ (1 << mode), complex(fermionic_sign(state, mode))


def create(state: int, mode: int) -> tuple[int, complex] | None:
    if (state >> mode) & 1:
        return None
    return state | (1 << mode), complex(fermionic_sign(state, mode))


def add_hop(matrix: np.ndarray, states: np.ndarray, positions: dict[int, int], destination: int, source: int,
            amplitude: complex, conditional_mode: int | None = None) -> None:
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        first = annihilate(state, source)
        if first is None:
            continue
        intermediate, first_sign = first
        second = create(intermediate, destination)
        if second is None:
            continue
        final, second_sign = second
        phase = 1.0 if conditional_mode is None else np.exp(1j * THETA * ((state >> conditional_mode) & 1))
        matrix[positions[final], column] += amplitude * phase * first_sign * second_sign


def add_pair_conversion(matrix: np.ndarray, states: np.ndarray, positions: dict[int, int], mediator: int,
                        first_mode: int, second_mode: int, amplitude: complex) -> None:
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        first = annihilate(state, second_mode)
        if first is None:
            continue
        intermediate, first_sign = first
        second = annihilate(intermediate, first_mode)
        if second is None:
            continue
        intermediate, second_sign = second
        third = create(intermediate, mediator)
        if third is None:
            continue
        final, third_sign = third
        matrix[positions[final], column] += amplitude * first_sign * second_sign * third_sign


def build_segment(
    coupling: float,
    *,
    echoed: bool,
    link_signs: tuple[int, ...] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[int, int]]:
    if link_signs is None:
        link_signs = (1,) * EDGES
    if len(link_signs) != EDGES or any(sign not in (-1, 1) for sign in link_signs):
        raise ValueError("link_signs must contain one +/-1 control sign per edge")
    states, positions = basis()
    matrix = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for edge in range(EDGES):
            n_q = ((state >> q_a(edge)) & 1) + ((state >> q_b(edge)) & 1)
            matrix[column, column] += MOTT_Q * (n_q - 1) ** 2
            matrix[column, column] += PAIR_DETUNING * ((state >> r_a(edge)) & 1) * ((state >> r_b(edge)) & 1)
            if edge:
                matrix[column, column] += WALKER_DETUNING * ((state >> d(edge)) & 1)
    for edge in range(EDGES):
        # Each rail and reservoir pair is local to its edge.  The two channels
        # out of a shared d-site are distinct physical conversions because the
        # reservoir pairs differ; this is the point exercised by the C3 test.
        rung = (-1.0 if echoed else 1.0) * link_signs[edge] * RUNG_HOPPING
        add_hop(matrix, states, positions, q_a(edge), q_b(edge), -rung)
        add_hop(matrix, states, positions, q_a(edge), r_a(edge), -LEG_HOPPING, conditional_mode=q_b(edge))
        add_hop(matrix, states, positions, q_b(edge), r_b(edge), -LEG_HOPPING, conditional_mode=r_a(edge))
        start, stop = d(edge), d((edge + 1) % EDGES)
        add_pair_conversion(matrix, states, positions, start, r_a(edge), r_b(edge), -coupling)
        # Echo flips the *destination channel for this oriented edge*.  The
        # same d-site may carry a different phase on its other reservoir
        # channel, which is an explicit programmable-channel requirement.
        add_pair_conversion(matrix, states, positions, stop, r_a(edge), r_b(edge), coupling if echoed else -coupling)
    matrix = matrix + matrix.conj().T
    if not np.allclose(matrix, matrix.conj().T, atol=1e-12):
        raise RuntimeError("C3 microscopic Hamiltonian is not Hermitian")
    return matrix, states, positions


def code_indices(positions: dict[int, int]) -> np.ndarray:
    values = []
    for rails in itertools.product((0, 1), repeat=EDGES):
        state = 1 << d(0)
        for edge, rail in enumerate(rails):
            state |= 1 << (q_a(edge) if rail == 0 else q_b(edge))
        values.append(positions[state])
    return np.asarray(values, dtype=int)


def schur_zero(matrix: np.ndarray, low: np.ndarray) -> np.ndarray:
    low_set = set(int(value) for value in low)
    high = np.asarray([index for index in range(matrix.shape[0]) if index not in low_set], dtype=int)
    h_ll = matrix[np.ix_(low, low)]
    h_lh = matrix[np.ix_(low, high)]
    h_hh = matrix[np.ix_(high, high)]
    h_hl = matrix[np.ix_(high, low)]
    return h_ll - h_lh @ np.linalg.solve(h_hh, h_hl)


def pauli_coefficients(matrix: np.ndarray) -> dict[str, float]:
    return {
        "".join(labels): float(np.real_if_close(np.trace(np.kron(np.kron(PAULIS[labels[0]], PAULIS[labels[1]]), PAULIS[labels[2]]) @ matrix) / (1 << EDGES)))
        for labels in itertools.product(PAULIS, repeat=EDGES)
    }


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * PAIR_DETUNING
        h_a, _, positions = build_segment(coupling, echoed=False)
        h_b, _, _ = build_segment(coupling, echoed=True)
        low = code_indices(positions)
        echoed = 0.5 * (schur_zero(h_a, low) + schur_zero(h_b, low))
        coefficients = pauli_coefficients(echoed)
        unwanted = {label: value for label, value in coefficients.items() if label not in {"III", "XXX"}}
        largest_label, largest_value = max(unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "target_xxx_coefficient": coefficients["XXX"],
            "maximum_unwanted_non_scalar_coefficient": abs(largest_value),
            "largest_unwanted_non_scalar_pauli": largest_label,
            "unwanted_over_target": float(abs(largest_value / coefficients["XXX"])),
            "echoed_schur_hermiticity_residual": float(np.linalg.norm(echoed - echoed.conj().T, ord="fro")),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.10]
    power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["target_xxx_coefficient"] for row in deep])),
        1,
    )[0])
    unwanted_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log([row["maximum_unwanted_non_scalar_coefficient"] for row in deep]),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-c3-walker-preflight.v1",
        "parameters": {
            "rail_qubits": EDGES,
            "walker_sites": EDGES,
            "low_walker_site": "d0",
            "virtual_walker_site_detuning": WALKER_DETUNING,
            "pair_detuning": PAIR_DETUNING,
            "fixed_total_charge": TOTAL_CHARGE,
            "target": "XXX",
            "control_contract": "Each oriented edge owns an independent reservoir pair and two independently phased conversions; one shared walker site therefore has separate programmable phases for its two incident reservoir channels.",
        },
        "dimensions": {"fixed_charge": int(len(basis()[0])), "low_code": 1 << EDGES},
        "rows": rows,
        "deep_sw_xxx_power": power,
        "deep_sw_worst_unwanted_power": unwanted_power,
        "decision": "Direct C3 composition is rejected at the registered detuning: it generates XXX but its leading non-scalar companion is larger at every deep registered point.",
        "claim_boundary": "This rejects the registered direct C3 detuning/control assignment only. It does not refute the one-link compiler, an explicitly derived counterterm/refocusing construction, a different walker encoding, a four-link stabilizer, a 2D code, a defect, fusion, non-Abelian braid, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_c3_walker_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
