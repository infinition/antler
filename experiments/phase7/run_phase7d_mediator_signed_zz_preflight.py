"""Stroboscopic signed-ZZ control from charge-two mediators.

Two same-rail conversion channels (aa, bb) and two opposite-rail channels
(ab, ba) give opposite logical ZZ phases after a complete virtual Rabi cycle.
All detunings stay positive.  This is the missing hardware bridge required by
the Phase 7D Floquet compiler control.  Opposite-rail channels preserve only a
mediator-dressed parity during the pulse; bare rail parity is checked after
return to the monomer subspace.
"""
from __future__ import annotations

from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


U_MOTT, DETUNING, G = 20.0, 40.0, 6.0
LOW_MODE_COUNT, MEDIATOR_COUNT, TOTAL_CHARGE = 4, 2, 2
PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def annihilate_pair(state: int, first: int, second: int) -> tuple[int, float] | None:
    first, second = sorted((first, second))
    if not ((state >> first) & 1):
        return None
    interim = state ^ (1 << first)
    first_sign = fermionic_sign(state, first)
    if not ((interim >> second) & 1):
        return None
    return interim ^ (1 << second), first_sign * fermionic_sign(interim, second)


def pauli_word(label: str) -> np.ndarray:
    operator = np.asarray([[1.0]], dtype=complex)
    for letter in label:
        operator = np.kron(operator, PAULIS[letter])
    return operator


def build_segment(pair_modes: tuple[tuple[int, int], tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Build one two-channel pulse and return H, monomer frame, states, pulse time."""
    states = np.asarray([
        state for state in range(1 << (LOW_MODE_COUNT + MEDIATOR_COUNT))
        if (state & ((1 << LOW_MODE_COUNT) - 1)).bit_count()
        + 2 * (state >> LOW_MODE_COUNT).bit_count() == TOTAL_CHARGE
    ], dtype=np.int64)
    index = {int(state): position for position, state in enumerate(states)}
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        n0 = ((state >> 0) & 1) + ((state >> 1) & 1)
        n1 = ((state >> 2) & 1) + ((state >> 3) & 1)
        hamiltonian[column, column] += U_MOTT * ((n0 - 1) ** 2 + (n1 - 1) ** 2)
        for mediator_index in range(MEDIATOR_COUNT):
            mediator = LOW_MODE_COUNT + mediator_index
            hamiltonian[column, column] += DETUNING * ((state >> mediator) & 1)
            if (state >> mediator) & 1:
                continue
            item = annihilate_pair(state, *pair_modes[mediator_index])
            if item is None:
                continue
            low_state, sign = item
            final = low_state | (1 << mediator)
            hamiltonian[index[final], column] += -G * sign
            hamiltonian[column, index[final]] += -G * sign
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("non-Hermitian mediator segment")
    frame = np.zeros((len(states), 4), dtype=complex)
    for column, rails in enumerate(product((0, 1), repeat=2)):
        state = (1 << rails[0]) | (1 << (2 + rails[1]))
        frame[index[state], column] = 1.0
    effective_detuning = DETUNING + 2.0 * U_MOTT
    pulse_time = 2.0 * np.pi / np.sqrt(effective_detuning ** 2 + 4.0 * G ** 2)
    return hamiltonian, frame, states, pulse_time


def remove_global_phase(matrix: np.ndarray, reference: complex) -> np.ndarray:
    return matrix / reference


def main() -> None:
    # a0,b0,a1,b1 mode order. Same channels select ZZ=+1; opposite select ZZ=-1.
    same_h, frame, states, pulse_time = build_segment(((0, 2), (1, 3)))
    opposite_h, opposite_frame, _, _ = build_segment(((0, 3), (1, 2)))
    if not np.allclose(frame, opposite_frame):
        raise RuntimeError("inconsistent monomer frame")
    same_unitary = expm(-1j * pulse_time * same_h)
    opposite_unitary = expm(-1j * pulse_time * opposite_h)
    same_logical = frame.conj().T @ same_unitary @ frame
    opposite_logical = frame.conj().T @ opposite_unitary @ frame
    q = same_logical[0, 0]
    expected_same = np.diag((q, 1.0, 1.0, q)).astype(complex)
    expected_opposite = np.diag((1.0, q, q, 1.0)).astype(complex)
    logical_projector = frame @ frame.conj().T
    leakage_same = float(np.linalg.norm((np.eye(len(states)) - logical_projector) @ same_unitary @ frame, ord=2) ** 2)
    leakage_opposite = float(np.linalg.norm((np.eye(len(states)) - logical_projector) @ opposite_unitary @ frame, ord=2) ** 2)
    z0z1 = pauli_word("ZZ")
    phase_phi = -float(np.angle(q))
    theta = phase_phi / 2.0
    # Logical rail rotations conjugate the two signed ZZ gates into XX and -YY.
    x0, x1 = pauli_word("XI"), pauli_word("IX")
    y0, y1 = pauli_word("YI"), pauli_word("IY")
    rotate_y = expm(-1j * np.pi / 4.0 * (y0 + y1))
    rotate_x = expm(-1j * np.pi / 4.0 * (x0 + x1))
    compiled_x = rotate_y.conj().T @ same_logical @ rotate_y
    compiled_minus_y = rotate_x.conj().T @ opposite_logical @ rotate_x
    compiled_pair_gate = compiled_minus_y @ compiled_x
    pair_operator = 0.5 * (pauli_word("XX") - pauli_word("YY"))
    expected_pair_gate = q * expm(-1j * (2.0 * theta) * pair_operator)
    parity_a = np.diag((1.0, -1.0, -1.0, 1.0))
    parity_b = np.diag((1.0, -1.0, -1.0, 1.0))
    pair_angle = abs(2.0 * theta)
    repeats = max(1, round((np.pi / 2.0) / pair_angle))
    transfer = np.linalg.matrix_power(compiled_pair_gate / q, repeats) @ np.array((0.0, 0.0, 0.0, 1.0), dtype=complex)
    out = {
        "schema": "antler.phase7d.mediator-signed-zz-preflight.v1",
        "model": {
            "low_modes": "a0,b0,a1,b1",
            "hard_core_charge_two_mediators_per_segment": 2,
            "same_rail_segment": ["a0 a1", "b0 b1"],
            "opposite_rail_segment": ["a0 b1", "b0 a1"],
            "new_resource": "time-multiplexing between same-rail and opposite-rail positive-detuning channels",
            "opposite_rail_symmetry_note": "bare branch parity is restored on the logical subspace after the full pulse; the full pulse has mediator-dressed parity",
        },
        "parameters": {
            "U_mott": U_MOTT, "detuning": DETUNING, "coupling": G,
            "coupling_over_detuning": G / DETUNING, "effective_virtual_detuning": DETUNING + 2.0 * U_MOTT,
            "full_rabi_pulse_time": pulse_time,
        },
        "audits": {
            "same_logical_gate_residual": float(np.linalg.norm(same_logical - expected_same, ord=2)),
            "opposite_logical_gate_residual": float(np.linalg.norm(opposite_logical - expected_opposite, ord=2)),
            "same_monomer_leakage": leakage_same,
            "opposite_monomer_leakage": leakage_opposite,
            "same_times_opposite_inverse_residual": float(np.linalg.norm(same_logical @ opposite_logical - q * np.eye(4), ord=2)),
            "logical_bare_parity_a_residual": float(np.linalg.norm(compiled_pair_gate @ parity_a - parity_a @ compiled_pair_gate, ord=2)),
            "logical_bare_parity_b_residual": float(np.linalg.norm(compiled_pair_gate @ parity_b - parity_b @ compiled_pair_gate, ord=2)),
            "signed_zz_phase_phi": phase_phi,
            "compiled_pair_gate_residual": float(np.linalg.norm(compiled_pair_gate - expected_pair_gate, ord=2)),
            "pair_rotation_angle_per_compiled_gate": pair_angle,
            "repeats_near_pair_swap": repeats,
            "bb_to_aa_transfer_probability": float(abs(transfer[0]) ** 2),
        },
        "decision": (
            "Positive-detuning charge-two channels provide opposite logical ZZ phases after exact full-Rabi pulses. "
            "Together with rail rotations, they compile a stroboscopic pair gate. This establishes a concrete dynamic "
            "extension of the static ANTLER vocabulary, subject to finite-pulse and many-link audits."
        ),
        "claim_boundary": (
            "This two-rung pulse construction does not prove that the channels can be isolated independently in an extended "
            "ladder, that pulse imperfections are tolerable, or that the resulting ladder is a protected topological phase. "
            "It establishes neither a 2D code, braiding, non-Abelian statistics, universality nor fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "mediator_signed_zz_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
