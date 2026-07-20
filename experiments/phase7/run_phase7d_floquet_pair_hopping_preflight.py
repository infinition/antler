"""Minimal exact Floquet pair-hopping control for the Phase 7D pivot.

This is a literature-inspired compiler control, not a reproduction of an
optical-lattice proposal.  Two rung qubits are embedded in their full fixed
two-particle Fock block.  Alternating signed, rail-rotated density interactions
compile the stroboscopic pair-hopping operator (XX - YY)/2 exactly.  The sign
modulation is an explicit *new control resource*, outside the frozen static
charge-two mediator grammar.
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


U_MOTT, V = 20.0, 1.0
TAU = np.pi / 40.0
PERIOD = 2.0 * TAU
MODE_COUNT, TOTAL_CHARGE = 4, 2
PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def annihilate(state: int, mode: int) -> tuple[int, float] | None:
    if not ((state >> mode) & 1):
        return None
    return state ^ (1 << mode), fermionic_sign(state, mode)


def create(state: int, mode: int) -> tuple[int, float] | None:
    if (state >> mode) & 1:
        return None
    return state | (1 << mode), fermionic_sign(state, mode)


def bilinear(states: np.ndarray, index: dict[int, int], destination: int, source: int) -> np.ndarray:
    operator = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        item = annihilate(int(raw_state), source)
        if item is None:
            continue
        intermediate, amplitude = item
        item = create(intermediate, destination)
        if item is None:
            continue
        final, sign = item
        operator[index[final], column] = amplitude * sign
    return operator


def pauli_word(label: str) -> np.ndarray:
    out = np.asarray([[1.0]], dtype=complex)
    for letter in label:
        out = np.kron(out, PAULIS[letter])
    return out


def main() -> None:
    states = np.asarray([state for state in range(1 << MODE_COUNT) if state.bit_count() == TOTAL_CHARGE], dtype=np.int64)
    index = {int(state): position for position, state in enumerate(states)}
    # Mode order: a0,b0,a1,b1.  The rail rotations are number conserving.
    x0 = bilinear(states, index, 0, 1) + bilinear(states, index, 1, 0)
    x1 = bilinear(states, index, 2, 3) + bilinear(states, index, 3, 2)
    y0 = -1j * bilinear(states, index, 0, 1) + 1j * bilinear(states, index, 1, 0)
    y1 = -1j * bilinear(states, index, 2, 3) + 1j * bilinear(states, index, 3, 2)
    z0 = np.diag([((int(state) >> 0) & 1) - ((int(state) >> 1) & 1) for state in states]).astype(complex)
    z1 = np.diag([((int(state) >> 2) & 1) - ((int(state) >> 3) & 1) for state in states]).astype(complex)
    mott = np.zeros((len(states), len(states)), dtype=complex)
    for position, state in enumerate(states):
        n0 = ((int(state) >> 0) & 1) + ((int(state) >> 1) & 1)
        n1 = ((int(state) >> 2) & 1) + ((int(state) >> 3) & 1)
        mott[position, position] = U_MOTT * ((n0 - 1) ** 2 + (n1 - 1) ** 2)

    # H_x and H_y are rail-rotated ZZ interactions.  The minus sign in H_y
    # requires a modulated interaction sign; it is deliberately explicit.
    h_x = mott + V * (x0 @ x1)
    h_y = mott - V * (y0 @ y1)
    h_target = mott + 0.5 * V * (x0 @ x1 - y0 @ y1)
    cycle = expm(-1j * TAU * h_y) @ expm(-1j * TAU * h_x)
    exact_target_cycle = expm(-1j * PERIOD * h_target)

    frame = np.zeros((len(states), 4), dtype=complex)
    for column, rails in enumerate(product((0, 1), repeat=2)):
        state = (1 << (2 * 0 + rails[0])) | (1 << (2 * 1 + rails[1]))
        frame[index[state], column] = 1.0
    effective = frame.conj().T @ h_target @ frame
    pair_hop = 0.5 * (pauli_word("XX") - pauli_word("YY"))
    pair_coefficient = float(np.trace(pair_hop @ effective).real / np.trace(pair_hop @ pair_hop).real)
    pauli_coefficients = {
        label: float(np.trace(pauli_word(label) @ effective).real / 4.0)
        for label in ("XX", "YY", "ZZ", "ZI", "IZ")
    }
    projector = frame @ frame.conj().T
    leakage = float(np.linalg.norm((np.eye(len(states)) - projector) @ cycle @ frame, ord=2) ** 2)
    parity_a = np.diag([(-1.0 if sum((int(state) >> mode) & 1 for mode in (0, 2)) & 1 else 1.0) for state in states])
    parity_b = np.diag([(-1.0 if sum((int(state) >> mode) & 1 for mode in (1, 3)) & 1 else 1.0) for state in states])
    cycles_to_swap = 10
    evolved = np.linalg.matrix_power(cycle, cycles_to_swap) @ frame[:, 3]  # |b0 b1>
    transfer_probability = float(abs(frame[:, 0].conj() @ evolved) ** 2)  # |a0 a1>
    out = {
        "schema": "antler.phase7d.floquet-pair-hopping-preflight.v1",
        "model": {
            "fixed_particle_number": TOTAL_CHARGE,
            "full_block_dimension": len(states),
            "logical_monomer_dimension": 4,
            "new_control_resource": "stroboscopic rail rotations plus sign-modulated ZZ interaction",
            "segments": ["H_x=U C + V X0 X1", "H_y=U C - V Y0 Y1"],
            "compiled_stroboscopic_term": "H_F=U C + V(XX-YY)/2",
        },
        "parameters": {"U_mott": U_MOTT, "V": V, "tau": TAU, "period": PERIOD, "cycles_to_pair_swap": cycles_to_swap},
        "audits": {
            "segment_commutator_frobenius": float(np.linalg.norm(h_x @ h_y - h_y @ h_x)),
            "one_period_compiler_residual": float(np.linalg.norm(cycle - exact_target_cycle, ord=2)),
            "stroboscopic_parity_a_residual": float(np.linalg.norm(cycle @ parity_a - parity_a @ cycle)),
            "stroboscopic_parity_b_residual": float(np.linalg.norm(cycle @ parity_b - parity_b @ cycle)),
            "monomer_leakage_per_period": leakage,
            "pair_hopping_coefficient": pair_coefficient,
            "logical_pauli_coefficients": pauli_coefficients,
            "bb_to_aa_transfer_probability_after_10_cycles": transfer_probability,
        },
        "decision": (
            "The exact two-rung control validates a distinct Floquet resource: stroboscopic, number-conserving pair hopping "
            "can be compiled without relying on a weak fourth-order static mediator process. It may proceed to a larger-ladder "
            "audit, not directly to a 2D code or non-Abelian braid claim."
        ),
        "claim_boundary": (
            "This is a deliberately designed two-rung compiler control inspired by number-conserving Floquet pair-hopping work. "
            "It does not reproduce a full experimental protocol, derive the sign-modulated interaction from frozen ANTLER hardware, "
            "or establish a protected phase, 2D topological order, braiding, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase7" / "floquet_pair_hopping_preflight.json"
    result.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
