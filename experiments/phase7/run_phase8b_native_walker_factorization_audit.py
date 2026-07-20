"""Minimal native-resource obstruction for a conditional neutral-walker link.

The frozen/registered ANTLER ingredients provide a rung hopping X on a
one-particle rail qubit and charge-two conversions between a reservoir pair
and detuned mediator modes.  This smallest fixed-charge block asks whether
coherent phases on those pair channels can already generate a *conditional*
walker transition W tensor X/Y/Z.

They cannot when the mediator reservoir and rail qubit are disjoint: the full
Hamiltonian factorizes exactly as H_q tensor I + I tensor H_w at every time.
Thus arbitrary piecewise/Floquet phase programming still factorizes its
propagator, and no conditional walker link can emerge.  This is a scoped
obstruction, not a no-go for a future non-separable microscopic coupling.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]

# q_a,q_b carry the rail qubit; r_0,r_1 form its separate charge-two
# reservoir; d_0,d_1 are hard-core charge-two mediator species.
CHARGES = (1, 1, 1, 1, 2, 2)
TOTAL_CHARGE = 3
Q_A, Q_B, R_0, R_1, D_0, D_1 = range(6)
DELTA = 10.0


def weighted_basis() -> tuple[np.ndarray, dict[int, int]]:
    states = np.asarray([
        state for state in range(1 << len(CHARGES))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(CHARGES)) == TOTAL_CHARGE
        and ((state >> Q_A) & 1) + ((state >> Q_B) & 1) == 1
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


def bilinear(states: np.ndarray, positions: dict[int, int], destination: int, source: int) -> np.ndarray:
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
        operator[positions[final], column] = amplitude * sign
    return operator


def convert_pair(states: np.ndarray, positions: dict[int, int], mediator: int) -> np.ndarray:
    """d^dag r0 r1 + h.c. on the fixed weighted-charge basis."""
    operator = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        first = annihilate(state, R_1)
        if first is None:
            continue
        intermediate, amplitude = first
        second = annihilate(intermediate, R_0)
        if second is None:
            continue
        intermediate, sign = second
        final = create(intermediate, mediator)
        if final is None:
            continue
        final_state, final_sign = final
        operator[positions[final_state], column] = amplitude * sign * final_sign
    return operator + operator.conj().T


def product_order(states: np.ndarray) -> tuple[list[int], dict[int, int]]:
    """Product basis q in {a,b} x walker in {pair,d0,d1}."""
    raw = {
        (0, 0): (1 << Q_A) | (1 << R_0) | (1 << R_1),
        (1, 0): (1 << Q_B) | (1 << R_0) | (1 << R_1),
        (0, 1): (1 << Q_A) | (1 << D_0),
        (1, 1): (1 << Q_B) | (1 << D_0),
        (0, 2): (1 << Q_A) | (1 << D_1),
        (1, 2): (1 << Q_B) | (1 << D_1),
    }
    ordered = [raw[q, walker] for q in range(2) for walker in range(3)]
    positions = {int(state): index for index, state in enumerate(states)}
    if set(ordered) != set(int(state) for state in states):
        raise RuntimeError("fixed-charge block does not factor into the declared q x walker basis")
    return ordered, positions


def build_hamiltonian(j_perp: float, g0: complex, g1: complex) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    states, positions = weighted_basis()
    order, old_positions = product_order(states)
    reorder = np.asarray([old_positions[state] for state in order], dtype=int)
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    hamiltonian += -j_perp * (bilinear(states, positions, Q_A, Q_B) + bilinear(states, positions, Q_B, Q_A))
    for mediator in (D_0, D_1):
        hamiltonian += np.diag([DELTA if (int(state) >> mediator) & 1 else 0.0 for state in states])
    conversion0 = convert_pair(states, positions, D_0)
    conversion1 = convert_pair(states, positions, D_1)
    # The real/imaginary phase is attached to the directed conversion before
    # adding its Hermitian conjugate; rebuild it from the upper directed block.
    directed0 = 0.5 * (conversion0 + conversion0.conj().T)
    directed1 = 0.5 * (conversion1 + conversion1.conj().T)
    # Conversion matrices are real symmetric in this convention. A coherent
    # channel phase is represented by its cosine component here; the separate
    # Floquet test below also uses a pi/2 phase through the explicit complex
    # directed operator returned by ``directed_conversion``.
    hamiltonian += g0.real * directed0 + g1.real * directed1
    return hamiltonian[np.ix_(reorder, reorder)], states, reorder


def directed_conversion(states: np.ndarray, positions: dict[int, int], mediator: int) -> np.ndarray:
    operator = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        first = annihilate(state, R_1)
        if first is None:
            continue
        intermediate, amplitude = first
        second = annihilate(intermediate, R_0)
        if second is None:
            continue
        intermediate, sign = second
        final = create(intermediate, mediator)
        if final is None:
            continue
        final_state, final_sign = final
        operator[positions[final_state], column] = amplitude * sign * final_sign
    return operator


def build_phased_hamiltonian(j_perp: float, channels: tuple[complex, complex]) -> np.ndarray:
    states, positions = weighted_basis()
    order, old_positions = product_order(states)
    reorder = np.asarray([old_positions[state] for state in order], dtype=int)
    hamiltonian = -j_perp * (bilinear(states, positions, Q_A, Q_B) + bilinear(states, positions, Q_B, Q_A))
    for mediator in (D_0, D_1):
        hamiltonian += np.diag([DELTA if (int(state) >> mediator) & 1 else 0.0 for state in states])
    for mediator, coupling in zip((D_0, D_1), channels):
        directed = directed_conversion(states, positions, mediator)
        hamiltonian += coupling * directed + np.conj(coupling) * directed.conj().T
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("non-Hermitian phased native block")
    return hamiltonian[np.ix_(reorder, reorder)]


def partial_components(hamiltonian: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, complex]:
    """Best trace decomposition H = H_q x I + I x H_w plus scalar gauge."""
    tensor = hamiltonian.reshape(2, 3, 2, 3)
    h_q = np.einsum("aibi->ab", tensor) / 3.0
    h_w = np.einsum("aiaj->ij", tensor) / 2.0
    scalar = np.trace(hamiltonian) / 6.0
    residual = hamiltonian - np.kron(h_q, np.eye(3)) - np.kron(np.eye(2), h_w) + scalar * np.eye(6)
    return h_q, h_w, float(np.linalg.norm(residual, ord="fro")), scalar


def json_complex_matrix(matrix: np.ndarray) -> list[list[list[float]]]:
    return [[[float(value.real), float(value.imag)] for value in row] for row in matrix]


def main() -> None:
    states, _ = weighted_basis()
    if len(states) != 6:
        raise RuntimeError(f"unexpected fixed-charge dimension {len(states)}")
    h0 = build_phased_hamiltonian(0.7, (0.5, 0.5j))
    hq, hw, static_residual, static_scalar = partial_components(h0)
    segments = (
        (0.7, (0.5, 0.5j), 0.11),
        (0.2, (0.35j, -0.4), 0.07),
        (-0.45, (0.2 * np.exp(0.37j), 0.45 * np.exp(-0.91j)), 0.13),
    )
    total = np.eye(6, dtype=complex)
    q_total = np.eye(2, dtype=complex)
    walker_total = np.eye(3, dtype=complex)
    segment_residuals = []
    global_phase_exponent = 0.0
    for j_perp, channels, duration in segments:
        hamiltonian = build_phased_hamiltonian(j_perp, channels)
        q_part, walker_part, residual, scalar = partial_components(hamiltonian)
        segment_residuals.append(residual)
        global_phase_exponent += duration * float(np.real_if_close(scalar))
        total = expm(-1j * duration * hamiltonian) @ total
        q_total = expm(-1j * duration * q_part) @ q_total
        walker_total = expm(-1j * duration * walker_part) @ walker_total
    factorized_unitary_residual = float(np.linalg.norm(
        total - np.exp(1j * global_phase_exponent) * np.kron(q_total, walker_total), ord="fro"
    ))
    # A conditional walker link would be a nonzero connected term, e.g.
    # |d0><d1|+h.c. tensored X_q.  All such terms are absent under exact
    # factorization, independently of the coherent channel phases.
    output = {
        "schema": "antler.phase8b.native-walker-factorization-audit.v1",
        "parameters": {
            "weighted_total_charge": TOTAL_CHARGE,
            "mode_order": ["q_a", "q_b", "r0", "r1", "d0", "d1"],
            "fixed_charge_dimension": len(states),
            "native_terms_tested": [
                "rung hopping q_a^dag q_b+h.c.",
                "two independently phase-programmed charge-two conversions d_i^dag r0 r1+h.c.",
                "positive mediator detunings",
            ],
            "excluded_new_term": "d_i^dag d_j times a rail-Pauli or density-conditioned factor",
        },
        "static_factorization": {
            "frobenius_residual_H_minus_Hq_tensor_I_minus_I_tensor_Hw": static_residual,
            "subtracted_scalar_energy": float(np.real_if_close(static_scalar)),
            "rail_factor": json_complex_matrix(hq),
            "walker_factor": json_complex_matrix(hw),
        },
        "complex_phase_floquet_factorization": {
            "segments": [
                {"j_perp": item[0], "g0": [float(item[1][0].real), float(item[1][0].imag)], "g1": [float(item[1][1].real), float(item[1][1].imag)], "duration": item[2]}
                for item in segments
            ],
            "maximum_segment_connected_residual": float(max(segment_residuals)),
            "accumulated_global_phase_exponent": global_phase_exponent,
            "full_piecewise_unitary_factorization_residual": factorized_unitary_residual,
        },
        "decision": (
            "In the smallest exact fixed-charge realization using only disjoint native rung hopping and phase-programmed "
            "charge-two conversions, the rail qubit and mediator walker factorize at every time. Complex channel phases and "
            "Floquet sequencing cannot produce a conditional neutral-walker Pauli link."
        ),
        "claim_boundary": (
            "This is a factorization obstruction for the disjoint-reservoir minimal block. It does not rule out a future "
            "non-separable coupling in which the mediator transition shares matter modes or a derived correlated-hopping path; "
            "such a construction must separately pass selectivity, leakage, symmetry and code-protection audits."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_native_walker_factorization_audit.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
