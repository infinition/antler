"""Four-rung explicit-mediator audit of the Phase 7D pulse compiler.

This is deliberately below a phase claim.  It evolves the complete fixed-charge
Hilbert space (eight rail modes plus twelve hard-core charge-two mediators),
including Mott costs, full Rabi pulses and rail rotations.  The inactive
channels can be given a controlled residual coupling epsilon*g to measure the
failure of the ideal closed-pulse/compiler contract.
"""
from __future__ import annotations

from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


L, LOW_MODES, LINKS = 4, 8, ((0, 1), (1, 2), (2, 3))
CHANNEL_KINDS = ("same_a", "same_b", "opposite_ab", "opposite_ba")
CHANNELS = tuple((link, kind) for link in range(len(LINKS)) for kind in CHANNEL_KINDS)
MEDIATORS, TOTAL_CHARGE = len(CHANNELS), 4
U_MOTT, DETUNING, G = 20.0, 40.0, 6.0
PULSE_TIME = 2.0 * np.pi / np.sqrt((DETUNING + 2.0 * U_MOTT) ** 2 + 4.0 * G ** 2)
EPSILONS = (0.0, 1e-4, 1e-3, 3e-3, 1e-2)
LEG_HOPPINGS = (0.0, 0.1, 0.3, 1.0)


def low_mode(rung: int, rail: int) -> int:
    return 2 * rung + rail


def sign_before(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def annihilate(state: int, mode: int) -> tuple[int, float] | None:
    if not ((state >> mode) & 1):
        return None
    return state ^ (1 << mode), sign_before(state, mode)


def create(state: int, mode: int) -> tuple[int, float] | None:
    if (state >> mode) & 1:
        return None
    return state | (1 << mode), sign_before(state, mode)


def pair_annihilation(state: int, first: int, second: int) -> tuple[int, float] | None:
    first, second = sorted((first, second))
    once = annihilate(state, first)
    if once is None:
        return None
    twice = annihilate(once[0], second)
    if twice is None:
        return None
    return twice[0], once[1] * twice[1]


def hopping(state: int, first: int, second: int) -> tuple[int, float] | None:
    """Apply c^dagger_second c_first."""
    once = annihilate(state, first)
    if once is None:
        return None
    twice = create(once[0], second)
    if twice is None:
        return None
    return twice[0], once[1] * twice[1]


def build_states() -> tuple[np.ndarray, dict[int, int]]:
    total_modes = LOW_MODES + MEDIATORS
    states = np.asarray([
        state for state in range(1 << total_modes)
        if (state & ((1 << LOW_MODES) - 1)).bit_count()
        + 2 * (state >> LOW_MODES).bit_count() == TOTAL_CHARGE
    ], dtype=np.int64)
    return states, {int(state): row for row, state in enumerate(states)}


STATES, INDEX = build_states()


def channel_pair(link: int, kind: str) -> tuple[int, int]:
    left, right = LINKS[link]
    if kind == "same_a":
        return low_mode(left, 0), low_mode(right, 0)
    if kind == "same_b":
        return low_mode(left, 1), low_mode(right, 1)
    if kind == "opposite_ab":
        return low_mode(left, 0), low_mode(right, 1)
    if kind == "opposite_ba":
        return low_mode(left, 1), low_mode(right, 0)
    raise ValueError(kind)


def sparse_from_entries(entries: dict[tuple[int, int], complex]) -> csr_matrix:
    rows, columns, values = zip(*((row, column, value) for (row, column), value in entries.items() if abs(value) > 0.0))
    return csr_matrix((values, (rows, columns)), shape=(len(STATES), len(STATES)), dtype=complex)


def pulse_hamiltonian(active_links: tuple[int, ...], kind_prefix: str, epsilon: float, leg_hopping: float = 0.0) -> csr_matrix:
    entries: dict[tuple[int, int], complex] = {}
    active = {(link, kind) for link in active_links for kind in CHANNEL_KINDS if kind.startswith(kind_prefix)}
    for column, raw_state in enumerate(STATES):
        state = int(raw_state)
        onsite = 0.0
        for rung in range(L):
            occupation = ((state >> low_mode(rung, 0)) & 1) + ((state >> low_mode(rung, 1)) & 1)
            onsite += U_MOTT * (occupation - 1) ** 2
        for mediator in range(MEDIATORS):
            onsite += DETUNING * ((state >> (LOW_MODES + mediator)) & 1)
        entries[(column, column)] = entries.get((column, column), 0.0) + onsite
        for mediator, channel in enumerate(CHANNELS):
            amplitude = G if channel in active else epsilon * G
            if amplitude == 0.0:
                continue
            mediator_mode = LOW_MODES + mediator
            if (state >> mediator_mode) & 1:
                continue
            item = pair_annihilation(state, *channel_pair(*channel))
            if item is None:
                continue
            low_state, sign = item
            final = low_state | (1 << mediator_mode)
            row = INDEX[final]
            value = -amplitude * sign
            entries[(row, column)] = entries.get((row, column), 0.0) + value
            entries[(column, row)] = entries.get((column, row), 0.0) + value
        if leg_hopping:
            for left_rung, right_rung in LINKS:
                for rail in (0, 1):
                    left, right = low_mode(left_rung, rail), low_mode(right_rung, rail)
                    for first, second in ((right, left), (left, right)):
                        item = hopping(state, first, second)
                        if item is not None:
                            final, fermion_sign = item
                            row = INDEX[final]
                            entries[(row, column)] = entries.get((row, column), 0.0) - leg_hopping * fermion_sign
    return sparse_from_entries(entries)


def rail_rotation(kind: str, rungs: tuple[int, ...], sign: float = 1.0) -> csr_matrix:
    entries: dict[tuple[int, int], complex] = {}
    for column, raw_state in enumerate(STATES):
        state = int(raw_state)
        for rung in rungs:
            a, b = low_mode(rung, 0), low_mode(rung, 1)
            if kind == "x":
                terms = ((b, a, 1.0), (a, b, 1.0))
            elif kind == "y":
                terms = ((b, a, -1j), (a, b, 1j))
            else:
                raise ValueError(kind)
            for first, second, amplitude in terms:
                item = hopping(state, first, second)
                if item is not None:
                    final, fermion_sign = item
                    row = INDEX[final]
                    entries[(row, column)] = entries.get((row, column), 0.0) + sign * amplitude * fermion_sign
    return sparse_from_entries(entries)


def evolve(state_vectors: np.ndarray, hamiltonian: csr_matrix, time: float) -> np.ndarray:
    return np.asarray(expm_multiply((-1j * time) * hamiltonian, state_vectors), dtype=complex)


def frame() -> np.ndarray:
    output = np.zeros((len(STATES), 1 << L), dtype=complex)
    for column, rails in enumerate(product((0, 1), repeat=L)):
        state = sum(1 << low_mode(rung, rail) for rung, rail in enumerate(rails))
        output[INDEX[state], column] = 1.0
    return output


FRAME = frame()
LOGICAL_PA = np.diag([(-1.0 if sum(rail == 0 for rail in rails) & 1 else 1.0) for rails in product((0, 1), repeat=L)]).astype(complex)
LOGICAL_PB = np.diag([(-1.0 if sum(rail == 1 for rail in rails) & 1 else 1.0) for rails in product((0, 1), repeat=L)]).astype(complex)


def pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], epsilon: float, leg_hopping: float = 0.0) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry = rail_rotation("y", rungs)
    rx = rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = evolve(vectors, pulse_hamiltonian(active_links, "same", epsilon, leg_hopping), PULSE_TIME)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = evolve(vectors, pulse_hamiltonian(active_links, "opposite", epsilon, leg_hopping), PULSE_TIME)
    return evolve(vectors, rx, -np.pi / 4.0)


def projected(vectors: np.ndarray) -> np.ndarray:
    return FRAME.conj().T @ vectors


def leakage(vectors: np.ndarray) -> float:
    return float(np.linalg.norm(vectors - FRAME @ projected(vectors), ord=2) ** 2)


def remove_global_phase(candidate: np.ndarray, reference: np.ndarray) -> np.ndarray:
    overlap = np.trace(reference.conj().T @ candidate)
    if abs(overlap) < 1e-15:
        return candidate
    return candidate * np.exp(-1j * np.angle(overlap))


def audit(epsilon: float) -> dict:
    even_state = pair_gate(FRAME.copy(), (0, 2), epsilon)
    complete_state = pair_gate(even_state, (1,), epsilon)
    g_even, g_odd, g_complete = projected(even_state), projected(pair_gate(FRAME.copy(), (1,), epsilon)), projected(complete_state)
    g01 = projected(pair_gate(FRAME.copy(), (0,), epsilon))
    g23 = projected(pair_gate(FRAME.copy(), (2,), epsilon))
    factorized_even = g23 @ g01
    product_logical = g_odd @ g_even
    return {
        "inactive_channel_coupling_over_g": epsilon,
        "even_layer_monomer_leakage": leakage(even_state),
        "complete_schedule_monomer_leakage": leakage(complete_state),
        "even_layer_factorization_residual": float(np.linalg.norm(g_even - remove_global_phase(factorized_even, g_even), ord=2)),
        "closed_pulse_composition_residual": float(np.linalg.norm(g_complete - remove_global_phase(product_logical, g_complete), ord=2)),
        "final_logical_parity_a_residual": float(np.linalg.norm(g_complete @ LOGICAL_PA - LOGICAL_PA @ g_complete, ord=2)),
        "final_logical_parity_b_residual": float(np.linalg.norm(g_complete @ LOGICAL_PB - LOGICAL_PB @ g_complete, ord=2)),
        "final_logical_singular_value_min": float(np.linalg.svd(g_complete, compute_uv=False)[-1]),
    }


def audit_leg_hopping(leg_hopping: float, reference: np.ndarray) -> dict:
    even_state = pair_gate(FRAME.copy(), (0, 2), 0.0, leg_hopping)
    complete_state = pair_gate(even_state, (1,), 0.0, leg_hopping)
    logical = projected(complete_state)
    return {
        "leg_hopping_during_pulses": leg_hopping,
        "complete_schedule_monomer_leakage": leakage(complete_state),
        "logical_deviation_from_zero_leg_schedule": float(np.linalg.norm(logical - remove_global_phase(reference, logical), ord=2)),
        "final_logical_parity_a_residual": float(np.linalg.norm(logical @ LOGICAL_PA - LOGICAL_PA @ logical, ord=2)),
        "final_logical_parity_b_residual": float(np.linalg.norm(logical @ LOGICAL_PB - LOGICAL_PB @ logical, ord=2)),
        "final_logical_singular_value_min": float(np.linalg.svd(logical, compute_uv=False)[-1]),
    }


def main() -> None:
    rows = [audit(epsilon) for epsilon in EPSILONS]
    zero_leg_even = pair_gate(FRAME.copy(), (0, 2), 0.0)
    zero_leg_reference = projected(pair_gate(zero_leg_even, (1,), 0.0))
    leg_hopping_rows = [audit_leg_hopping(value, zero_leg_reference) for value in LEG_HOPPINGS]
    out = {
        "schema": "antler.phase7d.four-rung-microscopic-pulse-audit.v1",
        "model": {
            "rungs": L, "low_fermion_modes": LOW_MODES,
            "hard_core_charge_two_mediators": MEDIATORS,
            "fixed_total_charge": TOTAL_CHARGE,
            "full_hilbert_dimension": int(len(STATES)),
            "schedule": "same-rail and opposite-rail closed pulses, with rail rotations; even links (0,1),(2,3) then odd link (1,2)",
        },
        "parameters": {
            "U_mott": U_MOTT, "detuning": DETUNING, "coupling": G,
            "coupling_over_detuning": G / DETUNING, "full_rabi_pulse_time": PULSE_TIME,
            "inactive_crosstalk_scan_over_g": list(EPSILONS),
            "leg_hopping_during_pulse_scan": list(LEG_HOPPINGS),
        },
        "rows": rows,
        "leg_hopping_rows": leg_hopping_rows,
        "decision": (
            "Microscopic four-rung pulse/crosstalk preflight only. The zero-crosstalk row tests closure and composition of the "
            "explicit mediator pulses; nonzero rows quantify this particular residual-channel error model. The separate leg-hopping "
            "scan shows that unrefocused hopping during opposite-rail pulses prevents exact stroboscopic closure. No many-body phase is inferred."
        ),
        "claim_boundary": (
            "The inactive-channel model is an idealized coherent coupling error, not a hardware noise model. This audit does not include "
            "a refocusing or recalibration solution for the observed leg-hopping failure, pulse calibration errors, dissipation, a thermodynamic gap, edge protection, 2D order, braiding, non-Abelian "
            "statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "four_rung_microscopic_pulse_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
