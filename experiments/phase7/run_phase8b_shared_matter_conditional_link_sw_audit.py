"""Shared-matter Schrieffer-Wolff route to a conditional walker link.

This is the first non-factorizing microscopic candidate after the disjoint
reservoir no-go.  A rail qubit q=(q_a,q_b), a neighboring physical pair
(r_a,r_b), and two charge-two mediator species d0,d1 share the same local
matter block.  Pair conversion produces the walker transition, while virtual
leg hopping dresses it by the rail X.

Two static microscopic segments are related by controls already registered as
dynamic resources: a Peierls sign flip of J_perp and a pi phase flip of one
pair channel.  Their *Schrieffer-Wolff effective-Hamiltonian average* cancels
the isolated rail and walker terms and retains W_X X_q.  This is a controlled
static/SW compiler audit, not yet a finite-pulse Floquet-gate claim.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

CHARGES = (1, 1, 1, 1, 2, 2)
TOTAL_CHARGE = 3
Q_A, Q_B, R_A, R_B, D_0, D_1 = range(6)
DELTA_PAIR = 10.0
MOTT_Q = 15.0
LEG_HOPPING = 0.7
RUNG_HOPPING = 0.2
THETA = np.pi
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)
PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.diag((1.0, -1.0)).astype(complex),
}


def weighted_basis() -> tuple[np.ndarray, dict[int, int]]:
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


def add_directed_hop(
    hamiltonian: np.ndarray,
    states: np.ndarray,
    positions: dict[int, int],
    destination: int,
    source: int,
    amplitude: complex,
    *,
    conditional_mode: int | None = None,
) -> None:
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        item = annihilate(state, source)
        if item is None:
            continue
        intermediate, sign = item
        item = create(intermediate, destination)
        if item is None:
            continue
        final, final_sign = item
        phase = 1.0 if conditional_mode is None else np.exp(1j * THETA * ((state >> conditional_mode) & 1))
        hamiltonian[positions[final], column] += amplitude * phase * sign * final_sign


def add_directed_pair_conversion(
    hamiltonian: np.ndarray,
    states: np.ndarray,
    positions: dict[int, int],
    mediator: int,
    amplitude: complex,
) -> None:
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        item = annihilate(state, R_B)
        if item is None:
            continue
        intermediate, sign = item
        item = annihilate(intermediate, R_A)
        if item is None:
            continue
        intermediate, second_sign = item
        item = create(intermediate, mediator)
        if item is None:
            continue
        final, final_sign = item
        hamiltonian[positions[final], column] += amplitude * sign * second_sign * final_sign


def build_segment(
    j_perp: complex,
    g0: complex,
    g1: complex,
    *,
    rail_bias: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, dict[int, int]]:
    states, positions = weighted_basis()
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    for position, raw_state in enumerate(states):
        state = int(raw_state)
        n_q = ((state >> Q_A) & 1) + ((state >> Q_B) & 1)
        hamiltonian[position, position] += MOTT_Q * (n_q - 1) ** 2
        hamiltonian[position, position] += DELTA_PAIR * ((state >> R_A) & 1) * ((state >> R_B) & 1)
        hamiltonian[position, position] += rail_bias * (((state >> Q_A) & 1) - ((state >> Q_B) & 1))
    # The rung hop is native and unconditioned in rung-major ordering.
    add_directed_hop(hamiltonian, states, positions, Q_A, Q_B, -j_perp)
    # These are the two correlated leg hops in the frozen rung-major convention.
    add_directed_hop(hamiltonian, states, positions, Q_A, R_A, -LEG_HOPPING, conditional_mode=Q_B)
    add_directed_hop(hamiltonian, states, positions, Q_B, R_B, -LEG_HOPPING, conditional_mode=R_A)
    add_directed_pair_conversion(hamiltonian, states, positions, D_0, -g0)
    add_directed_pair_conversion(hamiltonian, states, positions, D_1, -g1)
    hamiltonian = hamiltonian + hamiltonian.conj().T
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("shared-matter block is not Hermitian")
    return hamiltonian, states, positions


def code_indices(positions: dict[int, int]) -> np.ndarray:
    return np.asarray([
        positions[(1 << rail) | (1 << mediator)]
        for rail in (Q_A, Q_B) for mediator in (D_0, D_1)
    ], dtype=int)


def schur_effective(hamiltonian: np.ndarray, low: np.ndarray) -> tuple[np.ndarray, float, float]:
    high = np.asarray([index for index in range(hamiltonian.shape[0]) if index not in set(low)], dtype=int)
    h_ll = hamiltonian[np.ix_(low, low)]
    h_lh = hamiltonian[np.ix_(low, high)]
    h_hh = hamiltonian[np.ix_(high, high)]
    h_hl = hamiltonian[np.ix_(high, low)]
    effective = h_ll - h_lh @ np.linalg.solve(h_hh, h_hl)
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    low_vectors = eigenvectors[:, :4]
    frame = np.zeros((hamiltonian.shape[0], 4), dtype=complex)
    frame[low, np.arange(4)] = 1.0
    capture = float(np.min(np.linalg.svd(frame.conj().T @ low_vectors, compute_uv=False)) ** 2)
    gap = float(eigenvalues[4] - eigenvalues[3])
    return effective, capture, gap


def pauli_coefficients(matrix: np.ndarray) -> dict[str, float]:
    return {
        left + right: float(np.real_if_close(np.trace(np.kron(PAULIS[left], PAULIS[right]) @ matrix) / 4.0))
        for left, right in itertools.product(PAULIS, repeat=2)
    }


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        # A: p=+1, q=+1.  B: p=-1 (relative pair-channel phase pi), q=-1
        # (Peierls sign flip of the rung hop).  IX and XI are odd under one
        # flip; XX is even under both and survives their average.
        h_a, _, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        h_b, _, _ = build_segment(-RUNG_HOPPING, coupling, -coupling)
        low = code_indices(positions)
        effective_a, capture_a, gap_a = schur_effective(h_a, low)
        effective_b, capture_b, gap_b = schur_effective(h_b, low)
        echoed = 0.5 * (effective_a + effective_b)
        coefficients = pauli_coefficients(echoed)
        unwanted = {label: value for label, value in coefficients.items() if label not in {"II", "XX"}}
        largest_label, largest_value = max(unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "echoed_xx_coefficient": coefficients["XX"],
            "echoed_ix_coefficient": coefficients["IX"],
            "echoed_xi_coefficient": coefficients["XI"],
            "maximum_unwanted_non_scalar_coefficient": abs(largest_value),
            "largest_unwanted_non_scalar_pauli": largest_label,
            "minimum_low_frame_capture": min(capture_a, capture_b),
            "minimum_low_high_gap": min(gap_a, gap_b),
            "echoed_schur_hermiticity_residual": float(np.linalg.norm(echoed - echoed.conj().T, ord="fro")),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["echoed_xx_coefficient"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-conditional-link-sw-audit.v1",
        "parameters": {
            "mode_order": ["q_a", "q_b", "r_a", "r_b", "d0", "d1"],
            "weighted_total_charge": TOTAL_CHARGE,
            "theta": THETA,
            "pair_detuning": DELTA_PAIR,
            "mott_q": MOTT_Q,
            "leg_hopping": LEG_HOPPING,
            "rung_hopping_magnitude": RUNG_HOPPING,
            "segments": {
                "A": "J_perp=+J, g0=+g, g1=+g",
                "B": "J_perp=-J, g0=+g, g1=-g",
            },
            "target": "X_q tensor X_(d0,d1)",
            "control_type": "registered Peierls rung-sign flip plus coherent relative pair-channel phase pi",
        },
        "rows": rows,
        "deep_sw_echoed_xx_power": power,
        "decision": (
            "The non-separable shared-matter block produces a conditional X_q X_walker term. Averaging the two sign-correlated "
            "Schrieffer-Wolff segments cancels the isolated X_q and X_walker terms while retaining the conditional term."
        ),
        "claim_boundary": (
            "This is a static/downfolded effective-Hamiltonian compiler result only. It does not yet demonstrate a finite-pulse "
            "Floquet sequence with low leakage, derive a four-state neutral walker, establish a Y-conditioned link, integrate a "
            "full stabilizer patch, derive the new mediators from frozen ANTLER without extension, or establish twists, fusion, "
            "non-Abelian braiding, universality or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_shared_matter_conditional_link_sw_audit.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
