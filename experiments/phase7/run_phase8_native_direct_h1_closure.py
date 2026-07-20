"""Direct rotated-mediator implementation of the second Floquet segment.

Instead of implementing H1=P^dag H0 P by a sudden physical rail rotation,
factor its rotated two-body interaction into two coherent charge-two mediator
channels per link.  H0 and H1 then act on the same bare low-energy frame;
integer virtual-Rabi closure can be imposed on both free halves.

This is a new *controlled dynamic resource proposal*, not a feature of the
frozen static ANTLER Hamiltonian.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh, expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply
from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation
from run_phase8_native_micro_floquet_l3 import (
    ALPHA, LENGTH, PARTICLE_NUMBER, RATIO, TARGET_U0, phase_aligned_distance,
    polar_unitary,
)


def pair_state_indices(states: np.ndarray, length: int) -> list[int]:
    pairs = (
        (site_index(0, 0), site_index(1, 0)),  # aa
        (site_index(0, 0), site_index(1, 1)),  # ab
        (site_index(0, 1), site_index(1, 0)),  # ba
        (site_index(0, 1), site_index(1, 1)),  # bb
    )
    index = {int(state): position for position, state in enumerate(states)}
    return [index[(1 << first) | (1 << second)] for first, second in pairs]


def rotated_channels() -> tuple[np.ndarray, dict]:
    """Factor -P^dag H_int P into two normalized pair-annihilation rows."""
    h_density, rotation, states, _ = build_h0_and_rotation(2, 2, -1.0)
    h_hop, _, _, _ = build_h0_and_rotation(2, 2, 0.0)
    h1_interaction = rotation.conj().T @ (h_density - h_hop) @ rotation
    rows = pair_state_indices(states, 2)
    target = h1_interaction[np.ix_(rows, rows)]
    values, vectors = eigh(-target)
    selected = np.where(values > 1e-10)[0]
    if len(selected) != 2:
        raise RuntimeError(f"expected rank-two rotated density interaction, got eigenvalues {values}")
    channels = np.asarray([np.conj(vectors[:, col]) * np.sqrt(values[col]) for col in selected], dtype=complex)
    reconstructed = -channels.conj().T @ channels
    return channels, {
        "rotated_interaction_rank": int(len(selected)),
        "eigenvalues": [float(value) for value in values],
        "factorization_frobenius_residual": float(np.linalg.norm(reconstructed - target)),
    }


def mediator_mode(length: int, link: int, rail: int) -> int:
    return 2 * length + 2 * link + rail


def weighted_basis(length: int, particle_number: int) -> tuple[np.ndarray, dict[int, int]]:
    """Bare fermions carry charge one; every hard-core mediator carries two."""
    mode_charges = [1] * (2 * length) + [2] * (2 * (length - 1))
    states = np.asarray([
        state for state in range(1 << len(mode_charges))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(mode_charges)) == particle_number
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_micro(length: int, particle_number: int, channels: np.ndarray, g: float, detuning: float,
                nearest_link_crosstalk: complex = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the explicit-mediator block.

    ``nearest_link_crosstalk`` is a dimensionless residual conversion amplitude:
    a mediator assigned to link j also addresses pair configurations on j-1 and
    j+1 with amplitude epsilon*g.  It is zero for the registered bridge and
    exists solely for finite-block control audits.
    """
    if not np.isscalar(nearest_link_crosstalk):
        raise TypeError("nearest_link_crosstalk must be a scalar amplitude")
    states, index = weighted_basis(length, particle_number)
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for link in range(length - 1):
            for rail in (0, 1):
                first, second = site_index(link, rail), site_index(link + 1, rail)
                for operations in ((("ann", second), ("create", first)), (("ann", first), ("create", second))):
                    item = _apply(state, operations)
                    if item is not None:
                        new_state, amplitude = item
                        hamiltonian[index[new_state], column] += -1.0 * amplitude
            for mediator_slot in (0, 1):
                mediator = mediator_mode(length, link, mediator_slot)
                hamiltonian[column, column] += detuning * ((state >> mediator) & 1)
                if (state >> mediator) & 1:
                    continue
                addressed_links = [(link, 1.0)]
                if nearest_link_crosstalk:
                    addressed_links.extend(
                        (other_link, nearest_link_crosstalk)
                        for other_link in (link - 1, link + 1)
                        if 0 <= other_link < length - 1
                    )
                for pair_link, scale in addressed_links:
                    pair_modes = (
                        (site_index(pair_link, 0), site_index(pair_link + 1, 0)),
                        (site_index(pair_link, 0), site_index(pair_link + 1, 1)),
                        (site_index(pair_link, 1), site_index(pair_link + 1, 0)),
                        (site_index(pair_link, 1), site_index(pair_link + 1, 1)),
                    )
                    for pair_slot, (first, second) in enumerate(pair_modes):
                        if not ((state >> first) & 1 and (state >> second) & 1):
                            continue
                        item = _apply(state, (("ann", second), ("ann", first), ("create", mediator)))
                        if item is None:
                            raise RuntimeError("valid pair conversion vanished")
                        new_state, amplitude = item
                        value = -g * scale * channels[mediator_slot, pair_slot] * amplitude
                        hamiltonian[index[new_state], column] += value
                        hamiltonian[column, index[new_state]] += np.conj(value)
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("direct-channel microscopic Hamiltonian is not Hermitian")
    low = np.asarray([
        position for position, state in enumerate(states)
        if all(not ((int(state) >> mediator_mode(length, link, rail)) & 1) for link in range(length - 1) for rail in (0, 1))
    ], dtype=int)
    frame = np.zeros((len(states), len(low)), dtype=complex)
    frame[low, np.arange(len(low))] = 1.0
    return hamiltonian, states, frame


def metrics(unitary: np.ndarray, target: np.ndarray, frame: np.ndarray, projector: np.ndarray, identity: np.ndarray, pa: np.ndarray) -> dict:
    raw = frame.conj().T @ unitary @ frame
    logical = polar_unitary(raw)
    return {
        "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
        "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
        "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
        "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
    }


def main() -> None:
    channels_h1, factor = rotated_channels()
    channels_h0 = np.asarray(((1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0)), dtype=complex)
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 4.0 * np.pi / omega
    h0_micro, states, frame = build_micro(LENGTH, PARTICLE_NUMBER, channels_h0, g, detuning)
    h1_micro, states_h1, frame_h1 = build_micro(LENGTH, PARTICLE_NUMBER, channels_h1, g, detuning)
    if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
        raise RuntimeError("microsegment bases disagree")
    target_h0, target_p, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, col])))] for col in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("low physical frame disagrees with target basis")
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    # Check each low-space SW interaction factorization directly at L=2 via
    # the local construction, then audit the L=3 cycle itself.
    projector = frame @ frame.conj().T
    identity = np.eye(h0_micro.shape[0], dtype=complex)
    pa = np.diag([
        -1.0 if sum((int(state) >> site_index(rung, 0)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in low_states
    ])
    direct_cycle = expm(-1j * (1.0 - ALPHA) * period * h1_micro) @ expm(-1j * ALPHA * period * h0_micro)
    target_cycle = expm(-1j * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    rows = []
    for count in (1, 2, 4, 8):
        rows.append({"cycles": count, **metrics(
            np.linalg.matrix_power(direct_cycle, count),
            np.linalg.matrix_power(target_cycle, count),
            frame, projector, identity, pa,
        )})
    out = {
        "schema": "antler.phase8.native-direct-h1-closure.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "g_over_detuning": RATIO, "detuning": detuning, "g": g, "omega_pair_mediator": omega, "closure_period": period, "alpha": ALPHA, "segments": "direct H0 and direct P†H0P mediator channels"},
        "rotated_channel_factorization": {
            **factor,
            "channels": [[{"real": float(value.real), "imag": float(value.imag)} for value in row] for row in channels_h1],
        },
        "rows": rows,
        "decision": "Direct switched-mediator realization of the two leading Floquet segments at simultaneous virtual-Rabi closure.",
        "claim_boundary": "This adds coherent rotated pair-conversion channels as a dynamic resource. It is an L=3 ideal-pulse test, not a frozen-ANTLER, hardware, thermodynamic, braid, non-Abelian, universal or fault-tolerant result.",
    }
    path = ROOT / "results" / "phase7" / "native_direct_h1_closure.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
