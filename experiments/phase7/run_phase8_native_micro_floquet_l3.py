"""Three-rung microscopic Floquet bridge with explicit separate mediators.

This is the first composition test of the local Phase-8 SW bridge.  The model
has the theta=pi ANTLER leg hopping, two independent positive-detuned
charge-two mediators on every leg link, and an ideal switched global rung
pulse.  It compares the projected exact cycle to the published H_eff at the
same U0=-g^2/Delta.

The rung pulse is intentionally ideal (only J_perp acts during it).  Thus a
passing result is a microscopic finite-block bridge, not yet a hardware-pulse
derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm, svd


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply
from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation


LENGTH, PARTICLE_NUMBER, T_LEG = 3, 2, 1.0
ALPHA, ETA, RATIO, TARGET_U0 = 0.5, np.pi / 2.0, 0.05, -1.5
PERIODS = (0.20, 0.10, 0.05, 0.025)


def mediator_mode(link: int, rail: int) -> int:
    return 2 * LENGTH + 2 * link + rail


def weighted_basis() -> tuple[np.ndarray, dict[int, int]]:
    mode_charges = [1] * (2 * LENGTH) + [2] * (2 * (LENGTH - 1))
    states = np.asarray([
        state for state in range(1 << len(mode_charges))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(mode_charges)) == PARTICLE_NUMBER
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_micro_h0(g: float, detuning: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    states, index = weighted_basis()
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    jx = np.zeros_like(hamiltonian)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for link in range(LENGTH - 1):
            for rail in (0, 1):
                first, second = site_index(link, rail), site_index(link + 1, rail)
                mediator = mediator_mode(link, rail)
                hamiltonian[column, column] += detuning * ((state >> mediator) & 1)
                for operations in ((("ann", second), ("create", first)), (("ann", first), ("create", second))):
                    item = _apply(state, operations)
                    if item is not None:
                        new_state, amplitude = item
                        hamiltonian[index[new_state], column] += -T_LEG * amplitude
                if ((state >> first) & 1) and ((state >> second) & 1) and not ((state >> mediator) & 1):
                    item = _apply(state, (("ann", second), ("ann", first), ("create", mediator)))
                    if item is None:
                        raise RuntimeError("valid conversion vanished")
                    new_state, amplitude = item
                    hamiltonian[index[new_state], column] += -g * amplitude
                    hamiltonian[column, index[new_state]] += -g * amplitude
        for rung in range(LENGTH):
            a, b = site_index(rung, 0), site_index(rung, 1)
            for operations in ((("ann", b), ("create", a)), (("ann", a), ("create", b))):
                item = _apply(state, operations)
                if item is not None:
                    new_state, amplitude = item
                    jx[index[new_state], column] += 0.5 * amplitude
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("microscopic H0 is not Hermitian")
    if not np.allclose(jx, jx.conj().T, atol=1e-12):
        raise RuntimeError("microscopic Jx is not Hermitian")
    low_positions = np.asarray([
        position for position, state in enumerate(states)
        if all(not ((int(state) >> mediator_mode(link, rail)) & 1) for link in range(LENGTH - 1) for rail in (0, 1))
    ], dtype=int)
    frame = np.zeros((len(states), len(low_positions)), dtype=complex)
    frame[low_positions, np.arange(len(low_positions))] = 1.0
    return hamiltonian, jx, states, frame


def polar_unitary(matrix: np.ndarray) -> np.ndarray:
    left, _, right_h = svd(matrix)
    return left @ right_h


def phase_aligned_distance(actual: np.ndarray, target: np.ndarray) -> float:
    phase = np.angle(np.trace(target.conj().T @ actual))
    return float(np.linalg.norm(actual * np.exp(-1j * phase) - target) / np.sqrt(target.shape[0]))


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    h_micro, jx_micro, states, frame = build_micro_h0(g, detuning)
    h_target, p_target, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("low physical frame and target basis disagree")
    p_micro = expm(-1j * ETA * jx_micro)
    if np.linalg.norm(frame.conj().T @ p_micro @ frame - p_target) > 1e-12:
        raise RuntimeError("microscopic rail pulse failed exact low-frame check")
    # Align the exact low-energy eigenframe to the bare physical frame. This
    # is the finite-dimensional Schrieffer--Wolff dressing: virtual mediator
    # occupation is part of a low-energy state rather than a leakage event.
    _, eigenvectors = np.linalg.eigh(h_micro)
    low_eigenframe = eigenvectors[:, :frame.shape[1]]
    left, singular_values, right_h = svd(frame.conj().T @ low_eigenframe)
    dressed_frame = low_eigenframe @ right_h.conj().T @ left.conj().T
    if np.linalg.norm(dressed_frame.conj().T @ dressed_frame - np.eye(frame.shape[1])) > 1e-10:
        raise RuntimeError("dressed low frame is not orthonormal")
    full_pa = np.diag([
        -1.0 if sum((int(state) >> site_index(rung, 0)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in states
    ])
    pa = frame.conj().T @ full_pa @ frame
    identity = np.eye(h_micro.shape[0], dtype=complex)
    low_projector = frame @ frame.conj().T
    rows, dressed_rows = [], []
    for period in PERIODS:
        cycle = (
            p_micro.conj().T
            @ expm(-1j * (1.0 - ALPHA) * period * h_micro)
            @ p_micro
            @ expm(-1j * ALPHA * period * h_micro)
        )
        target = expm(-1j * period * (ALPHA * h_target + (1.0 - ALPHA) * (p_target.conj().T @ h_target @ p_target)))
        raw_logical = frame.conj().T @ cycle @ frame
        logical = polar_unitary(raw_logical)
        leakage = float(np.linalg.norm((identity - low_projector) @ cycle @ frame, ord=2) ** 2)
        rows.append({
            "period": period,
            "low_frame_leakage_worst": leakage,
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, target),
            "logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical @ pa - pa @ logical) / np.sqrt(pa.shape[0])),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw_logical, compute_uv=False))),
        })
        raw_dressed = dressed_frame.conj().T @ cycle @ dressed_frame
        logical_dressed = polar_unitary(raw_dressed)
        dressed_projector = dressed_frame @ dressed_frame.conj().T
        pa_dressed = dressed_frame.conj().T @ full_pa @ dressed_frame
        dressed_rows.append({
            "period": period,
            "dressed_frame_leakage_worst": float(np.linalg.norm((identity - dressed_projector) @ cycle @ dressed_frame, ord=2) ** 2),
            "dressed_polar_logical_vs_target_distance": phase_aligned_distance(logical_dressed, target),
            "dressed_logical_branch_parity_commutator_normalized": float(np.linalg.norm(logical_dressed @ pa_dressed - pa_dressed @ logical_dressed) / np.sqrt(pa.shape[0])),
            "dressed_raw_min_singular_value": float(np.min(np.linalg.svd(raw_dressed, compute_uv=False))),
        })
    out = {
        "schema": "antler.phase8.native-micro-floquet-l3.v1",
        "parameters": {
            "L": LENGTH,
            "N": PARTICLE_NUMBER,
            "theta": "pi",
            "t_leg": T_LEG,
            "target_u0": TARGET_U0,
            "g_over_detuning": RATIO,
            "detuning": detuning,
            "g": g,
            "alpha": ALPHA,
            "eta": "pi/2",
            "mediators": "one independent charge-two mediator per rail per link",
        },
        "physical_dimension": int(h_micro.shape[0]),
        "low_frame_dimension": int(frame.shape[1]),
        "rail_pulse_low_frame_residual": float(np.linalg.norm(frame.conj().T @ p_micro @ frame - p_target)),
        "minimum_bare_to_dressed_frame_overlap": float(min(singular_values)),
        "rows": rows,
        "dressed_rows": dressed_rows,
        "decision": "Exact L=3 microscopic composition audit with an ideal switched rung pulse.",
        "claim_boundary": "This small coherent block does not establish many-link convergence, physical pulse switching, noise robustness, a thermodynamic phase, braid, non-Abelian statistics, universality, or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "native_micro_floquet_l3.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
