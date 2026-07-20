"""Controlled Schrieffer--Wolff bridge from ANTLER resources to Floquet H0.

The Floquet reference needs two *decoupled* attractive density interactions
U0 n_{a,j}n_{a,j+1} and U0 n_{b,j}n_{b,j+1}.  A shared charge-two mediator
would additionally generate pair hopping at the same order.  This audit uses
one positive-detuned hard-core charge-two mediator per rail and per link,
which instead gives the desired density terms at O(g^2/Delta) and no
cross-rail pair transfer at that order.

The audit also checks that the frozen rung-major ANTLER hopping at theta=pi
equals the ordinary fermionic leg hopping used by the external Floquet model.
It is a local controlled derivation, not an extended-ladder phase claim.
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
from antler.model import build_hamiltonian
from antler.number_conserving_pairwire import _apply
from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation


LENGTH, PARTICLE_NUMBER, T_LEG, TARGET_U0 = 2, 2, 1.0, -1.5
RATIOS = (0.05, 0.10, 0.15)
MODE_CHARGES = (1, 1, 1, 1, 2, 2)  # a0,b0,a1,b1,d_a,d_b


def weighted_basis(total_charge: int) -> tuple[np.ndarray, dict[int, int]]:
    states = np.asarray([
        state for state in range(1 << len(MODE_CHARGES))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(MODE_CHARGES)) == total_charge
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_separate_mediator_block(g: float, detuning: float) -> tuple[np.ndarray, np.ndarray, list[int], list[int]]:
    """Two-rung physical ladder plus one charge-two mediator per rail."""
    states, index = weighted_basis(PARTICLE_NUMBER)
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    physical_pairs = ((site_index(0, 0), site_index(1, 0), 4), (site_index(0, 1), site_index(1, 1), 5))
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for mediator in (4, 5):
            hamiltonian[column, column] += detuning * ((state >> mediator) & 1)
        # Physical intraleg hopping on the rung-major fermionic convention.
        for left, right in ((site_index(0, 0), site_index(1, 0)), (site_index(0, 1), site_index(1, 1))):
            for operations in ((("ann", right), ("create", left)), (("ann", left), ("create", right))):
                item = _apply(state, operations)
                if item is not None:
                    new_state, amplitude = item
                    hamiltonian[index[new_state], column] += -T_LEG * amplitude
        # Pair-to-mediator conversion. Each pair has its own mediator: this is
        # the structural condition which removes the cross-rail pair transfer.
        for first, second, mediator in physical_pairs:
            if ((state >> first) & 1) and ((state >> second) & 1) and not ((state >> mediator) & 1):
                item = _apply(state, (("ann", second), ("ann", first), ("create", mediator)))
                if item is None:
                    raise RuntimeError("valid pair conversion unexpectedly vanished")
                new_state, amplitude = item
                hamiltonian[index[new_state], column] += -g * amplitude
                hamiltonian[column, index[new_state]] += -g * amplitude
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("separate-mediator block is not Hermitian")
    low = [position for position, state in enumerate(states) if not ((int(state) >> 4) & 1) and not ((int(state) >> 5) & 1)]
    high = [position for position in range(len(states)) if position not in low]
    return hamiltonian, states, low, high


def target_h0(interaction: float) -> tuple[np.ndarray, np.ndarray]:
    matrix, _, states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, interaction)
    return matrix, states


def analyze(ratio: float) -> dict:
    # Keep the desired U0 fixed while taking the controlled g/Delta limit.
    detuning = abs(TARGET_U0) / ratio**2
    g = ratio * detuning
    full, states, low, high = build_separate_mediator_block(g, detuning)
    target, target_states = target_h0(TARGET_U0)
    low_states = np.asarray([states[position] for position in low], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("low physical basis no longer agrees with target Fock ordering")
    h_pp = full[np.ix_(low, low)]
    h_pq = full[np.ix_(low, high)]
    h_qq = full[np.ix_(high, high)]
    h_sw = h_pp - h_pq @ np.linalg.inv(h_qq) @ h_pq.conj().T
    # |a0 a1> and |b0 b1> are the potential unwanted pair-transfer entries.
    low_index = {int(state): position for position, state in enumerate(low_states)}
    aa = low_index[(1 << site_index(0, 0)) | (1 << site_index(1, 0))]
    bb = low_index[(1 << site_index(0, 1)) | (1 << site_index(1, 1))]
    full_values, full_vectors = eigh(full)
    low_values = full_values[:len(low)]
    capture = [float(np.sum(np.abs(full_vectors[low, column]) ** 2)) for column in range(len(low))]
    target_values = np.linalg.eigvalsh(target)
    return {
        "g_over_detuning": ratio,
        "detuning": detuning,
        "g": g,
        "target_u0": TARGET_U0,
        "sw_target_frobenius_residual": float(np.linalg.norm(h_sw - target)),
        "sw_unwanted_aa_to_bb_abs": float(abs(h_sw[aa, bb])),
        "exact_low_spectrum_max_abs_error": float(np.max(np.abs(low_values - target_values))),
        "minimum_low_manifold_capture": float(min(capture)),
    }


def main() -> None:
    antler_hop, antler_states, _ = build_hamiltonian(LENGTH, PARTICLE_NUMBER, np.pi, T_LEG, T_LEG, 0.0)
    reference_hop, reference_states = target_h0(0.0)
    if not np.array_equal(antler_states, reference_states):
        raise RuntimeError("ANTLER and reference hopping bases disagree")
    # With a pi Peierls sign on the rung coupling, Jperp=-1/2 gives the
    # generator +Jx used by P=exp(-i eta Jx). The opposite sign gives P^dag.
    antler_jx, rotation_states, _ = build_hamiltonian(LENGTH, PARTICLE_NUMBER, np.pi, 0.0, 0.0, -0.5)
    _, reference_rotation, rotation_reference_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, 0.0)
    if not np.array_equal(rotation_states, rotation_reference_states):
        raise RuntimeError("ANTLER and reference rail-pulse bases disagree")
    native_rotation = expm(-1j * (np.pi / 2.0) * antler_jx)
    rows = [analyze(ratio) for ratio in RATIOS]
    exponent = float(np.polyfit(
        np.log([row["g_over_detuning"] for row in rows]),
        np.log([row["exact_low_spectrum_max_abs_error"] for row in rows]),
        1,
    )[0])
    out = {
        "schema": "antler.phase8.native-h0-sw-bridge.v1",
        "target": {
            "h0": "-t sum_leg(c†_j c_{j+1}+h.c.) + U0 sum_leg n_j n_{j+1}",
            "u0": TARGET_U0,
            "separate_charge_two_mediators_per_link": 2,
            "frozen_antler_hopping_setting": {"theta": "pi", "J1": T_LEG, "J2": T_LEG, "Jperp_free": 0.0},
        },
        "antler_theta_pi_hopping_frobenius_residual": float(np.linalg.norm(antler_hop - reference_hop)),
        "global_rail_rotation": {
            "native_pulse_setting": {"theta": "pi", "J1": 0.0, "J2": 0.0, "Jperp": -0.5, "duration": "pi/2"},
            "native_generator_vs_Jx_frobenius_residual": float(np.linalg.norm(antler_jx + 0.5 * build_hamiltonian(LENGTH, PARTICLE_NUMBER, np.pi, 0.0, 0.0, 1.0)[0])),
            "native_rotation_vs_reference_P_frobenius_residual": float(np.linalg.norm(native_rotation - reference_rotation)),
            "control_requirement": "a globally synchronous pi Peierls sign on Jperp, or an exactly equivalent reversed pulse, is required",
        },
        "rows": rows,
        "exact_spectrum_error_loglog_slope": exponent,
        "decision": (
            "A separate positive-detuned charge-two mediator on each rail gives the required attractive density interaction "
            "at second order without same-order cross-rail pair hopping. Together with theta=pi uniform ANTLER leg hopping, "
            "this is a controlled local bridge to H0."
        ),
        "claim_boundary": (
            "This is a two-rung Schrieffer--Wolff derivation. It has not yet established a many-link mediator schedule, "
            "global rail pulse isolation, finite-frequency full-ladder fidelity, a thermodynamic phase, braid, non-Abelian "
            "statistics, universality, or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "native_h0_sw_bridge.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
