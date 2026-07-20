"""Exact Schrieffer--Wolff audit of a local coarse-block pair-link mediator.

The Phase-8B two-block transport used an inserted pair-transfer coefficient.
Here a single positive-detuned charge-two mediator is coupled coherently to a
pair on either of two adjacent b=2 blocks:

    V = g m^dagger (a_0 a_1 + a_2 a_3) + h.c.

In the rung-major Fock convention the low branch contains the derived term
-(g^2/Delta)(P_L^dagger P_R + h.c.) at SW order two, plus equal pair-projector
shifts.  The audit checks the exact finite block, its parity algebra and the
full Rabi dynamics rather than inserting that coefficient.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply


DETUNING = 10.0
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)
TIME_SAMPLES = 401


def block_a_parity(state: int, block: int) -> float:
    count = sum((state >> site_index(2 * block + offset, 0)) & 1 for offset in (0, 1))
    return -1.0 if count % 2 else 1.0


def build_model(coupling: float) -> tuple[sparse.csr_matrix, sparse.csr_matrix, np.ndarray, np.ndarray, dict]:
    """Fixed total charge two: matter N=2,m=0 plus matter vacuum,m=1."""
    states, positions = build_basis(8, 2)
    matter_dimension = len(states)
    mediator = matter_dimension
    dimension = matter_dimension + 1
    hamiltonian = sparse.lil_matrix((dimension, dimension), dtype=complex)
    pair_bridge = sparse.lil_matrix((dimension, dimension), dtype=complex)
    p_left = np.ones(dimension, dtype=float)
    p_right = np.ones(dimension, dtype=float)
    hamiltonian[mediator, mediator] = DETUNING
    pairs = (
        (site_index(0, 0), site_index(1, 0)),
        (site_index(2, 0), site_index(3, 0)),
    )
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        p_left[column] = block_a_parity(state, 0)
        p_right[column] = block_a_parity(state, 1)
        for first, second in pairs:
            # Rightmost operation acts first; the two annihilations retain
            # the frozen rung-major fermionic sign convention.
            item = _apply(state, (("ann", second), ("ann", first)))
            if item is None:
                continue
            new_state, amplitude = item
            if new_state != 0:
                continue
            hamiltonian[mediator, column] += coupling * amplitude
            hamiltonian[column, mediator] += coupling * np.conj(amplitude)
            pair_bridge[mediator, column] += coupling * amplitude
            pair_bridge[column, mediator] += coupling * np.conj(amplitude)

    # Negative control: a single a-particle hop between blocks.  It is not in
    # H and must be sector changing under both coarse parities.
    single = sparse.lil_matrix((dimension, dimension), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for operations in (
            (("ann", site_index(1, 0)), ("create", site_index(2, 0))),
            (("ann", site_index(2, 0)), ("create", site_index(1, 0))),
        ):
            item = _apply(state, operations)
            if item is not None:
                new_state, amplitude = item
                single[positions[new_state], column] += amplitude
    return hamiltonian.tocsr(), pair_bridge.tocsr(), p_left, p_right, {
        "dimension": dimension,
        "matter_dimension": matter_dimension,
        "mediator_index": mediator,
        "source_index": positions[(1 << site_index(0, 0)) | (1 << site_index(1, 0))],
        "target_index": positions[(1 << site_index(2, 0)) | (1 << site_index(3, 0))],
        "single_boundary": single.tocsr(),
    }


def frobenius(matrix: sparse.spmatrix) -> float:
    return float(np.sqrt(np.sum(np.abs(matrix.data) ** 2)))


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DETUNING
        hamiltonian, bridge, p_left, p_right, metadata = build_model(coupling)
        source = metadata["source_index"]
        target = metadata["target_index"]
        mediator = metadata["mediator_index"]
        p_l_op = sparse.diags(p_left, format="csr")
        p_r_op = sparse.diags(p_right, format="csr")
        physical = (p_left == 1.0) & (p_right == 1.0)
        p_physical = sparse.diags(physical.astype(float), format="csr")
        p_nonphysical = sparse.eye(hamiltonian.shape[0], format="csr") - p_physical

        # The coupled bright pair has E_-=(Delta-sqrt(Delta^2+8g^2))/2.
        # Its coefficient on |L>+|R> is E_-/2; SW predicts -g^2/Delta.
        bright_energy = (DETUNING - np.sqrt(DETUNING ** 2 + 8.0 * coupling ** 2)) / 2.0
        exact_pair_coefficient = bright_energy / 2.0
        sw_pair_coefficient = -(coupling ** 2) / DETUNING
        expected_time = np.pi / (2.0 * abs(sw_pair_coefficient))
        initial = np.zeros(hamiltonian.shape[0], dtype=complex)
        initial[source] = 1.0
        states_time = expm_multiply(
            -1j * hamiltonian, initial, start=0.0, stop=2.0 * expected_time,
            num=TIME_SAMPLES, endpoint=True,
        )
        target_population = np.abs(states_time[:, target]) ** 2
        peak = int(np.argmax(target_population))
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": coupling,
            "exact_pair_transfer_coefficient": float(exact_pair_coefficient),
            "sw_second_order_pair_transfer_coefficient": float(sw_pair_coefficient),
            "relative_sw_coefficient_error": float(abs(exact_pair_coefficient - sw_pair_coefficient) / abs(sw_pair_coefficient)),
            "equal_local_pair_projector_coefficient": float(exact_pair_coefficient),
            "hamiltonian_left_parity_commutator_frobenius": frobenius(hamiltonian @ p_l_op - p_l_op @ hamiltonian),
            "hamiltonian_right_parity_commutator_frobenius": frobenius(hamiltonian @ p_r_op - p_r_op @ hamiltonian),
            "physical_pair_bridge_frobenius": frobenius(p_physical @ bridge @ p_physical),
            "physical_single_boundary_projection_frobenius": frobenius(p_physical @ metadata["single_boundary"] @ p_physical),
            "single_boundary_leaves_physical_sector_frobenius": frobenius(p_nonphysical @ metadata["single_boundary"] @ p_physical),
            "expected_sw_transfer_time": float(expected_time),
            "maximum_target_pair_population": float(target_population[peak]),
            "time_of_maximum_target_pair_population": float(2.0 * expected_time * peak / (TIME_SAMPLES - 1)),
            "maximum_virtual_mediator_population": float(np.max(np.abs(states_time[:, mediator]) ** 2)),
            "maximum_parity_sector_leakage": float(np.max(1.0 - np.sum(np.abs(states_time[:, physical]) ** 2, axis=1))),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    error_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log([row["relative_sw_coefficient_error"] for row in deep]),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-pair-link-sw-audit.v1",
        "parameters": {
            "matter_modes": "four a-rail modes on two adjacent b=2 blocks in rung-major order",
            "fixed_total_charge": 2,
            "mediator_charge": 2,
            "detuning": DETUNING,
            "ratios": list(RATIOS),
            "pair_channels": "m^dagger a_0 a_1 and m^dagger a_2 a_3 plus Hermitian conjugates",
            "derived_low_energy_term": "-(g^2/Delta)[P_L^dag P_L + P_R^dag P_R + P_L^dag P_R + P_R^dag P_L]",
        },
        "dimensions": {"total": 29, "low_mediator_empty": 28, "charge_two_mediator": 1},
        "rows": rows,
        "deep_sw_relative_coefficient_error_power": error_power,
        "decision": (
            "A single positive-detuned charge-two mediator coherently coupled to the two adjacent local pair channels derives "
            "a parity-preserving pair link at SW order two. Its equal pair-projector shift is scalar in the one-mobile-pair "
            "link subspace but must be retained in any many-pair embedding."
        ),
        "claim_boundary": (
            "This is a 29-dimensional local bridge with a new four-leg pair-channel topology. It is not yet integrated with "
            "the Lambda walkers, derived from the frozen two-mediator-per-link ladder grammar, or shown to form a many-pair "
            "phase, local code, T junction, fusion space or non-Abelian braid."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_pair_link_sw_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
