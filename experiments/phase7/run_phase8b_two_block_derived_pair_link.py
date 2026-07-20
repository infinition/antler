"""Integrate the derived charge-two pair link with two coarse-Gauss walkers.

This is the first Phase-8B two-block calculation that replaces the former
inserted P_L^dagger P_R coefficient by an explicit positive-detuned charge-two
mediator.  The neutral Lambda walkers remain the separately declared new
Gauss resource.  Total charge is fixed to two: the pair-link mediator carries
charge two when the matter block is empty.
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
RATIOS = (0.10, 0.05, 0.025)
TIME_SAMPLES = 301


def basis_index(matter: int, gauge: int, walker0: int, walker1: int, matter_dimension: int) -> int:
    return ((((walker1 * 4) + walker0) * 8 + gauge) * matter_dimension) + matter


def x_eigenvalue(gauge: int, link: int) -> float:
    return -1.0 if (gauge >> link) & 1 else 1.0


def block_parity(state: int, block: int) -> float:
    count = sum((state >> site_index(2 * block + offset, 0)) & 1 for offset in (0, 1))
    return -1.0 if count % 2 else 1.0


def frobenius(matrix: sparse.spmatrix) -> float:
    return float(np.sqrt(np.sum(np.abs(matrix.data) ** 2)))


def build_model(coupling: float) -> tuple[sparse.csr_matrix, sparse.csr_matrix, np.ndarray, np.ndarray, dict]:
    """Two b=2 blocks, two neutral walkers and one charge-two link mediator."""
    two_particle_states, positions = build_basis(8, 2)
    states = [int(state) for state in two_particle_states]
    vacuum_link_mediator = len(states)
    matter_dimension = len(states) + 1
    dimension = matter_dimension * 8 * 4 * 4
    hamiltonian = sparse.lil_matrix((dimension, dimension), dtype=complex)
    pair_link = sparse.lil_matrix((dimension, dimension), dtype=complex)
    single_boundary = sparse.lil_matrix((dimension, dimension), dtype=complex)
    g0_values = np.empty(dimension, dtype=float)
    g1_values = np.empty(dimension, dtype=float)
    source_matter = positions[(1 << site_index(0, 0)) | (1 << site_index(1, 0))]
    target_matter = positions[(1 << site_index(2, 0)) | (1 << site_index(3, 0))]

    def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
        matrix[row, column] += value
        matrix[column, row] += np.conj(value)

    for matter in range(matter_dimension):
        state = states[matter] if matter != vacuum_link_mediator else 0
        p0 = block_parity(state, 0)
        p1 = block_parity(state, 1)
        rung_parities = [
            -1.0 if ((state >> site_index(rung, 0)) & 1) else 1.0
            for rung in range(4)
        ]
        for gauge in range(8):
            gauss0 = p0 * x_eigenvalue(gauge, 0) * x_eigenvalue(gauge, 1)
            gauss1 = p1 * x_eigenvalue(gauge, 1) * x_eigenvalue(gauge, 2)
            for walker0 in range(4):
                for walker1 in range(4):
                    column = basis_index(matter, gauge, walker0, walker1, matter_dimension)
                    g0_values[column] = gauss0
                    g1_values[column] = gauss1
                    if matter == vacuum_link_mediator:
                        hamiltonian[column, column] += DETUNING
                    if walker0:
                        hamiltonian[column, column] += DETUNING
                    if walker1:
                        hamiltonian[column, column] += DETUNING
                    for first, second, amplitude in (
                        (0, 1, coupling * x_eigenvalue(gauge, 0)),
                        (1, 2, coupling * rung_parities[0]),
                        (2, 3, coupling * rung_parities[1]),
                        (3, 0, coupling * x_eigenvalue(gauge, 1)),
                    ):
                        if walker0 == first:
                            row = basis_index(matter, gauge, second, walker1, matter_dimension)
                            add_symmetric(hamiltonian, row, column, amplitude)
                    for first, second, amplitude in (
                        (0, 1, coupling * x_eigenvalue(gauge, 1)),
                        (1, 2, coupling * rung_parities[2]),
                        (2, 3, coupling * rung_parities[3]),
                        (3, 0, coupling * x_eigenvalue(gauge, 2)),
                    ):
                        if walker1 == first:
                            row = basis_index(matter, gauge, walker0, second, matter_dimension)
                            add_symmetric(hamiltonian, row, column, amplitude)

                    # Explicit link mediator: it can annihilate either local
                    # a-pair into its charged vacuum state.  The reverse term
                    # is added by Hermitian conjugation in add_symmetric.
                    if matter != vacuum_link_mediator:
                        for first, second in (
                            (site_index(0, 0), site_index(1, 0)),
                            (site_index(2, 0), site_index(3, 0)),
                        ):
                            item = _apply(state, (("ann", second), ("ann", first)))
                            if item is None or item[0] != 0:
                                continue
                            row = basis_index(vacuum_link_mediator, gauge, walker0, walker1, matter_dimension)
                            add_symmetric(hamiltonian, row, column, coupling * item[1])
                            add_symmetric(pair_link, row, column, coupling * item[1])

                    if matter != vacuum_link_mediator:
                        for operations in (
                            (("ann", site_index(2, 0)), ("create", site_index(1, 0))),
                            (("ann", site_index(1, 0)), ("create", site_index(2, 0))),
                        ):
                            item = _apply(state, operations)
                            if item is not None:
                                row = basis_index(positions[item[0]], gauge, walker0, walker1, matter_dimension)
                                single_boundary[row, column] += item[1]
    return hamiltonian.tocsr(), pair_link.tocsr(), g0_values, g1_values, {
        "dimension": dimension,
        "matter_dimension": matter_dimension,
        "source_matter": source_matter,
        "target_matter": target_matter,
        "link_mediator_matter": vacuum_link_mediator,
        "single_boundary": single_boundary.tocsr(),
    }


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DETUNING
        hamiltonian, pair_link, g0, g1, metadata = build_model(coupling)
        dimension = hamiltonian.shape[0]
        g0_operator = sparse.diags(g0, format="csr")
        g1_operator = sparse.diags(g1, format="csr")
        physical = (g0 == 1.0) & (g1 == 1.0)
        p_physical = sparse.diags(physical.astype(float), format="csr")
        p_nonphysical = sparse.eye(dimension, format="csr") - p_physical
        source = basis_index(metadata["source_matter"], 0, 0, 0, metadata["matter_dimension"])
        target_indices = np.asarray([
            basis_index(metadata["target_matter"], gauge, walker0, walker1, metadata["matter_dimension"])
            for gauge in range(8) for walker0 in range(4) for walker1 in range(4)
        ], dtype=int)
        link_mediator_indices = np.asarray([
            basis_index(metadata["link_mediator_matter"], gauge, walker0, walker1, metadata["matter_dimension"])
            for gauge in range(8) for walker0 in range(4) for walker1 in range(4)
        ], dtype=int)
        walker_excited = np.asarray([
            basis_index(matter, gauge, walker0, walker1, metadata["matter_dimension"])
            for matter in range(metadata["matter_dimension"])
            for gauge in range(8) for walker0 in range(4) for walker1 in range(4)
            if walker0 != 0 or walker1 != 0
        ], dtype=int)
        initial = np.zeros(dimension, dtype=complex)
        initial[source] = 1.0
        expected_time = np.pi * DETUNING / (2.0 * coupling ** 2)
        states_time = expm_multiply(
            -1j * hamiltonian, initial, start=0.0, stop=2.0 * expected_time,
            num=TIME_SAMPLES, endpoint=True,
        )
        target_population = np.sum(np.abs(states_time[:, target_indices]) ** 2, axis=1)
        peak = int(np.argmax(target_population))
        gauss_leakage = 1.0 - np.sum(np.abs(states_time[:, physical]) ** 2, axis=1)
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": coupling,
            "expected_sw_pair_transfer_time": float(expected_time),
            "maximum_target_pair_population": float(target_population[peak]),
            "time_of_maximum_target_pair_population": float(2.0 * expected_time * peak / (TIME_SAMPLES - 1)),
            "maximum_link_mediator_population": float(np.max(np.sum(np.abs(states_time[:, link_mediator_indices]) ** 2, axis=1))),
            "maximum_neutral_walker_population": float(np.max(np.sum(np.abs(states_time[:, walker_excited]) ** 2, axis=1))),
            "maximum_gauss_sector_leakage": float(np.max(np.abs(gauss_leakage))),
            "hamiltonian_g0_commutator_frobenius": frobenius(hamiltonian @ g0_operator - g0_operator @ hamiltonian),
            "hamiltonian_g1_commutator_frobenius": frobenius(hamiltonian @ g1_operator - g1_operator @ hamiltonian),
            "pair_link_g0_commutator_frobenius": frobenius(pair_link @ g0_operator - g0_operator @ pair_link),
            "pair_link_g1_commutator_frobenius": frobenius(pair_link @ g1_operator - g1_operator @ pair_link),
            "physical_pair_link_frobenius": frobenius(p_physical @ pair_link @ p_physical),
            "physical_single_boundary_projection_frobenius": frobenius(p_physical @ metadata["single_boundary"] @ p_physical),
            "single_boundary_leaves_physical_sector_frobenius": frobenius(p_nonphysical @ metadata["single_boundary"] @ p_physical),
        })
    output = {
        "schema": "antler.phase8b.two-block-derived-pair-link.v1",
        "parameters": {
            "blocks": 2,
            "block_size_b": 2,
            "fixed_total_charge": 2,
            "detuning": DETUNING,
            "ratios": list(RATIOS),
            "resources": "two neutral Lambda walkers plus one positive-detuned charge-two mediator shared by adjacent pair channels",
            "replaced_term": "the formerly inserted coherent P_L^dag P_R boundary transfer",
        },
        "dimensions": {"total": 3712, "matter_charge2_plus_link_vacuum": 29, "gauge_x_basis": 8, "walker_states": 16},
        "rows": rows,
        "decision": (
            "The explicit charge-two pair-link mediator can replace the inserted two-block pair-transfer term while preserving "
            "both coarse Gauss generators exactly in this fixed-charge Fock block."
        ),
        "claim_boundary": (
            "This establishes only a two-block, one-pair derived-link control. The neutral walker and four-leg pair-channel "
            "topology are declared new resources; no many-pair phase, thermodynamic gap, local indistinguishability, T "
            "junction, fusion or non-Abelian braid is established."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_two_block_derived_pair_link.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
