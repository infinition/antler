"""Two-block fixed-charge Gauss audit with an explicitly pair-only boundary.

The b=2 Lambda walker is placed on each of two blocks.  The only inter-block
matter term is a coherent transfer of two a-rail particles, which preserves
the a-parity of each block.  A one-particle boundary hop is retained solely as
a negative control.  The neutral-walker conditioned hop remains an inserted
new primitive.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply


DETUNING = 10.0
PARTICLE_NUMBER = 2
PAIR_TRANSFER = 0.05
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)


def basis_index(matter: int, gauge: int, walker0: int, walker1: int, matter_dimension: int) -> int:
    return ((((walker1 * 4) + walker0) * 8 + gauge) * matter_dimension) + matter


def x_eigenvalue(gauge: int, link: int) -> float:
    """Gauge labels are in the X eigenbasis, so X is diagonal."""
    return -1.0 if (gauge >> link) & 1 else 1.0


def block_parity(state: int, block: int) -> float:
    return -1.0 if sum((state >> site_index(2 * block + offset, 0)) & 1 for offset in (0, 1)) % 2 else 1.0


def build_model(coupling: float) -> tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix, np.ndarray, np.ndarray, dict]:
    states, positions = build_basis(8, PARTICLE_NUMBER)
    matter_dimension = len(states)
    dimension = matter_dimension * 8 * 4 * 4
    hamiltonian = sparse.lil_matrix((dimension, dimension), dtype=complex)
    pair_transfer = sparse.lil_matrix((dimension, dimension), dtype=complex)
    single_boundary = sparse.lil_matrix((dimension, dimension), dtype=complex)
    g0_values = np.empty(dimension, dtype=float)
    g1_values = np.empty(dimension, dtype=float)

    def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
        matrix[row, column] += value
        matrix[column, row] += np.conj(value)

    for matter, raw_state in enumerate(states):
        state = int(raw_state)
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
                    if walker0:
                        hamiltonian[column, column] += DETUNING
                    if walker1:
                        hamiltonian[column, column] += DETUNING
                    # Walker on block 0: links 0 and 1, rungs 0 and 1.
                    for first, second, amplitude in (
                        (0, 1, coupling * x_eigenvalue(gauge, 0)),
                        (1, 2, coupling * rung_parities[0]),
                        (2, 3, coupling * rung_parities[1]),
                        (3, 0, coupling * x_eigenvalue(gauge, 1)),
                    ):
                        if walker0 != first:
                            continue
                        row = basis_index(matter, gauge, second, walker1, matter_dimension)
                        add_symmetric(hamiltonian, row, column, amplitude)
                    # Walker on block 1: links 1 and 2, rungs 2 and 3.
                    for first, second, amplitude in (
                        (0, 1, coupling * x_eigenvalue(gauge, 1)),
                        (1, 2, coupling * rung_parities[2]),
                        (2, 3, coupling * rung_parities[3]),
                        (3, 0, coupling * x_eigenvalue(gauge, 2)),
                    ):
                        if walker1 != first:
                            continue
                        row = basis_index(matter, gauge, walker0, second, matter_dimension)
                        add_symmetric(hamiltonian, row, column, amplitude)
                    # Pair-only transfer across the two-block boundary:
                    # a_0^dag a_1^dag a_3 a_2 + h.c., hence Delta N_a per
                    # block is +/-2 and both Gauss parities are conserved.
                    for operations in (
                        (("ann", site_index(3, 0)), ("ann", site_index(2, 0)),
                         ("create", site_index(0, 0)), ("create", site_index(1, 0))),
                        (("ann", site_index(1, 0)), ("ann", site_index(0, 0)),
                         ("create", site_index(2, 0)), ("create", site_index(3, 0))),
                    ):
                        item = _apply(state, operations)
                        if item is not None:
                            new_state, amplitude = item
                            row = basis_index(positions[new_state], gauge, walker0, walker1, matter_dimension)
                            pair_transfer[row, column] += PAIR_TRANSFER * amplitude
                    # Negative control: a single a-rail hop across the same
                    # boundary flips both block parities.
                    for operations in (
                        (("ann", site_index(2, 0)), ("create", site_index(1, 0))),
                        (("ann", site_index(1, 0)), ("create", site_index(2, 0))),
                    ):
                        item = _apply(state, operations)
                        if item is not None:
                            new_state, amplitude = item
                            row = basis_index(positions[new_state], gauge, walker0, walker1, matter_dimension)
                            single_boundary[row, column] += amplitude
    hamiltonian = (hamiltonian.tocsr() + pair_transfer.tocsr())
    pair_transfer = pair_transfer.tocsr()
    single_boundary = single_boundary.tocsr()
    return hamiltonian, pair_transfer, single_boundary, g0_values, g1_values, {
        "dimension": dimension, "matter_dimension": matter_dimension,
    }


def frobenius(matrix: sparse.spmatrix) -> float:
    return float(np.sqrt(np.sum(np.abs(matrix.data) ** 2)))


def sector_ground_energy(hamiltonian: sparse.csr_matrix, g0: np.ndarray, g1: np.ndarray, sign0: int, sign1: int) -> float:
    indices = np.flatnonzero((g0 == sign0) & (g1 == sign1))
    block = hamiltonian[indices][:, indices]
    return float(sparse_linalg.eigsh(block, k=1, which="SA", return_eigenvectors=False, tol=1e-11)[0])


def main() -> None:
    rows = []
    for ratio in RATIOS:
        hamiltonian, pairs, single, g0, g1, metadata = build_model(ratio * DETUNING)
        g0_operator = sparse.diags(g0, format="csr")
        g1_operator = sparse.diags(g1, format="csr")
        p_physical = sparse.diags(((g0 == 1) & (g1 == 1)).astype(float), format="csr")
        p_nonphysical = sparse.eye(hamiltonian.shape[0], format="csr") - p_physical
        energies = {
            f"G0_{sign0:+d}_G1_{sign1:+d}": sector_ground_energy(hamiltonian, g0, g1, sign0, sign1)
            for sign0, sign1 in itertools.product((-1, 1), repeat=2)
        }
        physical_energy = energies["G0_+1_G1_+1"]
        excited = [value for key, value in energies.items() if key != "G0_+1_G1_+1"]
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": ratio * DETUNING,
            "sector_ground_energies": energies,
            "minimum_gauss_sector_gap_from_physical": float(min(excited) - physical_energy),
            "hamiltonian_g0_commutator_frobenius": frobenius(hamiltonian @ g0_operator - g0_operator @ hamiltonian),
            "hamiltonian_g1_commutator_frobenius": frobenius(hamiltonian @ g1_operator - g1_operator @ hamiltonian),
            "pair_transfer_g0_commutator_frobenius": frobenius(pairs @ g0_operator - g0_operator @ pairs),
            "pair_transfer_g1_commutator_frobenius": frobenius(pairs @ g1_operator - g1_operator @ pairs),
            "physical_pair_transfer_frobenius": frobenius(p_physical @ pairs @ p_physical),
            "single_boundary_g0_anticommutator_frobenius": frobenius(single @ g0_operator + g0_operator @ single),
            "single_boundary_g1_anticommutator_frobenius": frobenius(single @ g1_operator + g1_operator @ single),
            "physical_single_boundary_projection_frobenius": frobenius(p_physical @ single @ p_physical),
            "single_boundary_leaves_physical_sector_frobenius": frobenius(p_nonphysical @ single @ p_physical),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log([row["minimum_gauss_sector_gap_from_physical"] for row in deep]),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.two-block-pair-transport-audit.v1",
        "parameters": {
            "blocks": 2, "block_size_b": 2, "matter_modes": 8,
            "fixed_matter_particle_number": PARTICLE_NUMBER, "detuning": DETUNING,
            "pair_transfer": PAIR_TRANSFER, "ratios": list(RATIOS),
            "gauge_representation": "X eigenbasis on left, middle, right boundary qubits",
            "allowed_boundary_process": "coherent two-a-particle transfer only",
            "negative_control": "single a-particle boundary hop",
            "new_inserted_primitive": "neutral walker hop conditioned on (1-2 n_a,j)",
        },
        "dimensions": metadata,
        "rows": rows,
        "deep_sw_minimum_gauss_gap_power": power,
        "decision": (
            "The two-block inserted-primitive model has exact G0 and G1 conservation while its pair-only boundary transfer "
            "acts inside the physical sector. The negative-control single boundary hop is exactly sector-changing. This is "
            "a transport-algebra gate, not evidence for a paired phase or a native hardware implementation."
        ),
        "claim_boundary": (
            "The inter-block pair transfer is inserted as an effective coherent term. This audit does not derive it from the "
            "charge-two mediators in the gauged architecture, demonstrate phase stability, local indistinguishability, a "
            "thermodynamic gap, T-junction holonomies, fusion or non-Abelian braid statistics."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_two_block_pair_transport_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    import itertools
    main()
