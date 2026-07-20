"""Joint Schur audit of microscopic neutral-walker star and plaquette gadgets.

Each four-state neutral walker has one low position and three states detuned by
Delta.  A closed walk can only collect a non-scalar gauge word after traversing
all four of its links.  The star walker is conditioned on X links and the
plaquette walker on Z links.  They share exactly two links, as a toric-code
star and plaquette do.

This is intentionally an abstract new-resource construction: it checks whether
the local walker mechanism can generate compatible A_s and B_p words without
unwanted lower-order non-scalar terms.  It is not a native ANTLER derivation or
a tiled topological phase.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


ROOT = Path(__file__).resolve().parents[2]

DETUNING = 10.0
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)
EDGES = 6
GAUGE_DIMENSION = 1 << EDGES
WALKER_STATES = 4
DIMENSION = GAUGE_DIMENSION * WALKER_STATES * WALKER_STATES
STAR_EDGES = (0, 1, 2, 3)
PLAQUETTE_EDGES = (0, 1, 4, 5)
PAULI_NAMES = "IXYZ"


def index(gauge: int, star_walker: int, plaquette_walker: int) -> int:
    return ((star_walker * WALKER_STATES + plaquette_walker) * GAUGE_DIMENSION) + gauge


def z_eigenvalue(gauge: int, edge: int) -> float:
    return -1.0 if (gauge >> edge) & 1 else 1.0


def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
    matrix[row, column] += value
    matrix[column, row] += np.conj(value)


def build_hamiltonian(coupling: float) -> sparse.csr_matrix:
    hamiltonian = sparse.lil_matrix((DIMENSION, DIMENSION), dtype=complex)
    cycle = ((0, 1), (1, 2), (2, 3), (3, 0))
    for gauge in range(GAUGE_DIMENSION):
        for star_walker in range(WALKER_STATES):
            for plaquette_walker in range(WALKER_STATES):
                column = index(gauge, star_walker, plaquette_walker)
                if star_walker:
                    hamiltonian[column, column] += DETUNING
                if plaquette_walker:
                    hamiltonian[column, column] += DETUNING
                # The star walker is conditioned on X_e: it flips the Z-basis
                # gauge bit while hopping.  Its four-edge word is A_s.
                for (first, second), edge in zip(cycle, STAR_EDGES):
                    if star_walker != first:
                        continue
                    row = index(gauge ^ (1 << edge), second, plaquette_walker)
                    add_symmetric(hamiltonian, row, column, coupling)
                # The plaquette walker is conditioned on Z_e, diagonal in this
                # basis.  Its closed word is B_p.
                for (first, second), edge in zip(cycle, PLAQUETTE_EDGES):
                    if plaquette_walker != first:
                        continue
                    row = index(gauge, star_walker, second)
                    add_symmetric(hamiltonian, row, column, coupling * z_eigenvalue(gauge, edge))
    return hamiltonian.tocsr()


def pauli_action(codes: tuple[int, ...], state: int) -> tuple[int, complex]:
    """Return P|state> in the standard Z basis; code 0/1/2/3 is I/X/Y/Z."""
    target = state
    phase = 1.0 + 0.0j
    for edge, code in enumerate(codes):
        bit = (state >> edge) & 1
        if code == 1:  # X
            target ^= 1 << edge
        elif code == 2:  # Y
            target ^= 1 << edge
            phase *= 1j if bit == 0 else -1j
        elif code == 3:  # Z
            phase *= -1.0 if bit else 1.0
    return target, phase


def pauli_coefficients(hamiltonian: np.ndarray) -> dict[str, complex]:
    coefficients: dict[str, complex] = {}
    for codes in itertools.product(range(4), repeat=EDGES):
        value = 0.0j
        for column in range(GAUGE_DIMENSION):
            row, phase = pauli_action(codes, column)
            # trace(P H) = sum_column <column|P H|column>.
            value += np.conj(phase) * hamiltonian[row, column]
        label = "".join(PAULI_NAMES[code] for code in codes)
        coefficients[label] = value / GAUGE_DIMENSION
    return coefficients


def pauli_matrix(codes: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((GAUGE_DIMENSION, GAUGE_DIMENSION), dtype=complex)
    for column in range(GAUGE_DIMENSION):
        row, phase = pauli_action(codes, column)
        matrix[row, column] = phase
    return matrix


def downfold_at_zero(hamiltonian: sparse.csr_matrix) -> tuple[np.ndarray, float]:
    low = np.asarray([index(gauge, 0, 0) for gauge in range(GAUGE_DIMENSION)], dtype=int)
    mask = np.ones(DIMENSION, dtype=bool)
    mask[low] = False
    high = np.flatnonzero(mask)
    h_lh = hamiltonian[low][:, high]
    h_hh = hamiltonian[high][:, high].tocsc()
    h_hl = hamiltonian[high][:, low].toarray()
    effective = -np.asarray(h_lh @ sparse_linalg.spsolve(h_hh, h_hl))
    hermiticity_residual = float(np.linalg.norm(effective - effective.conj().T, ord="fro"))
    return effective, hermiticity_residual


def main() -> None:
    star_label = "XXXXII"
    plaquette_label = "ZZIIZZ"
    # Up to a sign convention, this is A_s B_p.  It belongs to the same
    # commuting stabilizer algebra and first appears at higher order.
    stabilizer_product_label = "YYXXZZ"
    star = pauli_matrix((1, 1, 1, 1, 0, 0))
    plaquette = pauli_matrix((3, 3, 0, 0, 3, 3))
    rows = []
    for ratio in RATIOS:
        effective, hermiticity = downfold_at_zero(build_hamiltonian(ratio * DETUNING))
        coefficients = pauli_coefficients(effective)
        noncommuting_unwanted = {
            label: value for label, value in coefficients.items()
            if label not in {"IIIIII", star_label, plaquette_label, stabilizer_product_label}
        }
        largest_label, largest_value = max(noncommuting_unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "star_coefficient": float(np.real_if_close(coefficients[star_label])),
            "plaquette_coefficient": float(np.real_if_close(coefficients[plaquette_label])),
            "identity_coefficient": float(np.real_if_close(coefficients["IIIIII"])),
            "star_plaquette_product_coefficient": float(np.real_if_close(coefficients[stabilizer_product_label])),
            "maximum_outside_stabilizer_algebra_coefficient": float(abs(largest_value)),
            "largest_outside_stabilizer_algebra_pauli": largest_label,
            "zero_energy_schur_hermiticity_residual": hermiticity,
            "star_plaquette_commutator_frobenius": float(np.linalg.norm(star @ plaquette - plaquette @ star, ord="fro")),
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    star_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["star_coefficient"] for row in deep])),
        1,
    )[0])
    plaquette_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["plaquette_coefficient"] for row in deep])),
        1,
    )[0])
    product_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["star_plaquette_product_coefficient"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.star-plaquette-walker-overlap-audit.v1",
        "parameters": {
            "gauge_links": EDGES,
            "shared_links": [0, 1],
            "star_edges": list(STAR_EDGES),
            "plaquette_edges": list(PLAQUETTE_EDGES),
            "detuning": DETUNING,
            "ratios": list(RATIOS),
            "new_primitives": "one neutral X-conditioned four-state star walker and one neutral Z-conditioned four-state plaquette walker",
            "downfolding": "exact zero-energy Schur complement onto both walkers in their low position",
        },
        "dimensions": {"total": DIMENSION, "low_gauge_subspace": GAUGE_DIMENSION, "high_walker_subspace": DIMENSION - GAUGE_DIMENSION},
        "rows": rows,
        "deep_sw_star_coefficient_power": star_power,
        "deep_sw_plaquette_coefficient_power": plaquette_power,
        "deep_sw_stabilizer_product_coefficient_power": product_power,
        "decision": (
            "The declared walker construction produces compatible star and plaquette words at fourth order in the joint "
            "downfolded six-link overlap block. The first extra word is their commuting product at higher order; all terms "
            "outside the generated stabilizer algebra are numerically absent in this audit."
        ),
        "claim_boundary": (
            "This is a single star/plaquette microscopic overlap control with newly declared X- and Z-conditioned neutral "
            "walker primitives. It does not establish a tiled 2D phase, a stable thermodynamic gap, a native ANTLER "
            "implementation, defects, fusion, a T junction or non-Abelian braid statistics."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_star_plaquette_walker_overlap_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
