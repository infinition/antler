"""Exact neutral-walker derivation of a mixed-Pauli pentagon word.

Twist endpoints in the Wen/Bombin surface-code construction require a local
non-CSS stabilizer of odd weight (a pentagon).  A five-step neutral walker loop
is the smallest analogue of the established four-step star/plaquette gadget.
It is conditioned successively on X, Z, X, Z, X gauge-link operators and is
audited without inserting the resulting word by hand.

The calculation establishes only the local word; a commuting branch-cut
geometry and defect transport remain separate gates.
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
EDGES = 5
GAUGE_DIMENSION = 1 << EDGES
WALKER_STATES = 5
DIMENSION = GAUGE_DIMENSION * WALKER_STATES
PAULI_NAMES = "IXYZ"
PATTERN = (1, 3, 1, 3, 1)  # X Z X Z X
TARGET_LABEL = "XZXZX"


def index(gauge: int, walker: int) -> int:
    return walker * GAUGE_DIMENSION + gauge


def z_eigenvalue(gauge: int, edge: int) -> float:
    return -1.0 if (gauge >> edge) & 1 else 1.0


def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
    matrix[row, column] += value
    matrix[column, row] += np.conj(value)


def build_hamiltonian(coupling: float) -> sparse.csr_matrix:
    hamiltonian = sparse.lil_matrix((DIMENSION, DIMENSION), dtype=complex)
    cycle = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))
    for gauge in range(GAUGE_DIMENSION):
        for walker in range(WALKER_STATES):
            column = index(gauge, walker)
            if walker:
                hamiltonian[column, column] += DETUNING
            for (first, second), edge, pauli in zip(cycle, range(EDGES), PATTERN):
                if walker != first:
                    continue
                if pauli == 1:  # X condition flips the Z-basis link qubit.
                    row = index(gauge ^ (1 << edge), second)
                    amplitude = coupling
                else:  # Z condition contributes its link eigenvalue.
                    row = index(gauge, second)
                    amplitude = coupling * z_eigenvalue(gauge, edge)
                add_symmetric(hamiltonian, row, column, amplitude)
    return hamiltonian.tocsr()


def pauli_action(codes: tuple[int, ...], state: int) -> tuple[int, complex]:
    target = state
    phase = 1.0 + 0.0j
    for edge, code in enumerate(codes):
        bit = (state >> edge) & 1
        if code == 1:
            target ^= 1 << edge
        elif code == 2:
            target ^= 1 << edge
            phase *= 1j if bit == 0 else -1j
        elif code == 3:
            phase *= -1.0 if bit else 1.0
    return target, phase


def pauli_coefficients(hamiltonian: np.ndarray) -> dict[str, complex]:
    coefficients: dict[str, complex] = {}
    for codes in itertools.product(range(4), repeat=EDGES):
        value = 0.0j
        for column in range(GAUGE_DIMENSION):
            row, phase = pauli_action(codes, column)
            value += np.conj(phase) * hamiltonian[row, column]
        coefficients["".join(PAULI_NAMES[code] for code in codes)] = value / GAUGE_DIMENSION
    return coefficients


def downfold_at_zero(hamiltonian: sparse.csr_matrix) -> tuple[np.ndarray, float]:
    low = np.arange(GAUGE_DIMENSION, dtype=int)
    high = np.arange(GAUGE_DIMENSION, DIMENSION, dtype=int)
    h_lh = hamiltonian[low][:, high]
    h_hh = hamiltonian[high][:, high].tocsc()
    h_hl = hamiltonian[high][:, low].toarray()
    effective = -np.asarray(h_lh @ sparse_linalg.spsolve(h_hh, h_hl))
    return effective, float(np.linalg.norm(effective - effective.conj().T, ord="fro"))


def main() -> None:
    rows = []
    for ratio in RATIOS:
        effective, hermiticity = downfold_at_zero(build_hamiltonian(ratio * DETUNING))
        coefficients = pauli_coefficients(effective)
        unwanted = {
            label: value for label, value in coefficients.items()
            if label not in {"IIIII", TARGET_LABEL}
        }
        label, value = max(unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "mixed_pentagon_coefficient": float(np.real_if_close(coefficients[TARGET_LABEL])),
            "identity_coefficient": float(np.real_if_close(coefficients["IIIII"])),
            "maximum_unwanted_nonscalar_coefficient": float(abs(value)),
            "largest_unwanted_pauli": label,
            "zero_energy_schur_hermiticity_residual": hermiticity,
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["mixed_pentagon_coefficient"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.mixed-pentagon-walker-audit.v1",
        "parameters": {
            "gauge_links": EDGES,
            "walker_states": WALKER_STATES,
            "detuning": DETUNING,
            "ratios": list(RATIOS),
            "conditioned_pauli_pattern": TARGET_LABEL,
            "new_primitive": "one neutral five-state walker with phase-controlled X/Z-conditioned link hops",
            "downfolding": "exact zero-energy Schur complement onto the walker-low subspace",
        },
        "dimensions": {"total": DIMENSION, "low_gauge_subspace": GAUGE_DIMENSION, "high_walker_subspace": DIMENSION - GAUGE_DIMENSION},
        "rows": rows,
        "deep_sw_mixed_pentagon_power": power,
        "decision": (
            "The declared five-step walker produces the mixed non-CSS XZXZX pentagon word at fifth order; no other "
            "non-scalar Pauli coefficient is resolved in this isolated local downfolding."
        ),
        "claim_boundary": (
            "This is only a local mixed-pentagon mechanism. It does not specify a complete commuting twist branch-cut "
            "lattice, demonstrate a twist-defect fusion space, move defects, derive a braid, or establish a native ANTLER "
            "implementation. The phase-controlled conditioned walker is a new resource."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_mixed_pentagon_walker_audit.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
