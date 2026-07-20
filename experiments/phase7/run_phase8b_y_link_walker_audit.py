"""Test whether a phase-conditioned walker hop derives a local Pauli-Y link.

The Phase-8B walker vocabulary previously audited real X-conditioned and
Z-conditioned links.  A twist geometry also needs mixed checks with a local
Y action.  This script does not postulate a Y term: it uses the Hermitian
walker hopping matrix element

    <1|H|0> = i g Z,

which is exactly g Y on the addressed gauge qubit.  A four-state neutral
walker then traverses the explicit word Y X Z X.  The full finite Hamiltonian
is downfolded exactly and tomographed in the Pauli basis.

The result is a local compiler audit only.  The phase-conditioned neutral
walker remains a declared, non-native Phase-8B resource.
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
EDGES = 4
GAUGE_DIMENSION = 1 << EDGES
WALKER_STATES = 4
DIMENSION = GAUGE_DIMENSION * WALKER_STATES
PAULI_NAMES = "IXYZ"
# I=0, X=1, Y=2, Z=3.  The first link applies Z to the source
# amplitude then flips the gauge bit, i.e. i X Z = Y exactly.
PATTERN = (2, 1, 3, 1)
TARGET_LABEL = "YXZX"


def index(gauge: int, walker: int) -> int:
    return walker * GAUGE_DIMENSION + gauge


def z_eigenvalue(gauge: int, edge: int) -> float:
    return -1.0 if (gauge >> edge) & 1 else 1.0


def add_symmetric(matrix: sparse.lil_matrix, row: int, column: int, value: complex) -> None:
    matrix[row, column] += value
    matrix[column, row] += np.conj(value)


def link_transition(gauge: int, edge: int, pauli: int, coupling: float) -> tuple[int, complex]:
    """Return the gauge target and forward matrix element of one Pauli link."""
    if pauli == 1:  # X
        return gauge ^ (1 << edge), coupling
    if pauli == 2:  # Y = i Z X in the computational Z convention.
        return gauge ^ (1 << edge), 1j * coupling * z_eigenvalue(gauge, edge)
    if pauli == 3:  # Z
        return gauge, coupling * z_eigenvalue(gauge, edge)
    raise ValueError(f"unsupported Pauli code {pauli}")


def build_hamiltonian(coupling: float) -> sparse.csr_matrix:
    hamiltonian = sparse.lil_matrix((DIMENSION, DIMENSION), dtype=complex)
    for gauge in range(GAUGE_DIMENSION):
        for walker in range(WALKER_STATES):
            column = index(gauge, walker)
            if walker:
                hamiltonian[column, column] += DETUNING
            next_walker = (walker + 1) % WALKER_STATES
            next_gauge, amplitude = link_transition(gauge, walker, PATTERN[walker], coupling)
            add_symmetric(hamiltonian, index(next_gauge, next_walker), column, amplitude)
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
    low = np.asarray([index(gauge, 0) for gauge in range(GAUGE_DIMENSION)], dtype=int)
    mask = np.ones(DIMENSION, dtype=bool)
    mask[low] = False
    high = np.flatnonzero(mask)
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
        unwanted = {label: value for label, value in coefficients.items() if label not in {"IIII", TARGET_LABEL}}
        largest_label, largest_value = max(unwanted.items(), key=lambda item: abs(item[1]))
        rows.append({
            "coupling_over_detuning": ratio,
            "target_yxzx_coefficient": float(np.real_if_close(coefficients[TARGET_LABEL])),
            "maximum_unwanted_non_scalar_coefficient": float(abs(largest_value)),
            "largest_unwanted_non_scalar_pauli": largest_label,
            "zero_energy_schur_hermiticity_residual": hermiticity,
        })
    deep = [row for row in rows if row["coupling_over_detuning"] <= 0.075]
    power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log(np.abs([row["target_yxzx_coefficient"] for row in deep])),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.y-link-walker-audit.v1",
        "parameters": {
            "gauge_links": EDGES,
            "walker_states": WALKER_STATES,
            "detuning": DETUNING,
            "ratios": list(RATIOS),
            "target_word": TARGET_LABEL,
            "y_link_realization": "forward hop i*g*Z followed by gauge-bit flip; this is exactly Pauli Y",
            "downfolding": "exact zero-energy Schur complement onto walker vacuum",
        },
        "dimensions": {
            "total": DIMENSION,
            "low_gauge_subspace": GAUGE_DIMENSION,
            "high_walker_subspace": DIMENSION - GAUGE_DIMENSION,
        },
        "rows": rows,
        "deep_sw_yxzx_power": power,
        "decision": (
            "A coherent pi/2 phase on an X-conditioned walker hop implements its local Pauli-Y action exactly. "
            "The explicit four-link loop isolates the mixed YXZX word at fourth order."
        ),
        "claim_boundary": (
            "This is a local phase-controlled-walker derivation, not a full twist geometry, a branch cut, a fusion space, "
            "defect motion, a non-Abelian braid, a native ANTLER realization or a fault-tolerance result. It requires that "
            "the declared neutral walker can maintain and calibrate a coherent pi/2 conditional phase."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_y_link_walker_audit.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
