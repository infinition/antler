"""Four-rung composition audit for the signed-mediator Floquet pair gate.

Each nearest-neighbour pair gate is the stroboscopic output of the Phase 7D
positive-detuning signed-ZZ construction.  Even links are applied in parallel,
then the odd link.  This script measures the digital Trotter error against the
target sum of number-conserving pair-hopping generators on a four-rung ladder.
It deliberately audits a compiler, not a topological phase.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


RUNG_COUNT = 4
ANGLES = (0.008, 0.012, 0.018, 0.027, 0.040)
PAULIS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "Y": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    "Z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


@lru_cache(maxsize=None)
def pauli_word(label: str) -> np.ndarray:
    output = np.asarray([[1.0]], dtype=complex)
    for letter in label:
        output = np.kron(output, PAULIS[letter])
    return output


def pair_generator(left: int, right: int) -> np.ndarray:
    """S+_left S+_right + h.c. = (XX-YY)/2 on one logical link."""
    xx = ["I"] * RUNG_COUNT
    yy = ["I"] * RUNG_COUNT
    xx[left] = xx[right] = "X"
    yy[left] = yy[right] = "Y"
    return 0.5 * (pauli_word("".join(xx)) - pauli_word("".join(yy)))


def parity_operator(rail: int) -> np.ndarray:
    diagonal = []
    for rails in product((0, 1), repeat=RUNG_COUNT):
        count = sum(value == rail for value in rails)
        diagonal.append(-1.0 if count & 1 else 1.0)
    return np.diag(diagonal).astype(complex)


def fit_power(rows: list[dict]) -> dict:
    x = np.log(np.asarray([row["pair_angle"] for row in rows]))
    y = np.log(np.asarray([row["trotter_residual"] for row in rows]))
    power, intercept = np.polyfit(x, y, 1)
    fitted = power * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "power": float(power), "prefactor": float(np.exp(intercept)),
        "r_squared_log": float(1.0 - np.sum((y - fitted) ** 2) / total),
    }


def main() -> None:
    q01, q12, q23 = (pair_generator(*link) for link in ((0, 1), (1, 2), (2, 3)))
    target_generator = q01 + q12 + q23
    parity_a, parity_b = parity_operator(0), parity_operator(1)
    rows = []
    for angle in ANGLES:
        even_layer = expm(-1j * angle * q23) @ expm(-1j * angle * q01)
        odd_layer = expm(-1j * angle * q12)
        compiled = odd_layer @ even_layer
        target = expm(-1j * angle * target_generator)
        rows.append({
            "pair_angle": angle,
            "trotter_residual": float(np.linalg.norm(compiled - target, ord=2)),
            "parity_a_residual": float(np.linalg.norm(compiled @ parity_a - parity_a @ compiled, ord=2)),
            "parity_b_residual": float(np.linalg.norm(compiled @ parity_b - parity_b @ compiled, ord=2)),
            "unitarity_residual": float(np.linalg.norm(compiled.conj().T @ compiled - np.eye(1 << RUNG_COUNT), ord=2)),
        })
    out = {
        "schema": "antler.phase7d.multilink-floquet-compiler-audit.v1",
        "model": {
            "rung_count": RUNG_COUNT,
            "logical_dimension": 1 << RUNG_COUNT,
            "links": [[0, 1], [1, 2], [2, 3]],
            "schedule": "parallel even links (0,1),(2,3), followed by odd link (1,2)",
            "target_generator": "sum_j (S+_j S+_(j+1) + h.c.)",
            "microscopic_contract": "each pair gate is the closed positive-detuning mediator pulse from Phase 7D signed-ZZ preflight",
        },
        "noncommuting_layer_commutator_frobenius": float(np.linalg.norm((q01 + q23) @ q12 - q12 @ (q01 + q23))),
        "rows": rows,
        "trotter_power_fit": fit_power(rows),
        "decision": (
            "The dynamic pair-hopping primitive composes into a parity-preserving four-rung digital compiler with a measured "
            "Trotter error. Passing this gate permits a full finite-pulse ladder Hamiltonian audit; it does not by itself create "
            "a gapped topological phase or a commuting 2D parent."
        ),
        "claim_boundary": (
            "This audit acts entirely within the monomer logical sector after each closed mediator pulse. It does not yet include "
            "simultaneous mediator crosstalk, finite rail-rotation errors, leg hopping, a many-body gap, edge protection, 2D order, "
            "braiding, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "multilink_floquet_compiler_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
