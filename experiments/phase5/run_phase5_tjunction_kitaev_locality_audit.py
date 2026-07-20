"""Phase 5L: local-indistinguishability audit of the four-defect BdG target.

The four zero Majoranas are converted to a fixed-even-parity two-dimensional
code.  For every physical site the script projects the local fermion-parity
bilinear ``i a_j b_j`` onto that code and reports its traceless norm.  A
nonzero value would mean that a strictly local, parity-even probe can directly
distinguish logical states at zero-mode order.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from run_phase5_tjunction_kitaev_preflight import bdg_tjunction


I2 = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Y = np.array([[0.0, -1j], [1j, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def zero_majorana_basis(bdg: np.ndarray, tolerance: float) -> np.ndarray:
    """Return real orthonormal Majorana null vectors in (a-sites, b-sites) order."""

    n2 = bdg.shape[0]
    n = n2 // 2
    transform = np.block([
        [np.eye(n), np.eye(n)],
        [-1j * np.eye(n), 1j * np.eye(n)],
    ])
    # In this convention R H_BdG R^dagger = 2 i A for H=(i/4) gamma A gamma.
    majorana_generator = (-0.5j * transform @ bdg @ transform.conj().T).real
    majorana_generator = 0.5 * (majorana_generator - majorana_generator.T)
    _, singular, vh = np.linalg.svd(majorana_generator)
    null = vh[singular < tolerance].T
    if null.shape[1] != 4:
        raise RuntimeError(f"expected exactly four Majorana null vectors, found {null.shape[1]}")
    return null


def zero_gamma_matrices() -> list[np.ndarray]:
    return [np.kron(X, I2), np.kron(Y, I2), np.kron(Z, X), np.kron(Z, Y)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm-length", type=int, default=6)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--junction-chiral", type=float, default=1.0)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/tjunction_kitaev_locality_audit.json"))
    args = parser.parse_args()
    bdg, endpoints = bdg_tjunction(args.arm_length, 0.0, args.t, args.delta,
                                    args.junction_chiral, "bb")
    zero = zero_majorana_basis(bdg, args.tolerance)
    gamma = zero_gamma_matrices()
    code = np.zeros((4, 2), complex)
    code[0, 0] = 1.0
    code[3, 1] = 1.0
    rows = []
    n = bdg.shape[0] // 2
    for site in range(n):
        a = sum(zero[site, alpha] * gamma[alpha] for alpha in range(4))
        b = sum(zero[n + site, alpha] * gamma[alpha] for alpha in range(4))
        # The antisymmetrisation is exact for full Majoranas and remains the
        # well-defined projected parity bilinear after dropping gapped modes.
        local_parity = 0.5j * (a @ b - b @ a)
        projected = code.conj().T @ local_parity @ code
        traceless = projected - np.trace(projected) * np.eye(2) / 2.0
        rows.append({
            "site": site,
            "zero_weight_a": float(np.sum(zero[site] ** 2)),
            "zero_weight_b": float(np.sum(zero[n + site] ** 2)),
            "logical_projection_norm": float(np.linalg.norm(projected)),
            "logical_traceless_norm": float(np.linalg.norm(traceless)),
        })
    output = {
        "schema": "antler.phase5.tjunction-kitaev-locality-audit.v1",
        "claim_boundary": (
            "This checks zero-mode-order local parity projections in the phase-biased "
            "BdG target.  It does not establish locality against arbitrary finite-energy "
            "perturbations or derive the pairing target from frozen ANTLER."
        ),
        "arm_length": args.arm_length, "endpoints": endpoints,
        "zero_mode_count": 4,
        "max_local_traceless_norm": max(row["logical_traceless_norm"] for row in rows),
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({
        "max_local_traceless_norm": output["max_local_traceless_norm"],
        "largest_rows": sorted(rows, key=lambda row: row["logical_traceless_norm"], reverse=True)[:4],
    }, indent=2))


if __name__ == "__main__":
    main()
