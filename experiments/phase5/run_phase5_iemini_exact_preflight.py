"""Phase 5I: exact-Iemini-model reproduction and locality preflight.

This is a reproduction check for the published number-conserving parent
Hamiltonian, never a claim that frozen ANTLER realizes it.  At lambda=1 it
tests the exact zero-energy doublet and its local transfer profile in the same
finite-size audit used to reject the preceding minimal approximation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import issparse
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.number_conserving_pairwire import (
    build_iemini_hamiltonian,
    local_leg_raising_matrix,
    wire_a_parity,
)


def sector_ground(H: np.ndarray, states: np.ndarray, L: int, parity: int):
    rows = np.array([i for i, state in enumerate(states) if wire_a_parity(int(state), L) == parity])
    block = H[rows, :][:, rows] if issparse(H) else H[np.ix_(rows, rows)]
    if issparse(block):
        values, vectors = eigsh(block, k=3, which="SA", tol=1e-11)
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
    else:
        values, vectors = eigh(block, subset_by_index=[0, 2], driver="evr")
    full = np.zeros(len(states), complex)
    full[rows] = vectors[:, 0]
    return values, full, rows


def analyse(L: int, N: int, sparse: bool = False) -> dict:
    H, states, index = build_iemini_hamiltonian(L, N, lam=1.0, sparse=sparse)
    ee_values, ee, ee_rows = sector_ground(H, states, L, 0)
    oo_values, oo, oo_rows = sector_ground(H, states, L, 1)
    profile = []
    for j in range(L):
        O = local_leg_raising_matrix(states, index, L, j)
        profile.append(float(abs(oo.conj() @ O @ ee)))
    edge = float(np.mean([profile[0], profile[-1]]))
    bulk = float(max(profile[1:-1])) if L > 2 else 0.0
    return {
        "L": L,
        "N": N,
        "filling": N / (2.0 * L),
        "lambda": 1.0,
        "sector_energies_lowest_three": {"ee": ee_values.tolist(), "oo": oo_values.tolist()},
        "ground_doublet_split": float(abs(ee_values[0] - oo_values[0])),
        "fixed_parity_collective_gap": float(min(ee_values[1] - ee_values[0], oo_values[1] - oo_values[0])),
        "ground_energy_abs_max": float(max(abs(ee_values[0]), abs(oo_values[0]))),
        "cross_parity_block_norm": float(
            np.linalg.norm(H[ee_rows, :][:, oo_rows].data)
            if issparse(H) else np.linalg.norm(H[np.ix_(ee_rows, oo_rows)])
        ),
        "solver": "sparse-eigsh" if sparse else "dense-eigh",
        "local_a_dag_b_profile": profile,
        "edge_transfer_mean": edge,
        "bulk_transfer_max": bulk,
        "bulk_to_edge_transfer_ratio": float(bulk / edge) if edge > 1e-13 else None,
    }


def main() -> None:
    rows = [analyse(4, 4), analyse(6, 6), analyse(8, 8, sparse=True)]
    reproduction_pass = all(
        row["ground_energy_abs_max"] < 1e-10
        and row["ground_doublet_split"] < 1e-10
        and row["cross_parity_block_norm"] < 1e-12
        for row in rows
    )
    localized_trend = (
        rows[-1]["bulk_to_edge_transfer_ratio"] is not None
        and rows[-1]["bulk_to_edge_transfer_ratio"] < rows[0]["bulk_to_edge_transfer_ratio"]
    )
    out = {
        "schema": "antler.phase5.iemini-exact-preflight.v1",
        "reference": {
            "citation": "F. Iemini et al., Phys. Rev. Lett. 115, 156402 (2015), Eq. (3)",
            "arxiv": "https://arxiv.org/abs/1504.04230",
            "boundary_conditions": "open",
            "exact_line": "lambda=1",
        },
        "claim_boundary": (
            "This independently reproduces a published number-conserving reference Hamiltonian. "
            "It does not derive that Hamiltonian from ANTLER, does not establish a gapped many-body topological phase, and does not execute a physical braid."
        ),
        "rows": rows,
        "reproduction_checks_pass": reproduction_pass,
        "finite_size_localization_trend": localized_trend,
        "decision": (
            "exact-parent reproduction passed; the proper next task is a finite-size scaling and explicit published edge-braid-operator audit, still outside frozen ANTLER"
            if reproduction_pass else
            "reproduction failed; do not use this implementation for a localization or braid audit"
        ),
    }
    path = Path("results/phase5/iemini_exact_preflight.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
