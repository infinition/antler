"""Phase 5H: number-conserving pair-hopping doublet preflight.

The goal is deliberately narrower than a braid calculation.  For the minimal
two-wire, number-conserving pair-hopping family, it asks whether the two
wire-parity sectors provide a spectrally isolated and *locally
indistinguishable* doublet.  Only a positive finite-size trend would justify
the subsequent construction of defects, fusion spaces, and transport.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.number_conserving_pairwire import (
    build_pairwire_hamiltonian,
    local_leg_raising_matrix,
    local_transfer_matrix,
    wire_a_parity,
)


def sector_data(H: np.ndarray, states: np.ndarray, L: int, parity: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.array([i for i, state in enumerate(states) if wire_a_parity(int(state), L) == parity])
    # The two lowest eigenpairs give both the doublet energy and its sector gap.
    values, vectors = eigh(H[np.ix_(rows, rows)], subset_by_index=[0, 1], driver="evr")
    full = np.zeros(len(states), dtype=complex)
    full[rows] = vectors[:, 0]
    return values, full, rows


def occupancy(vector: np.ndarray, states: np.ndarray, site: int) -> float:
    probabilities = np.abs(vector) ** 2
    return float(sum(probabilities[i] * ((int(state) >> site) & 1) for i, state in enumerate(states)))


def analyse(L: int, w: float, v: float) -> dict:
    N, t = L, 1.0
    H, states, index = build_pairwire_hamiltonian(L, N, t=t, w=w, v=v)
    even_values, even, even_rows = sector_data(H, states, L, 0)
    odd_values, odd, odd_rows = sector_data(H, states, L, 1)
    split = float(abs(even_values[0] - odd_values[0]))
    isolation_gap = float(min(even_values[1] - even_values[0], odd_values[1] - odd_values[0]))
    transfer = []
    symmetric_transfer = []
    density_differences = []
    for j in range(L):
        raising = local_leg_raising_matrix(states, index, L, j)
        O = local_transfer_matrix(states, index, L, j)
        # ``a^dag b`` is the local complex operator.  Its magnitude is the
        # maximum response over the two independent Hermitian quadratures;
        # recording the symmetric quadrature separately exposes cancellations
        # caused by wire-exchange symmetry.
        transfer.append(float(abs(even.conj() @ raising @ odd)))
        symmetric_transfer.append(float(abs(even.conj() @ O @ odd)))
        for leg in (0, 1):
            site = 2 * j + leg
            density_differences.append(abs(occupancy(even, states, site) - occupancy(odd, states, site)))
    edge = float(np.mean([transfer[0], transfer[-1]]))
    bulk = float(max(transfer[1:-1])) if L > 2 else 0.0
    ratio = float(bulk / edge) if edge > 1e-13 else None
    cross_block_norm = float(np.linalg.norm(H[np.ix_(even_rows, odd_rows)]))
    return {
        "L": L,
        "N": N,
        "parameters": {"t": t, "w_pair": w, "v_nn": v},
        "sector_ground_energies": {"even_even": float(even_values[0]), "odd_odd": float(odd_values[0])},
        "logical_split": split,
        "sector_isolation_gap": isolation_gap,
        "split_over_gap": float(split / isolation_gap) if isolation_gap > 0 else None,
        "wire_parity_cross_block_norm": cross_block_norm,
        "local_transfer_by_rung": transfer,
        "symmetric_transfer_quadrature_by_rung": symmetric_transfer,
        "edge_transfer_mean": edge,
        "bulk_transfer_max": bulk,
        "bulk_to_edge_transfer_ratio": ratio,
        "max_density_difference": float(max(density_differences)),
        "max_bulk_density_difference": float(max(density_differences[2:-2])) if L > 3 else None,
    }


def score(row: dict) -> tuple[float, float, float]:
    """Conservative ordering only; it is not a protection certificate."""
    ratio = row["bulk_to_edge_transfer_ratio"]
    if ratio is None:
        ratio = 1e9
    return (row["split_over_gap"], ratio, row["max_density_difference"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="only scan L=4; useful as a smoke test")
    args = parser.parse_args()

    # A deliberately modest, reproducible scan.  It is a preflight, not an
    # optimization claim: no hidden fine tuning is accepted as evidence.
    parameters = [(w, v) for w in (0.5, 1.0, 1.5, 2.0) for v in (-1.0, 0.0, 1.0)]
    rows = [analyse(4, w, v) for w, v in parameters]
    if not args.quick:
        # Carry the four most promising *small-system* points to L=6.  The
        # ranking is just a resource allocation rule; all raw L=4 rows remain.
        promoted = sorted(rows, key=score)[:4]
        rows.extend(analyse(6, row["parameters"]["w_pair"], row["parameters"]["v_nn"]) for row in promoted)

    l6 = [row for row in rows if row["L"] == 6]
    conservative_candidates = [
        row for row in l6
        if row["split_over_gap"] is not None and row["split_over_gap"] < 1e-2
        and row["bulk_to_edge_transfer_ratio"] is not None and row["bulk_to_edge_transfer_ratio"] < 0.1
        and row["max_density_difference"] < 0.1
    ]
    decision = (
        "no finite-size protected-doublet candidate in this minimal scan; do not construct braids from it"
        if not conservative_candidates else
        "finite-size candidate(s) found only; require L-scaling, disorder, and defect/fusion construction before any braid claim"
    )
    out = {
        "schema": "antler.phase5.number-conserving-pairwire-preflight.v1",
        "model": {
            "description": "two ordinary fermion wires with leg hopping, pair transfer, and intra-wire NN density interaction",
            "conserved_quantities": ["total particle number", "wire-a parity", "wire-b parity"],
            "fock_ordering": "rung-major (j,a),(j,b)",
        },
        "claim_boundary": (
            "This is an additive, Iemini-inspired number-conserving reference family, not the frozen ANTLER correlated-hopping Hamiltonian and not an implementation of the exact solvable model. "
            "It does not establish Majorana modes, a fusion space, non-Abelian exchange, or a topological qubit."
        ),
        "criteria": {
            "finite_size_candidate_only": {
                "split_over_gap_below": 1e-2,
                "bulk_to_edge_transfer_ratio_below": 0.1,
                "max_density_difference_below": 0.1,
            },
            "warning": "Passing these small-L filters would be a reason to scale the test, not a protection proof.",
        },
        "rows": rows,
        "l6_conservative_candidates": conservative_candidates,
        "decision": decision,
        "next_step": (
            "If and only if a stable L-scaling emerges, introduce explicit defects and audit local indistinguishability of their fusion doublet before gate or braid dynamics."
        ),
    }
    path = Path("results/phase5/pair_hopping_two_wire_preflight.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": decision, "l6_rows": l6, "candidates": conservative_candidates}, indent=2))


if __name__ == "__main__":
    main()
