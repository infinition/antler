"""Small-ED reproduction control for the full number-conserving Floquet ladder.

This implements the high-frequency effective Hamiltonian H_eff=alpha H0+
(1-alpha) P^dag H0 P of Defossez et al. (arXiv:2412.14886), at eta=pi/2.
H0 contains intraleg hopping and attractive intraleg density interactions;
the rail rotation P generates the additional interleg density, swapping and
pair-hopping terms.  It is an external Floquet benchmark, not a derivation
from the frozen ANTLER hardware.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh, expm


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply, local_leg_raising_matrix, wire_a_parity


ETA = np.pi / 2.0
T_HOP = 1.0
L4_PARAMETERS = tuple((u0, alpha) for u0 in (-0.5, -1.0, -1.5, -2.0) for alpha in (0.25, 0.5, 0.75))


def build_h0_and_rotation(L: int, N: int, u0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[int, int]]:
    states, index = build_basis(2 * L, N)
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    j_x = np.zeros_like(hamiltonian)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for rung in range(L):
            a, b = site_index(rung, 0), site_index(rung, 1)
            for coefficient, operations in (
                (0.5, (("ann", b), ("create", a))),
                (0.5, (("ann", a), ("create", b))),
            ):
                item = _apply(state, operations)
                if item is not None:
                    new, amplitude = item
                    j_x[index[new], column] += coefficient * amplitude
        for rung in range(L - 1):
            for rail in (0, 1):
                left, right = site_index(rung, rail), site_index(rung + 1, rail)
                for operations in (
                    (("ann", right), ("create", left)),
                    (("ann", left), ("create", right)),
                ):
                    item = _apply(state, operations)
                    if item is not None:
                        new, amplitude = item
                        hamiltonian[index[new], column] += -T_HOP * amplitude
                hamiltonian[column, column] += u0 * (((state >> left) & 1) * ((state >> right) & 1))
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("H0 is not Hermitian")
    if not np.allclose(j_x, j_x.conj().T, atol=1e-12):
        raise RuntimeError("Jx is not Hermitian")
    return hamiltonian, expm(-1j * ETA * j_x), states, index


def analyze(L: int, N: int, u0: float, alpha: float) -> dict:
    h0, rotation, states, index = build_h0_and_rotation(L, N, u0)
    effective = alpha * h0 + (1.0 - alpha) * (rotation.conj().T @ h0 @ rotation)
    parity_residual = 0.0
    sectors = {}
    for parity in (0, 1):
        rows = np.asarray([position for position, state in enumerate(states) if wire_a_parity(int(state), L) == parity], dtype=int)
        values, vectors = eigh(effective[np.ix_(rows, rows)], subset_by_index=[0, 1], driver="evr")
        full = np.zeros(len(states), dtype=complex)
        full[rows] = vectors[:, 0]
        sectors[parity] = {"values": values, "vector": full, "rows": rows}
    even, odd = sectors[0], sectors[1]
    transfer = []
    for rung in range(L):
        raising = local_leg_raising_matrix(states, index, L, rung)
        transfer.append(float(abs(even["vector"].conj() @ raising @ odd["vector"])))
    edge = float(np.mean((transfer[0], transfer[-1])))
    bulk = float(max(transfer[1:-1])) if L > 2 else 0.0
    p_a = np.diag([(-1.0 if wire_a_parity(int(state), L) else 1.0) for state in states])
    p_b = (-1.0 if N & 1 else 1.0) * p_a
    parity_residual = max(
        float(np.linalg.norm(effective @ p_a - p_a @ effective)),
        float(np.linalg.norm(effective @ p_b - p_b @ effective)),
    )
    split = float(abs(even["values"][0] - odd["values"][0]))
    gap = float(min(even["values"][1] - even["values"][0], odd["values"][1] - odd["values"][0]))
    return {
        "L": L, "N": N, "filling": N / (2.0 * L),
        "parameters": {"u0_attractive_nn": u0, "alpha": alpha, "eta": ETA, "t_leg": T_HOP},
        "logical_split": split, "sector_isolation_gap": gap,
        "split_over_gap": float(split / gap) if gap > 1e-12 else None,
        "parity_commutator_frobenius_max": parity_residual,
        "local_transfer_by_rung": transfer,
        "edge_transfer_mean": edge,
        "bulk_transfer_max": bulk,
        "bulk_to_edge_transfer_ratio": float(bulk / edge) if edge > 1e-12 else None,
    }


def score(row: dict) -> tuple[float, float]:
    return (row["split_over_gap"] if row["split_over_gap"] is not None else np.inf,
            row["bulk_to_edge_transfer_ratio"] if row["bulk_to_edge_transfer_ratio"] is not None else np.inf)


def main() -> None:
    rows = [analyze(4, 2, u0, alpha) for u0, alpha in L4_PARAMETERS]
    promoted = sorted(rows, key=score)[:4]
    rows.extend(analyze(6, 4, row["parameters"]["u0_attractive_nn"], row["parameters"]["alpha"]) for row in promoted)
    l6 = [row for row in rows if row["L"] == 6]
    candidate_filter = [
        row for row in l6
        if row["parity_commutator_frobenius_max"] < 1e-10
        and row["split_over_gap"] is not None and row["split_over_gap"] < 1e-2
        and row["bulk_to_edge_transfer_ratio"] is not None and row["bulk_to_edge_transfer_ratio"] < 0.3
    ]
    out = {
        "schema": "antler.phase7d.full-floquet-ladder-preflight.v1",
        "citation": "Defossez et al., arXiv:2412.14886 (2025), high-frequency H_eff=alpha H0+(1-alpha)P^dag H0P at eta=pi/2",
        "model": {
            "H0": "intraleg hopping plus attractive intraleg nearest-neighbor density interaction",
            "P": "global rail rotation exp(-i eta Jx)",
            "effective_model_content": "includes the pair hopping, interleg density and swapping processes induced by the rotation",
        },
        "rows": rows,
        "L6_candidate_filter": {
            "split_over_gap_below": 1e-2, "bulk_to_edge_transfer_ratio_below": 0.3,
            "candidates": candidate_filter,
        },
        "decision": (
            "Small-ED reproduction control only. A passing L=6 filter would justify larger-size tensor-network or matrix-free work, "
            "not a topological or hardware claim; a failing filter closes only this sampled external Floquet benchmark window."
        ),
        "claim_boundary": (
            "This uses the published high-frequency effective model, not the finite-pulse charge-two mediator implementation. "
            "It does not establish an ANTLER derivation, an asymptotic phase, a protected local qubit, 2D order, braiding, "
            "non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "full_floquet_ladder_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"L6_rows": l6, "L6_candidates": candidate_filter}, indent=2))


if __name__ == "__main__":
    main()
