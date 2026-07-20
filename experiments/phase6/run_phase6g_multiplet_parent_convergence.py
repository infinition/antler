"""Phase 6G: state-to-Hamiltonian bridge for an external U(1) parent model.

No parameter search is performed.  The three mediator channels are fixed by
the exact local Gram factorization in Phase 6F, and Delta is only a controlled
Schrieffer--Wolff convergence parameter.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.multiplet_mediator_parent import (
    build_multiplet_mediator_parent,
    generalized_branch_parities,
    mediator_number,
)
from antler.number_conserving_pairwire import build_iemini_hamiltonian, wire_a_parity


L = N = 4
DELTAS = (10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0)


def sector_ground(H: np.ndarray, states: np.ndarray, sectors: tuple[int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.asarray([
        row for row, state in enumerate(states)
        if generalized_branch_parities(int(state), L) == sectors
    ], dtype=int)
    values, vectors = eigh(H[np.ix_(rows, rows)], subset_by_index=[0, 2], driver="evr")
    return rows, values, vectors


def target_frame(states: np.ndarray, index: dict[int, int]) -> np.ndarray:
    H, low_states, _ = build_iemini_hamiltonian(L, N, lam=1.0)
    frame = np.zeros((len(states), 2), dtype=complex)
    for column, parity in enumerate((0, 1)):
        rows = np.asarray([row for row, state in enumerate(low_states) if wire_a_parity(int(state), L) == parity], dtype=int)
        _, vectors = eigh(H[np.ix_(rows, rows)], subset_by_index=[0, 0], driver="evr")
        for source_row, amplitude in zip(rows, vectors[:, 0]):
            frame[index[int(low_states[source_row])], column] = amplitude
    return frame


def local_schur_residual(Delta: float) -> float:
    H, states, _ = build_multiplet_mediator_parent(L=2, total_charge=2, Delta=Delta)
    low = np.asarray([row for row, state in enumerate(states) if mediator_number(int(state), 2) == 0], dtype=int)
    high = np.asarray([row for row, state in enumerate(states) if mediator_number(int(state), 2) != 0], dtype=int)
    H_eff = H[np.ix_(low, low)] - H[np.ix_(low, high)] @ np.linalg.solve(H[np.ix_(high, high)], H[np.ix_(high, low)])
    H_target, _, _ = build_iemini_hamiltonian(L=2, N=2, lam=1.0)
    return float(np.linalg.norm(H_eff - H_target, ord="fro"))


def analyse(Delta: float) -> dict:
    H, states, index = build_multiplet_mediator_parent(L=L, total_charge=N, Delta=Delta)
    target = target_frame(states, index)
    records = []
    code = []
    for sector in ((1, 1), (-1, -1)):
        rows, values, vectors = sector_ground(H, states, sector)
        full = np.zeros(len(states), dtype=complex)
        full[rows] = vectors[:, 0]
        code.append(full)
        records.append({
            "sector": list(sector),
            "lowest_energies": values.tolist(),
            "mediator_number_expectation": float(sum(
                abs(full[row]) ** 2 * mediator_number(int(state), L) for row, state in enumerate(states)
            )),
            "zero_mediator_weight": float(sum(
                abs(full[row]) ** 2 for row, state in enumerate(states) if mediator_number(int(state), L) == 0
            )),
        })
    code_matrix = np.column_stack(code)
    overlap = target.conj().T @ code_matrix
    split = abs(records[1]["lowest_energies"][0] - records[0]["lowest_energies"][0])
    isolation = min(records[0]["lowest_energies"][1], records[1]["lowest_energies"][1]) - max(
        records[0]["lowest_energies"][0], records[1]["lowest_energies"][0]
    )
    maximum_cross_sector_element = 0.0
    labels = [generalized_branch_parities(int(state), L) for state in states]
    for row in range(len(states)):
        for column in range(len(states)):
            if labels[row] != labels[column] and abs(H[row, column]) > maximum_cross_sector_element:
                maximum_cross_sector_element = float(abs(H[row, column]))
    return {
        "Delta": Delta,
        "local_L2_zero_energy_schur_residual": local_schur_residual(Delta),
        "sectors": records,
        "logical_split": float(split),
        "isolation_gap": float(isolation),
        "minimum_target_frame_overlap": float(min(abs(overlap[0, 0]), abs(overlap[1, 1]))),
        "target_frame_overlap_abs": abs(overlap).tolist(),
        "maximum_cross_generalized_parity_matrix_element": maximum_cross_sector_element,
    }


def asymptotic_split_fit(rows: list[dict]) -> dict:
    """Fit the controlled large-detuning tail to ``a/Delta + b``."""
    tail = [row for row in rows if row["Delta"] >= 80.0]
    x = np.asarray([1.0 / row["Delta"] for row in tail], dtype=float)
    y = np.asarray([row["logical_split"] for row in tail], dtype=float)
    coefficient, intercept = np.polyfit(x, y, 1)
    fitted = coefficient * x + intercept
    residual = float(np.sum((y - fitted) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "fit_domain_Delta_at_least": 80.0,
        "coefficient_over_Delta": float(coefficient),
        "infinite_detuning_intercept": float(intercept),
        "r_squared": float(1.0 - residual / total) if total else 1.0,
    }


def main() -> None:
    rows = [analyse(Delta) for Delta in DELTAS]
    split_fit = asymptotic_split_fit(rows)
    qualifies = (
        max(row["local_L2_zero_energy_schur_residual"] for row in rows) < 1e-10
        and max(row["maximum_cross_generalized_parity_matrix_element"] for row in rows) < 1e-12
        and rows[-1]["minimum_target_frame_overlap"] > rows[0]["minimum_target_frame_overlap"]
        and rows[-1]["minimum_target_frame_overlap"] > 0.99
        and abs(split_fit["infinite_detuning_intercept"]) < 1e-3
        and split_fit["r_squared"] > 0.999
    )
    out = {
        "schema": "antler.phase6.multiplet-mediator-parent-convergence.v1",
        "construction": {
            "type": "state-to-Hamiltonian microscopic bridge",
            "external_low_energy_target": "published Iemini lambda=1 parent Hamiltonian",
            "fixed_design": "three charge-two mediator channels per bond from the rank-three Gram factorization; no coupling scan",
            "Delta_role": "controlled detuning for the Schrieffer--Wolff expansion only",
        },
        "symmetry": "generalized branch parities include the two mixed-pair mediator occupations",
        "rows": rows,
        "large_detuning_logical_split_fit": split_fit,
        "qualification": bool(qualifies),
        "claim_boundary": (
            "A successful convergence only gives a microscopic realization candidate for an external published parent Hamiltonian. "
            "It is not a distinct ANTLER topological discovery, a physical edge-braid implementation, or a topological quantum computer."
        ),
        "decision": (
            "the fixed mediator construction converges to the external parent frame; next audit quasi-conserved edge operators and physical implementation constraints"
            if qualifies else
            "the mediator construction fails to converge to the external parent frame; do not promote it as a microscopic bridge"
        ),
    }
    path = ROOT / "results" / "phase6" / "multiplet_mediator_parent_convergence.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "qualification": qualifies, "decision": out["decision"],
        "rows": [{key: row[key] for key in ("Delta", "minimum_target_frame_overlap", "logical_split", "isolation_gap")} for row in rows],
        "large_detuning_logical_split_fit": split_fit,
    }, indent=2))


if __name__ == "__main__":
    main()
