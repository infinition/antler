"""Phase 6J: first-order Schrieffer--Wolff dressing of a finite edge operator.

The mediator bridge has virtual charge-two components.  An edge operator that
acts only in the zero-mediator sector is therefore not the physical candidate.
This audit tests the analytically fixed first-order dressing O + [S,O], where
S_high,low = -H_high,low/Delta.  It is a convergence test, not a parameter
search or a braid calculation.
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

from antler.basis import build_basis
from antler.multiplet_mediator_parent import (
    build_multiplet_mediator_parent,
    generalized_branch_parities,
    mediator_number,
)
from antler.number_conserving_pairwire import build_iemini_braid_z


L = N = 4
SUPPORT = 1
DELTAS = (80.0, 160.0, 320.0, 640.0)


def code_frame(H: np.ndarray, states: np.ndarray) -> np.ndarray:
    frame = []
    for sector in ((1, 1), (-1, -1)):
        rows = np.asarray([
            row for row, state in enumerate(states)
            if generalized_branch_parities(int(state), L) == sector
        ], dtype=int)
        _, vectors = eigh(H[np.ix_(rows, rows)], subset_by_index=[0, 0], driver="evr")
        full = np.zeros(len(states), dtype=complex)
        full[rows] = vectors[:, 0]
        frame.append(full)
    return np.column_stack(frame)


def embed_low_operator(states: np.ndarray, index: dict[int, int]) -> np.ndarray:
    low_states, low_index = build_basis(2 * L, N)
    Z, _ = build_iemini_braid_z(L, N, "aR_bR", SUPPORT, basis=(low_states, low_index))
    full = np.zeros((len(states), len(states)), dtype=complex)
    low_rows = [index[int(state)] for state in low_states]
    full[np.ix_(low_rows, low_rows)] = Z.toarray()
    return full


def schrieffer_wolff_generator(H: np.ndarray, states: np.ndarray, Delta: float) -> np.ndarray:
    """Leading anti-Hermitian block-off-diagonal generator for P=no mediators."""
    low = np.asarray([row for row, state in enumerate(states) if mediator_number(int(state), L) == 0], dtype=int)
    high = np.asarray([row for row, state in enumerate(states) if mediator_number(int(state), L) != 0], dtype=int)
    S = np.zeros_like(H)
    S[np.ix_(high, low)] = -H[np.ix_(high, low)] / Delta
    S[np.ix_(low, high)] = -S[np.ix_(high, low)].conj().T
    if not np.allclose(S, -S.conj().T, atol=1e-12):
        raise RuntimeError("Schrieffer--Wolff generator is not anti-Hermitian")
    return S


def metrics(H: np.ndarray, G: np.ndarray, operator: np.ndarray) -> dict:
    action = operator @ G
    logical = G.conj().T @ action
    leakage = action - G @ logical
    commutator_action = H @ action - operator @ (H @ G)
    return {
        "logical_action_frobenius": float(np.linalg.norm(logical)),
        "logical_leakage_amplitude_frobenius": float(np.linalg.norm(leakage)),
        "code_commutator_action_frobenius": float(np.linalg.norm(commutator_action)),
        "code_commutator_action_normalized": float(np.linalg.norm(commutator_action) / np.linalg.norm(action)),
        "antihermiticity_error": float(np.linalg.norm(operator + operator.conj().T)),
    }


def analyse(Delta: float) -> dict:
    H, states, index = build_multiplet_mediator_parent(L, N, Delta)
    G = code_frame(H, states)
    raw = embed_low_operator(states, index)
    S = schrieffer_wolff_generator(H, states, Delta)
    # Since S_high,low=-H_high,low/Delta, the low-state embedding is
    # |physical> = e^S |low> to first order. Hence O_phys=e^S O e^-S.
    dressed = raw + (S @ raw - raw @ S)
    wrong_sign = raw - (S @ raw - raw @ S)
    return {
        "Delta": Delta,
        "raw_zero_mediator_operator": metrics(H, G, raw),
        "first_order_SW_dressed_operator": metrics(H, G, dressed),
        "opposite_sign_convention_diagnostic": metrics(H, G, wrong_sign),
    }


def main() -> None:
    rows = [analyse(Delta) for Delta in DELTAS]
    improves = all(
        row["first_order_SW_dressed_operator"]["code_commutator_action_normalized"]
        < row["raw_zero_mediator_operator"]["code_commutator_action_normalized"]
        for row in rows
    )
    out = {
        "schema": "antler.phase6.sw-dressed-edge-operator-audit.v1",
        "construction": {
            "external_target": "Iemini finite-support aR-bR edge generator embedded in the multiplet-mediator bridge",
            "support_rungs": SUPPORT,
            "dressing": "O_phys = O + [S,O] + O(Delta^-1), S_high,low=-H_high,low/Delta",
        },
        "rows": rows,
        "first_order_improves_all_tested_Delta": bool(improves),
        "claim_boundary": (
            "This is only a first-order dressed finite-support edge diagnostic on a microscopic bridge to an external parent. "
            "It does not establish a quasi-zero mode, an adiabatic braid, or a native ANTLER topological qubit."
        ),
        "decision": (
            "first-order dressing removes the divergent virtual-mediator contribution at every tested detuning, but the residual stays finite; next test support and L scaling"
            if improves else
            "first-order dressing is insufficient or convention-sensitive; do not promote the bridge to an edge-protected braid"
        ),
    }
    path = ROOT / "results" / "phase6" / "sw_dressed_edge_operator_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
