"""Phase 6L: positive control for the projected edge-commutator harness.

The Phase 6H--K diagnostic is useful only if it can distinguish an exact
zero mode from a nearby, non-edge operator.  This script applies the identical
code-projected commutator definition to an open Kitaev chain at its exact
``mu=0, t=Delta`` point.  The left Majorana is an exact zero mode while the
next-site Majorana is a deliberately failing bulk control.

This is a calibration of the *diagnostic*, not a mapping from ANTLER or the
number-conserving Iemini parent to a Kitaev wire.
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


def annihilator(L: int, site: int) -> np.ndarray:
    """Dense Fock-space annihilator in ordinary increasing-site ordering."""
    dimension = 1 << L
    out = np.zeros((dimension, dimension), dtype=complex)
    lower_mask = (1 << site) - 1
    for state in range(dimension):
        if (state >> site) & 1:
            target = state ^ (1 << site)
            sign = -1.0 if ((state & lower_mask).bit_count() & 1) else 1.0
            out[target, state] = sign
    return out


def ideal_kitaev_chain(L: int) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    """Return H = i sum_j b_j a_{j+1} and its Majorana quadratures."""
    c = [annihilator(L, site) for site in range(L)]
    a = [item + item.conj().T for item in c]
    b = [-1j * (item - item.conj().T) for item in c]
    H = sum((1j * b[j] @ a[j + 1] for j in range(L - 1)), start=np.zeros_like(a[0]))
    if not np.allclose(H, H.conj().T, atol=1e-13):
        raise RuntimeError("ideal Kitaev Hamiltonian is not Hermitian")
    return H, a, b


def projected_code_diagnostics(H: np.ndarray, G: np.ndarray, operator: np.ndarray) -> dict:
    """Same (I-P)[H,O]P normalization used for the Iemini edge audit."""
    operated = operator @ G
    logical = G.conj().T @ operated
    leakage = operated - G @ logical
    commutator_action = H @ operated - operator @ (H @ G)
    projected_commutator = commutator_action - G @ (G.conj().T @ commutator_action)
    denominator = float(np.linalg.norm(operated))
    return {
        "operator_on_code_frobenius": denominator,
        "logical_action_frobenius": float(np.linalg.norm(logical)),
        "logical_leakage_amplitude_frobenius": float(np.linalg.norm(leakage)),
        "full_commutator_on_code_frobenius": float(np.linalg.norm(commutator_action)),
        "code_commutator_action_frobenius": float(np.linalg.norm(projected_commutator)),
        "code_commutator_action_normalized": float(np.linalg.norm(projected_commutator) / denominator),
    }


def main() -> None:
    L = 6
    H, a, _ = ideal_kitaev_chain(L)
    values, vectors = eigh(H)
    G = vectors[:, :2]
    edge = projected_code_diagnostics(H, G, a[0])
    bulk = projected_code_diagnostics(H, G, a[1])
    edge_passes = edge["code_commutator_action_normalized"] < 1e-12
    bulk_rejected = bulk["code_commutator_action_normalized"] > 1e-3
    if not (edge_passes and bulk_rejected):
        raise RuntimeError("Kitaev positive/negative harness controls did not separate")
    out = {
        "schema": "antler.phase6.kitaev-edge-harness-control.v1",
        "model": "open ideal Kitaev chain H = i sum_j b_j a_{j+1}, mu=0, t=Delta",
        "L": L,
        "fock_dimension": 1 << L,
        "two_lowest_energies": [float(value) for value in values[:2]],
        "gap_above_ground_doublet": float(values[2] - values[1]),
        "projected_commutator_definition": "||(I-P)[H,O]P||_F / ||OP||_F",
        "exact_left_edge_majorana_a0": edge,
        "bulk_negative_control_majorana_a1": bulk,
        "passes_exact_edge_control": bool(edge_passes),
        "rejects_bulk_negative_control": bool(bulk_rejected),
        "claim_boundary": (
            "This validates the numerical projected-commutator harness on an exact, number-nonconserving Kitaev zero mode. "
            "It does not validate the external Iemini generator, establish a native ANTLER mapping, or demonstrate a physical braid."
        ),
        "decision": (
            "The harness resolves an exact zero mode from a bulk operator; apply its negative result to the external parent literally."
        ),
    }
    path = ROOT / "results" / "phase6" / "kitaev_edge_harness_control.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
