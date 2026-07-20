"""Calibrate the Phase 7 edge test on the revised Ising benchmark.

This is deliberately not a topological-candidate simulation.  It checks the
truncated strong-zero-mode formula at the left boundary of the transverse-field
Ising reduction and compares it to the same one-sided recurrence started in an
equal-sized bulk window.  The comparison validates a localization diagnostic
while retaining the benchmark's SSB/symmetry-only classification.
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

from antler.phase7_parent_audit import projected_edge_metrics


def xz_string(L: int, x_sites: set[int], z_sites: set[int]) -> np.ndarray:
    """Hermitian Pauli string in a bit basis with rung 0 as least-significant."""
    if x_sites & z_sites:
        raise ValueError("a site cannot carry X and Z in this real-string basis")
    dimension = 1 << L
    out = np.zeros((dimension, dimension), dtype=complex)
    x_mask = sum(1 << site for site in x_sites)
    z_mask = sum(1 << site for site in z_sites)
    for state in range(dimension):
        sign = -1.0 if ((state & z_mask).bit_count() & 1) else 1.0
        out[state ^ x_mask, state] = sign
    return out


def ising_hamiltonian(L: int, J: float, h: float) -> np.ndarray:
    identity = np.eye(1 << L, dtype=complex)
    H = np.zeros_like(identity)
    for rung in range(L - 1):
        H += 0.5 * J * (identity - xz_string(L, {rung, rung + 1}, set()))
    for rung in range(L):
        H += h * xz_string(L, set(), {rung})
    return H


def one_sided_zero_mode(L: int, r: float, start: int, support: int) -> np.ndarray:
    """Strong-zero-mode recurrence started at ``start`` and extended rightward."""
    if not 1 <= support <= L - start:
        raise ValueError("support must be between 1 and L")
    normalizer = np.sqrt((1.0 - r * r) / (1.0 - r ** (2 * support)))
    operator = np.zeros((1 << L, 1 << L), dtype=complex)
    for rung in range(support):
        site = start + rung
        # For H=-J_I XX+h Z, the recurrence is (-h/J_I)^rung.
        operator += normalizer * (-r) ** rung * xz_string(L, {site}, set(range(start, site)))
    return operator


def main() -> None:
    L, J, h = 8, 1.0, 0.1
    J_ising = J / 2.0
    r = h / J_ising
    H = ising_hamiltonian(L, J, h)
    energies, vectors = eigh(H)
    G = vectors[:, :2]
    rows = []
    for support in range(1, 5):
        operator = one_sided_zero_mode(L, r, 0, support)
        measured = projected_edge_metrics(H, G, operator)["code_commutator_action_normalized"]
        predicted = 2.0 * J_ising * r ** support * np.sqrt((1.0 - r * r) / (1.0 - r ** (2 * support)))
        bulk_start = (L - support) // 2
        bulk_operator = one_sided_zero_mode(L, r, bulk_start, support)
        bulk_metrics = projected_edge_metrics(H, G, bulk_operator)
        rows.append({
            "support": support,
            "left_truncated_mode_epsilon": measured,
            "left_formula_epsilon": float(predicted),
            "relative_formula_error": float(abs(measured - predicted) / predicted),
            "one_sided_bulk_window_start": bulk_start,
            "one_sided_bulk_mode_epsilon": bulk_metrics["code_commutator_action_normalized"],
            "one_sided_bulk_logical_action_frobenius": bulk_metrics["logical_action_frobenius"],
        })
    edge_decreases = all(
        right["left_truncated_mode_epsilon"] < left["left_truncated_mode_epsilon"]
        for left, right in zip(rows, rows[1:])
    )
    bulk_not_edge_like = all(
        row["one_sided_bulk_mode_epsilon"] > 5.0 * row["left_truncated_mode_epsilon"]
        for row in rows[1:]
    )
    out = {
        "schema": "antler.phase7.ising-edge-bulk-calibration.v1",
        "model": "open transverse-field Ising reduction of the revised Phase 7 benchmark",
        "parameters": {"L": L, "J": J, "h": h, "J_ising": J_ising, "r": r},
        "two_lowest_energies": [float(value) for value in energies[:2]],
        "gap_above_doublet": float(energies[2] - energies[1]),
        "rows": rows,
        "left_epsilon_strictly_decreases_with_support": bool(edge_decreases),
        "bulk_windows_remain_worse_than_left_windows": bool(bulk_not_edge_like),
        "claim_boundary": (
            "This calibrates an edge-versus-bulk projected-commutator diagnostic on a symmetry-restricted Ising/SSB benchmark. "
            "It does not restore full local indistinguishability, topological protection, a native microscopic realization, or a braid claim."
        ),
        "decision": (
            "Use the calibrated localization ratio only as a Phase 7 diagnostic; retain the benchmark's Ising/SSB classification."
        ),
    }
    path = ROOT / "results" / "phase7" / "ising_edge_bulk_calibration.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
