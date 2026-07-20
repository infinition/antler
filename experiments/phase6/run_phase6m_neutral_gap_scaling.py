"""Phase 6M: L=10 neutral-gap check before any larger edge-support run.

The exact lambda=1 Iemini parent has a fixed-N parity doublet.  Whether the
first excitation *within* a parity sector remains separated is a distinct
finite-size question.  Here the L=10 calculation uses the compiled,
matrix-free parent action already cross-checked in Phase 6K.  An L=6 result is
recomputed with that same route before quoting L=10.

The script deliberately stops at a gap diagnosis: it does not infer a
thermodynamic gap from four sizes and does not launch an L=12 braid audit.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis
from antler.number_conserving_pairwire import wire_a_parity
from experiments.phase5.run_phase5_iemini_braid_scaling import exact_parity_frame
from experiments.phase6.run_phase6k_edge_support_scaling import _apply_iemini_parent_numba


KNOWN_GAPS = {
    4: 1.4621250708201343,
    6: 0.7174655047194777,
    8: 0.41811198021175106,
}


def make_lookup(states: np.ndarray, modes: int) -> np.ndarray:
    lookup = np.full(1 << modes, -1, dtype=np.int64)
    for row, state in enumerate(states):
        lookup[int(state)] = row
    return lookup


def neutral_gap_in_parity_sector(
    L: int,
    states: np.ndarray,
    lookup: np.ndarray,
    G: np.ndarray,
    parity: int,
) -> dict:
    """Deflate the known zero state by a rank-one positive shift."""
    rows = np.asarray(
        [row for row, state in enumerate(states) if wire_a_parity(int(state), L) == parity],
        dtype=np.int64,
    )
    ground = G[rows, parity].copy()
    if not np.isclose(np.linalg.norm(ground), 1.0, atol=1e-12):
        raise RuntimeError("exact parity ground state did not match its sector")
    shifted_ground_energy = 16.0
    states_i64 = states.astype(np.int64, copy=False)
    full_dimension = len(states)

    def parent_matvec(vector: np.ndarray) -> np.ndarray:
        full = np.zeros((full_dimension, 1), dtype=complex)
        full[rows, 0] = vector
        applied = _apply_iemini_parent_numba(L, states_i64, lookup, full)[:, 0]
        sector = applied[rows]
        return sector + shifted_ground_energy * ground * np.vdot(ground, vector)

    operator = LinearOperator(
        (len(rows), len(rows)), matvec=parent_matvec, dtype=np.complex128,
    )
    values, vectors = eigsh(
        operator,
        k=1,
        which="SA",
        tol=1e-10,
        maxiter=10000,
        ncv=min(40, len(rows) - 1),
    )
    vector = vectors[:, 0]
    residual = np.linalg.norm(parent_matvec(vector) - values[0] * vector)
    parent_on_ground = parent_matvec(ground) - shifted_ground_energy * ground
    return {
        "wire_a_parity": parity,
        "sector_dimension": int(len(rows)),
        "deflation_shift": shifted_ground_energy,
        "first_neutral_excitation": float(values[0]),
        "eigensolver_residual": float(residual),
        "parent_action_on_exact_ground": float(np.linalg.norm(parent_on_ground)),
    }


def matrix_free_neutral_gap(L: int) -> dict:
    N = L
    states, _ = build_basis(2 * L, N)
    lookup = make_lookup(states, 2 * L)
    G = exact_parity_frame(L, N, states)
    sectors = [neutral_gap_in_parity_sector(L, states, lookup, G, parity) for parity in (0, 1)]
    return {
        "L": L,
        "N": N,
        "basis_dimension": int(len(states)),
        "sectors": sectors,
        "fixed_N_neutral_code_gap": float(min(item["first_neutral_excitation"] for item in sectors)),
    }


def line_fit(x: np.ndarray, y: np.ndarray) -> dict:
    slope, intercept = np.polyfit(x, y, deg=1)
    predicted = slope * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    residual = float(np.sum((y - predicted) ** 2))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(1.0 - residual / total),
    }


def main() -> None:
    # A same-route L=6 reconstruction makes the L=10 matrix-free result auditable.
    check_l6 = matrix_free_neutral_gap(6)
    l6_error = abs(check_l6["fixed_N_neutral_code_gap"] - KNOWN_GAPS[6])
    if l6_error >= 1e-8:
        raise RuntimeError("matrix-free neutral gap did not reproduce the L=6 sparse-ED result")
    l10 = matrix_free_neutral_gap(10)
    gaps = {**KNOWN_GAPS, 10: l10["fixed_N_neutral_code_gap"]}
    lengths = np.asarray(sorted(gaps), dtype=float)
    values = np.asarray([gaps[int(L)] for L in lengths], dtype=float)
    inverse_length_fit = line_fit(1.0 / lengths, values)
    log_fit = line_fit(np.log(lengths), np.log(values))
    power_law = {
        "amplitude": float(np.exp(log_fit["intercept"])),
        "exponent_p_in_a_times_L_to_minus_p": float(-log_fit["slope"]),
        "log_space_r_squared": log_fit["r_squared"],
    }
    strictly_falling = bool(all(right < left for left, right in zip(values, values[1:])))
    decision = (
        "The neutral gap continues to fall through L=10. These finite sizes do not establish a nonzero thermodynamic neutral gap; do not spend L=12 edge-support effort as evidence for a protected braid."
        if strictly_falling else
        "The L=10 value does not continue the falling neutral-gap sequence; reassess the finite-size model before choosing the next support calculation."
    )
    out = {
        "schema": "antler.phase6.neutral-gap-scaling.v1",
        "reference": "external Iemini lambda=1 parent at fixed N=L; parity-sector excitation above its exact ground state",
        "matrix_free_L6_reconstruction": check_l6,
        "known_L6_gap_absolute_error": l6_error,
        "matrix_free_L10": l10,
        "neutral_gaps": [{"L": int(L), "gap": float(gaps[int(L)])} for L in lengths],
        "strictly_decreasing_over_L4_to_L10": strictly_falling,
        "fit_gap_vs_inverse_L": inverse_length_fit,
        "pure_power_law_descriptive_fit": power_law,
        "claim_boundary": (
            "The fits are descriptive four-size finite-size diagnostics, not a proof of a thermodynamic exponent or gap closure. "
            "This audit concerns the external Iemini parent, not a native ANTLER Hamiltonian or a physical braid."
        ),
        "decision": decision,
    }
    path = ROOT / "results" / "phase6" / "neutral_gap_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
