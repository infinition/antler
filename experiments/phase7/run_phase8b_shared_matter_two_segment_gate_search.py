"""Classical exact search for a nontrivial two-segment shared-matter gate.

Integer Rabi closures return the low frame but also erase the desired
conditional signal.  Before proposing a new analytic pulse construction, this
bounded exact search asks a minimal question: do arbitrary durations of the
two already registered physical segments A and B contain a low-leakage,
nontrivial XX rotation at all?

The search is deliberately modest and fully reproducible: a rectangular grid
of two durations, exact 12-state propagation from precomputed eigensystems,
and fidelity to the best nonzero ``exp(-i phi XX)`` target.  It is a negative
or candidate-finding gate, not a black-box optimizer and not a gate claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    DELTA_PAIR,
    RUNG_HOPPING,
    build_segment,
    code_indices,
)
from run_phase8b_shared_matter_pulse_closure_audit import polar_unitary


RATIOS = (0.10, 0.05, 0.025)
TIME_MAX = 8.0
TIME_POINTS = 121
PHI_GRID = np.linspace(0.05, np.pi / 2.0, 301)
LEAKAGE_LIMIT = 1e-4
FIDELITY_LIMIT = 0.999


def propagator(eigenvalues: np.ndarray, eigenvectors: np.ndarray, duration: float) -> np.ndarray:
    return (eigenvectors * np.exp(-1j * duration * eigenvalues)) @ eigenvectors.conj().T


def best_xx_rotation(unitary: np.ndarray, xx: np.ndarray) -> tuple[float, float]:
    """Return (phi, normalized trace fidelity) for phi in the registered grid."""
    scalar = np.trace(unitary)
    xx_component = 1j * np.trace(xx @ unitary)
    overlap = np.cos(PHI_GRID) * scalar + np.sin(PHI_GRID) * xx_component
    index = int(np.argmax(np.abs(overlap)))
    return float(PHI_GRID[index]), float(abs(overlap[index]) / unitary.shape[0])


def main() -> None:
    xx = np.kron(np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex), np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex))
    times = np.linspace(0.0, TIME_MAX, TIME_POINTS)
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        h_a, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        h_b, states_b, positions_b = build_segment(-RUNG_HOPPING, coupling, -coupling)
        if not np.array_equal(states, states_b) or positions != positions_b:
            raise RuntimeError("search segment bases disagree")
        low = code_indices(positions)
        frame = np.zeros((len(states), len(low)), dtype=complex)
        frame[low, np.arange(len(low))] = 1.0
        values_a, vectors_a = np.linalg.eigh(h_a)
        values_b, vectors_b = np.linalg.eigh(h_b)
        evolutions_a = [propagator(values_a, vectors_a, time) for time in times]
        evolutions_b = [propagator(values_b, vectors_b, time) for time in times]
        candidates = []
        for index_a, u_a in enumerate(evolutions_a):
            state_after_a = u_a @ frame
            for index_b, u_b in enumerate(evolutions_b):
                evolved = u_b @ state_after_a
                raw = frame.conj().T @ evolved
                logical = polar_unitary(raw)
                phi, fidelity = best_xx_rotation(logical, xx)
                leakage = float(np.linalg.norm(evolved - frame @ raw, ord=2) ** 2)
                candidates.append({
                    "coupling_over_detuning": ratio,
                    "duration_a": float(times[index_a]),
                    "duration_b": float(times[index_b]),
                    "best_nonzero_xx_angle": phi,
                    "best_nonzero_xx_trace_fidelity": fidelity,
                    "low_frame_leakage_worst": leakage,
                    "passes_screen": bool(leakage < LEAKAGE_LIMIT and fidelity > FIDELITY_LIMIT),
                })
        candidates.sort(key=lambda row: (row["passes_screen"], row["best_nonzero_xx_trace_fidelity"] - 10.0 * row["low_frame_leakage_worst"]), reverse=True)
        rows.extend(candidates[:10])
    passing = [row for row in rows if row["passes_screen"]]
    output = {
        "schema": "antler.phase8b.shared-matter-two-segment-gate-search.v1",
        "parameters": {
            "coupling_ratios": list(RATIOS),
            "duration_range": [0.0, TIME_MAX],
            "duration_grid_points": TIME_POINTS,
            "nonzero_xx_angle_range": [float(PHI_GRID[0]), float(PHI_GRID[-1])],
            "screen": {"low_frame_leakage_worst": LEAKAGE_LIMIT, "best_nonzero_xx_trace_fidelity": FIDELITY_LIMIT},
        },
        "top_candidates_per_ratio": rows,
        "passing_candidates": passing,
        "decision": "No nontrivial low-leakage XX candidate is present in the registered two-segment duration box: all top rows collapse to the minimum allowed target angle with fidelity below the screen, so the box is rejected as a gate search.",
        "claim_boundary": "A passing grid point would not prove an analytic control law, Schrieffer-Wolff convergence, robustness, a walker loop, a code, fusion, a non-Abelian braid, universality or fault tolerance. An empty result rejects only this two-segment duration box and target family.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_two_segment_gate_search.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
