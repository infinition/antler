"""Deterministic classical-optimizer baseline before any RL control search.

The registered shared-matter local block has only two free hold durations once
the physically specified smooth zero-to-zero ramps are fixed.  This is small
enough for exact propagation and a seeded global classical search; using an
adaptive optimizer first would add no exploration capability and would obscure a
failed control grammar.  The baseline asks for a nontrivial `exp(-i 0.1 XX)`
gate while explicitly penalizing low-frame leakage.

An empty baseline is a negative optimization result, not a mathematical no-go
over all possible pulse families.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.optimize import differential_evolution


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    DELTA_PAIR,
    PAULIS,
    RUNG_HOPPING,
    build_segment,
    code_indices,
)
from run_phase8b_shared_matter_pulse_closure_audit import polar_unitary, phase_aligned_distance


RATIO = 0.05
RAMP_DURATION = 4.0
RAMP_SLICES = 32
TARGET_XX_ANGLE = 0.10
HOLD_BOUNDS = (0.0, 400.0)
SEEDS = (13, 19, 31)
LEAKAGE_TARGET = 1e-4
RELATIVE_TARGET_ERROR = 0.10


def pulse_parts(
    j_perp: complex,
    g0_sign: float,
    g1_sign: float,
    coupling: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    amplitudes = [float(np.sin(0.5 * np.pi * (index + 0.5) / RAMP_SLICES) ** 2) for index in range(RAMP_SLICES)]
    up = None
    down = None
    for amplitude in amplitudes:
        hamiltonian, _, _ = build_segment(j_perp, g0_sign * amplitude * coupling, g1_sign * amplitude * coupling)
        step = expm(-1j * RAMP_DURATION / RAMP_SLICES * hamiltonian)
        up = step if up is None else step @ up
    for amplitude in reversed(amplitudes):
        hamiltonian, _, _ = build_segment(j_perp, g0_sign * amplitude * coupling, g1_sign * amplitude * coupling)
        step = expm(-1j * RAMP_DURATION / RAMP_SLICES * hamiltonian)
        down = step if down is None else step @ down
    full, _, _ = build_segment(j_perp, g0_sign * coupling, g1_sign * coupling)
    eigenvalues, eigenvectors = np.linalg.eigh(full)
    if up is None or down is None:
        raise RuntimeError("empty ramp")
    return up, down, eigenvalues, eigenvectors


def held_propagator(parts: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], duration: float) -> np.ndarray:
    up, down, eigenvalues, eigenvectors = parts
    hold = (eigenvectors * np.exp(-1j * duration * eigenvalues)) @ eigenvectors.conj().T
    return down @ hold @ up


def main() -> None:
    coupling = RATIO * DELTA_PAIR
    _, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
    low = code_indices(positions)
    frame = np.zeros((len(states), len(low)), dtype=complex)
    frame[low, np.arange(len(low))] = 1.0
    projector = frame @ frame.conj().T
    target = expm(-1j * TARGET_XX_ANGLE * np.kron(PAULIS["X"], PAULIS["X"]))
    target_signal = phase_aligned_distance(target, np.eye(4, dtype=complex))
    parts_a = pulse_parts(+RUNG_HOPPING, +1.0, +1.0, coupling)
    parts_b = pulse_parts(-RUNG_HOPPING, +1.0, -1.0, coupling)

    def evaluate(parameters: np.ndarray) -> dict[str, float]:
        duration_a, duration_b = map(float, parameters)
        physical = held_propagator(parts_b, duration_b) @ held_propagator(parts_a, duration_a)
        raw = frame.conj().T @ physical @ frame
        logical = polar_unitary(raw)
        relative_error = phase_aligned_distance(logical, target) / target_signal
        leakage = float(np.linalg.norm((np.eye(len(states)) - projector) @ physical @ frame, ord=2) ** 2)
        return {
            "duration_a": duration_a,
            "duration_b": duration_b,
            "relative_target_error": float(relative_error),
            "low_frame_leakage_worst": leakage,
            "physical_polar_signal_distance_from_identity": phase_aligned_distance(logical, np.eye(4, dtype=complex)),
        }

    def objective(parameters: np.ndarray) -> float:
        metrics = evaluate(parameters)
        return metrics["relative_target_error"] + 1.0e4 * max(0.0, metrics["low_frame_leakage_worst"] - 1.0e-5)

    rows = []
    for seed in SEEDS:
        result = differential_evolution(
            objective,
            bounds=(HOLD_BOUNDS, HOLD_BOUNDS),
            seed=seed,
            popsize=24,
            maxiter=100,
            tol=1e-9,
            polish=True,
            workers=1,
            updating="immediate",
        )
        row = {
            "seed": seed,
            "objective": float(result.fun),
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            **evaluate(result.x),
        }
        row["passes_local_screen"] = bool(
            row["relative_target_error"] < RELATIVE_TARGET_ERROR
            and row["low_frame_leakage_worst"] < LEAKAGE_TARGET
        )
        rows.append(row)
    output = {
        "schema": "antler.phase8b.shared-matter-classical-control-baseline.v1",
        "parameters": {
            "coupling_over_detuning": RATIO,
            "ramp_shape": "sin^2 zero-to-zero",
            "ramp_duration": RAMP_DURATION,
            "ramp_slices": RAMP_SLICES,
            "target": "exp(-i 0.1 XX)",
            "hold_bounds": list(HOLD_BOUNDS),
            "seeds": list(SEEDS),
            "screen": {"relative_target_error": RELATIVE_TARGET_ERROR, "low_frame_leakage_worst": LEAKAGE_TARGET},
        },
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_local_screen"]],
        "decision": (
            "No registered seeded differential-evolution run finds a nontrivial low-leakage XX gate in the two-duration smooth-pulse box."
            if not any(row["passes_local_screen"] for row in rows)
            else "At least one registered classical baseline row passes the local screen."
        ),
        "claim_boundary": "This is a reproducible classical-optimizer baseline in a fixed two-parameter control box. Its failure does not prove that all smooth/composite/dressed pulse families or new microscopic resources fail.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_classical_control_baseline.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
