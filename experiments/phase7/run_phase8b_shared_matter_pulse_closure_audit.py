"""Finite-pulse closure test for the shared-matter conditional-link bridge.

The preceding static audit downfolds two sign-correlated microscopic segments.
Here those same finite Hamiltonians are propagated exactly.  Each segment
duration is an integer virtual pair--mediator Rabi period in the isolated
two-level limit.  The low-frame polar unitary is compared both with the
piecewise Schrieffer--Wolff target and with its first-order averaged Floquet
Hamiltonian.  This is deliberately a closure *test*, not an assumed result.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    DELTA_PAIR,
    RATIOS,
    RUNG_HOPPING,
    build_segment,
    code_indices,
    schur_effective,
)


RABI_MULTIPLIERS = (1, 2, 4, 8)
REPETITIONS = (1, 4, 16)


def polar_unitary(matrix: np.ndarray) -> np.ndarray:
    left, _, right = np.linalg.svd(matrix)
    return left @ right


def phase_aligned_distance(left: np.ndarray, right: np.ndarray) -> float:
    overlap = np.trace(right.conj().T @ left)
    phase = 1.0 if abs(overlap) < 1e-15 else np.exp(-1j * np.angle(overlap))
    return float(np.linalg.norm(phase * left - right, ord="fro") / np.sqrt(left.shape[0]))


def metrics(
    unitary: np.ndarray,
    target_piecewise: np.ndarray,
    target_average: np.ndarray,
    frame: np.ndarray,
    projector: np.ndarray,
) -> dict[str, float]:
    raw = frame.conj().T @ unitary @ frame
    logical = polar_unitary(raw)
    target_signal = phase_aligned_distance(target_piecewise, np.eye(target_piecewise.shape[0], dtype=complex))
    realized_signal = phase_aligned_distance(logical, np.eye(logical.shape[0], dtype=complex))
    target_error = phase_aligned_distance(logical, target_piecewise)
    return {
        "low_frame_leakage_worst": float(np.linalg.norm((np.eye(unitary.shape[0]) - projector) @ unitary @ frame, ord=2) ** 2),
        "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        "polar_distance_to_piecewise_sw": target_error,
        "polar_distance_to_average_sw": phase_aligned_distance(logical, target_average),
        "piecewise_sw_signal_distance_from_identity": target_signal,
        "physical_polar_signal_distance_from_identity": realized_signal,
        "relative_piecewise_sw_error": float(target_error / target_signal) if target_signal > 1e-15 else None,
    }


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        h_a, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        h_b, states_b, positions_b = build_segment(-RUNG_HOPPING, coupling, -coupling)
        if not np.array_equal(states, states_b) or positions != positions_b:
            raise RuntimeError("the two echoed segments do not share a basis")
        low = code_indices(positions)
        h_eff_a, _, _ = schur_effective(h_a, low)
        h_eff_b, _, _ = schur_effective(h_b, low)
        h_eff_average = 0.5 * (h_eff_a + h_eff_b)
        frame = np.zeros((len(states), len(low)), dtype=complex)
        frame[low, np.arange(len(low))] = 1.0
        projector = frame @ frame.conj().T
        omega = float(np.sqrt(DELTA_PAIR**2 + 4.0 * coupling**2))
        for multiplier in RABI_MULTIPLIERS:
            segment_duration = float(2.0 * np.pi * multiplier / omega)
            physical_cycle = expm(-1j * segment_duration * h_b) @ expm(-1j * segment_duration * h_a)
            sw_piecewise_cycle = expm(-1j * segment_duration * h_eff_b) @ expm(-1j * segment_duration * h_eff_a)
            sw_average_cycle = expm(-1j * 2.0 * segment_duration * h_eff_average)
            for repetitions in REPETITIONS:
                rows.append({
                    "coupling_over_detuning": ratio,
                    "rabi_multiplier_per_segment": multiplier,
                    "segment_duration": segment_duration,
                    "total_duration": float(2.0 * segment_duration * repetitions),
                    "repetitions": repetitions,
                    **metrics(
                        np.linalg.matrix_power(physical_cycle, repetitions),
                        np.linalg.matrix_power(sw_piecewise_cycle, repetitions),
                        np.linalg.matrix_power(sw_average_cycle, repetitions),
                        frame,
                        projector,
                    ),
                })
    selected = [
        row for row in rows
        if row["rabi_multiplier_per_segment"] == 1 and row["repetitions"] == 1
    ]
    deep = [row for row in selected if row["coupling_over_detuning"] <= 0.075]
    leakage_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log([row["low_frame_leakage_worst"] for row in deep]),
        1,
    )[0])
    best = min(rows, key=lambda row: (row["low_frame_leakage_worst"], row["polar_distance_to_piecewise_sw"]))
    # The local Peierls phase phi=pi/2 is the physical Y-axis version of the
    # same bridge.  Keep this as a separately serialized pulse control rather
    # than inferring it only from the static phase-link audit.
    coupling = 0.025 * DELTA_PAIR
    h_y_a, states_y, positions_y = build_segment(1j * RUNG_HOPPING, coupling, coupling)
    h_y_b, states_y_b, positions_y_b = build_segment(-1j * RUNG_HOPPING, coupling, -coupling)
    if not np.array_equal(states_y, states_y_b) or positions_y != positions_y_b:
        raise RuntimeError("the pi/2 echoed segments do not share a basis")
    low_y = code_indices(positions_y)
    h_y_eff_a, _, _ = schur_effective(h_y_a, low_y)
    h_y_eff_b, _, _ = schur_effective(h_y_b, low_y)
    frame_y = np.zeros((len(states_y), len(low_y)), dtype=complex)
    frame_y[low_y, np.arange(len(low_y))] = 1.0
    projector_y = frame_y @ frame_y.conj().T
    duration_y = float(2.0 * np.pi / np.sqrt(DELTA_PAIR**2 + 4.0 * coupling**2))
    y_cycle = expm(-1j * duration_y * h_y_b) @ expm(-1j * duration_y * h_y_a)
    y_piecewise = expm(-1j * duration_y * h_y_eff_b) @ expm(-1j * duration_y * h_y_eff_a)
    y_average = expm(-1j * duration_y * (h_y_eff_a + h_y_eff_b))
    y_control = {
        "coupling_over_detuning": 0.025,
        "rung_peierls_phase": float(np.pi / 2.0),
        "segment_duration": duration_y,
        **metrics(y_cycle, y_piecewise, y_average, frame_y, projector_y),
    }
    output = {
        "schema": "antler.phase8b.shared-matter-pulse-closure-audit.v1",
        "parameters": {
            "pair_detuning": DELTA_PAIR,
            "rabi_duration": "2*pi*m/sqrt(Delta_pair^2+4g^2) per microscopic segment",
            "segments": "A=(+J,+g,+g), B=(-J,+g,-g)",
            "comparison": "exact physical AB pulse against piecewise and averaged Schrieffer-Wolff low-frame targets",
        },
        "rows": rows,
        "single_cycle_deep_leakage_power": leakage_power,
        "best_registered_point": best,
        "pi_over_2_phase_pulse_control": y_control,
        "decision": "The registered integer virtual-Rabi closure is rejected as a phase-accumulating gate: its absolute error is small only because the SW target is small, while its relative error to the target signal is approximately one and the physical polar operation is nearly scalar.",
        "claim_boundary": "This rejects the isolated integer-Rabi timing family as a realization of the shared-matter SW link. It does not refute the static/downfolded bridge, an off-resonant dressed or adiabatic pulse schedule, a complete neutral walker, a tiled code, defects, fusion or a non-Abelian braid.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_pulse_closure_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
