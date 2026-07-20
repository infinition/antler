"""Exact smooth-envelope preflight for the shared-matter SW echo.

Abrupt A/B switching and integer virtual-Rabi returns were both rejected:
they either populate the mediator or erase the tiny conditional signal.  This
script tests the genuinely different off-resonant hypothesis.  Every A and B
segment starts and ends at zero pair conversion with a smooth sin^2 envelope,
so the bare low frame is restored before the sign-correlated echo changes
branch.  The full finite-dimensional propagator is compared with the
time-ordered instantaneous Schrieffer--Wolff target and with its isolated XX
component.  The target is therefore not inserted by hand.

This is a local-control preflight.  A passing row would still require timing,
noise, multi-link, code, defect and braid audits before any promotion.
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
    PAULIS,
    RUNG_HOPPING,
    build_segment,
    code_indices,
    pauli_coefficients,
    schur_effective,
)
from run_phase8b_shared_matter_pulse_closure_audit import polar_unitary, phase_aligned_distance


RATIOS = (0.05, 0.025)
RAMP_DURATIONS = (0.25, 0.5, 1.0, 2.0, 4.0)
CORE_DURATIONS = (0.025, 0.05, 0.10)
RAMP_SLICES = 16
TARGET_XX_ANGLE = 0.10
RELATIVE_ERROR_TARGET = 0.10
LEAKAGE_TARGET = 1e-4


def smooth_segment(
    j_perp: complex,
    g0_sign: float,
    g1_sign: float,
    coupling: float,
    ramp_duration: float,
    core_duration: float,
    low: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return exact and instantaneous-SW propagators for one zero-to-zero pulse."""
    amplitudes = [
        float(np.sin(0.5 * np.pi * (index + 0.5) / RAMP_SLICES) ** 2)
        for index in range(RAMP_SLICES)
    ]
    amplitudes = amplitudes + [1.0] + list(reversed(amplitudes))
    durations = [ramp_duration / RAMP_SLICES] * RAMP_SLICES + [core_duration] + [ramp_duration / RAMP_SLICES] * RAMP_SLICES
    physical = None
    sw = np.eye(len(low), dtype=complex)
    xx_integral = 0.0
    for amplitude, duration in zip(amplitudes, durations, strict=True):
        hamiltonian, _, _ = build_segment(
            j_perp,
            g0_sign * amplitude * coupling,
            g1_sign * amplitude * coupling,
        )
        effective, _, _ = schur_effective(hamiltonian, low)
        step = expm(-1j * duration * hamiltonian)
        physical = step if physical is None else step @ physical
        sw = expm(-1j * duration * effective) @ sw
        xx_integral += duration * pauli_coefficients(effective)["XX"]
    if physical is None:
        raise RuntimeError("empty pulse")
    return physical, sw, float(xx_integral)


def main() -> None:
    rows: list[dict[str, float | int | bool]] = []
    xx = np.kron(PAULIS["X"], PAULIS["X"])
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        reference_hamiltonian, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        low = code_indices(positions)
        frame = np.zeros((len(states), len(low)), dtype=complex)
        frame[low, np.arange(len(low))] = 1.0
        projector = frame @ frame.conj().T
        for ramp_duration in RAMP_DURATIONS:
            for core_duration in CORE_DURATIONS:
                physical_a, sw_a, phase_a = smooth_segment(
                    +RUNG_HOPPING, +1.0, +1.0, coupling, ramp_duration, core_duration, low
                )
                physical_b, sw_b, phase_b = smooth_segment(
                    -RUNG_HOPPING, +1.0, -1.0, coupling, ramp_duration, core_duration, low
                )
                physical_cycle = physical_b @ physical_a
                sw_cycle = sw_b @ sw_a
                xx_phase_per_cycle = phase_a + phase_b
                repetitions = max(1, int(round(TARGET_XX_ANGLE / abs(xx_phase_per_cycle))))
                physical = np.linalg.matrix_power(physical_cycle, repetitions)
                sw_target = np.linalg.matrix_power(sw_cycle, repetitions)
                ideal_xx = expm(-1j * repetitions * xx_phase_per_cycle * xx)
                raw = frame.conj().T @ physical @ frame
                logical = polar_unitary(raw)
                target_signal = phase_aligned_distance(sw_target, np.eye(4, dtype=complex))
                ideal_signal = phase_aligned_distance(ideal_xx, np.eye(4, dtype=complex))
                physical_to_sw = phase_aligned_distance(logical, sw_target)
                physical_to_ideal = phase_aligned_distance(logical, ideal_xx)
                target_to_ideal = phase_aligned_distance(sw_target, ideal_xx)
                leakage = float(np.linalg.norm((np.eye(len(states)) - projector) @ physical @ frame, ord=2) ** 2)
                rows.append({
                    "coupling_over_detuning": ratio,
                    "ramp_shape": "sin^2 zero-to-zero",
                    "ramp_duration": ramp_duration,
                    "core_duration": core_duration,
                    "ramp_slices": RAMP_SLICES,
                    "repetitions": repetitions,
                    "total_duration": float(repetitions * 2.0 * (2.0 * ramp_duration + core_duration)),
                    "target_xx_phase": float(repetitions * xx_phase_per_cycle),
                    "time_ordered_sw_signal_distance_from_identity": target_signal,
                    "ideal_xx_signal_distance_from_identity": ideal_signal,
                    "time_ordered_sw_to_ideal_xx_distance": target_to_ideal,
                    "physical_polar_signal_distance_from_identity": phase_aligned_distance(logical, np.eye(4, dtype=complex)),
                    "relative_physical_to_sw_error": float(physical_to_sw / target_signal) if target_signal > 1e-15 else None,
                    "relative_physical_to_ideal_xx_error": float(physical_to_ideal / ideal_signal) if ideal_signal > 1e-15 else None,
                    "low_frame_leakage_worst": leakage,
                    "passes_local_screen": bool(
                        physical_to_sw / target_signal < RELATIVE_ERROR_TARGET
                        and physical_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                        and target_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                        and leakage < LEAKAGE_TARGET
                    ),
                })
    output = {
        "schema": "antler.phase8b.shared-matter-adiabatic-echo-preflight.v1",
        "parameters": {
            "ratios": list(RATIOS),
            "ramp_durations": list(RAMP_DURATIONS),
            "core_durations": list(CORE_DURATIONS),
            "ramp_slices": RAMP_SLICES,
            "target_xx_angle": TARGET_XX_ANGLE,
            "screen": {
                "relative_physical_to_sw_error": RELATIVE_ERROR_TARGET,
                "relative_physical_to_ideal_xx_error": RELATIVE_ERROR_TARGET,
                "relative_sw_to_ideal_xx_error": RELATIVE_ERROR_TARGET,
                "low_frame_leakage_worst": LEAKAGE_TARGET,
            },
        },
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_local_screen"]],
        "decision": "Pending execution.",
        "claim_boundary": "This local smooth-envelope test neither inserts an XX gate nor claims a complete walker, a tiled code, protected defects, fusion, a non-Abelian braid, universality or fault tolerance.",
    }
    output["decision"] = (
        "At least one smooth off-resonant row passes the stated local signal/leakage screen."
        if output["passing_rows"]
        else "No registered smooth off-resonant echo row passes the stated local signal/leakage screen."
    )
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_adiabatic_echo_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
