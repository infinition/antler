"""Four-sign shared-matter echo with zero microscopic control average.

The two-word echo leaves one pair channel in the microscopic average, which
produces the documented `IZ` obstruction.  This is the minimal fully balanced
extension.  The four sign words obey `s_J=s_0*s_1`; hence every word has the
same static Schrieffer--Wolff `XX` sign, while each microscopic control has
zero group average.  It is therefore the correct next test of whether a
slow/dressed group can realize the average-of-SW compiler without the
two-word fast-limit spectator.

No gate is assumed: exact finite-block propagation is compared to both the
ordered SW product and the derived `XX` target.
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

from run_phase8b_shared_matter_adiabatic_echo_preflight import smooth_segment
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


SIGN_WORDS = ((+1, +1, +1), (-1, +1, -1), (-1, -1, +1), (+1, -1, -1))
RATIOS = (0.05, 0.025)
RAMP_DURATIONS = (0.25, 0.5, 1.0, 2.0, 4.0)
CORE_DURATIONS = (0.025, 0.05, 0.10)
TARGET_XX_ANGLE = 0.10
RELATIVE_ERROR_TARGET = 0.10
LEAKAGE_TARGET = 1e-4


def main() -> None:
    rows = []
    xx = np.kron(PAULIS["X"], PAULIS["X"])
    static_rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA_PAIR
        hamiltonians = []
        effective_words = []
        positions = None
        for sign_j, sign_0, sign_1 in SIGN_WORDS:
            hamiltonian, _, positions = build_segment(
                sign_j * RUNG_HOPPING, sign_0 * coupling, sign_1 * coupling
            )
            low = code_indices(positions)
            effective, _, _ = schur_effective(hamiltonian, low)
            hamiltonians.append(hamiltonian)
            effective_words.append(effective)
        if positions is None:
            raise RuntimeError("missing fixed-charge basis")
        static_average = sum(effective_words) / len(effective_words)
        microscopic_average = sum(hamiltonians) / len(hamiltonians)
        sw_microscopic_average, _, _ = schur_effective(microscopic_average, low)
        static_coefficients = pauli_coefficients(static_average)
        microscopic_coefficients = pauli_coefficients(sw_microscopic_average)
        static_rows.append({
            "coupling_over_detuning": ratio,
            "static_group_xx": static_coefficients["XX"],
            "static_group_maximum_unwanted": float(max(
                abs(value) for label, value in static_coefficients.items() if label not in {"II", "XX"}
            )),
            "sw_of_microscopic_group_average_xx": microscopic_coefficients["XX"],
            "sw_of_microscopic_group_average_maximum_non_scalar": float(max(
                abs(value) for label, value in microscopic_coefficients.items() if label != "II"
            )),
        })
        reference, states, positions = build_segment(+RUNG_HOPPING, coupling, coupling)
        frame = np.zeros((len(states), 4), dtype=complex)
        frame[low, np.arange(4)] = 1.0
        projector = frame @ frame.conj().T
        for ramp_duration in RAMP_DURATIONS:
            for core_duration in CORE_DURATIONS:
                physical_words = []
                sw_words = []
                phase_per_group = 0.0
                for sign_j, sign_0, sign_1 in SIGN_WORDS:
                    physical, sw, phase = smooth_segment(
                        sign_j * RUNG_HOPPING,
                        float(sign_0),
                        float(sign_1),
                        coupling,
                        ramp_duration,
                        core_duration,
                        low,
                    )
                    physical_words.append(physical)
                    sw_words.append(sw)
                    phase_per_group += phase
                physical_group = np.eye(len(states), dtype=complex)
                sw_group = np.eye(4, dtype=complex)
                for physical, sw in zip(physical_words, sw_words, strict=True):
                    physical_group = physical @ physical_group
                    sw_group = sw @ sw_group
                repetitions = max(1, int(round(TARGET_XX_ANGLE / abs(phase_per_group))))
                physical = np.linalg.matrix_power(physical_group, repetitions)
                sw_target = np.linalg.matrix_power(sw_group, repetitions)
                ideal_xx = expm(-1j * repetitions * phase_per_group * xx)
                raw = frame.conj().T @ physical @ frame
                logical = polar_unitary(raw)
                sw_signal = phase_aligned_distance(sw_target, np.eye(4, dtype=complex))
                ideal_signal = phase_aligned_distance(ideal_xx, np.eye(4, dtype=complex))
                physical_to_sw = phase_aligned_distance(logical, sw_target)
                physical_to_ideal = phase_aligned_distance(logical, ideal_xx)
                sw_to_ideal = phase_aligned_distance(sw_target, ideal_xx)
                leakage = float(np.linalg.norm((np.eye(len(states)) - projector) @ physical @ frame, ord=2) ** 2)
                rows.append({
                    "coupling_over_detuning": ratio,
                    "ramp_shape": "sin^2 zero-to-zero",
                    "ramp_duration": ramp_duration,
                    "core_duration": core_duration,
                    "repetitions": repetitions,
                    "total_duration": float(repetitions * len(SIGN_WORDS) * (2.0 * ramp_duration + core_duration)),
                    "target_xx_phase": float(repetitions * phase_per_group),
                    "time_ordered_sw_signal_distance_from_identity": sw_signal,
                    "time_ordered_sw_to_ideal_xx_distance": sw_to_ideal,
                    "physical_polar_signal_distance_from_identity": phase_aligned_distance(logical, np.eye(4, dtype=complex)),
                    "relative_physical_to_sw_error": float(physical_to_sw / sw_signal) if sw_signal > 1e-15 else None,
                    "relative_physical_to_ideal_xx_error": float(physical_to_ideal / ideal_signal) if ideal_signal > 1e-15 else None,
                    "low_frame_leakage_worst": leakage,
                    "passes_local_screen": bool(
                        physical_to_sw / sw_signal < RELATIVE_ERROR_TARGET
                        and physical_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                        and sw_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                        and leakage < LEAKAGE_TARGET
                    ),
                })
    output = {
        "schema": "antler.phase8b.shared-matter-four-sign-group-preflight.v1",
        "parameters": {
            "sign_words": [list(word) for word in SIGN_WORDS],
            "group_property": "s_J=s_0*s_1, all one-control microscopic averages vanish, and every static XX sign agrees",
            "ratios": list(RATIOS),
            "ramp_durations": list(RAMP_DURATIONS),
            "core_durations": list(CORE_DURATIONS),
            "target_xx_angle": TARGET_XX_ANGLE,
            "screen": {
                "relative_physical_to_sw_error": RELATIVE_ERROR_TARGET,
                "relative_physical_to_ideal_xx_error": RELATIVE_ERROR_TARGET,
                "low_frame_leakage_worst": LEAKAGE_TARGET,
            },
        },
        "static_rows": static_rows,
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_local_screen"]],
        "decision": "Pending execution.",
        "claim_boundary": "This is a local four-sign group preflight. It does not insert a gate and cannot establish a neutral walker, a tiled code, protection, defects, fusion, a non-Abelian braid, universality or fault tolerance.",
    }
    output["decision"] = (
        "At least one registered four-sign smooth group row passes the local signal/leakage screen."
        if output["passing_rows"]
        else "No registered four-sign smooth group row passes the local signal/leakage screen."
    )
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_four_sign_group_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
