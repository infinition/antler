"""Audit the physical X/Y/Z conditional-link library in one shared-matter block.

All three rail axes are obtained from microscopic controls before downfolding:

* X: real rung hopping and a pair-channel phase echo;
* Y: the same rung hopping with a physical Peierls phase pi/2;
* Z: a signed rail-potential bias and the same pair-channel echo.

The target is always a transition between the two charge-two walker states.
No Pauli word is placed directly in the low-energy Hamiltonian.  Exact pulse
controls at the deepest registered SW point are included for every axis.
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
    pauli_coefficients,
    schur_effective,
)
from run_phase8b_shared_matter_pulse_closure_audit import phase_aligned_distance, polar_unitary


AXES = {
    "X": {
        "target": "XX",
        "description": "real rung hop sign echo",
        "segments": lambda g: ((+RUNG_HOPPING, g, g, 0.0), (-RUNG_HOPPING, g, -g, 0.0)),
    },
    "Y": {
        "target": "YX",
        "description": "pi/2 rung Peierls phase plus sign echo",
        "segments": lambda g: ((+1j * RUNG_HOPPING, g, g, 0.0), (-1j * RUNG_HOPPING, g, -g, 0.0)),
    },
    "Z": {
        "target": "ZX",
        "description": "signed rail-potential bias plus pair-channel echo",
        "segments": lambda g: ((0.0, g, g, +RUNG_HOPPING), (0.0, g, -g, -RUNG_HOPPING)),
    },
}


def make_segment(spec: tuple[complex, complex, complex, float]) -> tuple[np.ndarray, np.ndarray, dict[int, int]]:
    rung, g0, g1, bias = spec
    return build_segment(rung, g0, g1, rail_bias=bias)


def pulse_metrics(
    first: np.ndarray,
    second: np.ndarray,
    first_eff: np.ndarray,
    second_eff: np.ndarray,
    frame: np.ndarray,
    duration: float,
) -> dict[str, float]:
    cycle = expm(-1j * duration * second) @ expm(-1j * duration * first)
    target = expm(-1j * duration * second_eff) @ expm(-1j * duration * first_eff)
    raw = frame.conj().T @ cycle @ frame
    logical = polar_unitary(raw)
    projector = frame @ frame.conj().T
    target_signal = phase_aligned_distance(target, np.eye(target.shape[0], dtype=complex))
    realized_signal = phase_aligned_distance(logical, np.eye(logical.shape[0], dtype=complex))
    target_error = phase_aligned_distance(logical, target)
    return {
        "low_frame_leakage_worst": float(np.linalg.norm((np.eye(cycle.shape[0]) - projector) @ cycle @ frame, ord=2) ** 2),
        "polar_distance_to_piecewise_sw": target_error,
        "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        "piecewise_sw_signal_distance_from_identity": target_signal,
        "physical_polar_signal_distance_from_identity": realized_signal,
        "relative_piecewise_sw_error": float(target_error / target_signal) if target_signal > 1e-15 else None,
    }


def main() -> None:
    rows = []
    pulse_controls = {}
    for axis, configuration in AXES.items():
        axis_rows = []
        for ratio in RATIOS:
            coupling = ratio * DELTA_PAIR
            segment_a, segment_b = (make_segment(spec) for spec in configuration["segments"](coupling))
            h_a, _, positions = segment_a
            h_b, _, _ = segment_b
            low = code_indices(positions)
            h_eff_a, capture_a, gap_a = schur_effective(h_a, low)
            h_eff_b, capture_b, gap_b = schur_effective(h_b, low)
            coefficients = pauli_coefficients(0.5 * (h_eff_a + h_eff_b))
            target = configuration["target"]
            unwanted = {label: value for label, value in coefficients.items() if label not in {"II", target}}
            row = {
                "axis": axis,
                "coupling_over_detuning": ratio,
                "target_pauli": target,
                "target_coefficient": coefficients[target],
                "maximum_unwanted_non_scalar_coefficient": float(max(abs(value) for value in unwanted.values())),
                "minimum_low_frame_capture": min(capture_a, capture_b),
                "minimum_low_high_gap": min(gap_a, gap_b),
            }
            rows.append(row)
            axis_rows.append(row)
            if abs(ratio - 0.025) < 1e-12:
                frame = np.zeros((h_a.shape[0], len(low)), dtype=complex)
                frame[low, np.arange(len(low))] = 1.0
                duration = float(2.0 * np.pi / np.sqrt(DELTA_PAIR**2 + 4.0 * coupling**2))
                pulse_controls[axis] = {
                    "target_pauli": target,
                    "segment_duration": duration,
                    **pulse_metrics(h_a, h_b, h_eff_a, h_eff_b, frame, duration),
                }
        deep = [row for row in axis_rows if row["coupling_over_detuning"] <= 0.075]
        pulse_controls[axis]["deep_sw_target_power"] = float(np.polyfit(
            np.log([row["coupling_over_detuning"] for row in deep]),
            np.log(np.abs([row["target_coefficient"] for row in deep])),
            1,
        )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-xyz-link-audit.v1",
        "parameters": {
            "pair_detuning": DELTA_PAIR,
            "rung_hopping_magnitude": RUNG_HOPPING,
            "walker_transition": "X_(d0,d1)",
            "axes": {axis: item["description"] for axis, item in AXES.items()},
        },
        "rows": rows,
        "deep_point_pulse_controls": pulse_controls,
        "decision": "The shared-matter controls derive selective conditional X, Y and Z axes in the static/downfolded Hamiltonian. The registered integer-Rabi pulse controls are rejected as gates because their relative error to the non-scalar SW signal is approximately one; they return almost scalar low-frame operations.",
        "claim_boundary": "This is a static one-link compiler library, not a pulse-realized link library. It does not derive a common four-state walker with all three link types, a simultaneous multi-link Hamiltonian, a stabilizer patch, a defect/fusion space, non-Abelian braid, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_xyz_link_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
