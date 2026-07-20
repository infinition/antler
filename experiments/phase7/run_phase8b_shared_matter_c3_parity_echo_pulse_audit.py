"""Exact finite-pulse test of the C3 even-parity echo compiler.

The static C3 group average cancels the XIX companion exactly after independent
Schur downfolding.  This audit propagates the full 1488-state microscopic
blocks through the actual eight-segment schedule (A then B for each of four
sign words) and compares the polar low-frame operation with both the ordered
piecewise SW product and the averaged-SW exponential.

It is the gate that distinguishes a useful refocusing compiler from a merely
formal average of effective Hamiltonians.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_c3_parity_echo_audit import PARITY_EVEN_SIGNS
from run_phase8b_shared_matter_c3_walker_preflight import (
    PAIR_DETUNING,
    build_segment,
    code_indices,
    pauli_coefficients,
    schur_zero,
)
from run_phase8b_shared_matter_pulse_closure_audit import phase_aligned_distance, polar_unitary


RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)


def main() -> None:
    rows = []
    for ratio in RATIOS:
        coupling = ratio * PAIR_DETUNING
        duration = float(2.0 * np.pi / np.sqrt(PAIR_DETUNING**2 + 4.0 * coupling**2))
        physical_segments: list[np.ndarray] = []
        effective_segments: list[np.ndarray] = []
        frame: np.ndarray | None = None
        grouped_effective = np.zeros((8, 8), dtype=complex)
        for signs in PARITY_EVEN_SIGNS:
            h_a, _, positions = build_segment(coupling, echoed=False, link_signs=signs)
            h_b, _, _ = build_segment(coupling, echoed=True, link_signs=signs)
            low = code_indices(positions)
            if frame is None:
                frame = np.zeros((h_a.shape[0], len(low)), dtype=complex)
                frame[low, np.arange(len(low))] = 1.0
            h_eff_a = schur_zero(h_a, low)
            h_eff_b = schur_zero(h_b, low)
            physical_segments.extend((h_a, h_b))
            effective_segments.extend((h_eff_a, h_eff_b))
            grouped_effective += h_eff_a + h_eff_b
        if frame is None:
            raise RuntimeError("C3 frame was not initialized")
        grouped_effective /= len(effective_segments)
        evolved = frame.copy()
        target_piecewise = np.eye(frame.shape[1], dtype=complex)
        for physical, effective in zip(physical_segments, effective_segments):
            evolved = expm_multiply((-1j * duration) * csr_matrix(physical), evolved)
            target_piecewise = expm(-1j * duration * effective) @ target_piecewise
        raw = frame.conj().T @ evolved
        logical = polar_unitary(raw)
        projected = frame @ raw
        target_average = expm(-1j * duration * len(physical_segments) * grouped_effective)
        coefficients = pauli_coefficients(grouped_effective)
        phase_per_group = float(duration * len(physical_segments) * abs(coefficients["XXX"]))
        cycles_to_pi_over_4 = float(np.pi / 4.0 / phase_per_group)
        leakage = float(np.linalg.norm(evolved - projected, ord=2) ** 2)
        target_signal = phase_aligned_distance(target_piecewise, np.eye(target_piecewise.shape[0], dtype=complex))
        realized_signal = phase_aligned_distance(logical, np.eye(logical.shape[0], dtype=complex))
        target_error = phase_aligned_distance(logical, target_piecewise)
        rows.append({
            "coupling_over_detuning": ratio,
            "segment_duration": duration,
            "segment_count": len(physical_segments),
            "total_duration": duration * len(physical_segments),
            "grouped_xxx_coefficient": coefficients["XXX"],
            "maximum_grouped_unwanted_non_scalar_coefficient": float(max(
                abs(value) for label, value in coefficients.items() if label not in {"III", "XXX"}
            )),
            "low_frame_leakage_worst": leakage,
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
            "polar_distance_to_piecewise_sw": target_error,
            "polar_distance_to_group_average_sw": phase_aligned_distance(logical, target_average),
            "piecewise_sw_signal_distance_from_identity": target_signal,
            "physical_polar_signal_distance_from_identity": realized_signal,
            "relative_piecewise_sw_error": float(target_error / target_signal) if target_signal > 1e-15 else None,
            "xxx_phase_magnitude_per_group": phase_per_group,
            "ideal_groups_to_pi_over_4": cycles_to_pi_over_4,
            "linear_leakage_budget_estimate_to_pi_over_4": float(cycles_to_pi_over_4 * leakage),
        })
    deep = rows
    leakage_power = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in deep]),
        np.log([row["low_frame_leakage_worst"] for row in deep]),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.shared-matter-c3-parity-echo-pulse-audit.v1",
        "parameters": {
            "sign_group": [list(signs) for signs in PARITY_EVEN_SIGNS],
            "schedule": "for each sign word in listed order: microscopic A then B, equal virtual-Rabi durations",
            "pair_detuning": PAIR_DETUNING,
            "fixed_charge_dimension": 1488,
        },
        "rows": rows,
        "deep_sw_group_pulse_leakage_power": leakage_power,
        "decision": "The integer-Rabi eight-segment group is rejected as a C3 gate: static parity averaging cancels its effective-Hamiltonian parasites, but the exact physical group returns an almost scalar polar operation and has relative error approximately one to the intended non-scalar SW signal.",
        "claim_boundary": "This rejects the registered integer-Rabi realization of the C3 parity group, not the static group algebra, an off-resonant dressed/adiabatic implementation, a different control sequence, a four-link stabilizer, a protected code, a 2D phase, defects, fusion, a non-Abelian braid, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "phase8b_shared_matter_c3_parity_echo_pulse_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
