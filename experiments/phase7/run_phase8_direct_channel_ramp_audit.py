"""Finite-bandwidth audit for the direct switched-mediator Phase-8 segment.

The two rank-two pair-channel subspaces have a particularly simple continuous
path.  In the pair basis (aa, ab, ba, bb), retain
u1=(aa-bb)/sqrt(2) and rotate u0=(aa+bb)/sqrt(2) into
v1=(ab+ba)/sqrt(2).  Every point therefore uses two normalized, orthogonal
charge-two channels with the same coupling g and detuning Delta.

This tests an explicitly declared control interpolation.  It is not yet a
derivation of hardware waveform bandwidth or calibration noise.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation
from run_phase8_native_micro_floquet_l3 import ALPHA, TARGET_U0, phase_aligned_distance, polar_unitary
from run_phase8_native_direct_h1_closure import build_micro, rotated_channels


LENGTH, PARTICLE_NUMBER = 3, 2
RATIO = 0.0125
LOGICAL_CYCLES = 16
RAMP_FRACTIONS = (0.0, 0.001, 0.003, 0.01, 0.03, 0.05)
SHAPES = ("linear", "sin2")
SEGMENTS_PER_RAMP = 16
SEGMENT_CONVERGENCE = (8, 16, 32, 64)


def channel_path(phi: float) -> np.ndarray:
    """Two orthonormal mediator conversion rows for phi in [0, pi/2]."""
    u0 = np.asarray((1.0, 0.0, 0.0, 1.0), dtype=complex) / np.sqrt(2.0)
    u1 = np.asarray((1.0, 0.0, 0.0, -1.0), dtype=complex) / np.sqrt(2.0)
    v1 = np.asarray((0.0, 1.0, 1.0, 0.0), dtype=complex) / np.sqrt(2.0)
    rows = np.asarray((u1, np.cos(phi) * u0 + np.sin(phi) * v1), dtype=complex)
    if not np.allclose(rows @ rows.conj().T, np.eye(2), atol=1e-12):
        raise RuntimeError("non-orthonormal switched pair channels")
    return rows


def shape_value(shape: str, position: float) -> float:
    if shape == "linear":
        return position
    if shape == "sin2":
        return float(np.sin(np.pi * position / 2.0) ** 2)
    raise ValueError(shape)


def append_evolution(unitary: np.ndarray, hamiltonian: np.ndarray, duration: float) -> np.ndarray:
    return expm(-1j * duration * hamiltonian) @ unitary


def cycle_with_ramps(h0: np.ndarray, h1: np.ndarray, ramp_hamiltonians: list[np.ndarray], period: float,
                     ramp_fraction: float, shape: str) -> np.ndarray:
    """Symmetric H0 -> H1 -> H0 cycle with two finite interpolation ramps."""
    ramp_duration = ramp_fraction * period
    dwell = (period - 2.0 * ramp_duration) / 2.0
    if dwell < -1e-15:
        raise ValueError("ramp duration exceeds the available Floquet period")
    unitary = np.eye(h0.shape[0], dtype=complex)
    unitary = append_evolution(unitary, h0, max(dwell, 0.0))
    if ramp_fraction:
        step = ramp_duration / len(ramp_hamiltonians)
        for position, hamiltonian in zip(
            ((index + 0.5) / len(ramp_hamiltonians) for index in range(len(ramp_hamiltonians))),
            ramp_hamiltonians,
        ):
            # The prepared list is ordered H0->H1; position is consumed only
            # to make the chronological direction explicit in this audit.
            del position
            unitary = append_evolution(unitary, hamiltonian, step)
    unitary = append_evolution(unitary, h1, max(dwell, 0.0))
    if ramp_fraction:
        step = ramp_duration / len(ramp_hamiltonians)
        for hamiltonian in reversed(ramp_hamiltonians):
            unitary = append_evolution(unitary, hamiltonian, step)
    return unitary


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 4.0 * np.pi / omega
    h0_channels = channel_path(0.0)
    h1_channels = channel_path(np.pi / 2.0)
    _, factorization = rotated_channels()
    h0_micro, states, frame = build_micro(LENGTH, PARTICLE_NUMBER, h0_channels, g, detuning)
    h1_micro, states_h1, frame_h1 = build_micro(LENGTH, PARTICLE_NUMBER, h1_channels, g, detuning)
    if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
        raise RuntimeError("endpoint channel frames disagree")
    target_h0, target_p, target_states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, TARGET_U0)
    low_states = np.asarray([states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError("bare low frame and target basis disagree")
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    target = expm(-1j * LOGICAL_CYCLES * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    projector = frame @ frame.conj().T
    identity = np.eye(h0_micro.shape[0], dtype=complex)
    rows = []
    endpoint_residuals = {
        "h0_channel_projector_vs_separate_rail_projector": float(np.linalg.norm(
            h0_channels.conj().T @ h0_channels
            - np.diag((1.0, 0.0, 0.0, 1.0))
        )),
        "h1_channel_factorization_residual": factorization["factorization_frobenius_residual"],
    }
    for shape in SHAPES:
        ramp_hamiltonians = [
            build_micro(
                LENGTH, PARTICLE_NUMBER,
                channel_path(np.pi * shape_value(shape, (index + 0.5) / SEGMENTS_PER_RAMP) / 2.0),
                g, detuning,
            )[0]
            for index in range(SEGMENTS_PER_RAMP)
        ]
        for fraction in RAMP_FRACTIONS:
            unitary = np.linalg.matrix_power(
                cycle_with_ramps(h0_micro, h1_micro, ramp_hamiltonians, period, fraction, shape),
                LOGICAL_CYCLES,
            )
            raw = frame.conj().T @ unitary @ frame
            logical = polar_unitary(raw)
            rows.append({
                "shape": shape,
                "ramp_fraction_of_period_per_transition": fraction,
                "ramp_segments_per_transition": SEGMENTS_PER_RAMP,
                "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
                "polar_logical_vs_instantaneous_target_distance": phase_aligned_distance(logical, target),
                "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
            })
    convergence_rows = []
    for segments in SEGMENT_CONVERGENCE:
        ramp_hamiltonians = [
            build_micro(
                LENGTH, PARTICLE_NUMBER,
                channel_path(np.pi * shape_value("sin2", (index + 0.5) / segments) / 2.0),
                g, detuning,
            )[0]
            for index in range(segments)
        ]
        unitary = np.linalg.matrix_power(
            cycle_with_ramps(h0_micro, h1_micro, ramp_hamiltonians, period, 0.01, "sin2"),
            LOGICAL_CYCLES,
        )
        raw = frame.conj().T @ unitary @ frame
        convergence_rows.append({
            "shape": "sin2",
            "ramp_fraction_of_period_per_transition": 0.01,
            "ramp_segments_per_transition": segments,
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
            "polar_logical_vs_instantaneous_target_distance": phase_aligned_distance(polar_unitary(raw), target),
        })
    out = {
        "schema": "antler.phase8.direct-channel-ramp-audit.v1",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0,
            "g_over_detuning": RATIO, "detuning": detuning, "g": g,
            "closure_period": period, "logical_cycles": LOGICAL_CYCLES,
            "total_duration": LOGICAL_CYCLES * period,
            "ramp_fractions": list(RAMP_FRACTIONS), "shapes": list(SHAPES),
            "segments_per_transition": SEGMENTS_PER_RAMP,
        },
        "endpoint_channel_checks": endpoint_residuals,
        "rows": rows,
        "ramp_discretization_convergence": {
            "registered_shape": "sin2",
            "registered_ramp_fraction": 0.01,
            "rows": convergence_rows,
            "relative_change_16_to_32": {
                "leakage": float(abs(convergence_rows[1]["low_frame_leakage_worst"] - convergence_rows[2]["low_frame_leakage_worst"]) / convergence_rows[2]["low_frame_leakage_worst"]),
                "distance": float(abs(convergence_rows[1]["polar_logical_vs_instantaneous_target_distance"] - convergence_rows[2]["polar_logical_vs_instantaneous_target_distance"]) / convergence_rows[2]["polar_logical_vs_instantaneous_target_distance"]),
            },
        },
        "strict_local_target": {
            "leakage_below": 1e-4,
            "instantaneous_target_distance_below": 1e-4,
            "passing_rows": [
                {"shape": row["shape"], "ramp_fraction": row["ramp_fraction_of_period_per_transition"]}
                for row in rows
                if row["low_frame_leakage_worst"] < 1e-4
                and row["polar_logical_vs_instantaneous_target_distance"] < 1e-4
            ],
        },
        "decision": "Finite-width direct-channel switching control at a deep-SW Rabi-closed point.",
        "claim_boundary": "The interpolation is an explicitly imposed coherent waveform. This finite L=3 ideal-control test does not derive laboratory bandwidth, amplitude limits, calibration error, crosstalk, a thermodynamic phase, a braid, non-Abelian statistics, universality or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "direct_channel_ramp_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
