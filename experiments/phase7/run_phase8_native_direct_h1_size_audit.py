"""Size/filling audit of the direct switched-mediator Phase-8 bridge.

The direct H1 construction removes the *ideal* bare rail pulse, but it does
not by itself prove that virtual mediator occupation is controlled beyond one
three-rung block.  This exact-diagonalization audit keeps the SW ratio and
the analytic virtual-Rabi closure fixed while changing both system size and
particle number.  It is deliberately a control audit, not a phase diagnosis.
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
from run_phase8_native_micro_floquet_l3 import ALPHA, RATIO, TARGET_U0, phase_aligned_distance, polar_unitary
from run_phase8_native_direct_h1_closure import build_micro, rotated_channels


SYSTEMS = ((3, 2), (4, 2), (4, 4), (6, 3))
CYCLES = (1, 4)


def audit_one(length: int, particle_number: int, channels_h0: np.ndarray, channels_h1: np.ndarray,
              g: float, detuning: float, period: float) -> dict:
    h0_micro, states, frame = build_micro(length, particle_number, channels_h0, g, detuning)
    h1_micro, states_h1, frame_h1 = build_micro(length, particle_number, channels_h1, g, detuning)
    if not np.array_equal(states, states_h1) or not np.allclose(frame, frame_h1):
        raise RuntimeError("the two direct microscopic segments use different bases")
    target_h0, target_p, target_states, _ = build_h0_and_rotation(length, particle_number, TARGET_U0)
    low_states = np.asarray([
        states[int(np.argmax(np.abs(frame[:, column])))] for column in range(frame.shape[1])
    ], dtype=np.int64)
    if not np.array_equal(low_states, target_states):
        raise RuntimeError(f"low-frame / target-basis mismatch at L={length}, N={particle_number}")
    direct_cycle = expm(-1j * (1.0 - ALPHA) * period * h1_micro) @ expm(-1j * ALPHA * period * h0_micro)
    target_h1 = target_p.conj().T @ target_h0 @ target_p
    target_cycle = expm(-1j * period * (ALPHA * target_h0 + (1.0 - ALPHA) * target_h1))
    projector = frame @ frame.conj().T
    identity = np.eye(h0_micro.shape[0], dtype=complex)
    rows = []
    for count in CYCLES:
        unitary = np.linalg.matrix_power(direct_cycle, count)
        raw = frame.conj().T @ unitary @ frame
        logical = polar_unitary(raw)
        rows.append({
            "cycles": count,
            "low_frame_leakage_worst": float(np.linalg.norm((identity - projector) @ unitary @ frame, ord=2) ** 2),
            "polar_logical_vs_target_distance": phase_aligned_distance(logical, np.linalg.matrix_power(target_cycle, count)),
            "raw_logical_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
        })
    return {
        "L": length,
        "N": particle_number,
        "filling": particle_number / (2.0 * length),
        "microscopic_dimension": int(h0_micro.shape[0]),
        "bare_low_frame_dimension": int(frame.shape[1]),
        "rows": rows,
    }


def main() -> None:
    channels_h1, factorization = rotated_channels()
    channels_h0 = np.asarray(((1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0)), dtype=complex)
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
    period = 4.0 * np.pi / omega
    rows = [audit_one(length, number, channels_h0, channels_h1, g, detuning, period) for length, number in SYSTEMS]
    reference = rows[0]["rows"][0]["low_frame_leakage_worst"]
    for row in rows:
        row["one_cycle_leakage_over_l3n2"] = float(row["rows"][0]["low_frame_leakage_worst"] / reference)
    out = {
        "schema": "antler.phase8.native-direct-h1-size-audit.v1",
        "parameters": {
            "systems": [{"L": length, "N": number} for length, number in SYSTEMS],
            "target_u0": TARGET_U0,
            "g_over_detuning": RATIO,
            "detuning": detuning,
            "g": g,
            "alpha": ALPHA,
            "virtual_rabi_omega": omega,
            "closure_period": period,
            "direct_segments": "H0 and P^dag H0 P from coherent pair-mediator channels",
        },
        "rotated_channel_factorization_frobenius_residual": factorization["factorization_frobenius_residual"],
        "rows": rows,
        "decision": "Exact finite-size control of the direct switched-mediator stroboscopic block at fixed SW ratio.",
        "claim_boundary": "Finite L=3,4 controls test only local microscopic leakage and target agreement. They do not show asymptotic size convergence, a gapped topological phase, a protected logical qubit, physical pulse bandwidth tolerance, braiding, non-Abelian statistics, universality, or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "native_direct_h1_size_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
