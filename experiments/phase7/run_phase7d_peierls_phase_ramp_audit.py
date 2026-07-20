"""Finite Peierls-phase-ramp stress test for the direct Phase 7D echo.

The instantaneous 0 <-> pi hopping-sign switches are replaced by linear phase
ramps with the pair-conversion Hamiltonian kept active.  This is the direct
phase analogue of the failed finite onsite-kick audit, with no large onsite
potential assumed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import norm as sparse_norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import (
    FRAME, INDEX, L, LINKS, LOGICAL_PA, LOGICAL_PB, PULSE_TIME, STATES,
    evolve, hopping, leakage, pair_gate, projected, pulse_hamiltonian,
    rail_rotation, remove_global_phase,
)


LEG_HOPPING, EPSILON, SUBCYCLES = 1.0, 1e-2, 16
RAMP_FRACTIONS = (0.0, 0.005, 0.02)
RAMP_STEPS = 2


def low_mode(rung: int, rail: int) -> int:
    return 2 * rung + rail


def leg_with_phase(phi: float) -> csr_matrix:
    entries: dict[tuple[int, int], complex] = {}
    forward, reverse = np.exp(1j * phi), np.exp(-1j * phi)
    for column, raw_state in enumerate(STATES):
        state = int(raw_state)
        for left_rung, right_rung in LINKS:
            for rail in (0, 1):
                left, right = low_mode(left_rung, rail), low_mode(right_rung, rail)
                for first, second, phase in ((right, left, forward), (left, right, reverse)):
                    item = hopping(state, first, second)
                    if item is not None:
                        final, fermion_sign = item
                        row = INDEX[final]
                        entries[(row, column)] = entries.get((row, column), 0.0) - LEG_HOPPING * phase * fermion_sign
    rows, columns, values = zip(*((row, column, value) for (row, column), value in entries.items()))
    return csr_matrix((values, (rows, columns)), shape=(len(STATES), len(STATES)), dtype=complex)


LEG_ZERO, LEG_PI = leg_with_phase(0.0), leg_with_phase(np.pi)


def ramp_evolve(vectors: np.ndarray, h_pair: csr_matrix, direction: int, duration: float) -> np.ndarray:
    if duration == 0.0:
        return vectors
    for step in range(RAMP_STEPS):
        fraction = (step + 0.5) / RAMP_STEPS
        phi = np.pi * fraction if direction > 0 else np.pi * (1.0 - fraction)
        vectors = evolve(vectors, h_pair + leg_with_phase(phi), duration / RAMP_STEPS)
    return vectors


def ramped_echo_pulse(vectors: np.ndarray, active_links: tuple[int, ...], kind_prefix: str,
                      ramp_fraction: float) -> np.ndarray:
    h_pair = pulse_hamiltonian(active_links, kind_prefix, EPSILON, 0.0)
    h_plus, h_minus = h_pair + LEG_ZERO, h_pair + LEG_PI
    delta = PULSE_TIME / SUBCYCLES
    plus_duration = delta * (0.25 - ramp_fraction / 2.0)
    minus_duration = delta * (0.50 - ramp_fraction)
    ramp_duration = delta * ramp_fraction
    if plus_duration < 0.0 or minus_duration < 0.0:
        raise ValueError("ramp fraction too large")
    for _ in range(SUBCYCLES):
        vectors = evolve(vectors, h_plus, plus_duration)
        vectors = ramp_evolve(vectors, h_pair, +1, ramp_duration)
        vectors = evolve(vectors, h_minus, minus_duration)
        vectors = ramp_evolve(vectors, h_pair, -1, ramp_duration)
        vectors = evolve(vectors, h_plus, plus_duration)
    return vectors


def ramped_pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], ramp_fraction: float) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry, rx = rail_rotation("y", rungs), rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = ramped_echo_pulse(vectors, active_links, "same", ramp_fraction)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = ramped_echo_pulse(vectors, active_links, "opposite", ramp_fraction)
    return evolve(vectors, rx, -np.pi / 4.0)


def metrics(vectors: np.ndarray, reference: np.ndarray) -> dict:
    logical = projected(vectors)
    return {
        "monomer_leakage": leakage(vectors),
        "logical_deviation_from_zero_leg_schedule": float(np.linalg.norm(
            logical - remove_global_phase(reference, logical), ord=2
        )),
        "logical_parity_a_residual": float(np.linalg.norm(logical @ LOGICAL_PA - LOGICAL_PA @ logical, ord=2)),
        "logical_parity_b_residual": float(np.linalg.norm(logical @ LOGICAL_PB - LOGICAL_PB @ logical, ord=2)),
        "logical_singular_value_min": float(np.linalg.svd(logical, compute_uv=False)[-1]),
    }


def main() -> None:
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    checks = []
    for kind in ("same", "opposite"):
        h_pair = pulse_hamiltonian((1,), kind, EPSILON, 0.0)
        checks.append({
            "channel_kind": kind,
            "phase_zero_leg_match_frobenius": float(sparse_norm(h_pair + LEG_ZERO - pulse_hamiltonian((1,), kind, EPSILON, LEG_HOPPING))),
            "phase_pi_leg_match_frobenius": float(sparse_norm(h_pair + LEG_PI - pulse_hamiltonian((1,), kind, EPSILON, -LEG_HOPPING))),
        })
    rows = []
    for fraction in RAMP_FRACTIONS:
        even = ramped_pair_gate(FRAME.copy(), (0, 2), fraction)
        complete = ramped_pair_gate(even, (1,), fraction)
        rows.append({
            "ramp_fraction_of_each_subcycle": fraction,
            "ramp_duration": fraction * PULSE_TIME / SUBCYCLES,
            "ramp_steps": RAMP_STEPS,
            "total_phase_ramp_time_per_mediator_pulse": 2.0 * fraction * PULSE_TIME,
            **metrics(complete, reference),
        })
    out = {
        "schema": "antler.phase7d.peierls-phase-ramp-audit.v1",
        "control_contract": "All intraleg hoppings share a smooth Peierls phase ramp 0->pi->0; charge-two pair conversion remains active.",
        "parameters": {"leg_hopping": LEG_HOPPING, "inactive_channel_coupling_over_g": EPSILON, "subcycles": SUBCYCLES, "ramp_fraction_scan": list(RAMP_FRACTIONS)},
        "endpoint_identity_checks": checks,
        "rows": rows,
        "decision": "Finite Peierls-ramp preflight only; it determines whether the direct sign-echo remains viable once switching is not instantaneous.",
        "claim_boundary": "No switching hardware, calibration noise, many-body phase, 2D code, braid, non-Abelianity, universality or fault tolerance is established.",
    }
    path = ROOT / "results" / "phase7" / "peierls_phase_ramp_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
