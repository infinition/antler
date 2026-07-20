"""Narrow finite-frequency amplitude optimization for the Phase 7D CDT control.

The first Bessel zero is asymptotically optimal.  At finite drive frequency,
Magnus corrections can shift the optimum.  This registered one-dimensional
scan varies only xi=A/omega near the first zero at fixed four-cycle pulse,
fixed full 472-state dynamics and fixed error model.
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

from run_phase7d_four_rung_microscopic_pulse_audit import (
    FRAME, L, LINKS, LOGICAL_PA, LOGICAL_PB, PULSE_TIME, evolve, leakage,
    pair_gate, projected, pulse_hamiltonian, rail_rotation, remove_global_phase,
)
from run_phase7d_continuous_cdt_refocusing_audit import K_VALUES


LEG_HOPPING, EPSILON = 1.0, 1e-2
CYCLES, STEPS_PER_CYCLE = 4, 16
XI_VALUES = (2.10, 2.20, 2.30, 2.35, 2.4048255576957724, 2.46, 2.55, 2.65, 2.75)


def precompute(active_links: tuple[int, ...], kind_prefix: str) -> tuple[np.ndarray, np.ndarray]:
    h_total = pulse_hamiltonian(active_links, kind_prefix, EPSILON, LEG_HOPPING)
    h_pair = pulse_hamiltonian(active_links, kind_prefix, EPSILON, 0.0)
    h_leg = h_total - h_pair
    dt = PULSE_TIME / (CYCLES * STEPS_PER_CYCLE)
    return expm((-0.5j * dt) * h_pair.toarray()), expm((-1j * dt) * h_leg.toarray())


PROPAGATORS = {
    (active, kind): precompute(active, kind)
    for active in ((0, 2), (1,)) for kind in ("same", "opposite")
}


def cdt_pulse(vectors: np.ndarray, active_links: tuple[int, ...], kind_prefix: str, xi: float) -> np.ndarray:
    u_pair_half, u_leg = PROPAGATORS[(active_links, kind_prefix)]
    dt = PULSE_TIME / (CYCLES * STEPS_PER_CYCLE)
    omega = 2.0 * np.pi * CYCLES / PULSE_TIME
    for step in range(CYCLES * STEPS_PER_CYCLE):
        phase = xi * np.sin(omega * (step + 0.5) * dt)
        rotating = np.exp(-1j * phase * K_VALUES)
        vectors = u_pair_half @ vectors
        vectors = rotating[:, None] * vectors
        vectors = u_leg @ vectors
        vectors = rotating.conj()[:, None] * vectors
        vectors = u_pair_half @ vectors
    return vectors


def cdt_pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], xi: float) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry, rx = rail_rotation("y", rungs), rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = cdt_pulse(vectors, active_links, "same", xi)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = cdt_pulse(vectors, active_links, "opposite", xi)
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
    omega = 2.0 * np.pi * CYCLES / PULSE_TIME
    rows = []
    for xi in XI_VALUES:
        even = cdt_pair_gate(FRAME.copy(), (0, 2), xi)
        complete = cdt_pair_gate(even, (1,), xi)
        rows.append({
            "xi": xi,
            "drive_frequency": omega,
            "base_potential_amplitude": xi * omega,
            "max_rail_onsite_amplitude": (L - 1) * xi * omega,
            **metrics(complete, reference),
        })
    best = min(rows, key=lambda row: (row["monomer_leakage"], row["logical_parity_a_residual"]))
    out = {
        "schema": "antler.phase7d.cdt-finite-frequency-optimization.v1",
        "registered_window": {
            "cycles_per_rabi_pulse": CYCLES, "steps_per_cycle": STEPS_PER_CYCLE,
            "xi_values": list(XI_VALUES), "center": "first J0 zero", "leg_hopping": LEG_HOPPING,
            "inactive_channel_coupling_over_g": EPSILON,
        },
        "rows": rows,
        "best_by_leakage_then_parity": best,
        "decision": (
            "A finite-frequency local amplitude optimization only. Any improved point remains a 472-state pulse-control result and "
            "must pass a separate temporal-resolution and hardware-bandwidth gate before it is used as a control primitive."
        ),
        "claim_boundary": (
            "This does not alter the absence of a demonstrated protected many-body phase, edge mode, 2D code, braid, non-Abelian "
            "statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "cdt_finite_frequency_optimization.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
