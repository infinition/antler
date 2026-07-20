"""Direct Peierls-sign implementation of the Phase 7D leg-hopping echo.

Instead of realizing H_leg -> -H_leg by instantaneous staggered onsite phases,
this control assumes a synchronous pi Peierls phase on all leg hoppings while
leaving the charge-two pair conversions unchanged.  It tests the exact same
full 472-state pulse schedule and records the required switching overhead.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse.linalg import norm as sparse_norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import (
    FRAME, L, LINKS, LOGICAL_PA, LOGICAL_PB, PULSE_TIME, evolve, leakage,
    pair_gate, projected, pulse_hamiltonian, rail_rotation, remove_global_phase,
)
from run_phase7d_staggered_echo_refocusing_audit import TOGGLE_MATRIX


LEG_HOPPING, EPSILON = 1.0, 1e-2
SUBCYCLES = (16, 32)


def direct_echo_pulse(vectors: np.ndarray, active_links: tuple[int, ...], kind_prefix: str,
                      subcycles: int) -> np.ndarray:
    h_plus = pulse_hamiltonian(active_links, kind_prefix, EPSILON, LEG_HOPPING)
    h_minus = pulse_hamiltonian(active_links, kind_prefix, EPSILON, -LEG_HOPPING)
    delta = PULSE_TIME / subcycles
    for _ in range(subcycles):
        vectors = evolve(vectors, h_plus, delta / 4.0)
        vectors = evolve(vectors, h_minus, delta / 2.0)
        vectors = evolve(vectors, h_plus, delta / 4.0)
    return vectors


def direct_pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], subcycles: int) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry, rx = rail_rotation("y", rungs), rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = direct_echo_pulse(vectors, active_links, "same", subcycles)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = direct_echo_pulse(vectors, active_links, "opposite", subcycles)
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
    identity_checks = []
    for kind in ("same", "opposite"):
        h_plus = pulse_hamiltonian((1,), kind, EPSILON, LEG_HOPPING)
        h_minus = pulse_hamiltonian((1,), kind, EPSILON, -LEG_HOPPING)
        identity_checks.append({
            "channel_kind": kind,
            "Q_Hplus_Q_minus_direct_pi_Peierls_hminus_frobenius": float(
                sparse_norm(TOGGLE_MATRIX @ h_plus @ TOGGLE_MATRIX - h_minus)
            ),
        })
    rows = []
    for subcycles in SUBCYCLES:
        even = direct_pair_gate(FRAME.copy(), (0, 2), subcycles)
        complete = direct_pair_gate(even, (1,), subcycles)
        rows.append({
            "subcycles_per_mediator_pulse": subcycles,
            "peierls_pi_sign_switches_per_mediator_pulse": 2 * subcycles,
            "peierls_pi_sign_switches_full_four_pulse_schedule": 8 * subcycles,
            **metrics(complete, reference),
        })
    out = {
        "schema": "antler.phase7d.direct-peierls-echo-audit.v1",
        "control_contract": (
            "Apply phi_leg=0 or pi synchronously to every intraleg hopping while pair-conversion amplitudes are unchanged. "
            "This is a new dynamic hardware resource; static local pi-flux cancellation does not by itself establish its availability."
        ),
        "parameters": {"leg_hopping": LEG_HOPPING, "inactive_channel_coupling_over_g": EPSILON, "subcycles": list(SUBCYCLES)},
        "identity_checks": identity_checks,
        "rows": rows,
        "decision": (
            "Direct Peierls-sign modulation exactly realizes the ideal echo Hamiltonian in the registered microscopic model. It is a "
            "conditional control bridge pending a derivation of synchronous phase-switch bandwidth, errors and effects on pair conversion."
        ),
        "claim_boundary": (
            "This establishes no experimental implementation, protected phase, edge mode, 2D code, braid, non-Abelian statistics, "
            "universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "direct_peierls_echo_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
