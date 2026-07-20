"""Finite-pulse stroboscopic audit of the published Floquet-ladder bridge.

The static effective Hamiltonian is useful only if it is reproduced by the
actual pulse sequence P^dag exp[-i(1-alpha) T H0] P exp[-i alpha T H0].
This exact finite-Hilbert-space control measures the one-period Trotter error
and the branch-parity commutator, including a deliberately off-angle pulse.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm, norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import site_index
from antler.number_conserving_pairwire import _apply, wire_a_parity
from run_phase7d_floquet_full_ladder_preflight import build_h0_and_rotation


LENGTH, PARTICLE_NUMBER, U0, ALPHA = 6, 4, -1.5, 0.5
PERIODS = (0.4, 0.2, 0.1, 0.05)
IDEAL_ETA = np.pi / 2.0
ANGLE_CASES = (IDEAL_ETA, IDEAL_ETA + 0.1)


def analyze(period: float, eta: float, h0: np.ndarray, jx: np.ndarray, parity: np.ndarray) -> dict:
    rotation = expm(-1j * eta * jx)
    h_eff = ALPHA * h0 + (1.0 - ALPHA) * (rotation.conj().T @ h0 @ rotation)
    cycle = (
        rotation.conj().T
        @ expm(-1j * (1.0 - ALPHA) * period * h0)
        @ rotation
        @ expm(-1j * ALPHA * period * h0)
    )
    eff = expm(-1j * period * h_eff)
    scale = np.sqrt(h0.shape[0])
    return {
        "period": period,
        "eta": eta,
        "one_period_unitary_distance_normalized": float(norm(cycle - eff) / scale),
        "stroboscopic_branch_parity_commutator_normalized": float(norm(cycle @ parity - parity @ cycle) / scale),
        "effective_branch_parity_commutator_normalized": float(norm(h_eff @ parity - parity @ h_eff) / scale),
    }


def main() -> None:
    h0, _, states, _ = build_h0_and_rotation(LENGTH, PARTICLE_NUMBER, U0)
    index = {int(state): position for position, state in enumerate(states)}
    jx = np.zeros_like(h0)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for rung in range(LENGTH):
            a, b = site_index(rung, 0), site_index(rung, 1)
            for operations in ((("ann", b), ("create", a)), (("ann", a), ("create", b))):
                item = _apply(state, operations)
                if item is not None:
                    new_state, amplitude = item
                    jx[index[new_state], column] += 0.5 * amplitude
    if not np.allclose(jx, jx.conj().T, atol=1e-12):
        raise RuntimeError("Jx construction is not Hermitian")
    parity = np.diag([(-1.0 if wire_a_parity(int(state), LENGTH) else 1.0) for state in states])
    rows = [analyze(period, eta, h0, jx, parity) for eta in ANGLE_CASES for period in PERIODS]
    ideal = [row for row in rows if abs(row["eta"] - IDEAL_ETA) < 1e-12]
    fit_slope = float(np.polyfit(np.log([row["period"] for row in ideal]), np.log([row["one_period_unitary_distance_normalized"] for row in ideal]), 1)[0])
    out = {
        "schema": "antler.phase8.finite-pulse-stroboscopic-audit.v1",
        "citation": "Defossez et al., arXiv:2412.14886v2 (2025)",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "u0_attractive_nn": U0, "alpha": ALPHA},
        "rows": rows,
        "ideal_eta_unitary_error_loglog_slope": fit_slope,
        "decision": "Finite-block validation of the actual pulsed external protocol against its high-frequency effective Hamiltonian.",
        "claim_boundary": "This validates neither an ANTLER-native pulse resource nor an asymptotic topological phase, braid, non-Abelian operation, or hardware gate.",
    }
    path = ROOT / "results" / "phase7" / "finite_pulse_stroboscopic_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
