"""Exact staggered-phase echo test for the Phase 7D intraleg-hopping no-go.

For nearest-neighbour pair conversion, Q=(-1)^(sum_j j n_j + N_mediator)
leaves every charge-two conversion invariant and flips the sign of intraleg
hopping.  A symmetric H_plus/H_minus/H_plus sequence is therefore a
pre-registered first attempt to refocus the hopping that spoiled the bare
closed Rabi pulse.  It is a control proposal, not a hardware claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import norm as sparse_norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import (
    FRAME, L, LINKS, LOGICAL_PA, LOGICAL_PB, LOW_MODES, PULSE_TIME, STATES,
    evolve, leakage, pair_gate, projected, pulse_hamiltonian,
    rail_rotation, remove_global_phase,
)


LEG_HOPPINGS = (0.1, 0.3, 1.0)
SUBCYCLES = (1, 2, 4, 8, 16)
COMBINED_EPSILONS = (1e-3, 1e-2)
COMBINED_SUBCYCLES = (8, 16, 32)


def staggered_toggle() -> np.ndarray:
    """Q: staggered rail phase times a pi phase on every charge-two mediator."""
    values = []
    for raw_state in STATES:
        state = int(raw_state)
        rail_exponent = sum(
            rung * (((state >> (2 * rung)) & 1) + ((state >> (2 * rung + 1)) & 1))
            for rung in range(L)
        )
        mediator_exponent = (state >> LOW_MODES).bit_count()
        values.append(-1.0 if (rail_exponent + mediator_exponent) & 1 else 1.0)
    return np.asarray(values, dtype=complex)


TOGGLE = staggered_toggle()
TOGGLE_MATRIX = diags(TOGGLE, format="csr")


def refocused_pulse(vectors: np.ndarray, active_links: tuple[int, ...], kind_prefix: str,
                     leg_hopping: float, subcycles: int, epsilon: float = 0.0) -> np.ndarray:
    h_plus = pulse_hamiltonian(active_links, kind_prefix, epsilon, leg_hopping)
    h_minus = TOGGLE_MATRIX @ h_plus @ TOGGLE_MATRIX
    delta = PULSE_TIME / subcycles
    for _ in range(subcycles):
        vectors = evolve(vectors, h_plus, delta / 4.0)
        vectors = evolve(vectors, h_minus, delta / 2.0)
        vectors = evolve(vectors, h_plus, delta / 4.0)
    return vectors


def refocused_pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], leg_hopping: float,
                         subcycles: int, epsilon: float = 0.0) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry, rx = rail_rotation("y", rungs), rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = refocused_pulse(vectors, active_links, "same", leg_hopping, subcycles, epsilon)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = refocused_pulse(vectors, active_links, "opposite", leg_hopping, subcycles, epsilon)
    return evolve(vectors, rx, -np.pi / 4.0)


def schedule_refocused(leg_hopping: float, subcycles: int, epsilon: float = 0.0) -> np.ndarray:
    even = refocused_pair_gate(FRAME.copy(), (0, 2), leg_hopping, subcycles, epsilon)
    return refocused_pair_gate(even, (1,), leg_hopping, subcycles, epsilon)


def logical_metrics(vectors: np.ndarray, reference: np.ndarray) -> dict:
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


def fit_power(rows: list[dict], field: str) -> dict:
    selected = [row for row in rows if row["method"] == "staggered_symmetric_echo" and row["subcycles"] >= 2]
    x = np.log(np.asarray([row["subcycles"] for row in selected], dtype=float))
    y = np.log(np.asarray([row[field] for row in selected], dtype=float))
    power, intercept = np.polyfit(x, y, 1)
    fitted = power * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "field": field,
        "fit_subcycles": [row["subcycles"] for row in selected],
        "power": float(power),
        "prefactor": float(np.exp(intercept)),
        "r_squared_log": float(1.0 - np.sum((y - fitted) ** 2) / total),
    }


def main() -> None:
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    toggling_checks = []
    for kind in ("same", "opposite"):
        h_plus = pulse_hamiltonian((1,), kind, 0.0, 1.0)
        h_minus = TOGGLE_MATRIX @ h_plus @ TOGGLE_MATRIX
        expected = pulse_hamiltonian((1,), kind, 0.0, -1.0)
        toggling_checks.append({
            "channel_kind": kind,
            "QHplusQ_minus_Hmed_minus_Hleg_frobenius": float(sparse_norm(h_minus - expected)),
        })
    rows = []
    for leg_hopping in LEG_HOPPINGS:
        bare = pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0, leg_hopping), (1,), 0.0, leg_hopping)
        rows.append({"leg_hopping_during_pulses": leg_hopping, "method": "bare", "subcycles": 0,
                     **logical_metrics(bare, reference)})
        for subcycles in SUBCYCLES:
            echoed = schedule_refocused(leg_hopping, subcycles)
            rows.append({"leg_hopping_during_pulses": leg_hopping, "method": "staggered_symmetric_echo",
                         "subcycles": subcycles, **logical_metrics(echoed, reference)})
    power_fits = []
    for leg_hopping in LEG_HOPPINGS:
        local_rows = [row for row in rows if row["leg_hopping_during_pulses"] == leg_hopping]
        power_fits.append({
            "leg_hopping_during_pulses": leg_hopping,
            "monomer_leakage": fit_power(local_rows, "monomer_leakage"),
            "logical_deviation": fit_power(local_rows, "logical_deviation_from_zero_leg_schedule"),
        })
    combined_rows = []
    for epsilon in COMBINED_EPSILONS:
        for subcycles in COMBINED_SUBCYCLES:
            combined = schedule_refocused(1.0, subcycles, epsilon)
            combined_rows.append({
                "leg_hopping_during_pulses": 1.0,
                "inactive_channel_coupling_over_g": epsilon,
                "subcycles": subcycles,
                **logical_metrics(combined, reference),
            })
    out = {
        "schema": "antler.phase7d.staggered-echo-refocusing-audit.v1",
        "identity": (
            "Q=(-1)^(sum_j j n_j + N_mediator): for every nearest-neighbour charge-two conversion, "
            "Q H_mediator Q=H_mediator and Q H_leg Q=-H_leg"
        ),
        "model": {
            "same_472_state_four_rung_block": True,
            "refocused_segment": "exp(-i Hplus delta/4) exp(-i Hminus delta/2) exp(-i Hplus delta/4)",
            "Hminus_implementation": "staggered rail pi phase plus mediator pi phase",
            "leg_hopping_scan": list(LEG_HOPPINGS), "subcycle_scan": list(SUBCYCLES),
            "combined_crosstalk_scan_over_g": list(COMBINED_EPSILONS),
        },
        "toggling_identity_checks": toggling_checks,
        "rows": rows,
        "power_fits": power_fits,
        "combined_leg_hopping_crosstalk_rows": combined_rows,
        "decision": (
            "This exact audit determines whether a specific control echo repairs the bare leg-hopping failure. Passing a finite-block "
            "echo would qualify a proposed control resource only; it would still require pulse-error, many-link, and phase-protection audits."
        ),
        "claim_boundary": (
            "The staggered rail and mediator pi phases are additional ideal controls. This is not a derivation of their experimental "
            "availability, nor a topological-phase, 2D-code, braid, non-Abelian, universal or fault-tolerance result."
        ),
    }
    path = ROOT / "results" / "phase7" / "staggered_echo_refocusing_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
