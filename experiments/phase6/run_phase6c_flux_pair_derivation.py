"""Phase 6C: local microscopic pair-transfer derivation preflight.

This is a mechanism audit, not a many-body protected-code or braid audit.
It tests whether a two-mediator pi-flux block cancels single-particle transfer
while retaining a nonzero, fourth-order low-pair splitting.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.native_fusion import (
    build_flux_pair_mediator_block,
    low_pair_masks,
    one_particle_schur_cross_norm,
)


PARAMETERS = {"Delta": 5.0, "U_mediator": 4.0, "E_bind": 2.0}
T_VALUES = (0.05, 0.075, 0.10, 0.125, 0.15, 0.20)
PHASES = (("zero_flux", 0.0), ("pi_flux", float(np.pi)))


def pair_sector_row(t: float, phi: float) -> dict:
    H, states, index = build_flux_pair_mediator_block(N=2, t=t, phi=phi, **PARAMETERS)
    values, vectors = eigh(H, subset_by_index=[0, 2], driver="evr")
    pair_a, pair_b = low_pair_masks()
    captures = [
        float(abs(vectors[index[pair_a], column]) ** 2 + abs(vectors[index[pair_b], column]) ** 2)
        for column in range(2)
    ]
    return {
        "t": t,
        "lowest_energies": values.tolist(),
        "pair_splitting": float(values[1] - values[0]),
        "pair_effective_offdiagonal_magnitude": float((values[1] - values[0]) / 2.0),
        "pair_isolation_gap": float(values[2] - values[1]),
        "pair_subspace_captures": captures,
        "minimum_pair_subspace_capture": float(min(captures)),
    }


def log_power_fit(rows: list[dict]) -> dict:
    x = np.log(np.asarray([row["t"] for row in rows], dtype=float))
    y = np.log(np.asarray([row["pair_splitting"] for row in rows], dtype=float))
    exponent, intercept = np.polyfit(x, y, 1)
    fitted = intercept + exponent * x
    residual = float(np.sum((y - fitted) ** 2))
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "power": float(exponent),
        "prefactor": float(np.exp(intercept)),
        "r_squared_log": float(1.0 - residual / total) if total else 1.0,
    }


def n1_cross_norm(t: float, phi: float) -> float:
    H, states, _ = build_flux_pair_mediator_block(N=1, t=t, phi=phi, **PARAMETERS)
    return one_particle_schur_cross_norm(H, states)


def interaction_off_control() -> list[dict]:
    """At pi flux, remove only the mediator interaction from the mechanism."""
    rows = []
    for t in T_VALUES:
        H, _, _ = build_flux_pair_mediator_block(
            N=2, t=t, phi=float(np.pi), Delta=PARAMETERS["Delta"],
            U_mediator=0.0, E_bind=PARAMETERS["E_bind"],
        )
        values = eigh(H, eigvals_only=True, subset_by_index=[0, 2], driver="evr")
        rows.append({
            "t": t,
            "pair_splitting": float(values[1] - values[0]),
            "pair_isolation_gap": float(values[2] - values[1]),
        })
    return rows


def main() -> None:
    phase_rows: dict[str, dict] = {}
    for label, phi in PHASES:
        rows = [pair_sector_row(t, phi) for t in T_VALUES]
        phase_rows[label] = {
            "phi": phi,
            "one_particle_schur_cross_norm_at_t_0p1": n1_cross_norm(0.1, phi),
            "two_particle_rows": rows,
            "pair_splitting_power_fit": log_power_fit(rows),
        }
    pi_rows = phase_rows["pi_flux"]["two_particle_rows"]
    interaction_control = interaction_off_control()
    qualifies = (
        phase_rows["pi_flux"]["one_particle_schur_cross_norm_at_t_0p1"] < 1e-12
        and min(row["pair_splitting"] for row in pi_rows) > 1e-9
        and min(row["pair_isolation_gap"] for row in pi_rows) > 0.1
        and min(row["minimum_pair_subspace_capture"] for row in pi_rows) > 0.99
        and 3.5 < phase_rows["pi_flux"]["pair_splitting_power_fit"]["power"] < 4.5
        and max(row["pair_splitting"] for row in interaction_control) < 1e-12
    )
    out = {
        "schema": "antler.phase6.flux-pair-mediator-local-preflight.v1",
        "model": {
            "name": "two-rail detuned-mediator pi-flux interferometer",
            "modes": ["a0", "a1", "b0", "b1", "p0", "p1", "m0", "m1"],
            "parameters": PARAMETERS,
            "relative_branch_parity": "Each retained pair transfer changes N_a and N_b by two, preserving (-1)^N_a and (-1)^N_b.",
        },
        "mechanism_test": {
            "single_particle_requirement": "pi flux cancels the exact N=1 zero-energy Schur cross block below 1e-12",
            "two_particle_requirement": "a nonzero, isolated low-pair splitting has a fourth-order hopping power",
            "qualification": bool(qualifies),
        },
        "results_by_flux": phase_rows,
        "interaction_off_control": {
            "description": "same pi-flux block with U_mediator=0; the pair splitting must vanish to numerical precision",
            "rows": interaction_control,
            "maximum_pair_splitting": max(row["pair_splitting"] for row in interaction_control),
        },
        "claim_boundary": (
            "This validates or rejects a local microscopic mechanism only. It does not construct a full correlated-Jordan--Wigner ANTLER lattice, "
            "a locally indistinguishable code, a topological phase, an adiabatic exchange, or a non-Abelian braid. "
            "The original frozen ladder remains unchanged."
        ),
        "decision": (
            "local mechanism qualifies: derive and pre-register a tiled full-ladder embedding before any code or braid calculation"
            if qualifies else
            "local mechanism fails its selectivity/perturbative criteria; do not embed it as a native ANTLER candidate"
        ),
    }
    path = ROOT / "results" / "phase6" / "flux_pair_mediator_local_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    summary = {
        "qualification": qualifies,
        "decision": out["decision"],
        "pi_flux_n1_cross_norm": phase_rows["pi_flux"]["one_particle_schur_cross_norm_at_t_0p1"],
        "pi_flux_pair_power": phase_rows["pi_flux"]["pair_splitting_power_fit"],
        "pi_flux_min_capture": min(row["minimum_pair_subspace_capture"] for row in pi_rows),
        "interaction_off_max_pair_splitting": max(row["pair_splitting"] for row in interaction_control),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
