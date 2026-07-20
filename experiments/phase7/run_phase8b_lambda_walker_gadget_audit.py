"""Exact low-sector audit of the repaired Phase-8B Lambda-walker gadget.

The declared new primitive is a neutral walker with vacuum plus three
single-occupation internal states.  Its two endpoint conversions and two
parity-conditioned hops form a closed virtual loop.  The loop flux is the
coarse Gauss word G_B=X_L P_1 P_2 X_R.  This script tests that claim exactly
at the walker layer, independently of a microscopic implementation of the
new density-conditioned hopping primitive.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DETUNING = 10.0
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)
LABELS = ("X_L", "p_1", "p_2", "X_R")


def walker_hamiltonian(signs: tuple[int, int, int, int], coupling: float) -> np.ndarray:
    """Vacuum--mu0--mu1--mu2--vacuum ring, at most one walker occupied."""
    x_left, p1, p2, x_right = signs
    hamiltonian = np.diag((0.0, DETUNING, DETUNING, DETUNING)).astype(complex)
    for first, second, value in (
        (0, 1, coupling * x_left),
        (1, 2, coupling * p1),
        (2, 3, coupling * p2),
        (3, 0, coupling * x_right),
    ):
        hamiltonian[first, second] = value
        hamiltonian[second, first] = value
    return hamiltonian


def walsh_coefficients(energies: dict[tuple[int, int, int, int], float]) -> dict[str, float]:
    coefficients = {}
    for mask in range(1 << len(LABELS)):
        label = "I" if mask == 0 else "_".join(LABELS[index] for index in range(len(LABELS)) if mask & (1 << index))
        coefficients[label] = float(sum(
            energy * np.prod([signs[index] for index in range(len(LABELS)) if mask & (1 << index)])
            for signs, energy in energies.items()
        ) / len(energies))
    return coefficients


def main() -> None:
    sign_configurations = list(itertools.product((-1, 1), repeat=len(LABELS)))
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DETUNING
        energies = {
            signs: float(np.linalg.eigvalsh(walker_hamiltonian(signs, coupling))[0])
            for signs in sign_configurations
        }
        coefficients = walsh_coefficients(energies)
        target = coefficients["X_L_p_1_p_2_X_R"]
        unwanted = {
            label: value for label, value in coefficients.items()
            if label not in {"I", "X_L_p_1_p_2_X_R"}
        }
        positive = np.mean([energy for signs, energy in energies.items() if np.prod(signs) > 0])
        negative = np.mean([energy for signs, energy in energies.items() if np.prod(signs) < 0])
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": coupling,
            "low_branch_energy_by_gauss_eigenvalue": {"G_plus": float(positive), "G_minus": float(negative)},
            "gauss_coefficient": target,
            "gauss_sector_gap": float(abs(positive - negative)),
            "max_unwanted_walsh_coefficient": float(max(abs(value) for value in unwanted.values())),
            "walsh_coefficients": coefficients,
        })

    x = np.log([row["coupling_over_detuning"] for row in rows])
    y = np.log([abs(row["gauss_coefficient"]) for row in rows])
    power = float(np.polyfit(x, y, 1)[0])
    target_sign = np.sign(rows[-1]["gauss_coefficient"])
    output = {
        "schema": "antler.phase8b.lambda-walker-gadget-audit.v1",
        "parameters": {
            "block_size_b": 2,
            "detuning": DETUNING,
            "coupling_ratios": list(RATIOS),
            "walker_basis": "|vac>, |mu0>, |mu1>, |mu2>; at most one neutral walker",
            "vertices": "lambda X_L, w p1, w p2, lambda X_R; all equal to coupling in this symmetric audit",
            "target": "G_B=X_L p1 p2 X_R",
        },
        "rows": rows,
        "fits": {
            "gauss_coefficient_power_vs_coupling_over_detuning": power,
            "deepest_ratio_target_sign": int(target_sign),
        },
        "decision": (
            "The repaired Lambda-walker layer produces the coarse Gauss word at fourth order and no other non-scalar "
            "Walsh word in this ideal one-walker model. This validates the algebraic loop mechanism only; the required "
            "parity-conditioned neutral-walker hopping remains a newly declared microscopic primitive, not an ANTLER derivation."
        ),
        "claim_boundary": (
            "This test does not implement matter Fock states, charge-two mediators, Floquet channels, block-boundary pair "
            "transport, noise, local indistinguishability, a phase, fusion or a braid. It is only the exact internal-walker "
            "layer of the repaired Phase-8B proposal."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_lambda_walker_gadget_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
