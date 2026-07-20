"""Exact selection-rule audit of the stated Phase-8B single-mediator gadget.

For the proposed block A=X_L+X_R+eta P_B and
H=Delta_G n_mu + lambda(mu^dagger A + A mu), returning from the low sector
n_mu=0 requires an even number of conversion vertices.  This script computes
the exact low branch and its complete abelian Pauli expansion to test whether
the claimed odd product G_B=X_L P_B X_R can occur at order lambda^3.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ETA = 0.37
DELTA = 10.0
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)


def kron3(left: np.ndarray, center: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.kron(np.kron(left, center), right)


def matrix_function_hermitian(matrix: np.ndarray, function) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * function(values)) @ vectors.conj().T


def coefficient(operator: np.ndarray, basis_operator: np.ndarray) -> float:
    return float(np.trace(basis_operator.conj().T @ operator).real / operator.shape[0])


def main() -> None:
    identity = np.eye(2, dtype=complex)
    x = np.asarray(((0.0, 1.0), (1.0, 0.0)), dtype=complex)
    z = np.asarray(((1.0, 0.0), (0.0, -1.0)), dtype=complex)
    # Tensor order: left boundary gauge qubit, matter-parity qubit, right
    # boundary gauge qubit.  All three factors commute.
    xl = kron3(x, identity, identity)
    p_block = kron3(identity, z, identity)
    xr = kron3(identity, identity, x)
    target_gauss = xl @ p_block @ xr
    a_operator = xl + xr + ETA * p_block
    a_squared = a_operator @ a_operator

    basis = {
        "I": np.eye(8, dtype=complex),
        "X_L": xl,
        "P_B": p_block,
        "X_R": xr,
        "X_L_X_R": xl @ xr,
        "X_L_P_B": xl @ p_block,
        "P_B_X_R": p_block @ xr,
        "G_B_X_L_P_B_X_R": target_gauss,
    }
    rows = []
    for ratio in RATIOS:
        coupling = ratio * DELTA
        # Exact low branch of the 2x2 mediator problem in every A eigenvalue:
        # e_-(A)=(Delta-sqrt(Delta^2+4 lambda^2 A^2))/2.
        h_eff = matrix_function_hermitian(
            a_squared,
            lambda values: (DELTA - np.sqrt(DELTA**2 + 4.0 * coupling**2 * values)) / 2.0,
        )
        coefficients = {name: coefficient(h_eff, operator) for name, operator in basis.items()}
        second_order = -(coupling**2 / DELTA) * a_squared
        rows.append({
            "lambda_over_delta": ratio,
            "lambda": coupling,
            "exact_low_branch_coefficients": coefficients,
            "exact_vs_second_order_frobenius": float(np.linalg.norm(h_eff - second_order)),
        })

    target_max = max(abs(row["exact_low_branch_coefficients"]["G_B_X_L_P_B_X_R"]) for row in rows)
    lower_order_coefficients = [
        abs(row["exact_low_branch_coefficients"]["X_L_X_R"]) for row in rows
    ]
    # The exact symmetry is H_eff=f(A^2): it lies in the even subalgebra
    # span{I, X_L X_R, X_L P_B, P_B X_R}; the odd triple product cannot occur.
    even_basis = np.stack([
        operator.reshape(-1) for name, operator in basis.items()
        if name in {"I", "X_L_X_R", "X_L_P_B", "P_B_X_R"}
    ], axis=1)
    residuals = []
    for row in rows:
        ratio = row["lambda_over_delta"]
        coupling = ratio * DELTA
        h_eff = matrix_function_hermitian(
            a_squared,
            lambda values: (DELTA - np.sqrt(DELTA**2 + 4.0 * coupling**2 * values)) / 2.0,
        )
        fit, *_ = np.linalg.lstsq(even_basis, h_eff.reshape(-1), rcond=None)
        residuals.append(float(np.linalg.norm(h_eff.reshape(-1) - even_basis @ fit)))

    output = {
        "schema": "antler.phase8b.gauss-odd-order-selection-audit.v1",
        "parameters": {
            "eta": ETA,
            "detuning_delta_g": DELTA,
            "ratios_lambda_over_delta_g": list(RATIOS),
            "stated_gadget": "H=Delta_G n_mu + lambda(mu^dagger A + A mu), A=X_L+X_R+eta P_B",
            "target": "G_B=X_L P_B X_R",
        },
        "rows": rows,
        "exact_selection_rule": {
            "low_energy_form": "H_eff=f(A^2) because each low-to-low virtual process has an even number of mu conversion vertices",
            "even_commutant_basis": ["I", "X_L X_R", "X_L P_B", "P_B X_R"],
            "max_target_gauss_coefficient": target_max,
            "max_even_subalgebra_reconstruction_residual": max(residuals),
            "minimum_abs_second_order_xl_xr_coefficient": min(lower_order_coefficients),
        },
        "decision": (
            "The stated single-linear-mediator gadget cannot generate G_B=X_L P_B X_R at third order: its exact low branch "
            "is an even function of A and the target coefficient is zero to numerical precision. This refutes T3 as written, "
            "but not every coarse-grained gauge construction; a viable revision must explicitly break the mediator-conversion "
            "selection rule or introduce a different microscopic primitive and then re-audit all lower-order terms."
        ),
        "claim_boundary": (
            "This exact algebra audit addresses only the Hamiltonian displayed in PHASE8B_Z2_GAUSS_RESOURCE_THEOREM.md. "
            "It does not refute T1, T2, a modified coarse-grained gadget, a gauge theory, a phase, fusion or a braid."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_gauss_odd_order_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
