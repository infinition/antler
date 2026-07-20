"""Minimal U(1)-conserving Z2 Gauss-law algebra control for Phase 8.

This is deliberately an *abstract local reference*, not an ANTLER
microscopic derivation.  It fixes the exact algebra that any proposed neutral
gauge ancilla must reproduce before a T-junction or braid calculation can be
interpreted.  A one-rung fixed-charge rail doublet is coupled to a neutral
gauge qubit through G=(-1)^N_a X_g.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LAMBDA = 3.0


def restricted_matrix(operator: np.ndarray, frame: np.ndarray) -> np.ndarray:
    return frame.conj().T @ operator @ frame


def nonscalar_norm(matrix: np.ndarray) -> float:
    dimension = matrix.shape[0]
    return float(np.linalg.norm(matrix - np.trace(matrix) * np.eye(dimension) / dimension, ord=2))


def main() -> None:
    identity_rail = np.eye(2, dtype=complex)
    identity_gauge = np.eye(2, dtype=complex)
    x = np.asarray(((0.0, 1.0), (1.0, 0.0)), dtype=complex)
    z = np.asarray(((1.0, 0.0), (0.0, -1.0)), dtype=complex)
    # Fixed total physical charge N=1: basis |a>, |b>.  Rail parity is -1 on
    # |a> and +1 on |b>, i.e. (-1)^N_a = diag(-1,+1).
    rail_parity = np.diag((-1.0, 1.0)).astype(complex)
    bare_rail_tunnel = np.kron(x, identity_gauge)
    gauge_dressed_tunnel = np.kron(x, z)
    gauss = np.kron(rail_parity, x)
    identity = np.eye(4, dtype=complex)
    physical_projector = (identity + gauss) / 2.0
    h_gauss = LAMBDA * (identity - gauss) / 2.0

    values, vectors = np.linalg.eigh(h_gauss)
    frame = vectors[:, np.abs(values) < 1e-12]
    if frame.shape[1] != 2:
        raise RuntimeError("minimal Gauss control did not have a two-dimensional physical sector")
    projected_bare = restricted_matrix(bare_rail_tunnel, frame)
    projected_dressed = restricted_matrix(gauge_dressed_tunnel, frame)
    if np.linalg.norm(projected_bare) > 1e-12:
        raise RuntimeError("bare rail tunnelling survived the Gauss projection")
    if nonscalar_norm(projected_dressed) < 0.9:
        raise RuntimeError("gauge-invariant dressed local tunnelling failed to expose local code readability")

    # H_G = Lambda/2(I-X_g) + Lambda n_a X_g, with n_a=(I-rail_parity)/2.
    number_a = (identity_rail - rail_parity) / 2.0
    decomposed_h_gauss = (
        LAMBDA / 2.0 * np.kron(identity_rail, identity_gauge - x)
        + LAMBDA * np.kron(number_a, x)
    )
    if np.linalg.norm(h_gauss - decomposed_h_gauss) > 1e-12:
        raise RuntimeError("Gauss constraint decomposition failed")

    output = {
        "schema": "antler.phase8.local-gauge-constraint-preflight.v1",
        "parameters": {
            "physical_rail_charge_sector": "N=1, basis |a>, |b>",
            "new_resource": "one neutral two-level gauge ancilla at the vertex",
            "lambda": LAMBDA,
            "gauss_generator": "G=(-1)^N_a X_g",
            "constraint": "H_G=lambda(1-G)/2=lambda/2(1-X_g)+lambda n_a X_g",
        },
        "spectrum": {
            "eigenvalues": [float(value) for value in values],
            "physical_sector_dimension": int(frame.shape[1]),
            "gauss_constraint_gap": float(values[2] - values[1]),
        },
        "exact_algebra": {
            "gauss_squared_residual": float(np.linalg.norm(gauss @ gauss - identity)),
            "h_gauss_commutator_with_gauss": float(np.linalg.norm(h_gauss @ gauss - gauss @ h_gauss)),
            "bare_rail_tunnel_anticommutator_with_gauss": float(np.linalg.norm(bare_rail_tunnel @ gauss + gauss @ bare_rail_tunnel)),
            "dressed_rail_tunnel_commutator_with_gauss": float(np.linalg.norm(gauge_dressed_tunnel @ gauss - gauss @ gauge_dressed_tunnel)),
            "projector_idempotence_residual": float(np.linalg.norm(physical_projector @ physical_projector - physical_projector)),
            "projected_bare_rail_tunnel_norm": float(np.linalg.norm(projected_bare, ord=2)),
            "projected_dressed_rail_tunnel_norm": float(np.linalg.norm(projected_dressed, ord=2)),
            "projected_dressed_rail_tunnel_nonscalar_norm": nonscalar_norm(projected_dressed),
            "constraint_decomposition_residual": float(np.linalg.norm(h_gauss - decomposed_h_gauss)),
        },
        "decision": (
            "A local Gauss constraint can algebraically forbid bare U(1)-conserving rail tunnelling while allowing a "
            "gauge-dressed process. The same minimal block is not a topological code: a gauge-invariant local dressed "
            "operator remains non-scalar on its two-dimensional physical sector."
        ),
        "implementation_gate": (
            "Any claimed ANTLER gauge extension must derive the neutral ancilla and the n_a X_g constraint from a "
            "declared microscopic resource, preserve the algebra below 1e-12 on its explicit weighted Fock block, and "
            "then pass a multi-vertex local-indistinguishability audit."
        ),
        "claim_boundary": (
            "This is a four-dimensional algebra calibration with an inserted neutral gauge ancilla. It does not derive "
            "that ancilla from ANTLER charge-two mediators, establish a T-junction, a thermodynamic phase, fusion, a braid, "
            "non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "local_gauge_constraint_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
