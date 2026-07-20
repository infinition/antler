"""Phase 5G: defect/fusion-code target for a T-junction extension.

The scalar ladder and the synthetic spinor code both fail passive local
indistinguishability.  This exact four-defect model states the next required
encoding target: a fixed-parity two-dimensional fusion space in which every
single local defect operator has zero logical projection, while adjacent
exchanges are noncommuting and satisfy the braid relation.

It is an architecture benchmark for a future T-junction/defect realization,
not a claim that ANTLER's present hopping Hamiltonian has generated Majorana
zero modes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


I2 = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Y = np.array([[0.0, -1j], [1j, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def main() -> None:
    # Jordan--Wigner representation of four spatially separated Majoranas.
    gamma = [kron(X, I2), kron(Y, I2), kron(Z, X), kron(Z, Y)]
    parity = kron(Z, Z)
    # Even total-parity fusion space |00>, |11>.
    code = np.zeros((4, 2), complex)
    code[0, 0] = 1.0
    code[3, 1] = 1.0
    assert np.allclose(code.conj().T @ parity @ code, np.eye(2))
    B1 = expm(np.pi * gamma[1] @ gamma[0] / 4.0)
    B2 = expm(np.pi * gamma[2] @ gamma[1] / 4.0)
    U1, U2 = code.conj().T @ B1 @ code, code.conj().T @ B2 @ code
    local = []
    for j, operator in enumerate(gamma, 1):
        projected = code.conj().T @ operator @ code
        leakage = (np.eye(4) - code @ code.conj().T) @ operator @ code
        local.append({
            "defect": j,
            "logical_projection_norm": float(np.linalg.norm(projected)),
            "takes_code_out_of_sector_norm": float(np.linalg.norm(leakage)),
        })
    out = {
        "schema": "antler.phase5.tjunction-defect-code-target.v1",
        "claim_boundary": (
            "This is a defect/fusion target model for a future T-junction extension. "
            "It does not assert emergent Majoranas or topological protection in the "
            "present correlated-hopping ladder."
        ),
        "code_dimension": 2,
        "fixed_parity": "+1",
        "B1": arr(U1), "B2": arr(U2),
        "commutator_norm": float(np.linalg.norm(U1 @ U2 - U2 @ U1)),
        "braid_relation_residual": float(np.linalg.norm(U1 @ U2 @ U1 - U2 @ U1 @ U2)),
        "single_defect_locality_audit": local,
        "required_ANTLER_extension": [
            "T-junction or equivalent branch graph with independently movable defects",
            "gapped fixed-parity/fusion manifold with exponentially separated defects",
            "microscopic derivation of defect operators from the extended synthetic Hamiltonian",
            "full ED gap, leakage, noise and braid-convergence audits",
        ],
    }
    path = Path("results/phase5/tjunction_defect_code_target.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "commutator_norm": out["commutator_norm"],
        "braid_relation_residual": out["braid_relation_residual"],
        "local_projection_norms": [r["logical_projection_norm"] for r in local],
    }, indent=2))


if __name__ == "__main__":
    main()
