"""Phase 5C: matrix-valued digital-transfer design for a synthetic dimension.

The N=3 pinned-mediator extension remains Abelian.  This script records the
minimal *different Hamiltonian ingredient* needed to cross that wall: a
two-component internal state with SU(2)-valued hopping links.  The digital
handoff lemma then transports a spinor by the link matrix, so two oriented
handoffs can implement noncommuting braid generators on one protected doublet.

This is a design target and algebraic verification, not evidence that the
current scalar correlated-hopping ladder already realises this extension.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


I = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Y = np.array([[0.0, -1j], [1j, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def su2(axis: np.ndarray, angle: float) -> np.ndarray:
    return expm(-0.5j * angle * axis)


def no_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def parallel_transport_link(U: np.ndarray, n: int = 801) -> np.ndarray:
    """Numerically verify |source,chi> -> |target,U chi> for a matrix link.

    The four-level Hamiltonian is the two-site handoff tensored with a spinor:
    H=[[eps_1 I, -J U^dag],[-J U, eps_2 I]].  Its two-dimensional ground
    subspace is Kato-parallel transported along the handoff.
    """

    D, J = 20.0, 1.0
    source = np.vstack((I, np.zeros((2, 2), complex)))
    target_spin = np.vstack((np.zeros((2, 2), complex), I))
    transported = source.copy()
    for s in np.linspace(0.0, 1.0, n):
        q = np.sin(0.5 * np.pi * s) ** 2
        H = np.block([
            [-D * (1.0 - q) * I, -J * U.conj().T],
            [-J * U, -D * q * I],
        ])
        _, vectors = np.linalg.eigh(H)
        ground = vectors[:, :2]
        overlap = ground.conj().T @ transported
        left, _, right = np.linalg.svd(overlap)
        transported = ground @ (left @ right)
    # In the fixed target-orbital spin basis, the transported columns are U.
    overlap = target_spin.conj().T @ transported
    left, _, right = np.linalg.svd(overlap)
    return left @ right


def main() -> None:
    # Adjacent Ising/Majorana braid generators, expressed as physical SU(2)
    # hopping links for a synthetic spin/orbital dimension.
    B1 = su2(Z, np.pi / 2.0)
    B2 = su2(X, np.pi / 2.0)
    transfer_1 = parallel_transport_link(B1)
    transfer_2 = parallel_transport_link(B2)
    commutator = B1 @ B2 - B2 @ B1
    braid = B1 @ B2 @ B1 - B2 @ B1 @ B2
    out = {
        "schema": "antler.phase5.synthetic-su2-link-design.v1",
        "claim_boundary": (
            "This is an algebraic synthetic-dimension design target.  It requires a "
            "new matrix-valued hopping Hamiltonian and does not upgrade the present "
            "scalar rung-major correlated-hopping ladder to a non-Abelian system."
        ),
        "required_extension": {
            "internal_dimension": 2,
            "link_hamiltonian": "H_link=[[eps1 I,-J U^dag],[-J U,eps2 I]]",
            "implementation_requirement": (
                "independently addressable non-collinear SU(2) Peierls links, e.g. "
                "a pseudo-spin or orbital synthetic dimension"
            ),
        },
        "B1": arr(B1), "B2": arr(B2),
        "handoff_transport_B1": arr(transfer_1),
        "handoff_transport_B2": arr(transfer_2),
        "transport_error_B1": float(np.linalg.norm(no_global(transfer_1.conj().T @ B1) - I)),
        "transport_error_B2": float(np.linalg.norm(no_global(transfer_2.conj().T @ B2) - I)),
        "commutator_norm": float(np.linalg.norm(commutator)),
        "braid_relation_residual": float(np.linalg.norm(braid)),
        "conclusion": (
            "Matrix-valued links can supply two noncommuting digital transfers.  The "
            "next implementation task is to extend ANTLER's ED Hamiltonian with this "
            "synthetic two-component degree of freedom and repeat the leakage, gap, "
            "convergence, and noise audits."
        ),
    }
    path = Path("results/phase5/synthetic_su2_link_design.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "commutator_norm": out["commutator_norm"],
        "braid_relation_residual": out["braid_relation_residual"],
        "transport_error_B1": out["transport_error_B1"],
        "transport_error_B2": out["transport_error_B2"],
    }, indent=2))


if __name__ == "__main__":
    main()
