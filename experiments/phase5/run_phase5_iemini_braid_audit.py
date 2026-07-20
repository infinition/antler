"""Phase 5J: finite-support Iemini braid-operator audit at fixed U(1) charge.

The operators are the published finite-support constructions, projected onto
the exact lambda=1 parity doublet.  This is a calibration of the audit stack
against a known number-conserving reference, not a microscopic ANTLER braid.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import polar
from scipy.sparse import issparse
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.number_conserving_pairwire import (
    build_iemini_braid_z,
    build_iemini_hamiltonian,
    wire_a_parity,
)


def ground_state(H, states: np.ndarray, L: int, parity: int) -> np.ndarray:
    rows = np.array([i for i, state in enumerate(states) if wire_a_parity(int(state), L) == parity])
    block = H[rows, :][:, rows]
    values, vectors = eigsh(block, k=1, which="SA", tol=1e-12)
    full = np.zeros(len(states), dtype=complex)
    full[rows] = vectors[:, 0]
    return full


def phase_aligned_distance(A: np.ndarray, B: np.ndarray) -> float:
    W = B.conj().T @ A
    return float(np.linalg.norm(np.exp(-1j * np.angle(np.trace(W))) * W - np.eye(W.shape[0])))


def gate_metrics(B, G: np.ndarray) -> tuple[np.ndarray, dict]:
    projected = G.conj().T @ (B @ G)
    residual = B @ G - G @ projected
    singular_values = np.linalg.svd(projected, compute_uv=False)
    return projected, {
        "projected_singular_values": singular_values.tolist(),
        "projected_unitarity_frobenius": float(np.linalg.norm(projected.conj().T @ projected - np.eye(2))),
        "leakage_amplitude_frobenius": float(np.linalg.norm(residual)),
    }


def main() -> None:
    L = N = 8
    H, states, index = build_iemini_hamiltonian(L, N, lam=1.0, sparse=True)
    if not issparse(H):
        raise RuntimeError("L=8 audit requires the sparse reference Hamiltonian")
    # The two columns have definite wire parity.  Their arbitrary individual
    # phases do not affect all reported norms or the Yang--Baxter residual.
    G = np.column_stack((ground_state(H, states, L, 0), ground_state(H, states, L, 1)))
    if not np.allclose(G.conj().T @ G, np.eye(2), atol=1e-10):
        raise RuntimeError("ground frame is not orthonormal")

    rows = []
    commutator_threshold = 1e-3
    for support in (1, 2, 3):
        Z_ab, F_ab = build_iemini_braid_z(L, N, "aR_bR", support, basis=(states, index))
        Z_aa, F_aa = build_iemini_braid_z(L, N, "aR_aL", support, basis=(states, index))
        I = __import__("scipy").sparse.identity(len(states), dtype=complex, format="csr")
        R_ab = (I + Z_ab) / np.sqrt(2.0)
        R_aa = (I + Z_aa) / np.sqrt(2.0)
        r_ab, ab_metrics = gate_metrics(R_ab, G)
        r_aa, aa_metrics = gate_metrics(R_aa, G)
        commutator = r_aa @ r_ab - r_ab @ r_aa
        commutator_norm = float(np.linalg.norm(commutator))
        aba = r_aa @ r_ab @ r_aa
        bab = r_ab @ r_aa @ r_ab
        yb_raw = float(np.linalg.norm(aba - bab))
        # The polar matrices diagnose the intended logical representation but
        # are never substituted for the raw finite-support operators.
        u_aa = polar(r_aa)[0]
        u_ab = polar(r_ab)[0]
        yb_polar = phase_aligned_distance(u_aa @ u_ab @ u_aa, u_ab @ u_aa @ u_ab)
        rows.append({
            "support_rungs": support,
            "normalizers": {"aR_bR": F_ab, "aR_aL": F_aa},
            "support_geometry": {
                "aR_bR": f"rightmost {support} rungs only",
                "aR_aL": f"bilocal: leftmost and rightmost {support} rungs only",
            },
            "aR_bR": ab_metrics,
            "aR_aL": aa_metrics,
            "commutator_norm_raw_projected": commutator_norm,
            "yang_baxter_residual_raw_projected": (
                yb_raw if commutator_norm > commutator_threshold else None
            ),
            "yang_baxter_status": (
                "reported_with_nonzero_commutator" if commutator_norm > commutator_threshold
                else "not_interpretable_commutator_below_threshold"
            ),
            "yang_baxter_projective_residual_polar_diagnostic": yb_polar,
        })

    best = rows[-1]
    out = {
        "schema": "antler.phase5.iemini-braid-audit.v1",
        "reference": {
            "citation": "F. Iemini et al., Phys. Rev. Lett. 115, 156402 (2015), braid construction and supplement",
            "arxiv": "https://arxiv.org/abs/1504.04230",
            "model": "published number-conserving two-wire parent Hamiltonian at lambda=1",
            "L": L, "N": N, "filling": N / (2.0 * L),
        },
        "claim_boundary": (
            "This is a finite-size projected-operator audit of the published Iemini braid construction. "
            "It calibrates the non-Abelian audit infrastructure on an external U(1)-conserving reference. "
            "It is not an adiabatic physical exchange, not a universal gate-set proof, and not a derivation from frozen ANTLER."
        ),
        "gates": {
            "commutator_threshold": commutator_threshold,
            "yang_baxter_rule": "Report the raw residual only when the raw projected commutator exceeds threshold.",
        },
        "rows": rows,
        "decision": (
            "finite-support operators are edge/bilocal-edge supported and display a nonzero projected commutator; finite-support leakage and raw Yang--Baxter residual must be scaled in L and support before promoting a braid representation"
            if best["commutator_norm_raw_projected"] > commutator_threshold else
            "commutator gate failed; do not interpret any Yang--Baxter residual"
        ),
    }
    path = Path("results/phase5/iemini_braid_audit.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
