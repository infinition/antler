"""Phase 5F: L=14 synthetic-code spectrum and local-distinguishability audit.

The active pseudo-spin is encoded at the left edge and a spin-polarised anchor
particle occupies the right edge.  This is the smallest full-ladder candidate
that combines an isolated spinor code with matrix-valued hopping.  The audit
reports both spectral isolation and the action of a local spin measurement;
the latter is a required honesty test for passive topological-memory claims.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.synthetic import build_synthetic_hamiltonian, mode, su2


X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def main() -> None:
    L, N, theta = 14, 2, 0.3
    n_sites = 2 * L
    active_site, anchor_site = 0, n_sites - 1
    mu = np.zeros(n_sites)
    mu[active_site] = -4.0
    mu[anchor_site] = -6.0
    # Fix the anchor to spin 0 while retaining a degenerate active spinor.
    anchor_field = np.diag([-1.0, 1.0])
    B1 = su2(Z, np.pi / 2.0)
    B2 = su2(X, np.pi / 2.0)
    rows = []
    for label, link in (("B1", B1), ("B2", B2)):
        H, states, index = build_synthetic_hamiltonian(
            L, N, theta, .4, 1.0, .1, mu=mu,
            link_matrices={(0, 2): link}, onsite_spin_matrices={anchor_site: anchor_field},
        )
        energies, vectors = eigsh(H, k=12, which="SA", tol=1e-10)
        order = np.argsort(energies); energies, vectors = energies[order], vectors[:, order]
        target_indices = [
            index[(1 << mode(active_site, spin)) | (1 << mode(anchor_site, 0))]
            for spin in range(2)
        ]
        target = np.zeros((len(states), 2), complex)
        for spin, row in enumerate(target_indices):
            target[row, spin] = 1.0
        score = np.sum(abs(target.conj().T @ vectors) ** 2, axis=0)
        selected = np.sort(np.argsort(-score)[:2])
        frame = vectors[:, selected]
        # Fix the code gauge against the two bare spin states.
        U, _, Vh = np.linalg.svd(target.conj().T @ frame)
        frame = frame @ (U @ Vh).conj().T
        outside = np.setdiff1d(np.arange(len(energies)), selected)
        gap = float(np.min(abs(energies[selected, None] - energies[outside])))
        # Local Z measurement on the active site acts nontrivially if the
        # proposed code is a conventional local spin qubit rather than a
        # passive topological memory.
        local_z = np.zeros(len(states))
        for i, state in enumerate(states):
            raw = int(state)
            if (raw >> mode(active_site, 0)) & 1:
                local_z[i] = 1.0
            elif (raw >> mode(active_site, 1)) & 1:
                local_z[i] = -1.0
        projected_z = frame.conj().T @ (local_z[:, None] * frame)
        traceless = projected_z - .5 * np.trace(projected_z) * np.eye(2)
        rows.append({
            "link": label, "selected_energies": energies[selected].tolist(),
            "capture": float(score[selected].sum()), "isolation_gap": gap,
            "logical_split": float(abs(energies[selected[1]] - energies[selected[0]])),
            "projected_local_spin_Z": arr(projected_z),
            "local_distinguishability_norm": float(np.linalg.norm(traceless)),
        })
    out = {
        "schema": "antler.phase5.synthetic-l14-code-preflight.v1",
        "claim_boundary": (
            "This is a full-L=14 spectral preflight for the synthetic extension. "
            "A nonzero local-distinguishability norm explicitly rules out calling "
            "this spinor code a passive topological memory."
        ),
        "rows": rows,
        "decision": (
            "matrix-link code is spectrally testable, but local spin distinguishability "
            "requires a further defect/fusion encoding for topological protection"
        ),
    }
    path = Path("results/phase5/synthetic_l14_code_preflight.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
