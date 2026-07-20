"""Phase 6F: factor the exact Iemini bond interaction into charge-two channels.

This is an analytic/matrix audit.  It asks whether the *published external*
Iemini bond interaction can arise at second order from positive-detuning
charge-two mediators.  It does not claim that the resulting model is a new
ANTLER topological phase.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.number_conserving_pairwire import build_iemini_hamiltonian


PAIR_BASIS = (
    "a0 b0", "a0 a1", "b0 a1", "a0 b1", "b0 b1", "a1 b1",
)
# Rows are mediator-pair conversion channels.  For a detuning Delta, the
# microscopic conversion amplitudes are sqrt(Delta) times these coefficients.
COMPACT_CHANNELS = np.asarray([
    [0.0, 2.0, 0.0, 0.0, 2.0, 0.0],
    [0.0, 0.0, np.sqrt(2.0), np.sqrt(2.0), 0.0, 0.0],
    [np.sqrt(2.0), 0.0, 0.0, 0.0, 0.0, -np.sqrt(2.0)],
])


def main() -> None:
    H_free, states, _ = build_iemini_hamiltonian(L=2, N=2, lam=0.0)
    H_exact, _, _ = build_iemini_hamiltonian(L=2, N=2, lam=1.0)
    interaction = H_exact - H_free
    reconstructed = -COMPACT_CHANNELS.conj().T @ COMPACT_CHANNELS
    values, vectors = np.linalg.eigh(interaction)
    factor = np.diag(np.sqrt(np.maximum(-values, 0.0))) @ vectors.conj().T
    spectral_reconstruction = -factor.conj().T @ factor
    residual = float(np.linalg.norm(interaction - reconstructed, ord="fro"))
    spectral_residual = float(np.linalg.norm(interaction - spectral_reconstruction, ord="fro"))
    out = {
        "schema": "antler.phase6.iemini-bond-mediator-factorization.v1",
        "external_target": "Iemini et al. PRL 115, 156402 (2015), local lambda=1 interaction on one bond",
        "pair_basis_rung_major": list(PAIR_BASIS),
        "fixed_N2_masks": [int(state) for state in states],
        "interaction_matrix_lambda1_minus_lambda0": interaction.real.tolist(),
        "interaction_eigenvalues": values.tolist(),
        "negative_semidefinite": bool(float(np.max(values)) < 1e-12),
        "interaction_rank": int(np.count_nonzero(np.abs(values) > 1e-12)),
        "compact_channel_matrix_C": COMPACT_CHANNELS.tolist(),
        "channel_symmetry_charges": [
            {"channel": 0, "content": "a0a1 + b0b1", "Z2_a": 1, "Z2_b": 1},
            {"channel": 1, "content": "b0a1 + a0b1", "Z2_a": -1, "Z2_b": -1},
            {"channel": 2, "content": "a0b0 - a1b1", "Z2_a": -1, "Z2_b": -1},
        ],
        "mediator_rule": (
            "With three charge-two mediators of positive detuning Delta and conversion amplitudes G=sqrt(Delta) C, "
            "the second-order low-energy interaction is -G^dagger G / Delta = -C^dagger C. "
            "The two mixed-pair mediators carry generalized branch-parity charges (-1,-1)."
        ),
        "compact_factorization_frobenius_residual": residual,
        "spectral_factorization_frobenius_residual": spectral_residual,
        "qualification": bool(
            float(np.max(values)) < 1e-12 and int(np.count_nonzero(np.abs(values) > 1e-12)) == 3
            and residual < 1e-12 and spectral_residual < 1e-12
        ),
        "claim_boundary": (
            "This proves only a local operator factorization of an external published parent interaction. "
            "It neither derives that interaction from frozen ANTLER degrees of freedom nor proves a protected phase, braid, or quantum computer."
        ),
    }
    path = ROOT / "results" / "phase6" / "iemini_bond_mediator_factorization.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"], "eigenvalues": out["interaction_eigenvalues"],
        "rank": out["interaction_rank"], "compact_residual": residual,
    }, indent=2))


if __name__ == "__main__":
    main()
