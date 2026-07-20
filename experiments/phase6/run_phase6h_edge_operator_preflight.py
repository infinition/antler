"""Phase 6H: protection-first edge-operator audit on the external parent.

For a finite-support candidate O and exact code projector P, the relevant
finite-size diagnostic is ``||(I-P)[H,O]P||``.  The unprojected operator norm
of ``[H,O]`` probes arbitrary excited states and is not, by itself, a logical
edge-protection criterion.  This script uses the published truncated Iemini
edge generator as a calibration target and is expected to be able to fail.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis
from antler.number_conserving_pairwire import build_iemini_braid_z, build_iemini_hamiltonian
from experiments.phase5.run_phase5_iemini_braid_scaling import exact_parity_frame


def analyse(L: int, support: int) -> dict:
    N = L
    states, index = build_basis(2 * L, N)
    H, _, _ = build_iemini_hamiltonian(L=L, N=N, lam=1.0, basis=(states, index), sparse=True)
    G = exact_parity_frame(L, N, states)
    Z, normalizer = build_iemini_braid_z(L, N, "aR_bR", support, basis=(states, index))
    ZG = Z @ G
    logical = G.conj().T @ ZG
    leakage = ZG - G @ logical
    # H G=0 for this exact parent, therefore H Z G = [H,Z]G.
    commutator_action = H @ ZG
    tail = (N / (2.0 * L)) ** (2 * support) + (1.0 - N / (2.0 * L)) ** (2 * support)
    return {
        "L": L,
        "N": N,
        "support_rungs": support,
        "analytic_truncation_tail_probability": float(tail),
        "normalizer": normalizer,
        "logical_action_frobenius": float(np.linalg.norm(logical)),
        "logical_leakage_amplitude_frobenius": float(np.linalg.norm(leakage)),
        "code_commutator_action_frobenius": float(np.linalg.norm(commutator_action)),
        "code_commutator_action_normalized": float(np.linalg.norm(commutator_action) / np.linalg.norm(ZG)),
    }


def main() -> None:
    rows = [analyse(L, support) for L in (6, 8) for support in range(1, L // 2)]
    out = {
        "schema": "antler.phase6.edge-operator-protection-preflight.v1",
        "reference": "external Iemini lambda=1 parent; published finite-support aR-bR edge generator",
        "criterion": {
            "finite_size_quantity": "||(I-P)[H,O]P|| = ||[H,O]P|| because H P=0",
            "why_not_full_operator_norm": "the full commutator includes arbitrary bulk-excited states and is not a code-protection certificate",
            "promotion_requirement": "the code-commutator action, leakage, and support-tail must all decrease under a controlled L/support scaling",
        },
        "rows": rows,
        "claim_boundary": (
            "This calibrates a protection-first diagnostic on an external parent. It does not prove that a finite-support row is an exact zero mode or physical braid, "
            "and it says nothing yet about a native ANTLER edge operator."
        ),
        "decision": (
            "finite-support published operators remain visibly non-conserved at these sizes; retain them only as convergence diagnostics and do not promote a protected physical braid"
        ),
    }
    path = ROOT / "results" / "phase6" / "edge_operator_protection_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
