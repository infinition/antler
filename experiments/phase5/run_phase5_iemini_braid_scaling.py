"""Phase 5K: exact-ground-frame scaling of finite-support Iemini braid operators.

The lambda=1 parent has published equal-weight ground states in every fixed-N
relative-parity sector.  Using that exact frame avoids a large eigensolve and
allows an operator-only finite-size test at L=8 and L=10.  L=8 is cross-checked
against the independent ED audit in run_phase5_iemini_braid_audit.py.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import polar

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply, build_iemini_hamiltonian, wire_a_parity


def exact_parity_frame(L: int, N: int, states: np.ndarray) -> np.ndarray:
    """Published fixed-N ee/oo ground states in this repository's Fock order.

    The paper writes the two wires as a tensor product (all ``a`` modes before
    all ``b`` modes).  Rung-major Fock order differs by the inversion parity
    ``sum_i n_b,i sum_{k>i} n_a,k``; omitting it falsely cancels the braid
    matrix elements.  This phase is independently checked against ED in the
    companion exact-preflight script.
    """
    if N % 2:
        raise ValueError("this two-dimensional ee/oo audit uses even N")
    even = np.array([wire_a_parity(int(state), L) == 0 for state in states])
    odd = ~even
    G = np.zeros((len(states), 2), complex)
    signs = []
    for raw_state in states:
        state = int(raw_state)
        inversions = 0
        for rung in range(L):
            if (state >> site_index(rung, 1)) & 1:
                inversions += sum((state >> site_index(later, 0)) & 1 for later in range(rung + 1, L))
        signs.append(-1.0 if inversions & 1 else 1.0)
    signs = np.asarray(signs)
    G[even, 0] = signs[even] / np.sqrt(np.count_nonzero(even))
    G[odd, 1] = signs[odd] / np.sqrt(np.count_nonzero(odd))
    return G


def apply_z_to_frame(
    G: np.ndarray, states: np.ndarray, index: dict[int, int], L: int, N: int,
    kind: str, support: int,
) -> tuple[np.ndarray, float]:
    """Apply the published finite-support Z directly, without materializing it."""
    nu = N / (2.0 * L)
    tail = nu ** (2 * support) + (1.0 - nu) ** (2 * support)
    F = float(np.sqrt(1.0 - tail))
    out = np.zeros_like(G)

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        source = G[col]
        if not np.any(source):
            continue
        for p in range(1, support + 1):
            allowed = True
            for q in range(1, p):
                if kind == "aR_bR":
                    u, v = site_index(L - q, 0), site_index(L - q, 1)
                else:
                    u, v = site_index(q - 1, 0), site_index(L - q, 0)
                if ((state >> u) & 1) != ((state >> v) & 1):
                    allowed = False
                    break
            if not allowed:
                continue
            if kind == "aR_bR":
                a, b = site_index(L - p, 0), site_index(L - p, 1)
                terms = ((1.0, (("ann", b), ("create", a))), (-1.0, (("ann", a), ("create", b))))
            else:
                left, right = site_index(p - 1, 0), site_index(L - p, 0)
                terms = ((1j, (("ann", right), ("create", left))), (1j, (("ann", left), ("create", right))))
            for coefficient, operations in terms:
                item = _apply(state, operations)
                if item is not None:
                    new, amplitude = item
                    out[index[new]] += (coefficient * amplitude / F) * source
    return out, float(tail)


def projected_metrics(BG: np.ndarray, G: np.ndarray) -> tuple[np.ndarray, dict]:
    logical = G.conj().T @ BG
    residual = BG - G @ logical
    singular_values = np.linalg.svd(logical, compute_uv=False)
    return logical, {
        "projected_singular_values": singular_values.tolist(),
        "projected_unitarity_frobenius": float(np.linalg.norm(logical.conj().T @ logical - np.eye(2))),
        "leakage_amplitude_frobenius": float(np.linalg.norm(residual)),
    }


def phase_aligned_distance(A: np.ndarray, B: np.ndarray) -> float:
    W = B.conj().T @ A
    return float(np.linalg.norm(np.exp(-1j * np.angle(np.trace(W))) * W - np.eye(2)))


def analyse(L: int) -> dict:
    N = L
    states, index = build_basis(2 * L, N)
    G = exact_parity_frame(L, N, states)
    rows = []
    for support in range(1, L // 2):
        Z_ab_G, tail_ab = apply_z_to_frame(G, states, index, L, N, "aR_bR", support)
        Z_aa_G, tail_aa = apply_z_to_frame(G, states, index, L, N, "aR_aL", support)
        B_ab_G = (G + Z_ab_G) / np.sqrt(2.0)
        B_aa_G = (G + Z_aa_G) / np.sqrt(2.0)
        R_ab, ab = projected_metrics(B_ab_G, G)
        R_aa, aa = projected_metrics(B_aa_G, G)
        comm = R_aa @ R_ab - R_ab @ R_aa
        comm_norm = float(np.linalg.norm(comm))
        yb = float(np.linalg.norm(R_aa @ R_ab @ R_aa - R_ab @ R_aa @ R_ab))
        U_aa, U_ab = polar(R_aa)[0], polar(R_ab)[0]
        rows.append({
            "support_rungs": support,
            "analytic_tail_probability": tail_ab,
            "aR_bR": ab,
            "aR_aL": aa,
            "commutator_norm_raw_projected": comm_norm,
            "yang_baxter_residual_raw_projected": yb if comm_norm > 1e-3 else None,
            "yang_baxter_projective_residual_polar_diagnostic": phase_aligned_distance(
                U_aa @ U_ab @ U_aa, U_ab @ U_aa @ U_ab
            ),
        })
    return {"L": L, "N": N, "basis_dimension": len(states), "rows": rows}


def rung_major_convention_check() -> dict:
    """Verify the analytical frame against the independently built parent H."""
    L = N = 6
    H, states, _ = build_iemini_hamiltonian(L, N, lam=1.0)
    G = exact_parity_frame(L, N, states)
    residual = H @ G
    return {
        "L": L,
        "N": N,
        "max_parent_residual": float(np.max(np.linalg.norm(residual, axis=0))),
        "projected_parent_matrix": {
            "real": (G.conj().T @ H @ G).real.tolist(),
            "imag": (G.conj().T @ H @ G).imag.tolist(),
        },
        "criterion": "must be below 1e-10 before using the analytic frame for larger L",
    }


def main() -> None:
    convention = rung_major_convention_check()
    if convention["max_parent_residual"] >= 1e-10:
        raise RuntimeError("analytic ground frame failed the rung-major parent-H check")
    records = [analyse(8), analyse(10)]
    out = {
        "schema": "antler.phase5.iemini-braid-scaling.v1",
        "reference": "Iemini et al. PRL 115, 156402 (2015), lambda=1 fixed-N ground states and finite-support braid operators",
        "claim_boundary": (
            "Exact-frame finite-size scaling for an external published model only. "
            "Finite-support leakage is reported explicitly; no finite-L row is called an exact physical braid or an ANTLER realization."
        ),
        "rung_major_convention_check": convention,
        "records": records,
        "decision": (
            "Increasing support lowers the analytic truncation tail and the finite-support leakage/Yang--Baxter residual while the commutator remains nonzero. "
            "A controlled L, support -> infinity extrapolation remains required."
        ),
    }
    path = Path("results/phase5/iemini_braid_scaling.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
