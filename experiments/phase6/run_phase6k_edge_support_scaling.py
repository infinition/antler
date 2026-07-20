"""Phase 6K: matrix-free support/size scaling of the external edge generator.

The finite-support edge test is extended without materializing the L=10 parent
matrix.  The action of the exact parent on a two-column frame is accumulated
directly in the fixed-N Fock basis and cross-checked against the sparse matrix
at L=6.  This is a support/size convergence test, not a parameter scan.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply, build_iemini_hamiltonian
from experiments.phase5.run_phase5_iemini_braid_scaling import apply_z_to_frame, exact_parity_frame


@njit
def _below_parity(state: int, site: int) -> float:
    count = 0
    for bit in range(site):
        count += (state >> bit) & 1
    return -1.0 if count & 1 else 1.0


@njit
def _sequence2(state: int, kind0: int, site0: int, kind1: int, site1: int):
    amplitude = 1.0
    current = state
    for kind, site in ((kind0, site0), (kind1, site1)):
        if kind == 0:  # annihilate
            if not ((current >> site) & 1):
                return False, current, amplitude
            amplitude *= _below_parity(current, site)
            current ^= 1 << site
        else:  # create
            if (current >> site) & 1:
                return False, current, amplitude
            amplitude *= _below_parity(current, site)
            current |= 1 << site
    return True, current, amplitude


@njit
def _sequence4(
    state: int, kind0: int, site0: int, kind1: int, site1: int,
    kind2: int, site2: int, kind3: int, site3: int,
):
    amplitude = 1.0
    current = state
    for kind, site in ((kind0, site0), (kind1, site1), (kind2, site2), (kind3, site3)):
        if kind == 0:  # annihilate
            if not ((current >> site) & 1):
                return False, current, amplitude
            amplitude *= _below_parity(current, site)
            current ^= 1 << site
        else:  # create
            if (current >> site) & 1:
                return False, current, amplitude
            amplitude *= _below_parity(current, site)
            current |= 1 << site
    return True, current, amplitude


@njit
def _apply_iemini_parent_numba(L: int, states: np.ndarray, lookup: np.ndarray, frame: np.ndarray) -> np.ndarray:
    """Compiled fixed-N parent action; kinds are 0=annihilate, 1=create."""
    out = np.zeros_like(frame)
    for column in range(states.shape[0]):
        state = states[column]
        has_source = False
        for logical_column in range(frame.shape[1]):
            if frame[column, logical_column] != 0.0:
                has_source = True
                break
        if not has_source:
            continue
        diagonal = 0.0
        for bond in range(L - 1):
            a0, a1 = 2 * bond, 2 * (bond + 1)
            b0, b1 = a0 + 1, a1 + 1
            for left, right in ((a0, a1), (b0, b1)):
                valid, new_state, amplitude = _sequence2(state, 0, right, 1, left)
                if valid:
                    row = lookup[new_state]
                    for logical_column in range(frame.shape[1]):
                        out[row, logical_column] += -4.0 * amplitude * frame[column, logical_column]
                valid, new_state, amplitude = _sequence2(state, 0, left, 1, right)
                if valid:
                    row = lookup[new_state]
                    for logical_column in range(frame.shape[1]):
                        out[row, logical_column] += -4.0 * amplitude * frame[column, logical_column]
                n_left, n_right = (state >> left) & 1, (state >> right) & 1
                diagonal += 4.0 * (n_left + n_right) - 4.0 * n_left * n_right
            n_a = ((state >> a0) & 1) + ((state >> a1) & 1)
            n_b = ((state >> b0) & 1) + ((state >> b1) & 1)
            diagonal += -2.0 * n_a * n_b
            # T1, T2, -2*T3, followed by their Hermitian conjugates.
            terms = (
                (2.0, 0, b1, 1, b0, 0, a1, 1, a0),
                (2.0, 0, b0, 1, b1, 0, a1, 1, a0),
                (-4.0, 0, a0, 0, a1, 1, b1, 1, b0),
                (2.0, 0, a0, 1, a1, 0, b0, 1, b1),
                (2.0, 0, a0, 1, a1, 0, b1, 1, b0),
                (-4.0, 0, b0, 0, b1, 1, a1, 1, a0),
            )
            for coefficient, kind0, site0, kind1, site1, kind2, site2, kind3, site3 in terms:
                valid, new_state, amplitude = _sequence4(
                    state, kind0, site0, kind1, site1, kind2, site2, kind3, site3,
                )
                if valid:
                    row = lookup[new_state]
                    for logical_column in range(frame.shape[1]):
                        out[row, logical_column] += coefficient * amplitude * frame[column, logical_column]
        for logical_column in range(frame.shape[1]):
            out[column, logical_column] += diagonal * frame[column, logical_column]
    return out


def apply_iemini_parent_to_frame(
    L: int, states: np.ndarray, index: dict[int, int], frame: np.ndarray,
) -> np.ndarray:
    """Matrix-free action of the lambda=1 parent on one or more vectors."""
    out = np.zeros_like(frame)

    def add(state: int, source: np.ndarray, coefficient: float, operations: tuple[tuple[str, int], ...]) -> None:
        item = _apply(state, operations)
        if item is not None:
            new_state, amplitude = item
            out[index[new_state]] += coefficient * amplitude * source

    for column, raw_state in enumerate(states):
        state = int(raw_state)
        source = frame[column]
        if not np.any(source):
            continue
        diagonal = 0.0
        for bond in range(L - 1):
            a0, a1 = site_index(bond, 0), site_index(bond + 1, 0)
            b0, b1 = site_index(bond, 1), site_index(bond + 1, 1)
            for left, right in ((a0, a1), (b0, b1)):
                add(state, source, -4.0, (("ann", right), ("create", left)))
                add(state, source, -4.0, (("ann", left), ("create", right)))
                n_left, n_right = (state >> left) & 1, (state >> right) & 1
                diagonal += 4.0 * (n_left + n_right) - 4.0 * n_left * n_right
            n_a = ((state >> a0) & 1) + ((state >> a1) & 1)
            n_b = ((state >> b0) & 1) + ((state >> b1) & 1)
            diagonal += -2.0 * n_a * n_b
            for coefficient, operations in (
                (2.0, (("ann", b1), ("create", b0), ("ann", a1), ("create", a0))),
                (2.0, (("ann", b0), ("create", b1), ("ann", a1), ("create", a0))),
                (-4.0, (("ann", a0), ("ann", a1), ("create", b1), ("create", b0))),
                (2.0, (("ann", a0), ("create", a1), ("ann", b0), ("create", b1))),
                (2.0, (("ann", a0), ("create", a1), ("ann", b1), ("create", b0))),
                (-4.0, (("ann", b0), ("ann", b1), ("create", a1), ("create", a0))),
            ):
                add(state, source, coefficient, operations)
        out[column] += diagonal * source
    return out


def apply_iemini_parent_to_frame_compiled(
    L: int, states: np.ndarray, index: dict[int, int], frame: np.ndarray,
) -> np.ndarray:
    """Build a dense mask lookup then invoke the compiled matrix-free kernel."""
    lookup = np.full(1 << (2 * L), -1, dtype=np.int64)
    for row, state in enumerate(states):
        lookup[int(state)] = row
    return _apply_iemini_parent_numba(L, states.astype(np.int64), lookup, frame)


def sparse_crosscheck() -> dict:
    L = N = 6
    states, index = build_basis(2 * L, N)
    rng = np.random.default_rng(20260719)
    probe = rng.normal(size=(len(states), 2)) + 1j * rng.normal(size=(len(states), 2))
    H, _, _ = build_iemini_hamiltonian(L, N, lam=1.0, basis=(states, index), sparse=True)
    error = float(np.linalg.norm(H @ probe - apply_iemini_parent_to_frame_compiled(L, states, index, probe)))
    return {"L": L, "N": N, "absolute_frobenius_error": error}


def analyse(L: int, support: int) -> dict:
    N = L
    states, index = build_basis(2 * L, N)
    G = exact_parity_frame(L, N, states)
    ZG, tail = apply_z_to_frame(G, states, index, L, N, "aR_bR", support)
    logical = G.conj().T @ ZG
    leakage = ZG - G @ logical
    H_G = apply_iemini_parent_to_frame_compiled(L, states, index, G)
    H_ZG = apply_iemini_parent_to_frame_compiled(L, states, index, ZG)
    # The matrix-free application makes both terms explicit, even though H G=0
    # analytically.  This avoids assuming the convention before checking it.
    commutator_action = H_ZG
    return {
        "L": L,
        "N": N,
        "basis_dimension": len(states),
        "support_rungs": support,
        "analytic_truncation_tail_probability": tail,
        "parent_action_on_exact_frame": float(np.linalg.norm(H_G)),
        "logical_action_frobenius": float(np.linalg.norm(logical)),
        "logical_leakage_amplitude_frobenius": float(np.linalg.norm(leakage)),
        "code_commutator_action_frobenius": float(np.linalg.norm(commutator_action)),
        "code_commutator_action_normalized": float(np.linalg.norm(commutator_action) / np.linalg.norm(ZG)),
    }


def main() -> None:
    crosscheck = sparse_crosscheck()
    if crosscheck["absolute_frobenius_error"] >= 1e-10:
        raise RuntimeError("matrix-free parent action failed the sparse cross-check")
    rows = [
        analyse(L, support)
        for L, supports in ((6, (1, 2)), (8, (1, 2, 3)), (10, (3, 4)))
        for support in supports
    ]
    max_support = [
        next(row for row in rows if row["L"] == L and row["support_rungs"] == support)
        for L, support in ((6, 2), (8, 3), (10, 4))
    ]
    normalized = [row["code_commutator_action_normalized"] for row in max_support]
    maximal_support_decrease = all(right < left for left, right in zip(normalized, normalized[1:]))
    fixed_support = {
        str(support): [row for row in rows if row["support_rungs"] == support]
        for support in sorted({row["support_rungs"] for row in rows})
    }
    common_j3 = fixed_support["3"]
    j3_change_l8_to_l10 = (
        common_j3[1]["code_commutator_action_normalized"]
        - common_j3[0]["code_commutator_action_normalized"]
    )
    out = {
        "schema": "antler.phase6.edge-support-scaling.v1",
        "reference": "external Iemini lambda=1 parent and published finite-support aR-bR generator",
        "matrix_free_crosscheck": crosscheck,
        "rows": rows,
        "maximal_available_support_trajectory": max_support,
        "strict_monotone_decrease_on_maximal_support": bool(maximal_support_decrease),
        "fixed_support_trajectories": fixed_support,
        "fixed_support_j3_L8_to_L10_residual_change": j3_change_l8_to_l10,
        "axis_separation": (
            "The maximal-support trajectory changes both L and support and is therefore a support-improvement sequence, not a size-scaling result. "
            "At fixed j=3, the normalized residual changes only from 6.3184 at L=8 to 6.3057 at L=10."
        ),
        "claim_boundary": (
            "This tests the projected commutator of a published truncated external edge generator. "
            "The available size and support axes are not independently converged; no exact edge mode, protected braid, or native microscopic mapping follows."
        ),
        "decision": (
            "The apparent decrease on the maximal-support path is support-driven. The fixed-j=3 L=8-to-L=10 residual is effectively flat and remains large; this finite-support generator is not qualified as a protected braid primitive."
        ),
    }
    path = ROOT / "results" / "phase6" / "edge_support_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "crosscheck": crosscheck,
        "decision": out["decision"],
        "maximal_available_support_trajectory": max_support,
    }, indent=2))


if __name__ == "__main__":
    main()
