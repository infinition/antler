"""Independent audit of the OBC fixed-point parent supplied for Phase 7.

This run tests the analytic parent before attempting its noncommuting
microscopic realization.  In particular, it distinguishes exact symmetry-
restricted Ising protection from full local indistinguishability and checks
whether every advertised ``projector`` is actually idempotent.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.native_charge2_ladder import build_weighted_basis
from antler.phase7_ising_parent import (
    build_fixed_parent,
    build_weighted_basis_fast,
    code_frame,
    local_bond_projector_matrix,
    local_cell_constraint_matrix,
    local_mediator_number_matrix,
    local_na_matrix,
    local_x_matrix,
    local_z_matrix,
    parity_a_labels,
)
from antler.phase7_parent_audit import (
    local_indistinguishability,
    local_projector_algebra,
    projected_edge_metrics,
    symmetry_audit,
)


U, DELTA, J = 4.0, 2.0, 1.0


def low_spectrum(L: int, charge: int, count: int) -> list[float]:
    H, _, _ = build_fixed_parent(L, charge, U=U, Delta=DELTA, J=J, sparse=True)
    values = eigsh(H, k=count, which="SA", return_eigenvectors=False, tol=1e-11)
    return [float(value) for value in np.sort(values)]


def spectral_row(L: int) -> dict:
    neutral = low_spectrum(L, L, count=4)
    e_minus = low_spectrum(L, L - 1, count=1)[0]
    e_plus = low_spectrum(L, L + 1, count=1)[0]
    ground = 0.5 * (neutral[0] + neutral[1])
    return {
        "L": L,
        "Q": L,
        "lowest_neutral_energies": neutral,
        "logical_splitting": float(neutral[1] - neutral[0]),
        "neutral_gap_above_doublet": float(neutral[2] - neutral[1]),
        "ground_energy_average": float(ground),
        "addition_energy": float(e_plus - ground),
        "removal_energy": float(e_minus - ground),
    }


def encoded_matrix(matrix: np.ndarray) -> dict:
    return {"real": matrix.real.tolist(), "imag": matrix.imag.tolist()}


def main() -> None:
    L = 4
    states, index = build_weighted_basis_fast(L, L)
    reference_states, _ = build_weighted_basis(L, L)
    if not np.array_equal(states, reference_states):
        raise RuntimeError("fast weighted basis disagrees with the established Phase 6 basis")

    H, _, _ = build_fixed_parent(L, L, U=U, Delta=DELTA, J=J, basis=(states, index), sparse=False)
    G = code_frame(L, states, index)
    cell_terms = [local_cell_constraint_matrix(L, states, rung, sparse=False) for rung in range(L)]
    mediator_terms = [local_mediator_number_matrix(L, states, bond, sparse=False) for bond in range(L - 1)]
    bond_terms = [local_bond_projector_matrix(L, states, index, bond, sparse=False) for bond in range(L - 1)]
    terms = cell_terms + mediator_terms + bond_terms
    cell_supports = []
    for rung in range(L):
        support = {2 * rung, 2 * rung + 1}
        support.update(2 * L + bond for bond in (rung - 1, rung) if 0 <= bond < L - 1)
        cell_supports.append(tuple(sorted(support)))
    supports = cell_supports + [
        (2 * L + bond,) for bond in range(L - 1)
    ] + [
        (2 * bond, 2 * bond + 1, 2 * (bond + 1), 2 * (bond + 1) + 1)
        for bond in range(L - 1)
    ]
    names = [f"C_{rung}" for rung in range(L)] + [f"n_d_{bond}" for bond in range(L - 1)] + [f"Pi_B_{bond}" for bond in range(L - 1)]
    algebra = local_projector_algebra(terms, supports, names)
    pa = parity_a_labels(L, states)
    pb = np.asarray([
        -1.0 if sum((int(state) >> (2 * rung + 1)) & 1 for rung in range(L)) & 1 else 1.0
        for state in states
    ])
    symmetry = symmetry_audit(H, {"Q_fixed_sector": np.full(len(states), L), "P_a": pa, "P_b": pb})
    X0 = local_x_matrix(L, states, index, 0, sparse=False)
    Xbulk = local_x_matrix(L, states, index, 2, sparse=False)
    Z0 = local_z_matrix(L, states, 0, sparse=False)
    n_a0 = local_na_matrix(L, states, 0, sparse=False)
    edge = {
        "left_X_0": projected_edge_metrics(H, G, X0),
        "bulk_X_2": projected_edge_metrics(H, G, Xbulk),
    }
    all_local = local_indistinguishability(G, {"X_0": X0, "Z_0": Z0, "n_a_0": n_a0})
    symmetry_preserving = local_indistinguishability(G, {"Z_0": Z0, "n_a_0": n_a0})
    projected_x = {
        "P_X_0_P": encoded_matrix(G.conj().T @ X0 @ G),
        "P_X_2_P": encoded_matrix(G.conj().T @ Xbulk @ G),
        "P_Z_0_P": encoded_matrix(G.conj().T @ Z0 @ G),
    }
    rows = [spectral_row(length) for length in (4, 6, 8)]
    code_residual = float(np.linalg.norm(H @ G))
    classification = (
        "The fixed-point OBC parent has the advertised finite-size doublet and flat gaps, but fails full local indistinguishability: X_j is a local logical operator at both edge and bulk. "
        "It is therefore an Ising symmetry-restricted cat/SSB benchmark, not a localized-edge or intrinsically topological code."
    )
    out = {
        "schema": "antler.phase7.fixed-point-parent-independent-audit.v1",
        "source": "docs/PHASE7_STATE_TO_HAMILTONIAN_DERIVATION.md, OBC fixed-point H_fix only",
        "parameters": {"U": U, "Delta": DELTA, "J": J},
        "basis_crosscheck_against_phase6": {"L": L, "Q": L, "exact_match": True, "dimension": len(states)},
        "code_frame_parent_residual_frobenius": code_residual,
        "local_term_algebra": algebra,
        "symmetry_audit_in_fixed_Q_sector": symmetry,
        "edge_vs_bulk_X_metrics": edge,
        "full_local_probe_indistinguishability": all_local,
        "symmetry_preserving_local_probe_indistinguishability": symmetry_preserving,
        "projected_local_matrices": projected_x,
        "open_boundary_spectra": rows,
        "periodic_boundary_status": (
            "not audited: the supplied OBC mode definition has L-1 mediators and does not specify the extra wraparound mediator/cell constraint required for a PBC Hamiltonian"
        ),
        "decision": classification,
        "claim_boundary": (
            "This independently audits the exactly soluble OBC parent, not the Schrieffer-Wolff microscopic Hamiltonian. "
            "The run does not test a braid, a non-Abelian phase, or an experimentally derived native realization."
        ),
    }
    path = ROOT / "results" / "phase7" / "fixed_point_parent_independent_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
