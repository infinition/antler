"""Fixed-density, even-N sparse-Lanczos test of the Phase 8 Floquet candidate.

The previous dense preflight could not continue nu=1/4 with even N beyond
L=4: the next point is L=8,N=4.  This script constructs the published
high-frequency H_eff exactly in the fixed-N Fock sector, using a sparse global
rail rotation P, validates that representation against the existing dense
L=6,N=4 implementation, then measures the L=8 candidate U0=-2, alpha=1/2.

This is still finite-size external-model evidence, not an ANTLER-native or
thermodynamic claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh, norm as sparse_norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply, wire_a_parity
from run_phase7d_floquet_full_ladder_preflight import ETA, build_h0_and_rotation


T_HOP, U0, ALPHA = 1.0, -2.0, 0.5
VALIDATION_CASE, TARGET_CASES = (6, 4), ((4, 2), (8, 4))


def sparse_from_entries(entries: dict[tuple[int, int], complex], dimension: int) -> csr_matrix:
    rows, columns, values = zip(*((row, column, value) for (row, column), value in entries.items() if abs(value) > 1e-15))
    return csr_matrix((values, (rows, columns)), shape=(dimension, dimension), dtype=complex)


def build_h0_sparse(L: int, N: int, u0: float) -> tuple[csr_matrix, np.ndarray, dict[int, int]]:
    states, index = build_basis(2 * L, N)
    entries: dict[tuple[int, int], complex] = {}
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        diagonal = 0.0
        for rung in range(L - 1):
            for rail in (0, 1):
                left, right = site_index(rung, rail), site_index(rung + 1, rail)
                for operations in (
                    (("ann", right), ("create", left)),
                    (("ann", left), ("create", right)),
                ):
                    item = _apply(state, operations)
                    if item is not None:
                        new, amplitude = item
                        key = (index[new], column)
                        entries[key] = entries.get(key, 0.0) - T_HOP * amplitude
                diagonal += u0 * (((state >> left) & 1) * ((state >> right) & 1))
        entries[(column, column)] = entries.get((column, column), 0.0) + diagonal
    return sparse_from_entries(entries, len(states)), states, index


def global_rail_rotation_sparse(L: int, states: np.ndarray, index: dict[int, int], eta: float) -> csr_matrix:
    """P=prod_j exp[-i eta (a^dag b+b^dag a)/2] in rung-major convention."""
    cos, sin = np.cos(eta / 2.0), -1j * np.sin(eta / 2.0)
    entries: dict[tuple[int, int], complex] = {}
    for column, raw_state in enumerate(states):
        amplitudes = {int(raw_state): 1.0 + 0.0j}
        for rung in range(L):
            a, b = site_index(rung, 0), site_index(rung, 1)
            evolved: dict[int, complex] = {}
            for state, amplitude in amplitudes.items():
                occupied_a, occupied_b = (state >> a) & 1, (state >> b) & 1
                if occupied_a ^ occupied_b:
                    evolved[state] = evolved.get(state, 0.0) + cos * amplitude
                    flipped = state ^ (1 << a) ^ (1 << b)
                    evolved[flipped] = evolved.get(flipped, 0.0) + sin * amplitude
                else:
                    evolved[state] = evolved.get(state, 0.0) + amplitude
            amplitudes = evolved
        for state, amplitude in amplitudes.items():
            entries[(index[state], column)] = amplitude
    return sparse_from_entries(entries, len(states))


def effective_sparse(L: int, N: int, u0: float, alpha: float) -> tuple[csr_matrix, np.ndarray, dict[int, int], csr_matrix]:
    h0, states, index = build_h0_sparse(L, N, u0)
    rotation = global_rail_rotation_sparse(L, states, index, ETA)
    effective = (alpha * h0 + (1.0 - alpha) * (rotation.conj().T @ h0 @ rotation)).tocsr()
    return effective, states, index, rotation


def central_cut_spectrum(vector: np.ndarray, states: np.ndarray, L: int) -> dict:
    left_modes = 2 * (L // 2)
    amplitude = np.zeros((1 << left_modes, 1 << (2 * L - left_modes)), dtype=complex)
    for coefficient, raw_state in zip(vector, states):
        state = int(raw_state)
        amplitude[state & ((1 << left_modes) - 1), state >> left_modes] = coefficient
    weights = np.linalg.eigvalsh(amplitude @ amplitude.conj().T)
    weights = np.sort(weights[weights > 1e-13])[::-1]
    levels = -np.log(weights)
    return {
        "schmidt_weights_descending": [float(value) for value in weights],
        "leading_entanglement_levels": [float(value) for value in levels[:12]],
        "leading_adjacent_pair_splittings": [float(abs(levels[i] - levels[i + 1])) for i in range(0, min(len(levels) - 1, 12), 2)],
        "von_neumann_entropy": float(-np.sum(weights * np.log(weights))),
    }


def two_lowest(matrix: csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    if matrix.shape[0] <= 4:
        values, vectors = eigh(matrix.toarray())
        return values[:2], vectors[:, :2]
    values, vectors = eigsh(matrix, k=2, which="SA", tol=1e-11, maxiter=200_000)
    order = np.argsort(values)
    return values[order], vectors[:, order]


def analyze(L: int, N: int, u0: float = U0, alpha: float = ALPHA) -> dict:
    effective, states, _, rotation = effective_sparse(L, N, u0, alpha)
    parity_values = np.asarray([(-1.0 if wire_a_parity(int(state), L) else 1.0) for state in states])
    parity = diags(parity_values, format="csr")
    sectors = {}
    for sector in (0, 1):
        rows = np.asarray([i for i, state in enumerate(states) if wire_a_parity(int(state), L) == sector], dtype=int)
        values, vectors = two_lowest(effective[rows][:, rows].tocsr())
        full_ground = np.zeros(len(states), dtype=complex)
        full_ground[rows] = vectors[:, 0]
        sectors[str(sector)] = {
            "ground_energy": float(values[0]),
            "first_neutral_excitation": float(values[1]),
            "neutral_gap": float(values[1] - values[0]),
            "central_cut": central_cut_spectrum(full_ground, states, L),
        }
    split = abs(sectors["0"]["ground_energy"] - sectors["1"]["ground_energy"])
    gap = min(sectors["0"]["neutral_gap"], sectors["1"]["neutral_gap"])
    identity = csr_matrix(np.eye(len(states), dtype=complex))
    return {
        "L": L,
        "N": N,
        "filling_N_over_2L": N / (2.0 * L),
        "parameters": {"u0_attractive_nn": u0, "alpha": alpha, "eta": ETA, "t_leg": T_HOP},
        "hilbert_dimension": len(states),
        "rotation_unitarity_frobenius": float(sparse_norm(rotation.conj().T @ rotation - identity)),
        "parity_commutator_frobenius": float(sparse_norm(effective @ parity - parity @ effective)),
        "parity_sector_ground_split": float(split),
        "smallest_neutral_gap": float(gap),
        "split_over_gap": float(split / gap),
        "parity_sectors": sectors,
    }


def validate_sparse_rotation() -> dict:
    L, N = VALIDATION_CASE
    sparse_effective, _, _, _ = effective_sparse(L, N, U0, ALPHA)
    dense_h0, dense_rotation, _, _ = build_h0_and_rotation(L, N, U0)
    dense_effective = ALPHA * dense_h0 + (1.0 - ALPHA) * (dense_rotation.conj().T @ dense_h0 @ dense_rotation)
    return {
        "L": L,
        "N": N,
        "sparse_vs_existing_dense_frobenius": float(np.linalg.norm(sparse_effective.toarray() - dense_effective)),
        "threshold": 1e-11,
    }


def main() -> None:
    validation = validate_sparse_rotation()
    if validation["sparse_vs_existing_dense_frobenius"] > validation["threshold"]:
        raise RuntimeError(f"Sparse Floquet representation failed validation: {validation}")
    targets = [analyze(*case) for case in TARGET_CASES]
    controls = {
        "no_floquet_mixing_alpha_one": analyze(8, 4, U0, 1.0),
        "no_interaction_u0_zero": analyze(8, 4, 0.0, ALPHA),
    }
    out = {
        "schema": "antler.phase8.nc-majorana-l8-sparse-audit.v1",
        "citation": "Defossez et al., arXiv:2412.14886v2 (2025)",
        "parameters": {"u0_attractive_nn": U0, "alpha": ALPHA, "eta": ETA, "t_leg": T_HOP},
        "validation": validation,
        "constant_density_targets": targets,
        "registered_negative_controls": controls,
        "decision": (
            "L=8 fixed-density even-N finite-size audit only. It can support or reject this registered candidate at a larger "
            "exact size, but cannot establish a thermodynamic topological phase without subsequent MPS/DMRG scaling."
        ),
        "claim_boundary": (
            "This is an external Floquet effective model, not a finite-pulse ANTLER derivation. It establishes no ANTLER-native "
            "phase, edge qubit, braid, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "nc_majorana_l8_sparse_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
