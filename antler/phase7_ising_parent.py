"""Exact fixed-point parent proposed in the Phase 7 analytic derivation.

This module is an independent transcription for auditing, not an endorsement of
the construction as a topological code.  It implements the OBC Hamiltonian

    H = U sum_j (q_j - 1)^2 + Delta sum_j n_d,j
        + J/2 sum_j (N_j N_{j+1} - X_j X_{j+1}),

on the weighted-charge basis used by the charge-two mediator ladder.  The
result is useful for testing the distinction between a symmetry-protected
Ising/cat doublet and a locally indistinguishable topological code.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.sparse import csr_matrix, dok_matrix


def mode_a(rung: int) -> int:
    return 2 * rung


def mode_b(rung: int) -> int:
    return 2 * rung + 1


def mode_d(L: int, bond: int) -> int:
    return 2 * L + bond


def build_weighted_basis_fast(L: int, total_charge: int) -> tuple[np.ndarray, dict[int, int]]:
    """Generate the weighted hard-core basis without enumerating all masks."""
    if L < 2 or total_charge < 0:
        raise ValueError("invalid L or total charge")
    rail_modes, bonds = 2 * L, L - 1
    states: list[int] = []
    for mediator_count in range(bonds + 1):
        rail_count = total_charge - 2 * mediator_count
        if not 0 <= rail_count <= rail_modes:
            continue
        for mediator_sites in combinations(range(bonds), mediator_count):
            mediator_mask = sum(1 << mode_d(L, bond) for bond in mediator_sites)
            for rail_sites in combinations(range(rail_modes), rail_count):
                states.append(mediator_mask + sum(1 << site for site in rail_sites))
    states_array = np.asarray(sorted(states), dtype=np.int64)
    return states_array, {int(state): row for row, state in enumerate(states_array)}


def rung_single_occupancy(state: int, rung: int) -> int:
    a, b = (state >> mode_a(rung)) & 1, (state >> mode_b(rung)) & 1
    return a + b - 2 * a * b


def cell_charge(state: int, L: int, rung: int) -> int:
    charge = ((state >> mode_a(rung)) & 1) + ((state >> mode_b(rung)) & 1)
    if rung > 0:
        charge += (state >> mode_d(L, rung - 1)) & 1
    if rung < L - 1:
        charge += (state >> mode_d(L, rung)) & 1
    return charge


def flip_rung(state: int, rung: int) -> int | None:
    """Apply X_rung; it is zero unless exactly one rail mode is occupied."""
    a, b = mode_a(rung), mode_b(rung)
    occupied = ((state >> a) & 1) + ((state >> b) & 1)
    return state ^ (1 << a) ^ (1 << b) if occupied == 1 else None


def local_x_matrix(L: int, states: np.ndarray, index: dict[int, int], rung: int, sparse: bool = True):
    """Number-conserving rail swap X_j in the fixed weighted-charge basis."""
    matrix = dok_matrix((len(states), len(states)), dtype=complex) if sparse else np.zeros((len(states), len(states)), complex)
    for column, raw_state in enumerate(states):
        target = flip_rung(int(raw_state), rung)
        if target is not None:
            matrix[index[target], column] = 1.0
    return matrix.tocsr() if sparse else matrix


def local_z_matrix(L: int, states: np.ndarray, rung: int, sparse: bool = True):
    values = np.asarray([
        ((int(state) >> mode_a(rung)) & 1) - ((int(state) >> mode_b(rung)) & 1)
        for state in states
    ], dtype=float)
    return csr_matrix(np.diag(values)) if sparse else np.diag(values.astype(complex))


def local_na_matrix(L: int, states: np.ndarray, rung: int, sparse: bool = True):
    values = np.asarray([(int(state) >> mode_a(rung)) & 1 for state in states], dtype=float)
    return csr_matrix(np.diag(values)) if sparse else np.diag(values.astype(complex))


def local_cell_constraint_matrix(L: int, states: np.ndarray, rung: int, sparse: bool = True):
    values = np.asarray([(cell_charge(int(state), L, rung) - 1) ** 2 for state in states], dtype=float)
    return csr_matrix(np.diag(values)) if sparse else np.diag(values.astype(complex))


def local_mediator_number_matrix(L: int, states: np.ndarray, bond: int, sparse: bool = True):
    values = np.asarray([(int(state) >> mode_d(L, bond)) & 1 for state in states], dtype=float)
    return csr_matrix(np.diag(values)) if sparse else np.diag(values.astype(complex))


def local_bond_projector_matrix(
    L: int, states: np.ndarray, index: dict[int, int], bond: int, sparse: bool = True,
):
    """Pi^B = (N_j N_{j+1} - X_j X_{j+1}) / 2 on one OBC bond."""
    matrix = dok_matrix((len(states), len(states)), dtype=complex) if sparse else np.zeros((len(states), len(states)), complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        if rung_single_occupancy(state, bond) and rung_single_occupancy(state, bond + 1):
            matrix[column, column] += 0.5
            first = flip_rung(state, bond)
            target = flip_rung(first, bond + 1) if first is not None else None
            if target is None:
                raise RuntimeError("single-occupancy bond flip unexpectedly vanished")
            matrix[index[target], column] += -0.5
    return matrix.tocsr() if sparse else matrix


def build_fixed_parent(
    L: int, total_charge: int, U: float = 4.0, Delta: float = 2.0, J: float = 1.0,
    basis: tuple[np.ndarray, dict[int, int]] | None = None, sparse: bool = True,
):
    """Build the OBC fixed-point parent for an independently chosen Q sector."""
    states, index = build_weighted_basis_fast(L, total_charge) if basis is None else basis
    dimension = len(states)
    H = dok_matrix((dimension, dimension), dtype=complex) if sparse else np.zeros((dimension, dimension), complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        diagonal = U * sum((cell_charge(state, L, rung) - 1) ** 2 for rung in range(L))
        diagonal += Delta * sum((state >> mode_d(L, bond)) & 1 for bond in range(L - 1))
        for bond in range(L - 1):
            if rung_single_occupancy(state, bond) and rung_single_occupancy(state, bond + 1):
                diagonal += 0.5 * J
                first = flip_rung(state, bond)
                target = flip_rung(first, bond + 1) if first is not None else None
                if target is None:
                    raise RuntimeError("single-occupancy bond flip unexpectedly vanished")
                H[index[target], column] += -0.5 * J
        H[column, column] += diagonal
    return H.tocsr() if sparse else H, states, index


def code_frame(L: int, states: np.ndarray, index: dict[int, int]) -> np.ndarray:
    """Return the two product X-eigenstate ground vectors Omega_+, Omega_-."""
    frame = np.zeros((len(states), 2), dtype=complex)
    normalization = 1.0 / np.sqrt(2 ** L)
    for selection in range(1 << L):
        state = 0
        b_count = 0
        for rung in range(L):
            if (selection >> rung) & 1:
                state |= 1 << mode_b(rung)
                b_count += 1
            else:
                state |= 1 << mode_a(rung)
        frame[index[state], 0] = normalization
        frame[index[state], 1] = normalization * (-1.0 if b_count & 1 else 1.0)
    return frame


def parity_a_labels(L: int, states: np.ndarray) -> np.ndarray:
    return np.asarray([
        -1.0 if sum((int(state) >> mode_a(rung)) & 1 for rung in range(L)) & 1 else 1.0
        for state in states
    ])
