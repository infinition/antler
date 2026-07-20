"""Synthetic-dimension extension of the ANTLER correlated-hopping ladder.

This module is additive: it does not alter the frozen scalar model in
``antler.model``.  Each rung-major spatial site carries a two-component
internal state.  The Jordan--Wigner string counts *spatial* occupations, while
each directed hopping link can carry a 2x2 unitary matrix.
"""
from __future__ import annotations

from itertools import combinations, product
from typing import Mapping

import numpy as np
from scipy.sparse import csr_matrix, diags

from .model import hop_list


SPIN_DIM = 2


def mode(site: int, spin: int) -> int:
    return SPIN_DIM * site + spin


def site_of(mode_index: int) -> int:
    return mode_index // SPIN_DIM


def build_spin_hardcore_basis(n_sites: int, n_particles: int):
    """Basis with at most one particle per *spatial* site.

    A particle may occupy either internal component, but the hard-core
    constraint remains the one used by ANTLER: two particles cannot share a
    rung/leg spatial orbital even when their pseudo-spins differ.
    """

    states = []
    for spatial_sites in combinations(range(n_sites), n_particles):
        for spins in product(range(SPIN_DIM), repeat=n_particles):
            state = 0
            for site, spin in zip(spatial_sites, spins):
                state |= 1 << mode(site, spin)
            states.append(state)
    states.sort()
    array = np.asarray(states, dtype=np.int64)
    return array, {int(state): i for i, state in enumerate(array)}


def occupied_spin(state: int, site: int) -> int | None:
    for spin in range(SPIN_DIM):
        if (state >> mode(site, spin)) & 1:
            return spin
    return None


def spatial_occupation_between(state: int, left: int, right: int) -> int:
    return sum(occupied_spin(state, site) is not None for site in range(left + 1, right))


def su2(axis: np.ndarray, angle: float) -> np.ndarray:
    """Return exp(-i angle axis/2) for a Pauli axis matrix."""

    eigenvalues, eigenvectors = np.linalg.eigh(np.asarray(axis, complex))
    return (eigenvectors * np.exp(-0.5j * angle * eigenvalues)) @ eigenvectors.conj().T


def build_synthetic_hamiltonian(
    L: int,
    N: int,
    theta: float,
    J1: float,
    J2: float,
    Jperp: float,
    mu: np.ndarray | None = None,
    link_matrices: Mapping[tuple[int, int], np.ndarray] | None = None,
    onsite_spin_matrices: Mapping[int, np.ndarray] | None = None,
    basis=None,
) -> tuple[csr_matrix, np.ndarray, dict[int, int]]:
    """Build the matrix-link correlated-hopping Hamiltonian.

    For a canonical spatial link ``k < l``, ``link_matrices[(k,l)]`` acts on
    the spinor for the decreasing-index hop ``l -> k``.  The reverse hop is
    added by Hermitian conjugation.  With every matrix equal to identity this
    is the spin-degenerate version of the rung-major correlated-hopping rule.
    """

    n_sites = 2 * L
    if basis is None:
        states, index = build_spin_hardcore_basis(n_sites, N)
    else:
        states, index = basis
    if mu is None:
        mu_array = np.zeros(n_sites)
    else:
        mu_array = np.asarray(mu, dtype=float)
        if mu_array.shape != (n_sites,):
            raise ValueError("mu must contain one value per spatial ladder site")
    matrices = link_matrices or {}
    onsite = onsite_spin_matrices or {}
    rows: list[int] = []
    cols: list[int] = []
    values: list[complex] = []
    onsite_rows: list[int] = []
    onsite_cols: list[int] = []
    onsite_values: list[complex] = []
    diagonal = np.zeros(len(states), dtype=float)
    identity = np.eye(SPIN_DIM, dtype=complex)
    for col, raw_state in enumerate(states):
        state = int(raw_state)
        for site in range(n_sites):
            spin_here = occupied_spin(state, site)
            if spin_here is not None:
                diagonal[col] += mu_array[site]
                W = np.asarray(onsite.get(site, np.zeros((SPIN_DIM, SPIN_DIM))), dtype=complex)
                if W.shape != (SPIN_DIM, SPIN_DIM):
                    raise ValueError(f"onsite matrix at {site} is not a 2x2 matrix")
                for spin_target in range(SPIN_DIM):
                    amplitude = W[spin_target, spin_here]
                    if abs(amplitude) == 0:
                        continue
                    new_state = state ^ (1 << mode(site, spin_here)) ^ (1 << mode(site, spin_target))
                    onsite_rows.append(index[new_state])
                    onsite_cols.append(col)
                    onsite_values.append(amplitude)
        for k, l, J in hop_list(L, J1, J2, Jperp):
            spin_source = occupied_spin(state, l)
            if spin_source is None or occupied_spin(state, k) is not None:
                continue
            U = np.asarray(matrices.get((k, l), identity), dtype=complex)
            if U.shape != (SPIN_DIM, SPIN_DIM):
                raise ValueError(f"link ({k},{l}) is not a 2x2 matrix")
            n_mid = spatial_occupation_between(state, k, l)
            phase = np.exp(1j * theta * n_mid)
            removed = state ^ (1 << mode(l, spin_source))
            for spin_target in range(SPIN_DIM):
                amplitude = -J * phase * U[spin_target, spin_source]
                if abs(amplitude) == 0:
                    continue
                new_state = removed | (1 << mode(k, spin_target))
                rows.append(index[new_state])
                cols.append(col)
                values.append(amplitude)
    one = csr_matrix((values, (rows, cols)), shape=(len(states), len(states)))
    onsite_term = csr_matrix((onsite_values, (onsite_rows, onsite_cols)),
                             shape=(len(states), len(states)))
    return one + one.conj().T + diags(diagonal) + onsite_term, states, index
