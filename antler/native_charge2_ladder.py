"""Exact-branch-parity ladder with explicit charge-two bond mediators.

This is a Phase 6 candidate extension.  Each bond has a hard-core charge-two
mode ``d_j`` which converts to a nearest-neighbour pair on rail ``a`` or
``b``.  It conserves total weighted U(1) charge and both rail parities exactly.
It is not a modification of the frozen correlated-Jordan--Wigner ladder.
"""
from __future__ import annotations

import numpy as np


def mode_a(rung: int) -> int:
    return 2 * rung


def mode_b(rung: int) -> int:
    return 2 * rung + 1


def mode_d(L: int, bond: int) -> int:
    return 2 * L + bond


def build_weighted_basis(L: int, total_charge: int) -> tuple[np.ndarray, dict[int, int]]:
    """Basis with unit charge on `a,b` and charge two on each bond mediator."""
    if L < 2 or total_charge < 0:
        raise ValueError("invalid L or total_charge")
    n_modes = 3 * L - 1
    states = np.asarray([
        state for state in range(1 << n_modes)
        if sum((state >> site) & 1 for site in range(2 * L))
        + 2 * sum((state >> mode_d(L, bond)) & 1 for bond in range(L - 1)) == total_charge
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_charge2_mediator_ladder(
    L: int,
    total_charge: int,
    t_leg: float,
    g: float,
    Delta: float,
    V_rung: float = 0.0,
    V_leg: float = 0.0,
    basis=None,
):
    """Build a two-rail ladder with explicit charge-two bond mediators.

    The conversion on bond `j` is
    ``-g d_j^dagger (a_j a_{j+1} + b_j b_{j+1}) + h.c.``.  Ordinary hopping
    stays within a rail, so `(-1)^N_a` and `(-1)^N_b` commute with the full
    microscopic Hamiltonian.  `V_rung` and `V_leg` are scalar density terms.
    """
    states, index = build_weighted_basis(L, total_charge) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)

    def occ(state: int, site: int) -> int:
        return (state >> site) & 1

    def add_hop(column: int, state: int, left: int, right: int, amplitude: float) -> None:
        if occ(state, right) and not occ(state, left):
            H[index[state ^ (1 << left) ^ (1 << right)], column] += amplitude
        if occ(state, left) and not occ(state, right):
            H[index[state ^ (1 << left) ^ (1 << right)], column] += amplitude

    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for rung in range(L):
            H[column, column] += V_rung * occ(state, mode_a(rung)) * occ(state, mode_b(rung))
        for bond in range(L - 1):
            a0, a1 = mode_a(bond), mode_a(bond + 1)
            b0, b1 = mode_b(bond), mode_b(bond + 1)
            mediator = mode_d(L, bond)
            H[column, column] += Delta * occ(state, mediator)
            H[column, column] += V_leg * (
                occ(state, a0) * occ(state, a1) + occ(state, b0) * occ(state, b1)
            )
            add_hop(column, state, a0, a1, -t_leg)
            add_hop(column, state, b0, b1, -t_leg)
            for first, second in ((a0, a1), (b0, b1)):
                if occ(state, first) and occ(state, second) and not occ(state, mediator):
                    new_state = state ^ (1 << first) ^ (1 << second) ^ (1 << mediator)
                    H[index[new_state], column] += -g
                if occ(state, mediator) and not occ(state, first) and not occ(state, second):
                    new_state = state ^ (1 << first) ^ (1 << second) ^ (1 << mediator)
                    H[index[new_state], column] += -g
    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("charge-two mediator ladder is not Hermitian")
    return H, states, index


def branch_parities(state: int, L: int) -> tuple[int, int]:
    """Eigenvalues of the two exact rail-parity symmetries."""
    n_a = sum((state >> mode_a(rung)) & 1 for rung in range(L))
    n_b = sum((state >> mode_b(rung)) & 1 for rung in range(L))
    return (1 if n_a % 2 == 0 else -1, 1 if n_b % 2 == 0 else -1)


def local_density(states: np.ndarray, site: int) -> np.ndarray:
    """Diagonal hard-core density of any physical ladder mode."""
    return np.asarray([(int(state) >> site) & 1 for state in states], dtype=float)
