"""Microscopic charge-two mediator realization of an external parent target.

Three hard-core charge-two mediators live on every bond.  Their conversion
channels are the rank-three Gram factorization of the local interaction in the
published Iemini number-conserving parent Hamiltonian.  This is a bridge to an
external benchmark, not a claim of a new ANTLER phase or a change to frozen
ANTLER.
"""
from __future__ import annotations

import numpy as np

from .basis import popcount


def mode_a(rung: int) -> int:
    return 2 * rung


def mode_b(rung: int) -> int:
    return 2 * rung + 1


def mode_d(L: int, bond: int, channel: int) -> int:
    if channel not in (0, 1, 2):
        raise ValueError("channel must be 0, 1, or 2")
    return 2 * L + 3 * bond + channel


def build_weighted_basis(L: int, total_charge: int) -> tuple[np.ndarray, dict[int, int]]:
    """Fixed charge basis: low fermions have charge one, mediators charge two."""
    if L < 2 or total_charge < 0:
        raise ValueError("invalid L or total_charge")
    n_modes = 5 * L - 3
    states = np.asarray([
        state for state in range(1 << n_modes)
        if popcount(state & ((1 << (2 * L)) - 1))
        + 2 * popcount(state >> (2 * L)) == total_charge
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def _annihilate_low_fermion(state: int, site: int) -> tuple[int, float] | None:
    if not ((state >> site) & 1):
        return None
    return state ^ (1 << site), -1.0 if popcount(state & ((1 << site) - 1)) & 1 else 1.0


def _create_low_fermion(state: int, site: int) -> tuple[int, float] | None:
    if (state >> site) & 1:
        return None
    return state | (1 << site), -1.0 if popcount(state & ((1 << site) - 1)) & 1 else 1.0


def _pair_annihilation(state: int, first: int, second: int) -> tuple[int, float] | None:
    """Apply ``c_high c_low`` to the canonical pair with `first < second`."""
    first, second = sorted((first, second))
    item = _annihilate_low_fermion(state, first)
    if item is None:
        return None
    current, amplitude = item
    item = _annihilate_low_fermion(current, second)
    if item is None:
        return None
    current, sign = item
    return current, amplitude * sign


def bond_channels(L: int, bond: int) -> tuple[tuple[tuple[int, int, float], ...], ...]:
    """The compact rank-three factorization channels on one bond."""
    a0, b0 = mode_a(bond), mode_b(bond)
    a1, b1 = mode_a(bond + 1), mode_b(bond + 1)
    root2 = float(np.sqrt(2.0))
    return (
        ((a0, a1, 2.0), (b0, b1, 2.0)),
        ((b0, a1, root2), (a0, b1, root2)),
        ((a0, b0, root2), (a1, b1, -root2)),
    )


def build_multiplet_mediator_parent(
    L: int,
    total_charge: int,
    Delta: float,
    basis=None,
):
    """Build the explicit-mediator Hamiltonian for the external parent bridge.

    The low-mode one-body block is the ``lambda=0`` Iemini bond Hamiltonian.
    Mediator conversion amplitudes are ``sqrt(Delta) C``.  At second order in
    the low-energy subspace, each bond adds ``-C^dagger C``, exactly the
    ``lambda=1 - lambda=0`` interaction audited in Phase 6F.
    """
    if Delta <= 0:
        raise ValueError("Delta must be positive")
    states, index = build_weighted_basis(L, total_charge) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)

    def occupied(state: int, site: int) -> int:
        return (state >> site) & 1

    def add_fermion_hop(column: int, state: int, source: int, destination: int, coefficient: float) -> None:
        item = _annihilate_low_fermion(state, source)
        if item is None:
            return
        current, amplitude = item
        item = _create_low_fermion(current, destination)
        if item is None:
            return
        new_state, sign = item
        H[index[new_state], column] += coefficient * amplitude * sign

    conversion_scale = float(np.sqrt(Delta))
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for bond in range(L - 1):
            a0, b0 = mode_a(bond), mode_b(bond)
            a1, b1 = mode_a(bond + 1), mode_b(bond + 1)
            for left, right in ((a0, a1), (b0, b1)):
                # The exactly reproduced lambda=0 one-body bond block.
                add_fermion_hop(column, state, right, left, -4.0)
                add_fermion_hop(column, state, left, right, -4.0)
                H[column, column] += 4.0 * (occupied(state, left) + occupied(state, right))
            for channel, pairs in enumerate(bond_channels(L, bond)):
                mediator = mode_d(L, bond, channel)
                H[column, column] += Delta * occupied(state, mediator)
                for first, second, coefficient in pairs:
                    item = _pair_annihilation(state, first, second)
                    if item is not None and not occupied(state, mediator):
                        low_state, sign = item
                        new_state = low_state | (1 << mediator)
                        H[index[new_state], column] += -conversion_scale * coefficient * sign
                    if occupied(state, mediator) and not occupied(state, first) and not occupied(state, second):
                        low_state = state ^ (1 << mediator) | (1 << first) | (1 << second)
                        item = _pair_annihilation(low_state, first, second)
                        if item is None:
                            raise RuntimeError("pair creation sign lookup failed")
                        _, sign = item
                        H[index[low_state], column] += -conversion_scale * coefficient * sign
    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("multiplet mediator parent is not Hermitian")
    return H, states, index


def generalized_branch_parities(state: int, L: int) -> tuple[int, int]:
    """Exact Z2 charges, including the parity carried by mixed-pair mediators."""
    n_a = sum((state >> mode_a(rung)) & 1 for rung in range(L))
    n_b = sum((state >> mode_b(rung)) & 1 for rung in range(L))
    mixed_mediators = sum(
        ((state >> mode_d(L, bond, channel)) & 1)
        for bond in range(L - 1) for channel in (1, 2)
    )
    return (
        1 if (n_a + mixed_mediators) % 2 == 0 else -1,
        1 if (n_b + mixed_mediators) % 2 == 0 else -1,
    )


def mediator_number(state: int, L: int) -> int:
    """Number of occupied hard-core charge-two mediators."""
    return popcount(state >> (2 * L))
