"""Charge-two mediator candidate for an exactly parity-preserving extension.

The elementary local interaction is atom--molecule conversion,
``d^dagger a0 a1 + d^dagger b0 b1 + h.c.``, where the hard-core mediator
``d`` carries U(1) charge two.  Unlike a network of single-particle mediator
hops, this interaction preserves each low-rail fermion/boson parity exactly.
It is a candidate extension with an additional charge-two degree of freedom;
it is not part of the frozen ANTLER model.
"""
from __future__ import annotations

import numpy as np


MODE_NAMES = ("a0", "a1", "b0", "b1", "d")
MODE_CHARGES = (1, 1, 1, 1, 2)


def build_charge_basis(total_charge: int) -> tuple[np.ndarray, dict[int, int]]:
    """Fixed total-U(1)-charge basis; the molecule has charge two."""
    if total_charge < 0:
        raise ValueError("total_charge must be nonnegative")
    states = np.asarray([
        state for state in range(1 << len(MODE_NAMES))
        if sum(((state >> site) & 1) * MODE_CHARGES[site] for site in range(len(MODE_NAMES))) == total_charge
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_charge2_mediator_block(
    total_charge: int,
    g: float,
    Delta: float = 5.0,
    V_mixed: float = 3.0,
    basis=None,
):
    """Build the minimal atom--molecule conversion block.

    ``V_mixed N_a N_b`` makes the two same-rail pairs the low manifold in the
    charge-two sector.  Conversion through the high-energy charge-two mode
    produces an effective pair transfer but no single-particle transfer.
    """
    if Delta <= 0:
        raise ValueError("Delta must be positive")
    states, index = build_charge_basis(total_charge) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)

    def occ(state: int, site: int) -> int:
        return (state >> site) & 1

    for column, raw_state in enumerate(states):
        state = int(raw_state)
        n_a = occ(state, 0) + occ(state, 1)
        n_b = occ(state, 2) + occ(state, 3)
        H[column, column] += Delta * occ(state, 4) + V_mixed * n_a * n_b
        for first, second in ((0, 1), (2, 3)):
            pair_present = occ(state, first) and occ(state, second)
            molecule_present = occ(state, 4)
            if pair_present and not molecule_present:
                new_state = state ^ (1 << first) ^ (1 << second) ^ (1 << 4)
                H[index[new_state], column] += -g
            if molecule_present and not occ(state, first) and not occ(state, second):
                new_state = state ^ (1 << first) ^ (1 << second) ^ (1 << 4)
                H[index[new_state], column] += -g
    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("charge-two mediator Hamiltonian is not Hermitian")
    return H, states, index


def branch_parity_operator(states: np.ndarray, branch: str) -> np.ndarray:
    """Exact low-rail parity; the molecule is parity even on both branches."""
    if branch not in {"a", "b"}:
        raise ValueError("branch must be 'a' or 'b'")
    sites = (0, 1) if branch == "a" else (2, 3)
    return np.diag([
        -1.0 if sum((int(state) >> site) & 1 for site in sites) % 2 else 1.0
        for state in states
    ])


def low_pair_masks() -> tuple[int, int]:
    """Masks of the low `aa` and `bb` pair configurations."""
    return (1 << 0) | (1 << 1), (1 << 2) | (1 << 3)
