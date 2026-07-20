"""Local microscopic building blocks for a symmetry-first ANTLER extension.

This module is deliberately separate from the frozen two-leg ANTLER
Hamiltonian.  It implements a *local* hard-core mediator interferometer used
to test whether a parity-preserving pair transfer can be obtained from
ordinary one-particle hoppings, instead of inserting an effective pair term by
hand.  A successful local test is not a protected code, a lattice Hamiltonian,
or a braid.
"""
from __future__ import annotations

import numpy as np

from .basis import build_basis


MODE_NAMES = ("a0", "a1", "b0", "b1", "p0", "p1", "m0", "m1")
"""Two low rails ``a,b`` and two detuned mediator rails ``p,m``."""


def pair_mask(first: int, second: int) -> int:
    """Return the hard-core Fock mask with the two supplied modes occupied."""
    return (1 << first) | (1 << second)


def build_flux_pair_mediator_block(
    N: int,
    t: float,
    Delta: float = 5.0,
    U_mediator: float = 4.0,
    E_bind: float = 2.0,
    phi: float = np.pi,
    basis=None,
):
    """Build an eight-mode two-rail flux-interferometer block.

    The low modes are ``(a0,a1,b0,b1)``.  For each position ``r`` a particle
    can travel from ``a_r`` to ``b_r`` through a ``p_r`` or an ``m_r``
    mediator.  The latter path carries a Peierls phase ``phi``.  Both mediator
    rails have detuning ``Delta``; double occupancy of either mediator rail
    costs ``U_mediator``.  Low-rail pairs have an explicit local binding
    energy ``E_bind``.

    At ``phi=pi``, the second-order one-particle cotunnelling amplitudes of
    the two paths cancel.  In the two-particle sector, a residual fourth-order
    pair matrix element can be measured without putting a pair-transfer term
    directly in the Hamiltonian.

    This uses ordinary hard-core local hopping, not the frozen ANTLER
    correlated-Jordan--Wigner link.  Any future embedding in a full ANTLER
    ladder must derive the corresponding ordered-string convention again.
    """
    if N < 0 or N > len(MODE_NAMES):
        raise ValueError("invalid particle number")
    if Delta <= 0:
        raise ValueError("Delta must be positive for the detuned-mediator expansion")
    states, index = build_basis(len(MODE_NAMES), N) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)

    def occupied(state: int, site: int) -> int:
        return (state >> site) & 1

    def add_oriented_hop(col: int, state: int, source: int, destination: int, amplitude: complex) -> None:
        """Add ``amplitude c^dag_destination c_source + h.c.`` for hard cores."""
        if occupied(state, source) and not occupied(state, destination):
            new_state = state ^ (1 << source) ^ (1 << destination)
            H[index[new_state], col] += amplitude
        if occupied(state, destination) and not occupied(state, source):
            new_state = state ^ (1 << source) ^ (1 << destination)
            H[index[new_state], col] += np.conjugate(amplitude)

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        # Detuned mediators and their local two-particle interaction.
        H[col, col] += Delta * sum(occupied(state, site) for site in (4, 5, 6, 7))
        H[col, col] += U_mediator * (
            occupied(state, 4) * occupied(state, 5)
            + occupied(state, 6) * occupied(state, 7)
        )
        # A pair on either low rail is the retained low-energy manifold.
        H[col, col] -= E_bind * (
            occupied(state, 0) * occupied(state, 1)
            + occupied(state, 2) * occupied(state, 3)
        )
        for r in (0, 1):
            a, b, p, m = r, 2 + r, 4 + r, 6 + r
            add_oriented_hop(col, state, a, p, -t)
            add_oriented_hop(col, state, p, b, -t)
            add_oriented_hop(col, state, a, m, -t)
            add_oriented_hop(col, state, m, b, -t * np.exp(1j * phi))

    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("flux-pair mediator Hamiltonian is not Hermitian")
    return H, states, index


def low_pair_masks() -> tuple[int, int]:
    """Fock masks of the retained low pairs ``|aa>`` and ``|bb>``."""
    return pair_mask(0, 1), pair_mask(2, 3)


def one_particle_schur_cross_norm(H: np.ndarray, states: np.ndarray) -> float:
    """Second-order low-rail cross norm from the exact N=1 Schur complement.

    The retained subspace is all one-particle states on the two low rails.
    This diagnostic is intentionally evaluated at zero low-subspace energy;
    the low block has zero bare energy in the supplied construction.
    """
    if len(states) != 8 or any(int(state).bit_count() != 1 for state in states):
        raise ValueError("one_particle_schur_cross_norm requires the complete N=1 block")
    low = [i for i, state in enumerate(states) if int(state) & 0b1111]
    high = [i for i, state in enumerate(states) if not (int(state) & 0b1111)]
    Hpp = H[np.ix_(low, low)]
    Hpq = H[np.ix_(low, high)]
    Hqq = H[np.ix_(high, high)]
    effective = Hpp - Hpq @ np.linalg.solve(Hqq, Hpq.conj().T)
    # Ordering inside ``low`` is a0,a1,b0,b1 because the basis masks ascend.
    return float(np.linalg.norm(effective[:2, 2:], ord="fro"))
