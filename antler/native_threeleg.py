"""Minimal native three-leg extension of the frozen ANTLER ladder.

The extension retains the project's scalar, rung-major correlated-hopping
convention and U(1) number conservation.  It deliberately contains neither
matrix-valued links nor the Iemini parent-Hamiltonian terms.  The module is a
Phase 6 candidate generator, separate from the frozen two-leg model.
"""
from __future__ import annotations

import numpy as np

from .basis import between_mask, build_basis, popcount


def mode(rung: int, leg: int) -> int:
    return 3 * rung + leg


def build_native_threeleg_hamiltonian(
    L: int,
    N: int,
    theta: float,
    J1: float,
    J2: float,
    Jperp: float,
    U_rung: float = 0.0,
    V_leg: float = 0.0,
    K_ring: float = 0.0,
    phi_ring: float = 0.0,
    basis=None,
):
    """Exact fixed-N Hamiltonian for the initial native candidate family.

    The hopping phase is the same scalar correlated Jordan--Wigner string as
    frozen ANTLER, now in the ordering ``(rung, leg=0,1,2)``.  Added terms are
    only local scalar density interactions:

    ``U_rung sum_{i,a<b} n_ia n_ib + V_leg sum_{i,a} n_ia n_{i+1,a}``.

    ``K_ring`` adds an oriented scalar plaquette exchange on the adjacent
    leg pairs ``(0,1)`` and ``(1,2)``.  It is a native three-leg ansatz, not a
    term copied from the two-wire Iemini parent Hamiltonian.

    This is intentionally a broad *preflight* family rather than a claim that
    these interactions have already been microscopically derived.
    """
    if L < 2 or N < 0 or N > 3 * L:
        raise ValueError("invalid L or N")
    states, index = build_basis(3 * L, N) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)
    hops: list[tuple[int, int, float]] = []
    for rung in range(L):
        hops.extend([
            (mode(rung, 0), mode(rung, 1), Jperp),
            (mode(rung, 1), mode(rung, 2), Jperp),
        ])
    for rung in range(L - 1):
        J = J1 if rung % 2 == 0 else J2
        for leg in range(3):
            hops.append((mode(rung, leg), mode(rung + 1, leg), J))

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        for rung in range(L):
            occupations = [((state >> mode(rung, leg)) & 1) for leg in range(3)]
            H[col, col] += U_rung * sum(
                occupations[left] * occupations[right]
                for left in range(3) for right in range(left + 1, 3)
            )
        for rung in range(L - 1):
            for leg in range(3):
                left, right = mode(rung, leg), mode(rung + 1, leg)
                H[col, col] += V_leg * (((state >> left) & 1) * ((state >> right) & 1))
        for left, right, coupling in hops:
            if ((state >> right) & 1) and not ((state >> left) & 1):
                n_mid = popcount(state & between_mask(left, right))
                amplitude = -coupling * np.exp(1j * theta * n_mid)
                new_state = state ^ (1 << right) ^ (1 << left)
                H[index[new_state], col] += amplitude
        if K_ring:
            for rung in range(L - 1):
                for lower_leg, upper_leg in ((0, 1), (1, 2)):
                    # c†_{i,a} c†_{i+1,b} c_{i,b} c_{i+1,a} in the
                    # hard-core bit convention, with an explicitly oriented
                    # scalar plaquette phase.  The h.c. is added globally.
                    source_a, source_b = mode(rung + 1, lower_leg), mode(rung, upper_leg)
                    dest_a, dest_b = mode(rung, lower_leg), mode(rung + 1, upper_leg)
                    sources_full = ((state >> source_a) & 1) and ((state >> source_b) & 1)
                    destinations_empty = not ((state >> dest_a) & 1) and not ((state >> dest_b) & 1)
                    if sources_full and destinations_empty:
                        new_state = state ^ (1 << source_a) ^ (1 << source_b) ^ (1 << dest_a) ^ (1 << dest_b)
                        H[index[new_state], col] += -K_ring * np.exp(1j * phi_ring)
    offdiag = H - np.diag(np.diag(H))
    H = H + offdiag.conj().T
    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("native three-leg Hamiltonian is not Hermitian")
    return H, states, index


def local_density_operator(states: np.ndarray, site: int) -> np.ndarray:
    """Diagonal local density in a fixed-N basis."""
    return np.asarray([((int(state) >> site) & 1) for state in states], dtype=float)
