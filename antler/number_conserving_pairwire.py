"""Number-conserving two-wire pair-hopping reference model.

This module is deliberately *not* a modification of the frozen ANTLER
Hamiltonian.  It provides a small, explicit fermionic benchmark family for
the next research question: can a number-conserving extension support a
locally indistinguishable parity doublet before one attempts any braid?

The Fock ordering is the project's rung-major ordering ``(j, a), (j, b)``.
Unlike :mod:`antler.model`, the signs below are ordinary fermionic Fock signs;
this is a separate, explicitly labelled extension, not a claim about the
Jordan--Wigner correlated-hopping ladder.
"""
from __future__ import annotations

import numpy as np

from .basis import build_basis, popcount, site_index


def _below_parity(state: int, site: int) -> int:
    """Parity of occupied Fock modes strictly below ``site``."""
    return popcount(state & ((1 << site) - 1)) & 1


def _annihilate(state: int, site: int) -> tuple[int, complex] | None:
    if not ((state >> site) & 1):
        return None
    return state ^ (1 << site), -1.0 if _below_parity(state, site) else 1.0


def _create(state: int, site: int) -> tuple[int, complex] | None:
    if (state >> site) & 1:
        return None
    return state | (1 << site), -1.0 if _below_parity(state, site) else 1.0


def _apply(state: int, operations: tuple[tuple[str, int], ...]) -> tuple[int, complex] | None:
    """Apply operations in their physical time order (rightmost first)."""
    amplitude: complex = 1.0
    current = state
    for kind, site in operations:
        item = _annihilate(current, site) if kind == "ann" else _create(current, site)
        if item is None:
            return None
        current, sign = item
        amplitude *= sign
    return current, amplitude


def build_pairwire_hamiltonian(
    L: int,
    N: int,
    t: float = 1.0,
    w: float = 1.0,
    v: float = 0.0,
    basis=None,
):
    """Build the minimal number-conserving two-wire Hamiltonian.

    ``H = -t H_leg - w H_pair + v H_nn`` where ``H_pair`` transfers a
    nearest-neighbour pair from one wire to the other.  It conserves total
    number and the fermion parity of each individual wire.  ``H_nn`` is a
    conventional nearest-neighbour density interaction within each wire.

    This is Iemini-inspired but is *not* presented as their exactly solvable
    Hamiltonian; it is an additive preflight family only.
    """
    if L < 2:
        raise ValueError("L must be at least two")
    M = 2 * L
    states, index = build_basis(M, N) if basis is None else basis
    H = np.zeros((len(states), len(states)), dtype=complex)

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        # Leg hopping and density interaction.
        for j in range(L - 1):
            for leg in (0, 1):
                left, right = site_index(j, leg), site_index(j + 1, leg)
                forward = _apply(state, (("ann", right), ("create", left)))
                backward = _apply(state, (("ann", left), ("create", right)))
                if forward is not None:
                    new, amp = forward
                    H[index[new], col] += -t * amp
                if backward is not None:
                    new, amp = backward
                    H[index[new], col] += -t * amp
                H[col, col] += v * (((state >> left) & 1) * ((state >> right) & 1))

            # a_j^dag a_{j+1}^dag b_{j+1} b_j + h.c.
            a0, a1 = site_index(j, 0), site_index(j + 1, 0)
            b0, b1 = site_index(j, 1), site_index(j + 1, 1)
            b_to_a = _apply(
                state,
                (("ann", b0), ("ann", b1), ("create", a1), ("create", a0)),
            )
            a_to_b = _apply(
                state,
                (("ann", a0), ("ann", a1), ("create", b1), ("create", b0)),
            )
            if b_to_a is not None:
                new, amp = b_to_a
                H[index[new], col] += -w * amp
            if a_to_b is not None:
                new, amp = a_to_b
                H[index[new], col] += -w * amp

    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("pair-wire Hamiltonian is not Hermitian")
    return H, states, index


def build_iemini_hamiltonian(
    L: int, N: int, lam: float = 1.0, basis=None, sparse: bool = False,
):
    """Build the open-boundary two-wire Hamiltonian of Iemini *et al.*

    This is Eq. (3) of F. Iemini *et al.*, PRL 115, 156402 (2015), in the
    rung-major Fock ordering used by this repository.  At ``lam=1`` it is the
    positive-semidefinite parent Hamiltonian with the published exactly
    solvable line.  It is intentionally separate from the frozen ANTLER
    correlated-hopping ladder and from :func:`build_pairwire_hamiltonian`.
    """
    if L < 2:
        raise ValueError("L must be at least two")
    M = 2 * L
    states, index = build_basis(M, N) if basis is None else basis
    if sparse:
        from scipy.sparse import dok_matrix
        H = dok_matrix((len(states), len(states)), dtype=complex)
    else:
        H = np.zeros((len(states), len(states)), dtype=complex)

    def add(col: int, state: int, coefficient: float, operations: tuple[tuple[str, int], ...]) -> None:
        item = _apply(state, operations)
        if item is not None:
            new, amplitude = item
            H[index[new], col] += coefficient * amplitude

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        for j in range(L - 1):
            a0, a1 = site_index(j, 0), site_index(j + 1, 0)
            b0, b1 = site_index(j, 1), site_index(j + 1, 1)
            for left, right in ((a0, a1), (b0, b1)):
                # -4[(c^dag_j c_{j+1}+h.c.)-(n_j+n_{j+1})+lambda n_j n_{j+1}]
                add(col, state, -4.0, (("ann", right), ("create", left)))
                add(col, state, -4.0, (("ann", left), ("create", right)))
                n_left, n_right = (state >> left) & 1, (state >> right) & 1
                H[col, col] += 4.0 * (n_left + n_right) - 4.0 * lam * n_left * n_right

            # -2 lambda [(n_a0+n_a1)(n_b0+n_b1) - (T1+T2-2T3+h.c.)]
            n_a = ((state >> a0) & 1) + ((state >> a1) & 1)
            n_b = ((state >> b0) & 1) + ((state >> b1) & 1)
            H[col, col] += -2.0 * lam * n_a * n_b

            # T1 = a0^dag a1 b0^dag b1; T2 = a0^dag a1 b1^dag b0;
            # T3 = b0^dag b1^dag a1 a0.  The tuples are applied right-to-left.
            for coefficient, operations in (
                (2.0 * lam, (("ann", b1), ("create", b0), ("ann", a1), ("create", a0))),
                (2.0 * lam, (("ann", b0), ("create", b1), ("ann", a1), ("create", a0))),
                (-4.0 * lam, (("ann", a0), ("ann", a1), ("create", b1), ("create", b0))),
                # Hermitian conjugates of T1, T2 and T3.
                (2.0 * lam, (("ann", a0), ("create", a1), ("ann", b0), ("create", b1))),
                (2.0 * lam, (("ann", a0), ("create", a1), ("ann", b1), ("create", b0))),
                (-4.0 * lam, (("ann", b0), ("ann", b1), ("create", a1), ("create", a0))),
            ):
                add(col, state, coefficient, operations)

    if sparse:
        H = H.tocsr()
        hermiticity_error = float(np.max(np.abs((H - H.getH()).data))) if (H - H.getH()).nnz else 0.0
        if hermiticity_error > 1e-12:
            raise RuntimeError("Iemini Hamiltonian is not Hermitian")
    elif not np.allclose(H, H.conj().T, atol=1e-12):
        raise RuntimeError("Iemini Hamiltonian is not Hermitian")
    return H, states, index


def wire_a_parity(state: int, L: int) -> int:
    """Return the parity (0 even, 1 odd) on wire ``a``."""
    return sum((state >> site_index(j, 0)) & 1 for j in range(L)) & 1


def local_leg_raising_matrix(states: np.ndarray, index: dict[int, int], L: int, rung: int) -> np.ndarray:
    """Matrix of the local, number-conserving operator ``a_rung^dag b_rung``."""
    d = len(states)
    O = np.zeros((d, d), dtype=complex)
    a, b = site_index(rung, 0), site_index(rung, 1)
    for col, raw_state in enumerate(states):
        state = int(raw_state)
        item = _apply(state, (("ann", b), ("create", a)))
        if item is not None:
            new, amp = item
            O[index[new], col] += amp
    return O


def local_transfer_matrix(states: np.ndarray, index: dict[int, int], L: int, rung: int) -> np.ndarray:
    """Hermitian local transfer quadrature ``a^dag b + b^dag a``."""
    raising = local_leg_raising_matrix(states, index, L, rung)
    return raising + raising.conj().T


def build_iemini_braid_z(
    L: int,
    N: int,
    kind: str,
    support: int,
    basis=None,
):
    """Build a finite-support anti-Hermitian ``Z`` operator from Iemini et al.

    ``kind='aR_bR'`` is the right-edge interwire operator of the paper.  The
    ``kind='aR_aL'`` construction is its published edge transformation and is
    bilocal (the first and last ``support`` rungs), as number conservation
    requires for an exchange of two ends of the same wire.  Both are truncated
    finite-size versions: their braid action is only asymptotically unitary in
    the ground space.
    """
    if kind not in {"aR_bR", "aR_aL"}:
        raise ValueError("kind must be 'aR_bR' or 'aR_aL'")
    if support < 1 or 2 * support >= L:
        raise ValueError("support must satisfy 1 <= support < L/2")
    M = 2 * L
    states, index = build_basis(M, N) if basis is None else basis
    from scipy.sparse import dok_matrix

    nu = N / (2.0 * L)
    normalizer = float(np.sqrt(1.0 - (nu ** (2 * support) + (1.0 - nu) ** (2 * support))))
    if normalizer <= 0.0:
        raise ValueError("braid normalizer is zero at this filling/support")
    Z = dok_matrix((len(states), len(states)), dtype=complex)

    def product_y_is_one(state: int, p: int) -> bool:
        # q=1 is the outermost rung.  The q=0 factor in the paper is I.
        for q in range(1, p):
            if kind == "aR_bR":
                left, right = site_index(L - q, 0), site_index(L - q, 1)
            else:
                left, right = site_index(q - 1, 0), site_index(L - q, 0)
            if ((state >> left) & 1) != ((state >> right) & 1):
                return False
        return True

    for col, raw_state in enumerate(states):
        state = int(raw_state)
        for p in range(1, support + 1):
            if not product_y_is_one(state, p):
                continue
            if kind == "aR_bR":
                a, b = site_index(L - p, 0), site_index(L - p, 1)
                terms = (
                    (1.0, (("ann", b), ("create", a))),
                    (-1.0, (("ann", a), ("create", b))),
                )
            else:
                left, right = site_index(p - 1, 0), site_index(L - p, 0)
                terms = (
                    (1j, (("ann", right), ("create", left))),
                    (1j, (("ann", left), ("create", right))),
                )
            for coefficient, operations in terms:
                item = _apply(state, operations)
                if item is not None:
                    new, amplitude = item
                    Z[index[new], col] += coefficient * amplitude / normalizer

    Z = Z.tocsr()
    antihermiticity_error = float(np.max(np.abs((Z + Z.getH()).data))) if (Z + Z.getH()).nnz else 0.0
    if antihermiticity_error > 1e-12:
        raise RuntimeError("finite-support braid Z is not anti-Hermitian")
    return Z, normalizer
