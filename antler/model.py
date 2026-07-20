"""ANTLER 0.1, model.py
Hamiltonien SSH ladder anyonique en representation bosons hard-core.

Convention JW fractionnaire: a_k = exp(i theta sum_{m<k} n_m) b_k
sur l'ordre lineaire rung-by-rung k = 2*i + sigma.

Consequences exactes (k < l):
  a_k^dag a_l = b_k^dag exp(i theta sum_{m=k+1}^{l-1} n_m) b_l
  - saut de barreau (l = k+1): phase nulle
  - saut de jambe (l = k+2): phase exp(i theta n_{k+1})

ATTENTION (verrou physique): cet ordonnancement est une DEFINITION du
modele. Sur une echelle, la statistique fractionnaire via corde 1D est
dependante de la convention; H(theta) est un modele de saut correle,
pas un anyon universel. La symetrie de reflexion des jambes P est
brisee par la corde des que theta != 0 (voir README).

Dimerisation SSH (bond entre rungs i et i+1, i 0-based):
  i pair -> J1 (intra-cellule), i impair -> J2 (inter-cellule).
Phase topologique: |J1| < |J2|.

Version: 0.1.0
"""

import numpy as np

from .basis import between_mask, build_basis, popcount, site_index


def hop_list(L: int, J1: float, J2: float, Jperp: float):
    """Liste des sauts (k, l, J) avec k < l, sur la chaine lineaire."""
    hops = []
    for i in range(L):
        hops.append((site_index(i, 0), site_index(i, 1), Jperp))
    for i in range(L - 1):
        J = J1 if (i % 2 == 0) else J2
        for sigma in (0, 1):
            hops.append((site_index(i, sigma), site_index(i + 1, sigma), J))
    return hops


def single_particle_matrix(L: int, J1: float, J2: float, Jperp: float,
                           mu=None) -> np.ndarray:
    """Matrice tight-binding 2L x 2L (secteur N=1, phases JW inactives)."""
    M = 2 * L
    h = np.zeros((M, M), dtype=float)
    for k, l, J in hop_list(L, J1, J2, Jperp):
        h[k, l] += -J
        h[l, k] += -J
    if mu is not None:
        h[np.diag_indices(M)] += np.asarray(mu, dtype=float)
    return h


def build_hamiltonian(L: int, N: int, theta: float, J1: float, J2: float,
                      Jperp: float, mu=None, basis=None):
    """Hamiltonien exact dans le secteur a N particules.
    Retourne (H, states, index) avec H dense complexe (d petit en phase 0).
    """
    M = 2 * L
    if basis is None:
        states, index = build_basis(M, N)
    else:
        states, index = basis
    d = len(states)
    H = np.zeros((d, d), dtype=complex)
    hops = hop_list(L, J1, J2, Jperp)
    mu_arr = None if mu is None else np.asarray(mu, dtype=float)

    for col, s in enumerate(states):
        s = int(s)
        # potentiel chimique (diagonal)
        if mu_arr is not None:
            diag = 0.0
            for k in range(M):
                if (s >> k) & 1:
                    diag += mu_arr[k]
            H[col, col] += diag
        # termes de saut: a_k^dag a_l, k < l, particule en l vers k
        for k, l, J in hops:
            if ((s >> l) & 1) and not ((s >> k) & 1):
                n_mid = popcount(s & between_mask(k, l))
                amp = -J * np.exp(1j * theta * n_mid)
                new = s ^ (1 << l) ^ (1 << k)
                H[index[new], col] += amp
    # h.c. (les sauts sont purement hors diagonale)
    H_offdiag = H - np.diag(np.diag(H))
    H = H + H_offdiag.conj().T
    assert np.allclose(H, H.conj().T), "H non hermitien: convention JW violee"
    return H, states, index
