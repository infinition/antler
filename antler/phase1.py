"""ANTLER 0.1, phase1.py
Phase 1 Hypothesis Test.

Question unique: existe-t-il une region (theta, delta) ou la variete de
controle possede une courbure non abelienne non nulle dans le sous-espace
logique tout en gardant P_leak ~ 0 ?

Contenu:
  - table de sauts precalculee (le scan ne varie que theta et le diagonal)
  - identification du doublet instantane par recouvrement maximal avec
    P_L fige a lambda_ref
  - courbure non abelienne Wilczek-Zee par boucles de Wilson de plaquette
    (produit de liens unitarises, invariant de jauge)
  - test de robustesse de la degenerescence sous desordre diagonal

Version: 0.1.0
"""

import numpy as np

from .basis import between_mask, popcount
from .model import hop_list


def hop_table(L, J1, J2, Jperp, states, index):
    """Precalcul des elements de saut (direction k < l uniquement).
    Retourne (rows, cols, mJ, nmid): H_hop(theta)[rows, cols] += -J e^{i theta nmid}.
    """
    hops = hop_list(L, J1, J2, Jperp)
    rows, cols, mJ, nmid = [], [], [], []
    for col, s in enumerate(states):
        s = int(s)
        for k, l, J in hops:
            if ((s >> l) & 1) and not ((s >> k) & 1):
                new = s ^ (1 << l) ^ (1 << k)
                rows.append(index[new])
                cols.append(col)
                mJ.append(-J)
                nmid.append(popcount(s & between_mask(k, l)))
    return (np.array(rows), np.array(cols),
            np.array(mJ, dtype=float), np.array(nmid, dtype=float))


def build_h(theta, table, d, mu_diag):
    """H(theta) + diag(mu_diag) dans le secteur, dense complexe."""
    rows, cols, mJ, nmid = table
    H = np.zeros((d, d), dtype=complex)
    np.add.at(H, (rows, cols), mJ * np.exp(1j * theta * nmid))
    H = H + H.conj().T
    H[np.diag_indices(d)] += mu_diag
    return H


def mu_diagonal(states, mu_sites):
    """Diagonale du potentiel chimique dans la base du secteur."""
    d = len(states)
    out = np.zeros(d)
    for p, s in enumerate(states):
        s = int(s)
        acc = 0.0
        for k in range(len(mu_sites)):
            if (s >> k) & 1:
                acc += mu_sites[k]
        out[p] = acc
    return out


def instantaneous_doublet(H, P_ref):
    """Diagonalise H, retourne (frame d x 2, capture, split, E_pair).
    Le doublet instantane = les 2 etats propres de recouvrement maximal
    avec le code fige P_ref.
    """
    E, V = np.linalg.eigh(H)
    ov = np.real(np.einsum("in,ij,jn->n", V.conj(), P_ref, V))
    b = np.argsort(-ov)[:2]
    b = b[np.argsort(E[b])]
    frame = V[:, b]
    capture = float(ov[b].sum())
    split = float(abs(E[b[1]] - E[b[0]]))
    return frame, capture, split, E[b]


def _polar_unitary(M):
    """Partie unitaire de M (2x2) par SVD."""
    u, _, vh = np.linalg.svd(M)
    return u @ vh


def plaquette_curvature(frames, dtheta, ddelta):
    """Courbure non abelienne par plaquette de Wilson.
    frames: tableau objet [n_theta, n_delta] de matrices d x 2 (ou None).
    Retourne ||F|| (Frobenius) sur la grille des plaquettes, nan si frame
    manquante. W = U01 U12 U23 U30 (liens unitarises), F = -i log(W)/aire.
    """
    from scipy.linalg import logm
    nt, nd = frames.shape
    Fmap = np.full((nt - 1, nd - 1), np.nan)
    for i in range(nt - 1):
        for j in range(nd - 1):
            fs = (frames[i, j], frames[i + 1, j],
                  frames[i + 1, j + 1], frames[i, j + 1])
            if any(f is None for f in fs):
                continue
            W = np.eye(2, dtype=complex)
            for a in range(4):
                A, B = fs[a], fs[(a + 1) % 4]
                W = W @ _polar_unitary(A.conj().T @ B)
            F = -1j * logm(W) / (dtheta * ddelta)
            Fmap[i, j] = float(np.linalg.norm(F))
    return Fmap


def disorder_test(table, d, states, mu_base, P_ref, widths, n_real, rng):
    """Split et capture du doublet sous desordre diagonal uniforme [-W, W].
    Retourne dict W -> (split_mean, split_max, capture_min).
    """
    M = int(np.log2(int(states[-1])) + 1)
    out = {}
    for W in widths:
        splits, caps = [], []
        for _ in range(n_real):
            w = rng.uniform(-W, W, size=M)
            mu_d = mu_base + mu_diagonal(states, w)
            H = build_h(0.0, table, d, mu_d)
            _, cap, split, _ = instantaneous_doublet(H, P_ref)
            splits.append(split)
            caps.append(cap)
        out[W] = (float(np.mean(splits)), float(np.max(splits)),
                  float(np.min(caps)))
    return out
