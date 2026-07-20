"""ANTLER 0.1, edge.py
Extraction des 4 orbitales de bord monoparticulaires au point de
reference lambda_ref (theta = 0) et rotation vers la base localisee
{L1, L2, R1, R2}.

Procedure (deterministe, sans arbitraire de phase):
  1. Diagonaliser h (2L x 2L), prendre les 4 modes de plus petit |E|.
  2. Dans ce sous-espace, diagonaliser le projecteur moitie gauche:
     valeurs propres proches de 1 -> orbitales gauches, proches de 0 -> droites.
  3. Dans chaque doublet, diagonaliser le projecteur leg 1 (sigma = 0)
     pour separer les jambes.
  4. Fixer la jauge: premiere composante significative reelle positive.

Version: 0.1.0
"""

import numpy as np


def _fix_gauge(v: np.ndarray) -> np.ndarray:
    j = int(np.argmax(np.abs(v)))
    ph = v[j] / abs(v[j])
    return v / ph


def edge_orbitals(h: np.ndarray, L: int, e_ref: float = 0.0):
    """Retourne (orbs, E_edge, gap_bulk).
    orbs: dict {"L1","L2","R1","R2"} -> vecteur 2L (reel apres jauge).
    E_edge: energies des 4 modes de bord (les plus proches de e_ref).
    gap_bulk: distance a e_ref du premier mode hors quartet.
    """
    M = 2 * L
    E, V = np.linalg.eigh(h)
    order = np.argsort(np.abs(E - e_ref))
    edge_idx = order[:4]
    bulk_idx = order[4]
    W = V[:, edge_idx]                      # 2L x 4
    E_edge = E[edge_idx]
    gap_bulk = abs(E[bulk_idx] - e_ref)

    # projecteur moitie gauche: sites k < L (rungs 0..L//2-1)
    left = np.zeros(M)
    left[: 2 * (L // 2)] = 1.0
    A = W.conj().T @ (left[:, None] * W)    # 4 x 4
    a, U = np.linalg.eigh(A)
    Wrot = W @ U                            # colonnes triees par poids gauche croissant
    right_pair = Wrot[:, :2]
    left_pair = Wrot[:, 2:]

    # separation des jambes: projecteur leg 1 (sigma = 0, sites pairs)
    leg1 = np.zeros(M)
    leg1[0::2] = 1.0

    def split_legs(pair):
        B = pair.conj().T @ (leg1[:, None] * pair)
        b, Q = np.linalg.eigh(B)
        p = pair @ Q                        # tri: poids leg1 croissant
        return _fix_gauge(p[:, 1]), _fix_gauge(p[:, 0])   # (leg1, leg2)

    L1, L2 = split_legs(left_pair)
    R1, R2 = split_legs(right_pair)
    orbs = {"L1": L1, "L2": L2, "R1": R1, "R2": R2}
    return orbs, E_edge, gap_bulk


def edge_weight(v: np.ndarray, L: int, n_rungs: int = 1) -> float:
    """Poids de l'orbitale sur les n_rungs extremes (deux bouts)."""
    M = 2 * L
    w = np.abs(v) ** 2
    ends = list(range(0, 2 * n_rungs)) + list(range(M - 2 * n_rungs, M))
    return float(np.sum(w[ends]))


def side_weight(v: np.ndarray, L: int, side: str) -> float:
    """Poids sur la moitie gauche ('L') ou droite ('R') du ladder."""
    prof = (np.abs(v) ** 2).reshape(L, 2).sum(axis=1)
    half = L // 2
    return float(prof[:half].sum() if side == "L" else prof[L - half:].sum())
