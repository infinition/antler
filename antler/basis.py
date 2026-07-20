"""ANTLER 0.1, basis.py
Base de Hilbert hard-core, secteur a N particules sur M sites.
Convention: site (rung i, leg sigma), i dans 0..L-1, sigma dans {0,1}.
Index lineaire k = 2*i + sigma (ordonnancement rung-by-rung, 0-based).
Etats: masques binaires, tri lexicographique strict (entiers croissants).
Version: 0.1.0
"""

from itertools import combinations
from math import comb

import numpy as np


def site_index(i: int, sigma: int) -> int:
    """Indice lineaire k du site (rung i, leg sigma)."""
    return 2 * i + sigma


def build_basis(n_sites: int, n_particles: int):
    """Retourne (states, index).
    states: np.ndarray de masques binaires tries, dimension comb(M, N).
    index: dict masque -> position dans la base.
    """
    states = []
    for occ in combinations(range(n_sites), n_particles):
        mask = 0
        for s in occ:
            mask |= 1 << s
        states.append(mask)
    states.sort()
    assert len(states) == comb(n_sites, n_particles)
    index = {s: p for p, s in enumerate(states)}
    return np.array(states, dtype=np.int64), index


def between_mask(k: int, l: int) -> int:
    """Masque des sites strictement entre k et l (k < l)."""
    if l - k < 2:
        return 0
    return ((1 << l) - 1) ^ ((1 << (k + 1)) - 1)


def popcount(x: int) -> int:
    return bin(x).count("1")
