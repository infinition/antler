"""ANTLER 0.1, logical.py
Construction du sous-espace logique a 2 particules et du projecteur P_L.

Definition figee (correctif du draft ANTLER-0.1):
Les DEUX etats logiques vivent dans le secteur "une particule par bord"
(bords opposes), pour rester quasi degeneres. Les etats de type
|L1 L2> ou |R1 R2> (deux particules meme bord) sont exclus: leur energie
est decalee par J_perp et par le blocage hard-core, ils cassent la
quasi-degenerescence.

Base logique (P-symetrique au point de reference theta = 0):
  |0_L> = (|L1 R1> + |L2 R2>) / sqrt(2)
  |1_L> = (|L1 R2> + |L2 R1>) / sqrt(2)
puis orthonormalisation QR (les recouvrements residuels de taille finie
sont elimines). P_L est FIGE au point lambda_ref: il ne doit jamais etre
rediagonalise pendant l'optimisation (sinon cible mouvante, triche IA).

Version: 0.1.0
"""

import numpy as np

from .basis import popcount


def two_particle_state(phiA: np.ndarray, phiB: np.ndarray,
                       states: np.ndarray, index: dict) -> np.ndarray:
    """Etat bosonique symetrise |A B> pour 2 particules hard-core.
    |AB> = Norm * sum_{k<l} (phiA[k] phiB[l] + phiA[l] phiB[k]) |k, l>
    Les composantes doublement occupees (k = l) sont projetees hors
    (contrainte hard-core).
    """
    d = len(states)
    v = np.zeros(d, dtype=complex)
    M = len(phiA)
    for p, s in enumerate(states):
        s = int(s)
        if popcount(s) != 2:
            continue
        occ = [k for k in range(M) if (s >> k) & 1]
        k, l = occ
        v[p] = phiA[k] * phiB[l] + phiA[l] * phiB[k]
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("Etat a deux particules nul: orbitales colineaires")
    return v / n


def logical_basis(orbs: dict, states: np.ndarray, index: dict):
    """Retourne (U, info): U de taille d x 2, colonnes orthonormees
    |0_L>, |1_L>; info contient les recouvrements bruts.
    """
    s_11 = two_particle_state(orbs["L1"], orbs["R1"], states, index)
    s_22 = two_particle_state(orbs["L2"], orbs["R2"], states, index)
    s_12 = two_particle_state(orbs["L1"], orbs["R2"], states, index)
    s_21 = two_particle_state(orbs["L2"], orbs["R1"], states, index)
    v0 = (s_11 + s_22) / np.sqrt(2)
    v1 = (s_12 + s_21) / np.sqrt(2)
    raw_overlap = complex(np.vdot(v0, v1))
    Q, _ = np.linalg.qr(np.column_stack([v0, v1]))
    info = {"raw_overlap_01": raw_overlap,
            "norm_v0": float(np.linalg.norm(v0)),
            "norm_v1": float(np.linalg.norm(v1))}
    return Q, info


def projector(U: np.ndarray) -> np.ndarray:
    """P_L = U U^dag, avec verifications."""
    P = U @ U.conj().T
    assert np.allclose(P, P.conj().T), "P_L non hermitien"
    assert np.allclose(P @ P, P), "P_L non idempotent"
    assert abs(np.trace(P).real - U.shape[1]) < 1e-10, "rang(P_L) incorrect"
    return P


def leakage(psi: np.ndarray, P: np.ndarray) -> float:
    """P_leak = 1 - <psi| P_L |psi>."""
    return float(1.0 - np.real(np.vdot(psi, P @ psi)))
