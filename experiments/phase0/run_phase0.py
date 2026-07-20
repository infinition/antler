"""ANTLER 0.1, run_phase0.py
Phase 0: validation ED du ladder SSH anyonique et definition figee du qubit.

Definition du qubit (correctif issu des runs de validation):
  - orbitales de bord L1,L2,R1,R2 par rotation projecteur (edge.py)
  - modes rung bonding/antibonding: B = (X1+X2)/sqrt2, A = (X1-X2)/sqrt2
  - doublet chiral: |0_L> = |B_L A_R>, |1_L> = |A_L B_R>
    (energies e_B + e_A identiques par construction: degenerescence exacte)
  - biais de bord mu_edge = delta sur les 4 sites extremes: extrait le
    doublet (E ~ 2*delta) du shell chiral E=0 des paires de volume
  - P_L fige a lambda_ref (theta = 0), jamais rediagonalise ensuite

Criteres:
  C1. dim Hilbert conforme
  C2. quartet 1p converge vers delta +/- Jperp (hybridation L-R < seuil)
  C3. orbitales de bord localisees par cote
  C4. P_L rang 2, hermitien, idempotent
  C5. doublet 2p exactement degenere, capture logique > 0.999
  C6. robustesse statique en theta: split induit et fuite du code fige

Usage: python run_phase0.py
"""

from math import comb

import numpy as np

from antler.basis import build_basis
from antler.edge import edge_orbitals, side_weight
from antler.logical import projector, two_particle_state
from antler.model import build_hamiltonian, single_particle_matrix

# parametres MVP (lambda_ref)
L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA = -0.35          # biais de bord mu_edge
THETA_REF = 0.0


def main():
    M = 2 * L
    mu = np.zeros(M)
    mu[[0, 1, M - 2, M - 1]] = DELTA
    print("=" * 68)
    print(f"ANTLER Phase 0  |  L={L} rungs, M={M} sites, N={N}")
    print(f"lambda_ref: J1={J1}, J2={J2}, Jperp={JPERP}, "
          f"mu_edge={DELTA}, theta={THETA_REF}")
    print("=" * 68)

    # C1
    states, index = build_basis(M, N)
    d = len(states)
    print(f"[C1] dim H(N={N}) = {d} (attendu {comb(M, N)})")
    assert d == comb(M, N)

    # C2 + C3
    h1 = single_particle_matrix(L, J1, J2, JPERP, mu=mu)
    E1 = np.linalg.eigvalsh(h1)
    quartet = np.sort(E1[np.argsort(np.abs(E1 - DELTA))[:4]])
    t_lr = (quartet[1] - quartet[0]) / 2
    print(f"[C2] quartet 1p: {quartet.round(5)}")
    print(f"[C2] hybridation L-R t_LR = {t_lr:.2e} "
          f"(cible << Jperp = {JPERP})")
    assert t_lr < 0.05 * JPERP
    orbs, _, _ = edge_orbitals(h1, L)
    for name, v in orbs.items():
        ws = side_weight(v, L, name[0])
        print(f"[C3] {name}: poids demi-ladder {name[0]} = {ws:.5f}")
        assert ws > 0.99

    # C4: doublet chiral et projecteur
    s2 = np.sqrt(2)
    BL = (orbs["L1"] + orbs["L2"]) / s2
    AL = (orbs["L1"] - orbs["L2"]) / s2
    BR = (orbs["R1"] + orbs["R2"]) / s2
    AR = (orbs["R1"] - orbs["R2"]) / s2
    for nm, v in (("B_L", BL), ("A_L", AL), ("B_R", BR), ("A_R", AR)):
        print(f"[C4] {nm}: <h1> = {float(v @ h1 @ v):+.5f}")
    v0 = two_particle_state(BL, AR, states, index)
    v1 = two_particle_state(AL, BR, states, index)
    U, _ = np.linalg.qr(np.column_stack([v0, v1]))
    P = projector(U)
    print(f"[C4] P_L: rang = {np.trace(P).real:.6f}, fige a lambda_ref")

    # C5: spectre 2p a lambda_ref
    H0, _, _ = build_hamiltonian(L, N, THETA_REF, J1, J2, JPERP,
                                 mu=mu, basis=(states, index))
    E, V = np.linalg.eigh(H0)
    ov = np.array([np.real(np.vdot(V[:, n], P @ V[:, n])) for n in range(d)])
    band = np.sort(np.argsort(-ov)[:2])
    others = np.delete(np.arange(d), band)
    gap_local = min(np.min(np.abs(E[others] - E[band[0]])),
                    np.min(np.abs(E[others] - E[band[1]])))
    capture = ov[band].sum()
    print(f"[C5] doublet: E = {E[band].round(6)}, "
          f"split interne = {abs(E[band[0]] - E[band[1]]):.2e}")
    print(f"[C5] capture logique = {capture:.5f}/2, "
          f"distance spectrale locale = {gap_local:.4f}")
    print(f"[C5] gap interne bande de bord (vers |B_L B_R>, |A_L A_R>) "
          f"~ 2*Jperp = {2 * JPERP}")
    assert capture > 1.998
    assert abs(E[band[0]] - E[band[1]]) < 1e-10

    # C6: robustesse statique en theta (P_L fige)
    print("[C6] scan theta (P_L fige a theta=0):")
    for th in (0.25, 0.5, 1.0, np.pi / 2):
        Hth, _, _ = build_hamiltonian(L, N, th, J1, J2, JPERP,
                                      mu=mu, basis=(states, index))
        Eth, Vth = np.linalg.eigh(Hth)
        ovth = np.array([np.real(np.vdot(Vth[:, n], P @ Vth[:, n]))
                         for n in range(d)])
        b = np.argsort(-ovth)[:2]
        split = abs(Eth[b[0]] - Eth[b[1]])
        cap = ovth[b].sum()
        print(f"     theta={th:.4f}: split doublet = {split:.3e}, "
              f"capture code fige = {cap / 2:.4f}, "
              f"fuite statique = {1 - cap / 2:.4f}")

    print("=" * 68)
    print("Phase 0 PASSED. Verrous et findings: voir README.md")


if __name__ == "__main__":
    main()
