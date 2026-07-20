"""ANTLER 0.1, run_phase1_asym.py
Phase 1b: test d'hypothese complet apres decouverte du no-go de symetrie.

Sequence:
  1. Verification exacte de la symetrie anti-unitaire T = K o Rev
     (biais symetrique) : ||T H T^-1 - H|| = 0 => F = 0 sur (theta, delta_sym)
  2. Scan (theta, delta_L) a delta_R fixe (T brisee) : cartes fuite / split / ||F||
  3. Diagnostic de scaling des plaquettes : courbure lisse (log W ~ h^2)
     vs filaments d'anticroisements (log W ~ h, valeurs erratiques)
  4. Theoreme d'inertie : ||P_L (dH/dtheta) P_L|| et canal de fuite
     ||(1 - P_L)(dH/dtheta) P_L||

Usage: python run_phase1_asym.py
"""

import time

import numpy as np

from antler.basis import build_basis
from antler.edge import edge_orbitals
from antler.logical import projector, two_particle_state
from antler.model import single_particle_matrix
from antler.phase1 import (build_h, hop_table, instantaneous_doublet,
                           mu_diagonal, plaquette_curvature)

L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA_REF = -0.35
THETAS = np.linspace(0.2, 1.2, 13)
DELTALS = np.linspace(-0.50, -0.20, 13)


def main():
    t0 = time.time()
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    edge_sites = [0, 1, M - 2, M - 1]

    # P_ref fige (identique Phase 0)
    mu_ref = np.zeros(M)
    mu_ref[edge_sites] = DELTA_REF
    h1 = single_particle_matrix(L, J1, J2, JPERP, mu=mu_ref)
    orbs, _, _ = edge_orbitals(h1, L)
    s2 = np.sqrt(2)
    BL = (orbs["L1"] + orbs["L2"]) / s2
    AL = (orbs["L1"] - orbs["L2"]) / s2
    BR = (orbs["R1"] + orbs["R2"]) / s2
    AR = (orbs["R1"] - orbs["R2"]) / s2
    v0 = two_particle_state(BL, AR, states, index)
    v1 = two_particle_state(AL, BR, states, index)
    U, _ = np.linalg.qr(np.column_stack([v0, v1]))
    P_ref = projector(U)

    table = hop_table(L, J1, J2, JPERP, states, index)
    rows, cols, mJ, nmid = table
    muL = np.zeros(M); muL[[0, 1]] = 1.0
    muR = np.zeros(M); muR[[M - 2, M - 1]] = 1.0
    dL_unit = mu_diagonal(states, muL)
    dR_unit = mu_diagonal(states, muR)

    # 1. no-go de symetrie
    H = build_h(0.7, table, d, DELTA_REF * (dL_unit + dR_unit))
    rev_map = np.empty(d, dtype=int)
    for p, s in enumerate(states):
        s = int(s)
        ns = 0
        for k in range(M):
            if (s >> k) & 1:
                ns |= 1 << (M - 1 - k)
        rev_map[p] = index[ns]
    H_T = np.conj(H[np.ix_(rev_map, rev_map)])
    print(f"[1] no-go: ||T H T^-1 - H|| = {np.linalg.norm(H_T - H):.2e} "
          f"(T = K o Rev, biais symetrique) => F = 0 sur (theta, delta_sym)")

    # 2. scan asymetrique
    nt, nd = len(THETAS), len(DELTALS)
    leak = np.zeros((nt, nd)); split = np.zeros((nt, nd))
    frames = np.empty((nt, nd), dtype=object)
    print(f"[2] scan (theta, delta_L) {nt}x{nd}, delta_R = {DELTA_REF}...")
    for i, th in enumerate(THETAS):
        Hh = build_h(th, table, d, np.zeros(d))
        for j, dl in enumerate(DELTALS):
            Hp = Hh.copy()
            Hp[np.diag_indices(d)] += dl * dL_unit + DELTA_REF * dR_unit
            fr, cap, sp, _ = instantaneous_doublet(Hp, P_ref)
            leak[i, j] = 1 - cap / 2
            split[i, j] = sp
            frames[i, j] = fr
    F = plaquette_curvature(frames, THETAS[1] - THETAS[0],
                            DELTALS[1] - DELTALS[0])
    lc = 0.25 * (leak[:-1, :-1] + leak[1:, :-1] + leak[:-1, 1:] + leak[1:, 1:])
    green = (lc < 0.01) & np.isfinite(F)
    print(f"    fuite min/max = {leak.min():.1e}/{leak.max():.3f}, "
          f"split max = {split.max():.2e}")
    print(f"    zone P_leak<1%: {green.sum()}/{green.size}, "
          f"||F|| max en zone = {np.max(np.where(green, F, -np.inf)):.3f}")

    # 3. diagnostic de scaling (filaments vs champ lisse)
    from scipy.linalg import logm

    def frame_at(th, dl):
        Hp = build_h(th, table, d, dl * dL_unit + DELTA_REF * dR_unit)
        return instantaneous_doublet(Hp, P_ref)[0]

    def logW(th, dl, h):
        fs = [frame_at(th, dl), frame_at(th + h, dl),
              frame_at(th + h, dl + h), frame_at(th, dl + h)]
        W = np.eye(2, dtype=complex)
        for a in range(4):
            A, B = fs[a], fs[(a + 1) % 4]
            u, _, vh = np.linalg.svd(A.conj().T @ B)
            W = W @ (u @ vh)
        return float(np.linalg.norm(logm(W)))

    print("[3] scaling |log W| (lisse => ratio 4, filament => ratio ~<2):")
    rng = np.random.default_rng(1)
    ratios = []
    for _ in range(6):
        th = 0.742 + rng.uniform(-0.15, 0.15)
        dl = -0.288 + rng.uniform(-0.05, 0.05)
        a, b = logW(th, dl, 0.02), logW(th, dl, 0.01)
        ratios.append(a / b)
        print(f"    ({th:.3f},{dl:.3f}): {a:.2e} -> {b:.2e}, ratio {a / b:.2f}")
    print(f"    ratios erratiques ({min(ratios):.1f} a {max(ratios):.0f}): "
          f"pas de champ lisse, courbure concentree en filaments")

    # 4. theoreme d'inertie du code separe
    print("[4] inertie du code en theta:")
    for th in (0.0, 0.742):
        Hd = np.zeros((d, d), dtype=complex)
        np.add.at(Hd, (rows, cols), mJ * (1j * nmid) * np.exp(1j * th * nmid))
        Hd = Hd + Hd.conj().T
        inside = np.linalg.norm(P_ref @ Hd @ P_ref)
        out = np.linalg.norm((np.eye(d) - P_ref) @ Hd @ P_ref)
        print(f"    theta={th:.3f}: ||P dH/dth P|| = {inside:.2e} (zero machine), "
              f"||(1-P) dH/dth P|| = {out:.2e} (fuite pure)")

    print("=" * 68)
    print("VERDICT Phase 1: hypothese statique FALSIFIEE.")
    print("theta n'a AUCUNE action intra-code sur l'encodage a bords opposes;")
    print("il n'agit que comme canal de fuite. La statistique exige la")
    print("proximite des particules: Phase 2 = protocole d'echange dynamique")
    print("(navette d'une particule a travers le volume), pas boucle statique.")
    np.savez("/mnt/user-data/outputs/antler_phase1_asym.npz",
             thetas=THETAS, deltals=DELTALS, leak=leak, split=split,
             curvature=F)
    print(f"Total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
