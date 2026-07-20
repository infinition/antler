"""ANTLER 0.1, run_phase3_loop.py
Phase 3: les trois tests critiques du reviewer.

  TEST 1 (habillage): P_L habille = span des 2 etats propres exacts de
    H(lambda_ref) de recouvrement maximal avec le chat nu {|L1L2>,|R1R2>}.
    Capture 100% par construction au point de reference; on verifie le
    contenu physique (poids chat, localisation, gap local) et la capture
    aux coins du domaine de controle.

  TEST 2 (courbure logique): F Wilczek-Zee du doublet suivi, plaquettes
    de Wilson, scaling h^2, decomposition Pauli dans la base habillee
    (f0 abelien, fz differentiel Z, |fxy| non diagonal).

  TEST 3 (boucle fermee): rectangle (theta, delta_asym):
    (0,0) -> (0.8,0) -> (0.8,0.08) -> (0,0.08) -> (0,0).
    a. transport de Wilson (holonomie geometrique ideale U_geo)
    b. evolution temporelle reelle aller et retour (schema segmentaire
       exact par diagonalisation), extraction:
         P_leak max, U_logical, phase relative geometrique
         phi_geo = (phi_aller - phi_retour)/2 (la partie dynamique
         s'annule dans la demi-difference).

Usage: python run_phase3_loop.py
"""

import time

import numpy as np
from scipy.linalg import logm

from antler.basis import build_basis
from antler.edge import edge_orbitals
from antler.logical import two_particle_state
from antler.model import single_particle_matrix
from antler.phase1 import build_h, hop_table, mu_diagonal

L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA_REF = -2.0
TH_MAX, DA_MAX = 0.8, 0.30
NSEG_SIDE = 60
T_TOTAL = 800.0


def main():
    t0 = time.time()
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    table = hop_table(L, J1, J2, JPERP, states, index)
    mu_ref = np.zeros(M)
    mu_ref[[0, 1, M - 2, M - 1]] = DELTA_REF
    h1 = single_particle_matrix(L, J1, J2, JPERP, mu=mu_ref)
    orbs, _, _ = edge_orbitals(h1, L, e_ref=DELTA_REF)
    vLL = two_particle_state(orbs["L1"], orbs["L2"], states, index)
    vRR = two_particle_state(orbs["R1"], orbs["R2"], states, index)
    muL = np.zeros(M); muL[[0, 1]] = 1.0
    muR = np.zeros(M); muR[[M - 2, M - 1]] = 1.0
    dL = mu_diagonal(states, muL)
    dR = mu_diagonal(states, muR)

    def H_at(th, da):
        return build_h(th, table, d,
                       (DELTA_REF + da / 2) * dL + (DELTA_REF - da / 2) * dR)

    # ---------- TEST 1: habillage
    print("=" * 68)
    print("TEST 1: habillage du code")
    E0, V0 = np.linalg.eigh(H_at(0.0, 0.0))
    # a da=0 les etats propres sont les chats de Rev (LL +/- RR)/sqrt2:
    # on prend les 2 meilleurs par score combine puis on tourne le cadre
    # vers la base localisee {LL, RR}
    wL = np.abs(V0.conj().T @ vLL) ** 2
    wR = np.abs(V0.conj().T @ vRR) ** 2
    score0 = wL + wR
    pair = np.argsort(-score0)[:2]
    pair = pair[np.argsort(E0[pair])]
    iL, iR = int(pair[0]), int(pair[1])
    Upair = V0[:, [iL, iR]]
    D = Upair.conj().T @ (np.outer(vLL, vLL.conj())
                          - np.outer(vRR, vRR.conj())) @ Upair
    _, Q = np.linalg.eigh(D)
    U_ref = Upair @ Q[:, ::-1]   # colonne 0 -> type LL, colonne 1 -> type RR
    # phase de jauge: composante dominante reelle positive
    for c in range(2):
        j = int(np.argmax(np.abs(U_ref[:, c])))
        U_ref[:, c] /= U_ref[j, c] / abs(U_ref[j, c])
    w0L = float(abs(np.vdot(vLL, U_ref[:, 0])) ** 2)
    w1R = float(abs(np.vdot(vRR, U_ref[:, 1])) ** 2)
    print(f"  paire habillee: E = {E0[iL]:+.6f}, {E0[iR]:+.6f}, "
          f"split = {abs(E0[iR] - E0[iL]):.3e}")
    print(f"  cadre localise: |<LL|0>|^2 = {w0L:.4f}, |<RR|1>|^2 = {w1R:.4f}")
    others = np.delete(np.arange(d), [iL, iR])
    gap_loc = min(np.min(np.abs(E0[others] - E0[iL])),
                  np.min(np.abs(E0[others] - E0[iR])))
    print(f"  gap local vers le reste = {gap_loc:.4f}")

    def tracked(H, ref_frame):
        E, V = np.linalg.eigh(H)
        ov = np.abs(V.conj().T @ ref_frame) ** 2
        score = ov.sum(axis=1)
        b = np.argsort(-score)[:2]
        b = b[np.argsort(E[b])]
        cap = float(score[b].sum())
        return V[:, b], cap, E[b], (E, V)

    for th, da in ((TH_MAX, 0.0), (TH_MAX, DA_MAX), (0.0, DA_MAX)):
        _, cap, _, _ = tracked(H_at(th, da), U_ref)
        print(f"  capture code habille en ({th:.2f},{da:.2f}) = {cap / 2:.4f}")

    # ---------- TEST 2: courbure WZ du doublet habille
    print("=" * 68)
    print("TEST 2: courbure Wilczek-Zee (base habillee)")

    def polar(Mx):
        u, _, vh = np.linalg.svd(Mx)
        return u @ vh

    def F_matrix(th, da, h):
        pts = [(th, da), (th + h, da), (th + h, da + h), (th, da + h)]
        frames = []
        prev = U_ref
        for p in pts:
            fr, _, _, _ = tracked(H_at(*p), prev)
            frames.append(fr)
            prev = fr
        W = np.eye(2, dtype=complex)
        for a in range(4):
            A, B = frames[a], frames[(a + 1) % 4]
            W = W @ polar(A.conj().T @ B)
        # exprime dans la base habillee du coin de base
        G0 = polar(U_ref.conj().T @ frames[0])
        W = G0 @ W @ G0.conj().T
        return -1j * logm(W) / (h * h)

    for (th, da) in ((0.4, 0.04), (0.7, 0.06)):
        prev_n = None
        for h in (0.04, 0.02):
            F = F_matrix(th, da, h)
            f0 = float(np.real(np.trace(F)) / 2)
            fz = float(np.real(F[0, 0] - F[1, 1]) / 2)
            fxy = float(abs(F[0, 1]))
            n = float(np.linalg.norm(F))
            r = "" if prev_n is None else f", ratio ||F|| = {prev_n / n:.2f}"
            print(f"  ({th:.2f},{da:.2f}) h={h:.2f}: ||F||={n:.2f}, "
                  f"f0={f0:+.2f}, fz={fz:+.2f}, |fxy|={fxy:.2e}{r}")
            prev_n = n

    # ---------- TEST 3: boucle fermee
    print("=" * 68)
    print(f"TEST 3: boucle rectangle theta<= {TH_MAX}, da <= {DA_MAX}, "
          f"T = {T_TOTAL}, {4 * NSEG_SIDE} segments")

    def perimeter(n_side):
        pts = []
        for s in np.linspace(0, 1, n_side, endpoint=False):
            pts.append((s * TH_MAX, 0.0))
        for s in np.linspace(0, 1, n_side, endpoint=False):
            pts.append((TH_MAX, s * DA_MAX))
        for s in np.linspace(0, 1, n_side, endpoint=False):
            pts.append(((1 - s) * TH_MAX, DA_MAX))
        for s in np.linspace(0, 1, n_side, endpoint=False):
            pts.append((0.0, (1 - s) * DA_MAX))
        pts.append((0.0, 0.0))
        return pts

    pts = perimeter(NSEG_SIDE)

    # a. transport de Wilson
    prev = U_ref
    W = np.eye(2, dtype=complex)
    for p in pts[1:]:
        fr, _, _, _ = tracked(H_at(*p), prev)
        W = W @ polar(prev.conj().T @ fr)
        prev = fr
    W = W @ polar(prev.conj().T @ U_ref)
    phi_W = float(np.angle(W[0, 0]) - np.angle(W[1, 1]))
    print(f"  [a] Wilson: |offdiag| = {abs(W[0, 1]):.2e}, "
          f"phase relative geometrique = {phi_W:+.4f} rad")

    # b. evolution temporelle reelle
    def evolve(path):
        nseg = len(path) - 1
        dt = T_TOTAL / nseg
        psi = U_ref.copy().astype(complex)
        leak_max = 0.0
        ref_frame = U_ref
        for a in range(nseg):
            th = 0.5 * (path[a][0] + path[a + 1][0])
            da = 0.5 * (path[a][1] + path[a + 1][1])
            E, V = np.linalg.eigh(H_at(th, da))
            psi = V @ (np.exp(-1j * E * dt)[:, None] * (V.conj().T @ psi))
            if a % 6 == 0:
                fr, _, _, _ = tracked_from(E, V, ref_frame)
                ref_frame = fr
                P2 = fr @ fr.conj().T
                for c in range(2):
                    lk = 1 - float(np.real(np.vdot(psi[:, c], P2 @ psi[:, c])))
                    leak_max = max(leak_max, lk)
        Mfin = U_ref.conj().T @ psi
        leak_end = [1 - float(np.linalg.norm(Mfin[:, c]) ** 2)
                    for c in range(2)]
        return Mfin, leak_max, leak_end

    def tracked_from(E, V, ref_frame):
        ov = np.abs(V.conj().T @ ref_frame) ** 2
        score = ov.sum(axis=1)
        b = np.argsort(-score)[:2]
        b = b[np.argsort(E[b])]
        return V[:, b], float(score[b].sum()), E[b], None

    Mf, leak_f, leak_end_f = evolve(pts)
    Mr, leak_r, leak_end_r = evolve(pts[::-1])
    phi_f = float(np.angle(Mf[0, 0] * np.conj(Mf[1, 1])))
    phi_r = float(np.angle(Mr[0, 0] * np.conj(Mr[1, 1])))
    phi_geo = 0.5 * (phi_f - phi_r)
    phi_dyn = 0.5 * (phi_f + phi_r)
    print(f"  [b] aller : P_leak max = {leak_f:.4f}, "
          f"fin = {max(leak_end_f):.4f}, phi_rel = {phi_f:+.4f}")
    print(f"      retour: P_leak max = {leak_r:.4f}, "
          f"fin = {max(leak_end_r):.4f}, phi_rel = {phi_r:+.4f}")
    print(f"      phi_geometrique = {phi_geo:+.4f} rad, "
          f"phi_dynamique = {phi_dyn:+.4f} rad")
    print(f"      reference Wilson = {phi_W:+.4f} rad, "
          f"ecart = {abs(phi_geo - phi_W):.4f}")
    print(f"      |offdiag U_logical| aller = "
          f"{abs(Mf[0, 1]) + abs(Mf[1, 0]):.2e}")
    print("=" * 68)
    print(f"Total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
