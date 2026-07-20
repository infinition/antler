"""ANTLER 0.1, run_lemma_antler1.py
Phase 2 kickoff.

Partie 1, Lemme ANTLER-1 (regle de selection de parite):
  Soit Pi l'operateur de parite de rung: Pi |config> = (-1)^(somme des
  indices de rung occupes) |config>.
  (i)   les sauts de rung commutent avec Pi (Pi-pairs)
  (ii)  les sauts de jambe anticommutent avec Pi (Pi-impairs)
  (iii) dH/dtheta ne contient QUE des sauts de jambe (les sauts de rung
        ont n_mid = 0), donc Pi (dH/dtheta) Pi = -(dH/dtheta)
  (iv)  si P_L est contenu dans UN secteur propre de Pi, alors pour tout
        operateur Pi-impair O: P_L O P_L = 0 EXACTEMENT.
        Preuve: X = P O P vit dans le secteur, donc Pi X Pi = X;
        mais Pi X Pi = P (Pi O Pi) P = -X. Donc X = 0.
  Corollaire: le zero de V10 est un zero de symetrie exact, pas un zero
  exponentiel. L'action de theta dans tout code Pi-pur est au mieux du
  SECOND ordre, d'amplitude controlee par la proximite des particules.

Partie 2, modele constructif (encodage meme bord, type chat):
  |0_L> = |B_L A_L>  (deux particules bord gauche, secteur Pi = +1)
  |1_L> = |B_R A_R>  (deux particules bord droit,  secteur Pi = +1)
  Premier ordre nul (lemme), mais action de theta au second ordre
  ~ J^2/Delta car les particules sont adjacentes. Sous delta_asym,
  l'action devient differentielle (porte Z geometrique candidate).
  Tests: degenerescence, action differentielle dE/dtheta, courbure
  F(theta, delta_asym) avec verification de scaling h^2.

Usage: python run_lemma_antler1.py
"""

import numpy as np
from scipy.linalg import logm

from antler.basis import build_basis
from antler.edge import edge_orbitals
from antler.logical import projector, two_particle_state
from antler.model import single_particle_matrix
from antler.phase1 import (build_h, hop_table, instantaneous_doublet,
                           mu_diagonal)

L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA_REF = -0.35


def main():
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    edge_sites = [0, 1, M - 2, M - 1]
    table = hop_table(L, J1, J2, JPERP, states, index)
    rows, cols, mJ, nmid = table

    # orbitales de bord au point de reference
    mu_ref = np.zeros(M)
    mu_ref[edge_sites] = DELTA_REF
    h1 = single_particle_matrix(L, J1, J2, JPERP, mu=mu_ref)
    orbs, _, _ = edge_orbitals(h1, L)
    s2 = np.sqrt(2)
    BL = (orbs["L1"] + orbs["L2"]) / s2
    AL = (orbs["L1"] - orbs["L2"]) / s2
    BR = (orbs["R1"] + orbs["R2"]) / s2
    AR = (orbs["R1"] - orbs["R2"]) / s2

    # operateur Pi (diagonal)
    pi_diag = np.empty(d)
    for p, s in enumerate(states):
        s = int(s)
        tot = sum((k // 2) for k in range(M) if (s >> k) & 1)
        pi_diag[p] = -1.0 if tot % 2 else 1.0

    def dH_dtheta(theta):
        Hd = np.zeros((d, d), dtype=complex)
        np.add.at(Hd, (rows, cols),
                  mJ * (1j * nmid) * np.exp(1j * theta * nmid))
        return Hd + Hd.conj().T

    print("=" * 68)
    print("PARTIE 1: Lemme ANTLER-1")
    # (ii)+(iii): Pi-imparite de dH/dtheta
    Hd = dH_dtheta(0.742)
    odd = np.linalg.norm(pi_diag[:, None] * Hd * pi_diag[None, :] + Hd)
    print(f"[L1.a] ||Pi dH/dth Pi + dH/dth|| = {odd:.2e} (Pi-impair exact)")
    # (i): H a sauts de rung seuls est Pi-pair
    table_rung = hop_table(L, 0.0, 0.0, JPERP, states, index)
    Hr = build_h(0.3, table_rung, d, np.zeros(d))
    even = np.linalg.norm(pi_diag[:, None] * Hr * pi_diag[None, :] - Hr)
    print(f"[L1.b] ||Pi H_rung Pi - H_rung|| = {even:.2e} (Pi-pair exact)")

    # (iv): purete de secteur des deux encodages
    v0 = two_particle_state(BL, AR, states, index)
    v1 = two_particle_state(AL, BR, states, index)
    U_sep, _ = np.linalg.qr(np.column_stack([v0, v1]))
    P_sep = projector(U_sep)
    vLL = two_particle_state(BL, AL, states, index)
    vRR = two_particle_state(BR, AR, states, index)
    U_cat, _ = np.linalg.qr(np.column_stack([vLL, vRR]))
    P_cat = projector(U_cat)
    for name, Ucode, expect in (("bords opposes", U_sep, -1.0),
                                ("meme bord (chat)", U_cat, +1.0)):
        for col in range(2):
            v = Ucode[:, col]
            mean = float(np.real(np.vdot(v, pi_diag * v)))
            var = float(np.real(np.vdot(v, (pi_diag ** 2) * v)) - mean ** 2)
            print(f"[L1.c] code {name}, etat {col}: <Pi> = {mean:+.6f}, "
                  f"var(Pi) = {var:.2e} (attendu {expect:+.0f}, 0)")
        G = Ucode.conj().T @ Hd @ Ucode
        print(f"[L1.d] code {name}: ||P dH/dth P|| = "
              f"{np.linalg.norm(G):.2e} (zero exact par le lemme)")

    # mecanisme: une composante Pi opposee d'amplitude eps retablit
    # un element intra-code lineaire en eps
    chi = (np.eye(d) - P_sep) @ (Hd @ U_sep[:, 0])
    chi /= np.linalg.norm(chi)
    print("[L1.e] restauration lineaire: eps -> |<chi + code| dH |code>|")
    for eps in (1e-3, 1e-2, 1e-1):
        v_pert = U_sep[:, 0] + eps * chi
        v_pert /= np.linalg.norm(v_pert)
        el = abs(np.vdot(v_pert, Hd @ U_sep[:, 1])) + \
            abs(np.vdot(v_pert, Hd @ v_pert))
        print(f"       eps = {eps:.0e}: element = {el:.3e}")

    print("=" * 68)
    print("PARTIE 2: encodage meme bord {|B_L A_L>, |B_R A_R>}")
    mu_edge = np.zeros(M)
    mu_edge[edge_sites] = 1.0
    diag_sym = mu_diagonal(states, mu_edge)
    muL = np.zeros(M); muL[[0, 1]] = 1.0
    muR = np.zeros(M); muR[[M - 2, M - 1]] = 1.0
    dL_unit = mu_diagonal(states, muL)
    dR_unit = mu_diagonal(states, muR)

    def H_at(th, dasym):
        dl = DELTA_REF + dasym / 2
        dr = DELTA_REF - dasym / 2
        return build_h(th, table, d, dl * dL_unit + dr * dR_unit)

    # degenerescence et isolation a lambda_ref
    H0 = H_at(0.0, 0.0)
    fr, cap, sp, Epair = instantaneous_doublet(H0, P_cat)
    print(f"[P2.a] doublet chat: E = {Epair.round(6)}, "
          f"split = {sp:.2e}, capture = {cap:.5f}/2")

    # action de theta au SECOND ordre: E(theta) du doublet suivi
    print("[P2.b] E_doublet(theta), delta_asym = 0 puis 0.10:")
    for da in (0.0, 0.10):
        Es = []
        ths = (0.0, 0.3, 0.6, 0.9)
        for th in ths:
            _, capt, spt, Ep = instantaneous_doublet(H_at(th, da), P_cat)
            Es.append(Ep)
        Es = np.array(Es)
        dE = Es[-1] - Es[0]
        print(f"   da = {da:.2f}: E(th=0) = {Es[0].round(5)}, "
              f"E(th=0.9) = {Es[-1].round(5)}")
        print(f"            shift des 2 branches = {dE.round(6)}, "
              f"action differentielle = {dE[1] - dE[0]:+.3e}")

    # courbure lisse F(theta, delta_asym) avec test de scaling
    def frame_at(th, da):
        return instantaneous_doublet(H_at(th, da), P_cat)[0]

    def logW(th, da, h):
        fs = [frame_at(th, da), frame_at(th + h, da),
              frame_at(th + h, da + h), frame_at(th, da + h)]
        W = np.eye(2, dtype=complex)
        for a in range(4):
            A, B = fs[a], fs[(a + 1) % 4]
            u, _, vh = np.linalg.svd(A.conj().T @ B)
            W = W @ (u @ vh)
        return float(np.linalg.norm(logm(W)))

    print("[P2.c] courbure F(theta, delta_asym) au point (0.6, 0.05):")
    prev = None
    for h in (0.08, 0.04, 0.02):
        lw = logW(0.6, 0.05, h)
        F = lw / h ** 2
        r = "" if prev is None else f", ratio |logW| = {prev / lw:.2f}"
        print(f"   h = {h:.2f}: |log W| = {lw:.3e}, ||F|| = {F:.4f}{r}")
        prev = lw
    print("   (scaling lisse attendu: ratio ~ 4, ||F|| stable)")

    print("=" * 68)
    print("Voir README, section Lemme ANTLER-1 et Phase 2.")


if __name__ == "__main__":
    main()
