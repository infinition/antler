"""ANTLER 0.1, run_phase4_shuttle.py
Phase 4.0: navette manuelle, signature de tressage differentielle.

Choreographie (echange des deux particules du bord gauche):
  particule 1 piegee dans le puits mobile (0, jambe 1)
  particule 2 piegee statiquement en (0, jambe 2)
  EXCHANGE : p1 sort le long de la jambe 1 jusqu'au rung R, pendant le
    croisement p2 est transferee (0,j2) -> (0,j1) par bascule de pieges
    (anticroisement de gap 2*Jperp), p1 traverse le barreau R puis
    revient par la jambe 2 jusqu'en (0, j2). Lignes d'univers croisees.
  ROUND-TRIP : p1 sort et revient par la meme jambe, p2 immobile.
    Topologie triviale.

Comptage de cordes (prediction exacte):
  exchange   : un seul passage au-dessus de p2 -> phase e^{i theta}
  round-trip : passage aller e^{i theta}, retour e^{-i theta} -> nul
  donc Delta_phi_topo = theta, independant de T, A, w.

Mesure differentielle a double soustraction:
  phi_P(theta) = arg <psi_T^P(0) | psi_T^P(theta)>   (P = ex ou rt)
  Delta_phi_topo(theta) = phi_ex(theta) - phi_rt(theta)
  Toutes les phases dynamiques et geometriques non topologiques
  s'annulent (memes plannings, meme etat initial par theta).
  Test nul integre: Delta_phi_topo(0) = 0 par construction.

Plannings strictement lineaires (aucune IA, Phase 4.0).

Usage: python run_phase4_shuttle.py
"""

import time

import numpy as np

from antler.basis import build_basis
from antler.phase1 import build_h, hop_table

L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA = -4.0        # profondeur des pieges statiques
A_WELL = 2.6        # profondeur du puits mobile
W_WELL = 1.0       # largeur du puits mobile (en rungs)
R_LOOP = 6          # rung de retournement de la navette
T_TOTAL = 20000.0
NSEG = 2500
THETAS = (0.0, 0.3, -0.3, 0.6, -0.6, 0.9, -0.9)


def schedules(u, exchange):
    """Retourne (c, lam, trap_leg1_r0, trap_leg2_r0) au temps reduit u.
    c: centre du puits mobile (rung), lam: poids de jambe (0 = j1, 1 = j2).
    Plannings lineaires par morceaux."""
    def lin(a, b, x0, x1):
        s = np.clip((u - x0) / (x1 - x0), 0.0, 1.0)
        return a + (b - a) * s
    c = lin(0.0, R_LOOP, 0.00, 0.35) if u < 0.55 else lin(R_LOOP, 0.0, 0.55, 0.90)
    if exchange:
        lam = lin(0.0, 1.0, 0.35, 0.55)
        t1 = lin(0.0, DELTA, 0.35, 0.55)     # piege (0, j1) s'allume
        t2 = lin(DELTA, 0.0, 0.35, 0.55)     # piege (0, j2) s'eteint
    else:
        lam = 0.0
        t1 = 0.0
        t2 = DELTA
    return c, lam, t1, t2


def mu_sites(u, exchange, M):
    c, lam, t1, t2 = schedules(u, exchange)
    mu = np.zeros(M)
    rungs = np.arange(L)
    g = -A_WELL * np.exp(-(rungs - c) ** 2 / (2 * W_WELL ** 2))
    mu[0::2] += (1 - lam) * g            # puits mobile, jambe 1
    mu[1::2] += lam * g                  # puits mobile, jambe 2
    mu[0] += t1                          # piege statique (0, j1)
    mu[1] += t2                          # piege statique (0, j2)
    # etage final commun: les deux protocoles convergent vers le MEME H_f
    # (sites 0 et 1 a DELTA, puits eteint) pour que les termes de bord
    # geometriques soient identiques et se soustraient exactement
    if u > 0.88:
        r = np.sin(0.5 * np.pi * min((u - 0.88) / 0.12, 1.0)) ** 2
        mu_common = np.zeros(M)
        mu_common[0] = DELTA
        mu_common[1] = DELTA
        mu = (1 - r) * mu + r * mu_common
    return mu


def main():
    t0 = time.time()
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    table = hop_table(L, J1, J2, JPERP, states, index)
    # matrice d'occupation pour vectoriser mu_diagonal
    OCC = np.zeros((d, M))
    for p, s in enumerate(states):
        s = int(s)
        for k in range(M):
            if (s >> k) & 1:
                OCC[p, k] = 1.0

    from scipy.sparse import csr_matrix, diags
    from scipy.sparse.linalg import expm_multiply, eigsh, splu

    def sparse_hop(theta):
        rows, cols, mJ, nmid = table
        amp = mJ * np.exp(1j * theta * nmid)
        Hs = csr_matrix((amp, (rows, cols)), shape=(d, d))
        return Hs + Hs.conj().T

    def run(theta, exchange):
        Hhop = sparse_hop(theta)
        dt = T_TOTAL / NSEG
        H0 = Hhop + diags(OCC @ mu_sites(0.0, exchange, M))
        E0, V0 = eigsh(H0, k=2, which="SA")
        psi = V0[:, np.argmin(E0)].astype(complex)
        j = int(np.argmax(np.abs(psi)))
        psi *= np.conj(psi[j]) / abs(psi[j])   # jauge deterministe
        fid_inst_min = 1.0
        gap_min = np.inf
        for a in range(NSEG):
            u = (a + 0.5) / NSEG
            H = Hhop + diags(OCC @ mu_sites(u, exchange, M))
            psi = expm_multiply(-1j * dt * H.tocsc(), psi)
            if a % 250 == 0:
                e_loc = float(np.real(np.vdot(psi, H @ psi)))
                Ek, Vk = eigsh(H, k=5, sigma=e_loc, which="LM")
                ov = np.abs(Vk.conj().T @ psi) ** 2
                j = int(np.argmax(ov))
                fid_inst_min = min(fid_inst_min, float(ov[j]))
                others = np.delete(Ek, j)
                if len(others):
                    gap_min = min(gap_min, float(np.min(np.abs(others - Ek[j]))))
        psi /= np.linalg.norm(psi)
        # Adiabaticity at the common final Hamiltonian.
        Hf = Hhop + diags(OCC @ mu_sites(1.0, exchange, M))
        Ef, Vf = eigsh(Hf, k=2, which="SA")
        gs_f = Vf[:, np.argmin(Ef)]
        fid_final = float(abs(np.vdot(gs_f, psi)))
        return psi, fid_inst_min, gap_min, fid_final

    print("=" * 68)
    print(f"Phase 4.0 navette | L={L}, R={R_LOOP}, A={A_WELL}, w={W_WELL}, "
          f"T={T_TOTAL}, {NSEG} segments")
    print("=" * 68)
    ref = {}
    rows = []
    for proc in ("rt", "ex"):
        exchange = proc == "ex"
        for th in THETAS:
            psi, fmin, gmin, ffinal = run(th, exchange)
            if th == 0.0:
                ref[proc] = psi
                phi = 0.0
                fid = 1.0
            else:
                z = np.vdot(ref[proc], psi)
                phi = float(np.angle(z))
                fid = float(abs(z))
            rows.append((proc, th, phi, fid, fmin, gmin, ffinal))
            print(f"  {proc}  theta={th:.1f}: phi={phi:+.4f}, "
                  f"|<T(0)|T(th)>|={fid:.4f}, fid_inst_min={fmin:.4f}, "
                  f"gap_min={gmin:.4f}, fid_final={ffinal:.4f}  ({time.time() - t0:.0f} s)")

    print("-" * 68)
    print("Extraction impaire (le fond dynamique est PAIR en theta):")
    print("theta   odd_ex     odd_rt (test nul)   Delta_phi_topo   pred")
    def phi_of(proc, th):
        return next(r[2] for r in rows if r[0] == proc and abs(r[1] - th) < 1e-9)
    slopes = []
    for th in (0.3, 0.6, 0.9):
        oe = 0.5 * float(np.angle(np.exp(1j * (phi_of("ex", th) - phi_of("ex", -th)))))
        orr = 0.5 * float(np.angle(np.exp(1j * (phi_of("rt", th) - phi_of("rt", -th)))))
        dphi = float(np.angle(np.exp(1j * (oe - orr))))
        print(f"{th:.1f}   {oe:+.4f}    {orr:+.4f}            {dphi:+.4f}        +/-{th:.1f}")
        slopes.append(dphi / th)
    print("-" * 68)
    print(f"pente Delta_phi_topo / theta = {np.mean(slopes):+.4f} "
          f"(prediction de tressage: +/-1)")
    print(f"Total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
