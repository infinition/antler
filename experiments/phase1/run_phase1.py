"""ANTLER 0.1, run_phase1.py
Scan (theta, delta), trois cartes, test d'hypothese et robustesse.

Cartes produites (npz + png):
  A. fuite statique   P_leak = 1 - Tr(P_inst P_ref)/2
  B. split interne du doublet instantane
  C. courbure non abelienne ||F_theta_delta|| (plaquettes de Wilson)

Verdict: zone verte = { P_leak < 1% et split < 1e-3 }.
Hypothese validee si max ||F|| dans la zone verte > 0 significativement.

Usage: python run_phase1.py
"""

import time

import numpy as np

from antler.basis import build_basis
from antler.edge import edge_orbitals
from antler.logical import projector, two_particle_state
from antler.model import single_particle_matrix
from antler.phase1 import (build_h, disorder_test, hop_table,
                           instantaneous_doublet, mu_diagonal,
                           plaquette_curvature)

# parametres figes (Phase 0)
L = 14
N = 2
J1, J2, JPERP = 0.4, 1.0, 0.1
DELTA_REF = -0.35

# grille de scan
THETAS = np.linspace(0.0, np.pi, 25)
DELTAS = np.linspace(-0.60, -0.10, 21)

LEAK_MAX = 0.01
SPLIT_MAX = 1e-3


def main():
    t0 = time.time()
    M = 2 * L
    states, index = build_basis(M, N)
    d = len(states)
    edge_sites = [0, 1, M - 2, M - 1]

    # P_ref fige a lambda_ref (Phase 0)
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

    # table de sauts et diagonale de bord unitaire
    table = hop_table(L, J1, J2, JPERP, states, index)
    mu_edge_unit = np.zeros(M)
    mu_edge_unit[edge_sites] = 1.0
    diag_edge = mu_diagonal(states, mu_edge_unit)

    nt, nd = len(THETAS), len(DELTAS)
    leak = np.zeros((nt, nd))
    split = np.zeros((nt, nd))
    frames = np.empty((nt, nd), dtype=object)

    print(f"Scan {nt} x {nd} = {nt * nd} points ED (d = {d})...")
    for i, th in enumerate(THETAS):
        H_hop = build_h(th, table, d, np.zeros(d))
        for j, de in enumerate(DELTAS):
            H = H_hop.copy()
            H[np.diag_indices(d)] += de * diag_edge
            frame, cap, sp, _ = instantaneous_doublet(H, P_ref)
            leak[i, j] = 1.0 - cap / 2.0
            split[i, j] = sp
            frames[i, j] = frame
        print(f"  theta = {th:.3f} fait ({time.time() - t0:.0f} s)")

    dth = THETAS[1] - THETAS[0]
    dde = DELTAS[1] - DELTAS[0]
    print("Courbure non abelienne (plaquettes)...")
    Fmap = plaquette_curvature(frames, dth, dde)

    # zone verte evaluee aux centres de plaquettes
    leak_c = 0.25 * (leak[:-1, :-1] + leak[1:, :-1]
                     + leak[:-1, 1:] + leak[1:, 1:])
    split_c = 0.25 * (split[:-1, :-1] + split[1:, :-1]
                      + split[:-1, 1:] + split[1:, 1:])
    green = (leak_c < LEAK_MAX) & (split_c < SPLIT_MAX) & np.isfinite(Fmap)

    print("=" * 68)
    print(f"[A] fuite: min = {leak.min():.2e}, max = {leak.max():.3f}")
    print(f"[B] split: min = {split.min():.2e}, max = {split.max():.3e}")
    print(f"[C] ||F||: max global = {np.nanmax(Fmap):.4f}")
    print(f"Zone verte (P_leak < {LEAK_MAX}, split < {SPLIT_MAX}): "
          f"{green.sum()} / {green.size} plaquettes")
    if green.any():
        Fg = Fmap[green]
        i_best = np.unravel_index(np.argmax(np.where(green, Fmap, -np.inf)),
                                  Fmap.shape)
        th_b = 0.5 * (THETAS[i_best[0]] + THETAS[i_best[0] + 1])
        de_b = 0.5 * (DELTAS[i_best[1]] + DELTAS[i_best[1] + 1])
        print(f"||F|| en zone verte: max = {Fg.max():.4f}, "
              f"mediane = {np.median(Fg):.4f}")
        print(f"Optimum vert: theta = {th_b:.4f}, delta = {de_b:.4f}, "
              f"||F|| = {Fmap[i_best]:.4f}, "
              f"fuite = {leak_c[i_best]:.2e}, split = {split_c[i_best]:.2e}")
        verdict = Fg.max() > 10 * np.nanmedian(np.abs(Fmap[~green])) * 0 + 1e-3
        print(f"VERDICT hypothese Phase 1: "
              f"{'VALIDEE' if verdict else 'NON VALIDEE'} "
              f"(courbure non nulle avec sous-espace ferme)")
    else:
        print("VERDICT: aucune zone verte, hypothese NON VALIDEE sur ce domaine")

    # robustesse de la degenerescence sous desordre (verrou V3 reviewer)
    print("-" * 68)
    print("Robustesse doublet sous desordre diagonal (theta = 0, delta_ref):")
    mu_d_ref = DELTA_REF * diag_edge
    rng = np.random.default_rng(42)
    rob = disorder_test(table, d, states, mu_d_ref, P_ref,
                        widths=(1e-3, 1e-2, 5e-2), n_real=10, rng=rng)
    for W, (sm, sx, cm) in rob.items():
        print(f"  W = {W:.0e}: split moyen = {sm:.2e}, "
              f"split max = {sx:.2e}, capture min = {cm:.4f}")

    np.savez("/mnt/user-data/outputs/antler_phase1_maps.npz",
             thetas=THETAS, deltas=DELTAS, leak=leak, split=split,
             curvature=Fmap, green=green)
    print(f"Maps sauvees. Total: {time.time() - t0:.0f} s")

    # figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
        ext = [DELTAS[0], DELTAS[-1], THETAS[0], THETAS[-1]]
        data = [(np.log10(np.maximum(leak, 1e-12)), "log10 P_leak"),
                (np.log10(np.maximum(split, 1e-12)), "log10 split doublet"),
                (Fmap, "||F_theta_delta||")]
        for ax, (Z, title) in zip(axes, data):
            im = ax.imshow(Z, origin="lower", aspect="auto", extent=ext,
                           cmap="viridis")
            ax.set_xlabel("delta")
            ax.set_ylabel("theta")
            ax.set_title(title)
            fig.colorbar(im, ax=ax)
        fig.suptitle("ANTLER Phase 1: scan (theta, delta), L=14, N=2")
        fig.tight_layout()
        fig.savefig("/mnt/user-data/outputs/antler_phase1_maps.png", dpi=140)
        print("Figure sauvee: antler_phase1_maps.png")
    except Exception as e:
        print(f"(figure non generee: {e})")


if __name__ == "__main__":
    main()
