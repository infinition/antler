"""ANTLER Phase 4.1a: isolated recapture harness.

Purpose
-------
Determine whether the large final leakage of the active |LL> branch is created
before recapture or by the final mobile-well -> static-edge-trap handoff.

The script:
  1. builds an exactly orthonormal dressed cat frame with dense eigh;
  2. propagates only the active left branch through the original protocol up to
     u=0.88 (start of recapture);
  3. measures its instantaneous eigenstate purity at that point;
  4. reuses the same incoming state for a scan of recapture ramps;
  5. reports final logical leakage and left-state capture.

Profiles
--------
linear      : mobile 1-s, static s
sin2        : mobile cos^2(pi s/2), static sin^2(pi s/2)
overlap(eta): static starts immediately and finishes at 1-eta; mobile remains
              full until eta then turns off. eta>0 creates a deep overlap.

The right edge remains statically trapped throughout.  No ghost shuttle is
introduced here: this experiment isolates the active-arm recapture only.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags

from antler.basis import build_basis
from antler.phase1 import hop_table


@dataclass
class Config:
    L: int = 14
    N: int = 2
    J1: float = 0.4
    J2: float = 1.0
    JPERP: float = 0.1
    DELTA: float = -4.0
    A_WELL: float = 2.6
    W_WELL: float = 1.0
    R_LOOP: float = 4.0
    THETA: float = 0.3
    T_TOTAL: float = 20000.0
    NSEG_TOTAL: int = 2500
    U_RECAP: float = 0.88


def clip01(x: float) -> float:
    return float(np.clip(x, 0.0, 1.0))


def smooth_sin2(x: float) -> float:
    x = clip01(x)
    return float(np.sin(0.5 * np.pi * x) ** 2)


def lin(u: float, a: float, b: float, x0: float, x1: float) -> float:
    s = clip01((u - x0) / (x1 - x0))
    return float(a + (b - a) * s)


def original_mu(u: float, exchange: bool, cfg: Config) -> np.ndarray:
    """Original common-H protocol used by the logical-gate probe."""
    L, M = cfg.L, 2 * cfg.L
    mu = np.zeros(M)
    mu[M - 2] = cfg.DELTA
    mu[M - 1] = cfg.DELTA
    if u < 0.10:
        c = 0.0
        amp = lin(u, 0.0, cfg.A_WELL, 0.00, 0.10)
        lam = 0.0
        t1 = lin(u, cfg.DELTA, 0.0, 0.00, 0.10)
        t2 = cfg.DELTA
    elif u < 0.38:
        c = lin(u, 0.0, cfg.R_LOOP, 0.10, 0.38)
        amp = cfg.A_WELL
        lam = 0.0
        t1, t2 = 0.0, cfg.DELTA
    elif u < 0.58:
        c = cfg.R_LOOP
        amp = cfg.A_WELL
        if exchange:
            lam = lin(u, 0.0, 1.0, 0.38, 0.58)
            t1 = lin(u, 0.0, cfg.DELTA, 0.38, 0.58)
            t2 = lin(u, cfg.DELTA, 0.0, 0.38, 0.58)
        else:
            lam = 0.0
            t1, t2 = 0.0, cfg.DELTA
    elif u < cfg.U_RECAP:
        c = lin(u, cfg.R_LOOP, 0.0, 0.58, cfg.U_RECAP)
        amp = cfg.A_WELL
        if exchange:
            lam = 1.0
            t1, t2 = cfg.DELTA, 0.0
        else:
            lam = 0.0
            t1, t2 = 0.0, cfg.DELTA
    else:
        # Baseline final ramp, included for closure/reference only.
        s = (u - cfg.U_RECAP) / (1.0 - cfg.U_RECAP)
        c = 0.0
        amp = cfg.A_WELL * (1.0 - clip01(s))
        if exchange:
            lam = 1.0
            t1, t2 = cfg.DELTA, cfg.DELTA * clip01(s)
        else:
            lam = 0.0
            t1, t2 = cfg.DELTA * clip01(s), cfg.DELTA
    mu[0] += t1
    mu[1] += t2
    rungs = np.arange(L, dtype=float)
    g = -amp * np.exp(-((rungs - c) ** 2) / (2.0 * cfg.W_WELL**2))
    mu[0::2] += (1.0 - lam) * g
    mu[1::2] += lam * g
    return mu


def ramp_weights(s: float, profile: str, eta: float) -> tuple[float, float]:
    """Return (mobile_fraction, static_fraction)."""
    s = clip01(s)
    if profile == "linear":
        return 1.0 - s, s
    if profile == "sin2":
        q = smooth_sin2(s)
        return 1.0 - q, q
    if profile == "overlap":
        if not (0.0 <= eta < 0.5):
            raise ValueError("eta must satisfy 0 <= eta < 0.5")
        span = 1.0 - eta
        static = smooth_sin2(s / span)
        mobile = 1.0 - smooth_sin2((s - eta) / span)
        return mobile, static
    raise ValueError(f"unknown profile: {profile}")


def recapture_mu(s: float, exchange: bool, cfg: Config,
                 profile: str, eta: float) -> np.ndarray:
    """Recapture-only potential, start H equals original H(u=0.88)."""
    L, M = cfg.L, 2 * cfg.L
    mobile, static = ramp_weights(s, profile, eta)
    mu = np.zeros(M)
    # Right cat arm is held fixed.
    mu[M - 2] = cfg.DELTA
    mu[M - 1] = cfg.DELTA
    if exchange:
        lam = 1.0
        mu[0] = cfg.DELTA
        mu[1] = cfg.DELTA * static
    else:
        lam = 0.0
        mu[0] = cfg.DELTA * static
        mu[1] = cfg.DELTA
    rungs = np.arange(L, dtype=float)
    g = -cfg.A_WELL * mobile * np.exp(-(rungs**2) / (2.0 * cfg.W_WELL**2))
    mu[0::2] += (1.0 - lam) * g
    mu[1::2] += lam * g
    return mu


def build_occ(states, M: int) -> np.ndarray:
    occ = np.zeros((len(states), M), dtype=float)
    for p, s0 in enumerate(states):
        s = int(s0)
        for k in range(M):
            occ[p, k] = (s >> k) & 1
    return occ


def bare_index(index, *sites: int) -> int:
    mask = 0
    for k in sites:
        mask |= 1 << k
    return index[mask]


def exact_logical_frame(H0, index, M: int):
    """Dense, exactly orthonormal localized dressed {|LL>,|RR>} frame."""
    H = H0.toarray() if hasattr(H0, "toarray") else np.asarray(H0)
    E, V = np.linalg.eigh(H)
    iLL = bare_index(index, 0, 1)
    iRR = bare_index(index, M - 2, M - 1)
    score = np.abs(V[iLL])**2 + np.abs(V[iRR])**2
    pair = np.argsort(-score)[:2]
    W, _ = np.linalg.qr(V[:, pair])
    zdiag = np.zeros(H.shape[0])
    zdiag[iLL], zdiag[iRR] = 1.0, -1.0
    Z = W.conj().T @ (zdiag[:, None] * W)
    _, Q = np.linalg.eigh(Z)
    U = W @ Q[:, ::-1]
    for j in range(2):
        k = int(np.argmax(np.abs(U[:, j])))
        U[:, j] *= np.conj(U[k, j]) / abs(U[k, j])
    return U, {
        "orth_error": float(np.linalg.norm(U.conj().T @ U - np.eye(2))),
        "energies": E[pair].tolist(),
        "bare_LL": float(abs(U[iLL, 0])**2),
        "bare_RR": float(abs(U[iRR, 1])**2),
    }


def sparse_hop(theta, table, d):
    rows, cols, mJ, nmid = table
    amp = mJ * np.exp(1j * theta * nmid)
    one = csr_matrix((amp, (rows, cols)), shape=(d, d))
    return one + one.conj().T


def hop_unitary(Hhop, dt: float):
    """Exact exp(-i H_hop dt), built once for Strang splitting."""
    Hd = Hhop.toarray() if hasattr(Hhop, "toarray") else np.asarray(Hhop)
    E, V = np.linalg.eigh(Hd)
    return (V * np.exp(-1j * E * dt)) @ V.conj().T


def strang_step(psi, Uhop, vdiag, dt: float):
    phase = np.exp(-0.5j * dt * vdiag)
    return phase * (Uhop @ (phase * psi))


def propagate(psi, Uhop, OCC, mu_fn, T: float, nseg: int):
    dt = T / nseg
    for a in range(nseg):
        s = (a + 0.5) / nseg
        psi = strang_step(psi, Uhop, OCC @ mu_fn(s), dt)
    return psi / np.linalg.norm(psi)


def eigen_purity(H, psi):
    """Dense diagnostic: largest instantaneous eigenstate population."""
    Hd = H.toarray() if hasattr(H, "toarray") else np.asarray(H)
    E, V = np.linalg.eigh(Hd)
    ov = np.abs(V.conj().T @ psi)**2
    j = int(np.argmax(ov))
    order = np.argsort(-ov)
    return {
        "max_population": float(ov[j]),
        "energy": float(E[j]),
        "top5_populations": [float(ov[k]) for k in order[:5]],
        "nearest_gap": float(np.min(np.abs(np.delete(E, j) - E[j]))),
    }


def metrics(psi, Ufinal):
    amps = Ufinal.conj().T @ psi
    pops = np.abs(amps)**2
    return {
        "left_capture": float(pops[0]),
        "right_admixture": float(pops[1]),
        "logical_capture": float(pops.sum()),
        "leakage": float(1.0 - pops.sum()),
        "left_phase": float(np.angle(amps[0])),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", type=float, default=0.3)
    ap.add_argument("--exchange", action="store_true")
    ap.add_argument("--dt", type=float, default=8.0)
    ap.add_argument("--durations", default="1200,2400,4800")
    ap.add_argument("--etas", default="0.0,0.2,0.35,0.45")
    ap.add_argument("--outdir", type=Path, default=Path("results/recapture"))
    args = ap.parse_args()
    cfg = Config(THETA=args.theta)
    cfg.NSEG_TOTAL = int(round(cfg.T_TOTAL / args.dt))
    M = 2 * cfg.L
    states, index = build_basis(M, cfg.N)
    d = len(states)
    OCC = build_occ(states, M)
    table = hop_table(cfg.L, cfg.J1, cfg.J2, cfg.JPERP, states, index)
    Hhop = sparse_hop(cfg.THETA, table, d)

    # Common final logical frame.
    mu_final = original_mu(1.0, args.exchange, cfg)
    Hfinal = Hhop + diags(OCC @ mu_final)
    Ufinal, frame_info = exact_logical_frame(Hfinal, index, M)

    # Active left dressed state at common initial H.
    Hinitial = Hhop + diags(OCC @ original_mu(0.0, args.exchange, cfg))
    Uinitial, _ = exact_logical_frame(Hinitial, index, M)
    psi = Uinitial[:, 0].copy()

    # Actual incoming state generated by the original protocol up to u=0.88.
    npre = int(round(cfg.U_RECAP * cfg.NSEG_TOTAL))
    dt = cfg.T_TOTAL / cfg.NSEG_TOTAL
    Uhop = hop_unitary(Hhop, dt)
    t0 = time.time()
    for a in range(npre):
        u = (a + 0.5) / cfg.NSEG_TOTAL
        psi = strang_step(psi, Uhop, OCC @ original_mu(u, args.exchange, cfg), dt)
    psi /= np.linalg.norm(psi)
    mu_start = original_mu(cfg.U_RECAP, args.exchange, cfg)
    Hstart = Hhop + diags(OCC @ mu_start)
    pre_diag = eigen_purity(Hstart, psi)

    args.outdir.mkdir(parents=True, exist_ok=True)
    durations = [float(x) for x in args.durations.split(",")]
    etas = [float(x) for x in args.etas.split(",")]
    rows = []
    tests = [("linear", 0.0), ("sin2", 0.0)] + [("overlap", e) for e in etas if e > 0]
    for Trec in durations:
        nseg = max(10, int(round(Trec / args.dt)))
        for profile, eta in tests:
            # Verify exact endpoint matching before propagation.
            err0 = np.linalg.norm(recapture_mu(0, args.exchange, cfg, profile, eta) - mu_start)
            err1 = np.linalg.norm(recapture_mu(1, args.exchange, cfg, profile, eta) - mu_final)
            if err0 > 1e-10 or err1 > 1e-10:
                raise RuntimeError(f"endpoint mismatch {profile}/{eta}: {err0}, {err1}")
            dt_rec = Trec / nseg
            Uhop_rec = Uhop if abs(dt_rec - dt) < 1e-12 else hop_unitary(Hhop, dt_rec)
            psif = propagate(
                psi.copy(), Uhop_rec, OCC,
                lambda s, p=profile, e=eta: recapture_mu(s, args.exchange, cfg, p, e),
                Trec, nseg,
            )
            met = metrics(psif, Ufinal)
            row = {
                "process": "exchange" if args.exchange else "roundtrip",
                "theta": cfg.THETA,
                "profile": profile,
                "eta": eta,
                "Trec": Trec,
                "nseg": nseg,
                **met,
            }
            rows.append(row)
            print(
                f"{row['process']:9s} {profile:7s} eta={eta:.2f} "
                f"T={Trec:6.0f}: leak={met['leakage']:.5f}, "
                f"left={met['left_capture']:.5f}, right={met['right_admixture']:.2e}",
                flush=True,
            )

    rows.sort(key=lambda r: r["leakage"])
    payload = {
        "config": asdict(cfg),
        "process": "exchange" if args.exchange else "roundtrip",
        "pre_recapture": pre_diag,
        "frame": frame_info,
        "best": rows[0],
        "rows": rows,
        "runtime_s": time.time() - t0,
    }
    stem = "ex" if args.exchange else "rt"
    (args.outdir / f"recapture_{stem}_{cfg.THETA:+.3f}.json").write_text(json.dumps(payload, indent=2))
    with (args.outdir / f"recapture_{stem}_{cfg.THETA:+.3f}.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("pre-recapture eigenstate diagnostic:", pre_diag)
    print("BEST:", rows[0])
    print(f"runtime={payload['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
