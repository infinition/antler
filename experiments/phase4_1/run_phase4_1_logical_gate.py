"""ANTLER Phase 4.1: full logical-subspace gate test.

A single common Hamiltonian acts on the cat code {|LL>, |RR>}.
The left branch experiences either an exchange shuttle or a trivial round trip;
the right branch remains trapped and is the spectator arm of the interferometer.
The script propagates the two dressed logical basis vectors simultaneously and
reconstructs the projected 2x2 logical operation.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from scipy.linalg import polar
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh, expm_multiply

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
    T_TOTAL: float = 20000.0
    NSEG: int = 5000


def lin(u: float, a: float, b: float, x0: float, x1: float) -> float:
    s = np.clip((u - x0) / (x1 - x0), 0.0, 1.0)
    return float(a + (b - a) * s)


def protocol_mu(u: float, exchange: bool, cfg: Config) -> np.ndarray:
    """Common start/end Hamiltonian with four static edge traps.

    Stages:
      0.00-0.10 load left-leg-1 particle into mobile well
      0.10-0.38 shuttle outward on leg 1
      0.38-0.58 exchange leg at R and swap the stationary left trap
      0.58-0.88 return (leg 2 for exchange, leg 1 for round trip)
      0.88-1.00 unload into the common four-edge-trap Hamiltonian
    Right-edge traps stay on throughout, providing the spectator logical arm.
    """
    L, M = cfg.L, 2 * cfg.L
    mu = np.zeros(M)
    # Right logical arm is always trapped.
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
    elif u < 0.88:
        c = lin(u, cfg.R_LOOP, 0.0, 0.58, 0.88)
        amp = cfg.A_WELL
        if exchange:
            lam = 1.0
            t1, t2 = cfg.DELTA, 0.0
        else:
            lam = 0.0
            t1, t2 = 0.0, cfg.DELTA
    else:
        c = 0.0
        amp = lin(u, cfg.A_WELL, 0.0, 0.88, 1.00)
        if exchange:
            lam = 1.0
            t1 = cfg.DELTA
            t2 = lin(u, 0.0, cfg.DELTA, 0.88, 1.00)
        else:
            lam = 0.0
            t1 = lin(u, 0.0, cfg.DELTA, 0.88, 1.00)
            t2 = cfg.DELTA

    mu[0] += t1
    mu[1] += t2
    rungs = np.arange(L, dtype=float)
    g = -amp * np.exp(-((rungs - c) ** 2) / (2.0 * cfg.W_WELL**2))
    mu[0::2] += (1.0 - lam) * g
    mu[1::2] += lam * g
    return mu


def build_occ(states: np.ndarray, M: int) -> np.ndarray:
    occ = np.zeros((len(states), M), dtype=float)
    for p, s0 in enumerate(states):
        s = int(s0)
        for k in range(M):
            occ[p, k] = (s >> k) & 1
    return occ


def bare_index(index: dict[int, int], *sites: int) -> int:
    mask = 0
    for k in sites:
        mask |= 1 << k
    return index[mask]


def logical_frame(H0, index: dict[int, int], M: int) -> tuple[np.ndarray, dict]:
    """Extract and localize the two dressed cat states around |LL>, |RR>."""
    # A small low-energy window is enough because deep edge traps isolate the cats.
    nev = min(12, H0.shape[0] - 2)
    E, V = eigsh(H0, k=nev, which="SA")
    order = np.argsort(E)
    E, V = E[order], V[:, order]
    iLL = bare_index(index, 0, 1)
    iRR = bare_index(index, M - 2, M - 1)
    score = np.abs(V[iLL, :]) ** 2 + np.abs(V[iRR, :]) ** 2
    pair = np.argsort(-score)[:2]
    W = V[:, pair]
    # Rotate the exact eigen-pair into maximally left/right localized dressed states.
    zdiag = np.zeros(H0.shape[0])
    zdiag[iLL] = 1.0
    zdiag[iRR] = -1.0
    Z = W.conj().T @ (zdiag[:, None] * W)
    zval, Q = np.linalg.eigh(Z)
    U = W @ Q[:, ::-1]  # left first, right second
    for j in range(2):
        k = int(np.argmax(np.abs(U[:, j])))
        U[:, j] *= np.conj(U[k, j]) / abs(U[k, j])
    info = {
        "energies": E[pair].tolist(),
        "bare_LL_weights": [float(abs(U[iLL, j]) ** 2) for j in range(2)],
        "bare_RR_weights": [float(abs(U[iRR, j]) ** 2) for j in range(2)],
        "selected_score": score[pair].tolist(),
        "orth_error": float(np.linalg.norm(U.conj().T @ U - np.eye(2))),
    }
    return U, info


def unitary_part(S: np.ndarray) -> np.ndarray:
    U, _ = polar(S)
    return U


def remove_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def gate_fidelity(U: np.ndarray, target: np.ndarray) -> float:
    d = U.shape[0]
    return float((abs(np.trace(target.conj().T @ U)) ** 2 + d) / (d * (d + 1)))


def simulate(theta: float, exchange: bool, cfg: Config, states, index, table, OCC):
    d, M = len(states), 2 * cfg.L
    rows, cols, mJ, nmid = table
    amp = mJ * np.exp(1j * theta * nmid)
    one_way = csr_matrix((amp, (rows, cols)), shape=(d, d))
    Hhop = one_way + one_way.conj().T
    mu0 = protocol_mu(0.0, exchange, cfg)
    muf = protocol_mu(1.0, exchange, cfg)
    if not np.allclose(mu0, muf, atol=1e-12):
        raise RuntimeError("Protocol does not close to the common Hamiltonian")
    H0 = Hhop + diags(OCC @ mu0)
    U0, frame_info = logical_frame(H0, index, M)
    Psi = U0.copy()
    dt = cfg.T_TOTAL / cfg.NSEG
    max_outside = np.zeros(2)
    # Sparse midpoint propagation; both logical states are propagated together.
    for a in range(cfg.NSEG):
        u = (a + 0.5) / cfg.NSEG
        H = Hhop + diags(OCC @ protocol_mu(u, exchange, cfg))
        Psi = expm_multiply((-1j * dt) * H.tocsc(), Psi)
        if a % max(1, cfg.NSEG // 20) == 0:
            S_now = U0.conj().T @ Psi
            max_outside = np.maximum(max_outside, 1.0 - np.sum(np.abs(S_now) ** 2, axis=0))
    # Projected logical map and exact final leakage.
    S = U0.conj().T @ Psi
    final_leak = 1.0 - np.sum(np.abs(S) ** 2, axis=0)
    Up = unitary_part(S)
    return {
        "theta": theta,
        "exchange": exchange,
        "S": S,
        "U": Up,
        "final_leak": final_leak,
        "sampled_max_leak": max_outside,
        "frame": frame_info,
        "closure_error": float(np.linalg.norm(protocol_mu(0, exchange, cfg) - protocol_mu(1, exchange, cfg))),
    }


def analyze_pair(ex, rt, theta: float):
    # Round-trip is the matched dynamical/geometric reference.
    D = remove_global(rt["U"].conj().T @ ex["U"])
    rel_phase = float(np.angle(D[0, 0]) - np.angle(D[1, 1]))
    rel_phase = float(np.angle(np.exp(1j * rel_phase)))
    offdiag = float(np.linalg.norm(D - np.diag(np.diag(D))))
    targets = {
        "minus": remove_global(np.diag([np.exp(-1j * theta), 1.0])),
        "plus": remove_global(np.diag([np.exp(+1j * theta), 1.0])),
    }
    fids = {k: gate_fidelity(D, v) for k, v in targets.items()}
    return D, rel_phase, offdiag, fids


def serializable_run(r):
    return {
        "theta": r["theta"],
        "exchange": r["exchange"],
        "S_real": r["S"].real.tolist(), "S_imag": r["S"].imag.tolist(),
        "U_real": r["U"].real.tolist(), "U_imag": r["U"].imag.tolist(),
        "final_leak": r["final_leak"].tolist(),
        "sampled_max_leak": r["sampled_max_leak"].tolist(),
        "frame": r["frame"], "closure_error": r["closure_error"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", type=float, default=0.3)
    ap.add_argument("--T", type=float, default=20000.0)
    ap.add_argument("--nseg", type=int, default=5000)
    ap.add_argument("--R", type=float, default=4.0)
    ap.add_argument("--out", type=Path, default=Path("results/phase4_1_logical_gate.json"))
    args = ap.parse_args()
    cfg = Config(T_TOTAL=args.T, NSEG=args.nseg, R_LOOP=args.R)
    M = 2 * cfg.L
    states, index = build_basis(M, cfg.N)
    table = hop_table(cfg.L, cfg.J1, cfg.J2, cfg.JPERP, states, index)
    OCC = build_occ(states, M)
    t0 = time.time()
    results = {}
    for th in (args.theta, -args.theta):
        for name, exchange in (("rt", False), ("ex", True)):
            print(f"run {name}, theta={th:+.3f}", flush=True)
            results[(name, th)] = simulate(th, exchange, cfg, states, index, table, OCC)
            r = results[(name, th)]
            print(f"  final leakage={r['final_leak']}, sampled max={r['sampled_max_leak']}", flush=True)
    analyses = {}
    for th in (args.theta, -args.theta):
        D, phase, offdiag, fids = analyze_pair(results[("ex", th)], results[("rt", th)], th)
        analyses[str(th)] = {
            "D_real": D.real.tolist(), "D_imag": D.imag.tolist(),
            "relative_phase": phase, "offdiag_norm": offdiag, "fidelities": fids,
        }
        print(f"theta={th:+.3f}: rel phase={phase:+.6f}, offdiag={offdiag:.3e}, fids={fids}")
    odd_phase = 0.5 * np.angle(np.exp(1j * (analyses[str(args.theta)]["relative_phase"] - analyses[str(-args.theta)]["relative_phase"])))
    print(f"odd logical relative phase={odd_phase:+.6f}, slope={odd_phase/args.theta:+.6f}")
    payload = {
        "config": asdict(cfg),
        "runs": {f"{k[0]}_{k[1]:+.3f}": serializable_run(v) for k, v in results.items()},
        "analysis": analyses,
        "odd_relative_phase": float(odd_phase),
        "odd_slope": float(odd_phase / args.theta),
        "runtime_s": time.time() - t0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"saved {args.out} in {payload['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
