"""Phase 5B: projected N=3 exchange-holonomy algebra.

This is the first dynamic-algebra gatekeeper after the N=3 static scan.  It
tracks the same sparse-ED logical doublet through two local digital exchanges:
one at each code edge, with the third particle pinned as a mediator.  Parallel
transport removes dynamical phases and returns the geometric matrices U_A and
U_B, their commutator, and the minimum branch-isolation gaps.

The test is intentionally capable of returning a no-go: a pinned mediator may
leave both local exchanges diagonal and hence commuting.  That result directs
the project to a *mobile* mediator, synthetic dimension, or T-junction rather
than being mislabelled non-Abelian.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import polar
from scipy.optimize import linear_sum_assignment
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
for path in (ROOT, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from run_phase5_n3_mediator_preflight import N3Preflight, PreflightConfig


def sin2(s: float) -> float:
    return float(np.sin(0.5 * np.pi * np.clip(s, 0.0, 1.0)) ** 2)


def add_trap(mu: np.ndarray, leg: int, x: float, depth: float, L: int) -> None:
    x = float(np.clip(x, 0.0, L - 1.0))
    rung = int(np.floor(x))
    if rung >= L - 1:
        mu[2 * (L - 1) + leg] += depth
        return
    q = sin2(x - rung)
    mu[2 * rung + leg] += depth * (1.0 - q)
    mu[2 * (rung + 1) + leg] += depth * q


def left_exchange_mu(system: N3Preflight, u: float, loop_rung: int) -> np.ndarray:
    """One sequential digital exchange of the left code pair, mediator fixed."""

    cfg = system.cfg
    D = cfg.code_depth
    mu = system.mu().copy()
    # Replace only the left code wells by the compact Phase-4.3 sequence.
    mu[0] = 0.0
    mu[1] = 0.0
    if u < 0.35:
        add_trap(mu, 0, loop_rung * (u / 0.35), D, cfg.L)
        mu[1] += D
    elif u < 0.45:
        q = sin2((u - 0.35) / 0.10)
        mu[2 * loop_rung] += D
        mu[1] += D * (1.0 - q)
        mu[0] += D * q
    elif u < 0.55:
        q = sin2((u - 0.45) / 0.10)
        mu[0] += D
        mu[2 * loop_rung] += D * (1.0 - q)
        mu[2 * loop_rung + 1] += D * q
    else:
        add_trap(mu, 1, loop_rung * (1.0 - (u - 0.55) / 0.45), D, cfg.L)
        mu[0] += D
    return mu


def exchange_mu(system: N3Preflight, u: float, loop_rung: int, side: str) -> np.ndarray:
    left = left_exchange_mu(system, u, loop_rung)
    if side == "left":
        return left
    if side != "right":
        raise ValueError("side must be left or right")
    # Reflect only the change to the code wells.  The asymmetric mediator trap
    # remains fixed, so the two controls act on one common N=3 system.
    base = system.mu()
    return base + (left - base)[::-1]


def remove_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def sparse_frame(system: N3Preflight, H, n_eigs: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    energies, vectors = eigsh(H, k=n_eigs, which="SA", tol=1e-10)
    order = np.argsort(energies)
    energies, vectors = energies[order], vectors[:, order]
    score = abs(vectors[system.i_left, :]) ** 2 + abs(vectors[system.i_right, :]) ** 2
    pair = np.sort(np.argsort(-score)[:2])
    return energies, vectors, pair


def holonomy(system: N3Preflight, theta: float, loop_rung: int, side: str,
             steps: int, n_eigs: int) -> tuple[np.ndarray, dict]:
    """Discrete Kato parallel transport of the tracked two-dimensional branch."""

    if steps < 3:
        raise ValueError("at least three path samples are required")
    H0 = system.hamiltonian(theta)
    e0, v0, pair0 = sparse_frame(system, H0, n_eigs)
    frame0 = v0[:, pair0]
    previous = frame0
    transported = frame0.copy()
    min_gap = float("inf")
    min_link_sv = 1.0
    trace = []
    for u in np.linspace(0.0, 1.0, steps):
        mu = exchange_mu(system, float(u), loop_rung, side)
        H = system.hamiltonian(theta) + diags(system.occ @ (mu - system.mu()))
        energies, vectors = eigsh(H, k=n_eigs, which="SA", tol=1e-10)
        order = np.argsort(energies)
        energies, vectors = energies[order], vectors[:, order]
        score = abs(previous.conj().T @ vectors) ** 2
        source_rows, selected = linear_sum_assignment(-score)
        selected = selected[np.argsort(source_rows)]
        frame = vectors[:, selected]
        outside = np.setdiff1d(np.arange(n_eigs), selected)
        gap = float(np.min(abs(energies[selected, None] - energies[outside])))
        overlap = frame.conj().T @ transported
        singular = np.linalg.svd(overlap, compute_uv=False)
        transported = frame @ polar(overlap)[0]
        previous = frame
        min_gap = min(min_gap, gap)
        min_link_sv = min(min_link_sv, float(min(singular)))
        trace.append({
            "u": float(u), "selected_indices": selected.tolist(),
            "selected_energies": energies[selected].tolist(), "isolation_gap": gap,
            "link_sigma_min": float(min(singular)),
        })
    U = remove_global(frame0.conj().T @ transported)
    return U, {
        "side": side,
        "loop_rung": loop_rung,
        "theta": theta,
        "steps": steps,
        "n_eigs": n_eigs,
        "minimum_isolation_gap": min_gap,
        "minimum_parallel_transport_link_singular_value": min_link_sv,
        "trace": trace,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--code-depth", type=float, default=-2.5)
    parser.add_argument("--mediator-rung", type=int, default=7)
    parser.add_argument("--mediator-leg", type=int, default=0)
    parser.add_argument("--mediator-depth-scale", type=float, default=1.5)
    parser.add_argument("--loop-rung", type=int, default=3)
    parser.add_argument("--steps", type=int, default=41)
    parser.add_argument("--n-eigs", type=int, default=16)
    parser.add_argument("--noncommuting-threshold", type=float, default=1e-5,
                        help="minimum ||[U_A,U_B]|| before a braid residual is meaningful")
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/n3_local_exchange_holonomy.json"))
    args = parser.parse_args()
    if args.code_depth >= 0 or args.mediator_depth_scale <= 0:
        raise ValueError("depths must define attractive (negative) code wells and positive scale")
    if not 1 <= args.loop_rung < args.mediator_rung:
        raise ValueError("loop-rung must be between the left edge and the mediator")

    cfg = PreflightConfig(code_depth=args.code_depth, mediator_rung=args.mediator_rung,
                          mediator_leg=args.mediator_leg,
                          mediator_depth_scale=args.mediator_depth_scale)
    system = N3Preflight(cfg)
    started = time.time()
    UA, audit_a = holonomy(system, args.theta, args.loop_rung, "left",
                            args.steps, args.n_eigs)
    UB, audit_b = holonomy(system, args.theta, args.loop_rung, "right",
                            args.steps, args.n_eigs)
    commutator = UA @ UB - UB @ UA
    braid_residual = UA @ UB @ UA - UB @ UA @ UB
    commutator_norm = float(np.linalg.norm(commutator))
    raw_braid_residual = float(np.linalg.norm(braid_residual))
    braid_is_interpretable = commutator_norm >= args.noncommuting_threshold
    out = {
        "schema": "antler.phase5.n3-local-exchange-holonomy.v2",
        "claim_boundary": (
            "This is a projected adiabatic-holonomy diagnostic, not a finite-time "
            "gate fidelity calculation.  It tests whether pinned-mediator local "
            "exchanges can even generate a matrix-valued noncommuting action."
        ),
        "config": asdict(cfg), "theta": args.theta,
        "UA": arr(UA), "UB": arr(UB),
        "commutator_norm": commutator_norm,
        "noncommuting_threshold": args.noncommuting_threshold,
        "raw_braid_relation_residual": raw_braid_residual,
        "braid_relation_interpretable": braid_is_interpretable,
        "braid_relation_residual": raw_braid_residual if braid_is_interpretable else None,
        "braid_relation_status": (
            "reported_with_nonzero_commutator" if braid_is_interpretable else
            "not_interpretable_commutator_below_threshold"
        ),
        "UA_offdiag_norm": float(np.linalg.norm(UA - np.diag(np.diag(UA)))),
        "UB_offdiag_norm": float(np.linalg.norm(UB - np.diag(np.diag(UB)))),
        "audit_A": audit_a, "audit_B": audit_b,
        "runtime_s": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "commutator_norm": out["commutator_norm"],
        "raw_braid_relation_residual": out["raw_braid_relation_residual"],
        "braid_relation_status": out["braid_relation_status"],
        "UA_offdiag_norm": out["UA_offdiag_norm"],
        "UB_offdiag_norm": out["UB_offdiag_norm"],
        "gap_A": audit_a["minimum_isolation_gap"],
        "gap_B": audit_b["minimum_isolation_gap"],
    }, indent=2), flush=True)
    print(f"saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
