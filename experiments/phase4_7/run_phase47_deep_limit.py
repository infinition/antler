"""Phase 4.7A -- close the digital deep-trap quantisation test.

This campaign deliberately separates two limits:

* physical localisation: |D| -> infinity while T(D) = T_ref (|D|/D_ref)^2;
* numerical propagation: dt -> 0 at fixed D=6 and D=8.

The JSON summary records the full logical reconstruction, leakage, singular
values, off-diagonal mixing, and an independently tracked handoff gap for
every completed point.  It does not promote a result automatically.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from phase47_common import DigitalConfig, require_json, run_gate


def point_name(depth: float, dt: float) -> str:
    return f"D{depth:g}_dt{dt:g}".replace(".", "p")


def compatible(cached: dict, cfg: DigitalConfig, theta: float, dt: float) -> bool:
    return (
        cached.get("schema") == "antler.phase47.gate.v1"
        and cached.get("config") == cfg.__dict__
        and float(cached.get("theta", np.nan)) == theta
        and float(cached.get("dt_requested", np.nan)) == dt
        and int(cached.get("cycles", 0)) == 1
    )


def run_or_resume(path: Path, cfg: DigitalConfig, theta: float, dt: float,
                  gap_samples: int, resume: bool) -> dict:
    cached = require_json(path) if resume else None
    if cached is not None and compatible(cached, cfg, theta, dt):
        print(f"resume {path.name}", flush=True)
        return cached
    print(
        f"run depth={abs(cfg.DEPTH):g}, T={cfg.T_TOTAL:g}, dt={dt:g}, "
        f"theta={theta:g}",
        flush=True,
    )
    started = time.time()
    result = run_gate(cfg, theta=theta, dt=dt, include_gap=True,
                      gap_samples=gap_samples)
    result["runtime_s"] = time.time() - started
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["metrics"], indent=2), flush=True)
    return result


def fit_deep_limit(rows: list[dict]) -> dict:
    """Fit slope(D)=-1+a D^-p; report residual diagnostics, not a claim."""

    if len(rows) < 3:
        return {"status": "need_at_least_three_completed_depths"}
    x = np.asarray([row["depth"] for row in rows], float)
    y = np.asarray([row["odd_slope"] + 1.0 for row in rows], float)
    if np.any(y <= 0.0):
        return {
            "status": "fixed_limit_fit_not_applicable",
            "reason": "at least one slope is at or beyond -1; inspect data without enforcing sign",
        }

    def model(depth: np.ndarray, a: float, p: float) -> np.ndarray:
        return a * depth ** (-p)

    params, covariance = curve_fit(model, x, y, p0=(0.8, 2.0), maxfev=100_000)
    predicted = model(x, *params)
    residual = y - predicted
    ss_total = float(np.sum((y - y.mean()) ** 2))
    r2 = float(1.0 - np.sum(residual ** 2) / ss_total) if ss_total else 1.0
    return {
        "model": "odd_slope(D) = -1 + a * D**(-p)",
        "a": float(params[0]),
        "p": float(params[1]),
        "stderr": np.sqrt(np.diag(covariance)).tolist(),
        "r2": r2,
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "observed_minus_one": y.tolist(),
        "predicted_minus_one": predicted.tolist(),
        "residual": residual.tolist(),
    }


def convergence_table(rows: list[dict], finest_dt: float) -> list[dict]:
    """Compare dt values at fixed physical D, always against the finest one."""

    by_depth: dict[float, list[dict]] = {}
    for row in rows:
        by_depth.setdefault(row["depth"], []).append(row)
    output = []
    for depth, group in sorted(by_depth.items()):
        ref = next((r for r in group if abs(r["dt"] - finest_dt) < 1e-12), None)
        if ref is None:
            output.append({"depth": depth, "status": "missing_finest_dt"})
            continue
        for row in sorted(group, key=lambda r: r["dt"], reverse=True):
            output.append({
                "depth": depth,
                "dt": row["dt"],
                "reference_dt": finest_dt,
                "delta_slope_to_finest": row["odd_slope"] - ref["odd_slope"],
                "delta_phase_to_finest": row["odd_phase"] - ref["odd_phase"],
                "delta_leak_to_finest": row["leak_worst"] - ref["leak_worst"],
                "delta_sigma_min_to_finest": row["sigma_min"] - ref["sigma_min"],
                "delta_offdiag_to_finest": row["odd_offdiag_norm"] - ref["odd_offdiag_norm"],
            })
    return output


def compact_row(result: dict, depth: float, dt: float, role: str, path: Path) -> dict:
    metrics = result["metrics"]
    gap = result.get("gap") or {}
    return {
        "role": role,
        "depth": depth,
        "T": result["config"]["T_TOTAL"],
        "dt": dt,
        "odd_slope": metrics["odd_slope"],
        "odd_phase": metrics["odd_phase"],
        "leak_worst": metrics["leak_worst"],
        "sigma_min": metrics["sigma_min"],
        "sigma_max": metrics["sigma_max"],
        "odd_offdiag_norm": metrics["odd_offdiag_norm"],
        "unitarity_frob_max": metrics["unitarity_frob_max"],
        "favg_target": metrics["favg_target"],
        "minimum_handoff_isolation_gap": gap.get("minimum_handoff_isolation_gap"),
        "gap_at_u": gap.get("handoff_at_u"),
        "file": path.as_posix(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depths", type=float, nargs="+", default=[4, 6, 8, 10, 12])
    parser.add_argument("--convergence-depths", type=float, nargs="+", default=[6, 8])
    parser.add_argument("--timesteps", type=float, nargs="+", default=[0.5, 0.25, 0.125])
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--reference-depth", type=float, default=4.0)
    parser.add_argument("--reference-time", type=float, default=20_000.0)
    parser.add_argument("--time-power", type=float, default=2.0)
    parser.add_argument("--primary-dt", type=float, default=0.25)
    parser.add_argument("--gap-samples", type=int, default=61)
    parser.add_argument("--outdir", type=Path, default=Path("results/phase4_7/deep_limit"))
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.primary_dt not in args.timesteps:
        raise ValueError("primary-dt must be listed in --timesteps")
    if args.reference_depth <= 0 or args.reference_time <= 0:
        raise ValueError("reference depth and time must be positive")

    started = time.time()
    primary_rows: list[dict] = []
    primary_results: dict[float, tuple[dict, Path]] = {}
    convergence_rows: list[dict] = []
    for depth in args.depths:
        if depth <= 0:
            raise ValueError("depths are positive magnitudes; the trap sign is applied internally")
        T = args.reference_time * (depth / args.reference_depth) ** args.time_power
        cfg = DigitalConfig(DEPTH=-depth, T_TOTAL=T)
        path = args.outdir / f"{point_name(depth, args.primary_dt)}.json"
        result = run_or_resume(path, cfg, args.theta, args.primary_dt, args.gap_samples,
                               resume=not args.no_resume)
        primary_results[depth] = (result, path)
        primary_rows.append(compact_row(result, depth, args.primary_dt, "deep_limit", path))

    for depth in args.convergence_depths:
        if depth not in args.depths:
            raise ValueError("every convergence depth must also be included in --depths")
        T = args.reference_time * (depth / args.reference_depth) ** args.time_power
        cfg = DigitalConfig(DEPTH=-depth, T_TOTAL=T)
        for dt in args.timesteps:
            if abs(dt - args.primary_dt) < 1e-12:
                # This is exactly the primary physical point, not an
                # approximation.  Reusing it avoids an unnecessary long run.
                result, path = primary_results[depth]
            else:
                path = args.outdir / "dt_convergence" / f"{point_name(depth, dt)}.json"
                result = run_or_resume(path, cfg, args.theta, dt, args.gap_samples,
                                       resume=not args.no_resume)
            convergence_rows.append(compact_row(result, depth, dt, "dt_convergence", path))

    summary = {
        "schema": "antler.phase47.deep-limit.v1",
        "claim_boundary": (
            "A completed fit supports only the correlated-hopping ladder with the "
            "frozen rung-major Jordan-Wigner convention.  It is not a claim of "
            "universal anyonic or non-Abelian computation."
        ),
        "protocol": {
            "theta": args.theta,
            "T_of_D": "T_ref * (D / D_ref)**time_power",
            "reference_depth": args.reference_depth,
            "reference_time": args.reference_time,
            "time_power": args.time_power,
            "primary_dt": args.primary_dt,
            "gap_definition": "overlap-tracked logical-branch isolation gap",
        },
        "deep_limit_rows": primary_rows,
        "fixed_limit_fit": fit_deep_limit(primary_rows),
        "dt_convergence_rows": convergence_rows,
        "dt_convergence_against_finest": convergence_table(
            convergence_rows, min(args.timesteps)
        ),
        "runtime_s": time.time() - started,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"saved {summary_path}", flush=True)


if __name__ == "__main__":
    main()
