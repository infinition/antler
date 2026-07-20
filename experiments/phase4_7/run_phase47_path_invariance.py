"""Phase 4.7B -- topology-preserving digital path-deformation campaign.

Every variant keeps a sequential (never simultaneous) exchange and the same
initial/final Hamiltonian.  It changes only the ramp profile, parametrisation,
parking wells, distance, or temporal ordering of the two independent
handoffs.  The output is a falsifiable invariance table, not a pass/fail
assertion hidden in the code.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from phase47_common import (DigitalConfig, matrix_from_json, remove_global,
                            require_json, run_gate)


def wrapped(angle: float) -> float:
    return float(np.angle(np.exp(1j * angle)))


def variants(base: DigitalConfig) -> dict[str, DigitalConfig]:
    """One-control-at-a-time deformations around the Phase 4.7 baseline."""

    return {
        "baseline": base,
        "ramp_linear": replace(base, ramp="linear"),
        "ramp_smoothstep": replace(base, ramp="smoothstep"),
        "handoff_short": replace(base, handoff_fraction=0.075),
        "handoff_long": replace(base, handoff_fraction=0.125),
        "intermediate_pauses": replace(base, pause_fraction=0.03),
        "parking_shallower": replace(base, parking_depth_scale=0.90),
        "parking_deeper": replace(base, parking_depth_scale=1.10),
        "spectator_shallower": replace(base, spectator_depth_scale=0.90),
        "distance_R3": replace(base, R_LOOP=3),
        "distance_R5": replace(base, R_LOOP=5),
        "rung_then_left": replace(base, handoff_order="rung_then_left"),
    }


def compatible(cached: dict, cfg: DigitalConfig, theta: float, dt: float) -> bool:
    return (
        cached.get("schema") == "antler.phase47.gate.v1"
        and cached.get("config") == asdict(cfg)
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
    print(f"run {path.stem}: {asdict(cfg)}", flush=True)
    started = time.time()
    result = run_gate(cfg, theta=theta, dt=dt, include_gap=True,
                      gap_samples=gap_samples)
    result["runtime_s"] = time.time() - started
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["metrics"], indent=2), flush=True)
    return result


def compact(name: str, result: dict, theta: float, reference: dict | None) -> dict:
    metrics = result["metrics"]
    gap = result.get("gap") or {}
    row = {
        "variant": name,
        "odd_phase": metrics["odd_phase"],
        "odd_slope": metrics["odd_slope"],
        "phase_error_to_minus_theta": wrapped(metrics["odd_phase"] + theta),
        "leak_worst": metrics["leak_worst"],
        "sigma_min": metrics["sigma_min"],
        "odd_offdiag_norm": metrics["odd_offdiag_norm"],
        "favg_target": metrics["favg_target"],
        "minimum_handoff_isolation_gap": gap.get("minimum_handoff_isolation_gap"),
        "gap_stage": gap.get("handoff_stage"),
        "config": result["config"],
    }
    if reference is not None:
        U = matrix_from_json(result["Uodd"])
        U0 = matrix_from_json(reference["Uodd"])
        comparison = remove_global(U0.conj().T @ U)
        row.update({
            "delta_phase_to_baseline": wrapped(metrics["odd_phase"] - reference["metrics"]["odd_phase"]),
            "unitary_distance_to_baseline": float(np.linalg.norm(comparison - np.eye(2))),
            "offdiag_of_relative_unitary": float(
                np.linalg.norm(comparison - np.diag(np.diag(comparison)))
            ),
        })
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--depth", type=float, default=8.0,
                        help="positive trap-depth magnitude")
    parser.add_argument("--T", type=float, default=80_000.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--R", type=int, default=4)
    parser.add_argument("--gap-samples", type=int, default=61)
    parser.add_argument("--variants", nargs="+", default=["all"],
                        help="variant names, or the single token 'all'")
    parser.add_argument("--outdir", type=Path,
                        default=Path("results/phase4_7/path_invariance"))
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.depth <= 0 or args.T <= 0 or args.dt <= 0:
        raise ValueError("depth, T, and dt must be positive")
    if not 1 <= args.R < 14:
        raise ValueError("R must be an accessible ladder rung")
    base = DigitalConfig(DEPTH=-args.depth, R_LOOP=args.R, T_TOTAL=args.T)
    all_variants = variants(base)
    selected = list(all_variants) if args.variants == ["all"] else args.variants
    unknown = sorted(set(selected) - set(all_variants))
    if unknown:
        raise ValueError(f"unknown variants: {', '.join(unknown)}")
    # Baseline is always calculated first, because every comparison is relative
    # to the same fully reconstructed odd logical gate.
    ordered = ["baseline"] + [name for name in selected if name != "baseline"]
    started = time.time()
    results: dict[str, dict] = {}
    for name in ordered:
        path = args.outdir / f"{name}.json"
        results[name] = run_or_resume(path, all_variants[name], args.theta, args.dt,
                                      args.gap_samples, resume=not args.no_resume)
    baseline = results["baseline"]
    rows = [compact(name, results[name], args.theta,
                    None if name == "baseline" else baseline) for name in ordered]
    summary = {
        "schema": "antler.phase47.path-invariance.v1",
        "claim_boundary": (
            "These are digital-protocol deformations for the frozen correlated-hopping "
            "ladder.  Agreement would support parametrisation invariance at finite "
            "depth; it does not by itself establish topological protection or "
            "non-Abelian statistics."
        ),
        "theta": args.theta,
        "baseline_config": asdict(base),
        "rows": rows,
        "max_abs_phase_error_to_minus_theta": float(
            max(abs(row["phase_error_to_minus_theta"]) for row in rows)
        ),
        "max_abs_phase_shift_from_baseline": float(max(
            abs(row.get("delta_phase_to_baseline", 0.0)) for row in rows
        )),
        "max_unitary_distance_from_baseline": float(max(
            row.get("unitary_distance_to_baseline", 0.0) for row in rows
        )),
        "runtime_s": time.time() - started,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    path = args.outdir / "summary.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"saved {path}", flush=True)


if __name__ == "__main__":
    main()
