"""Phase 5A.1: optimise the N=3 mediator candidate before dynamics.

The first N=3 scan can find a doublet with excellent capture but an isolation
gap too small for a useful adiabatic braid.  This scan maps that trade-off over
code-well depth, mediator position, and mediator depth.  It is deliberately a
static gatekeeper: no time-dependent braid is claimed until both localisation
and an absolute isolation gap are adequate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
for path in (ROOT, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from run_phase5_n3_mediator_preflight import N3Preflight, PreflightConfig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--code-depths", type=float, nargs="+",
                        default=[-2.5, -3.0, -3.5, -4.0, -5.0])
    parser.add_argument("--rungs", type=int, nargs="+", default=[5, 7, 9])
    parser.add_argument("--depth-scales", type=float, nargs="+", default=[1.0, 1.25, 1.5])
    parser.add_argument("--n-eigs", type=int, default=12)
    parser.add_argument("--minimum-gap", type=float, default=1e-3,
                        help="minimum absolute isolation gap for dynamic compilation")
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/n3_mediator_gap_optimization.json"))
    args = parser.parse_args()
    if any(depth >= 0 for depth in args.code_depths):
        raise ValueError("all code depths must be negative")
    if any(scale <= 0 for scale in args.depth_scales):
        raise ValueError("all mediator depth scales must be positive")

    started = time.time()
    rows = []
    total = len(args.code_depths) * len(args.rungs) * 2 * len(args.depth_scales)
    count = 0
    for code_depth in args.code_depths:
        for rung in args.rungs:
            for leg in (0, 1):
                for scale in args.depth_scales:
                    count += 1
                    cfg = PreflightConfig(code_depth=code_depth, mediator_rung=rung,
                                          mediator_leg=leg, mediator_depth_scale=scale)
                    print(f"scan {count}/{total}: D={code_depth:g}, r={rung}, leg={leg}, scale={scale:g}",
                          flush=True)
                    system = N3Preflight(cfg)
                    spectra = [system.analyse(theta, args.n_eigs) for theta in (0.0, args.theta)]
                    capture = min(scan["capture"] for scan in spectra)
                    mediator_occ = min(min(scan["mediator_occupation_selected"]) for scan in spectra)
                    gap = min(scan["isolation_gap_to_low_energy_complement"] for scan in spectra)
                    split = max(scan["logical_split"] for scan in spectra)
                    rows.append({
                        "config": cfg.__dict__, "spectra": spectra,
                        "capture_min": capture,
                        "mediator_occupation_min": mediator_occ,
                        "isolation_gap_min": gap,
                        "logical_split_max": split,
                        "gap_over_split": gap / split if split else float("inf"),
                        "dynamic_candidate": (
                            capture >= 0.80 and mediator_occ >= 0.80 and gap >= args.minimum_gap
                        ),
                    })
    rows.sort(key=lambda row: (
        not row["dynamic_candidate"], -row["isolation_gap_min"],
        -row["capture_min"], -row["mediator_occupation_min"],
    ))
    passing = [row for row in rows if row["dynamic_candidate"]]
    out = {
        "schema": "antler.phase5.n3-gap-optimization.v1",
        "criterion": {
            "capture_min": 0.80,
            "mediator_occupation_min": 0.80,
            "minimum_absolute_isolation_gap": args.minimum_gap,
        },
        "decision": (
            "compile_two_exchange_dynamics" if passing else
            "N3_pinned_mediator_has_no_adiabatically_practical_candidate_in_this_scan"
        ),
        "best": rows[0],
        "passing_rows": passing,
        "rows": rows,
        "runtime_s": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "best": out["best"]}, indent=2), flush=True)
    print(f"saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
