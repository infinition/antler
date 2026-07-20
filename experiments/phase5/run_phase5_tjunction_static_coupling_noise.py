"""Phase 5J: ensemble robustness of the finite-time T-junction exchange.

Each realization applies independent, static multiplicative calibration errors
to the three central-to-arm Majorana couplings.  The endpoints and exchange
order are retained, while the path through coupling space and its gap change.
This is an effective-Majorana-network audit, deliberately separate from the
number-conserving ANTLER ladder and from a microscopic noise model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from run_phase5_tjunction_braid_dynamics import braid_ensemble


def positive_scales(rng: np.random.Generator, sigma: float, samples: int) -> np.ndarray:
    """Draw physical positive coupling scales from the stated rms error."""

    scales = 1.0 + sigma * rng.standard_normal((samples, 3))
    while np.any(scales <= 0.0):
        bad = np.any(scales <= 0.0, axis=1)
        scales[bad] = 1.0 + sigma * rng.standard_normal((int(bad.sum()), 3))
    return scales


def summarize(rows: list[dict], sigma: float) -> dict:
    values = lambda key: np.asarray([row[key] for row in rows], float)
    return {
        "sigma": sigma,
        "samples": len(rows),
        "favg_mean": float(values("favg_target").mean()),
        "favg_variance": float(values("favg_target").var(ddof=1)),
        "leak_worst_mean": float(values("leak_worst").mean()),
        "leak_worst_variance": float(values("leak_worst").var(ddof=1)),
        "phase_matrix_error_mean": float(values("matrix_error").mean()),
        "phase_matrix_error_variance": float(values("matrix_error").var(ddof=1)),
        "minimum_gap_mean": float(values("minimum_even_parity_gap").mean()),
        "minimum_gap_worst": float(values("minimum_even_parity_gap").min()),
        "favg_min": float(values("favg_target").min()),
        "leak_worst_max": float(values("leak_worst").max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sigmas", type=float, nargs="+",
                        default=[0.01, 0.03, 0.05, 0.10, 0.20])
    parser.add_argument("--samples", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--T", type=float, default=80.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/tjunction_static_coupling_noise.json"))
    args = parser.parse_args()
    if args.samples < 2 or any(sigma < 0.0 for sigma in args.sigmas):
        raise ValueError("need at least two samples and non-negative sigmas")

    rng = np.random.default_rng(args.seed)
    ensembles = []
    for sigma in args.sigmas:
        rows = []
        results = braid_ensemble(args.T, args.dt, positive_scales(rng, sigma, args.samples))
        for sample, result in enumerate(results):
            rows.append({
                "sample": sample,
                "arm_scales": result["arm_scales"],
                "favg_target": result["favg_target"],
                "leak_worst": result["leak_worst"],
                "matrix_error": result["matrix_error"],
                "minimum_even_parity_gap": result["minimum_even_parity_gap"],
            })
        ensembles.append({"summary": summarize(rows, sigma), "realizations": rows})
        print(json.dumps(ensembles[-1]["summary"], indent=2), flush=True)

    output = {
        "schema": "antler.phase5.tjunction-static-coupling-noise.v1",
        "claim_boundary": (
            "This is a static calibration-noise ensemble for the effective Majorana "
            "T-junction exchange.  It neither derives pairing from the ANTLER ladder "
            "nor establishes a device-level topological noise threshold."
        ),
        "noise_model": "independent positive scales 1 + sigma * Normal(0,1) on the three arms",
        "T": args.T, "dt": args.dt, "seed": args.seed, "ensembles": ensembles,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
