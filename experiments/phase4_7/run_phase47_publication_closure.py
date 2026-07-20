"""Assemble the all-fine-step Phase 4.7 deep-limit publication table."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


ROOT = Path(__file__).resolve().parents[2]
POINTS = (
    (4, ROOT / "results" / "phase4_7" / "deep_limit_refined_d4" / "D4_dt0p125.json"),
    (6, ROOT / "results" / "phase4_7" / "deep_limit" / "dt_convergence" / "D6_dt0p125.json"),
    (8, ROOT / "results" / "phase4_7" / "deep_limit" / "dt_convergence" / "D8_dt0p125.json"),
    (10, ROOT / "results" / "phase4_7" / "deep_limit_refined_d10" / "D10_dt0p125.json"),
    (12, ROOT / "results" / "phase4_7" / "deep_limit_refined_d12" / "D12_dt0p125.json"),
)


def main() -> None:
    missing = [str(path) for _, path in POINTS if not path.exists()]
    if missing:
        raise RuntimeError("publication closure awaits: " + ", ".join(missing))
    rows = []
    for depth, path in POINTS:
        payload = json.loads(path.read_text())
        metrics = payload["metrics"]
        rows.append({"depth": depth, "dt": payload["dt_requested"], "T": payload["config"]["T_TOTAL"], "odd_slope": metrics["odd_slope"], "leak_worst": metrics["leak_worst"], "sigma_min": metrics["sigma_min"], "offdiag": metrics["odd_offdiag_norm"], "minimum_handoff_isolation_gap": payload["gap"]["minimum_handoff_isolation_gap"], "file": path.relative_to(ROOT).as_posix()})
    x = np.asarray([row["depth"] for row in rows], float)
    y = np.asarray([row["odd_slope"] + 1.0 for row in rows], float)
    def model(depth, a, p): return a * depth ** (-p)
    (a, p), covariance = curve_fit(model, x, y, p0=(0.8, 2.0), maxfev=100000)
    predicted = model(x, a, p)
    residual = y - predicted
    r2 = 1.0 - float(np.sum(residual**2) / np.sum((y - y.mean())**2))
    output = {"schema": "antler.phase47.publication-closure.v1", "rows": rows, "fit": {"model": "odd_slope(D)=-1+a D^-p", "a": float(a), "p": float(p), "stderr": np.sqrt(np.diag(covariance)).tolist(), "r2": r2, "rmse": float(np.sqrt(np.mean(residual**2))), "residual": residual.tolist()}, "gates": {"all_points_dt_0p125": all(abs(row["dt"] - 0.125) < 1e-15 for row in rows), "leakage_below_1e-4": all(row["leak_worst"] < 1e-4 for row in rows), "sigma_min_above_0p9999": all(row["sigma_min"] > 0.9999 for row in rows), "positive_handoff_gap": all(row["minimum_handoff_isolation_gap"] > 0.0 for row in rows)}, "claim_boundary": "This closes only the all-fine-step numerical deep-limit table for the frozen correlated-hopping ladder. It does not establish non-Abelian statistics, universality, fault tolerance or a topological quantum computer."}
    target = ROOT / "results" / "phase4_7" / "publication_closure.json"
    target.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
