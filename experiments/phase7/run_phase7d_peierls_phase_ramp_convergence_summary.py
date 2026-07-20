"""Summarize phase-segmentation convergence of the CUDA Peierls-ramp audit."""
from __future__ import annotations

import json
from math import log2
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "phase7"
FILES = {
    2: RESULTS / "peierls_phase_ramp_gpu_audit.json",
    4: RESULTS / "peierls_phase_ramp_gpu_4step_refinement.json",
    8: RESULTS / "peierls_phase_ramp_gpu_8step_refinement.json",
}
FRACTIONS = (0.005, 0.02)


def row_for(data: dict, fraction: float) -> dict:
    for row in data["rows"]:
        if abs(row["ramp_fraction_of_each_subcycle"] - fraction) < 1e-15:
            return row
    raise RuntimeError(f"missing fraction {fraction}")


def order(values: list[float]) -> dict:
    first, second = abs(values[0] - values[1]), abs(values[1] - values[2])
    return {
        "difference_2_to_4": first,
        "difference_4_to_8": second,
        "ratio": first / second if second else None,
        "observed_order": log2(first / second) if first and second else None,
    }


def main() -> None:
    data = {steps: json.loads(path.read_text(encoding="utf-8")) for steps, path in FILES.items()}
    rows = []
    for fraction in FRACTIONS:
        selected = {steps: row_for(dataset, fraction) for steps, dataset in data.items()}
        leakage = [selected[steps]["monomer_leakage"] for steps in (2, 4, 8)]
        parity = [selected[steps]["logical_parity_a_residual"] for steps in (2, 4, 8)]
        rows.append({
            "ramp_fraction_of_each_subcycle": fraction,
            "values_by_phase_segments": {str(steps): selected[steps] for steps in (2, 4, 8)},
            "leakage_convergence": order(leakage),
            "parity_convergence": order(parity),
        })
    out = {
        "schema": "antler.phase7d.peierls-phase-ramp-segmentation-convergence.v1",
        "source_files": {str(steps): str(path.relative_to(ROOT)).replace("\\", "/") for steps, path in FILES.items()},
        "rows": rows,
        "decision": (
            "The 2/4/8 segment sequence is a numerical convergence audit of the ramp representation. It does not validate physical "
            "waveform bandwidth, calibration, noise, or any many-body topological property."
        ),
    }
    path = RESULTS / "peierls_phase_ramp_segmentation_convergence.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
