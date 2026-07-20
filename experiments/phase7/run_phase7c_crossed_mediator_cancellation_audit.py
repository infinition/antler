"""Audit the submitted crossed charge-two mediator with calibrated static counterterms.

One mediator couples two disjoint same-branch pairs on a four-edge plaquette.
The second mediator supplies the reversed configuration.  Static Z and ZZ
counterterms are fitted from the uncorrected finite-block Pauli expansion at
each coupling value, then reapplied to test whether the four-body target is
both enhanced and selective.  This is a local calibration, not an analytic SW
proof or a global tiling claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    EDGE_COUNT,
    MicroscopicCandidate2D,
    PairConversionChannel,
    ZZCoupling,
    evaluate_local_candidate,
    interaction_connectivity,
    rail_mode,
)


U, DELTA = 20.0, 10.0
G_VALUES = (0.30, 0.40, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50)


def crossed_channels(g: float) -> tuple[PairConversionChannel, ...]:
    return (
        PairConversionChannel(
            "cross_aa01_bb23", DELTA, g, 0.0,
            ((rail_mode(0, 0), rail_mode(1, 0), 1.0), (rail_mode(2, 1), rail_mode(3, 1), 1.0)),
        ),
        PairConversionChannel(
            "cross_bb01_aa23", DELTA, g, 0.0,
            ((rail_mode(0, 1), rail_mode(1, 1), 1.0), (rail_mode(2, 0), rail_mode(3, 0), 1.0)),
        ),
    )


def base_candidate(g: float) -> MicroscopicCandidate2D:
    return MicroscopicCandidate2D(mott_u=U, maximum_sw_ratio=0.15, channels=crossed_channels(g))


def static_counterterms(coefficients: dict[str, float]) -> tuple[tuple[float, ...], tuple[ZZCoupling, ...]]:
    biases = [0.0] * EDGE_COUNT
    couplings = []
    for label, coefficient in coefficients.items():
        z_sites = [index for index, letter in enumerate(label) if letter == "Z"]
        if any(letter not in {"I", "Z"} for letter in label):
            continue
        if len(z_sites) == 1:
            biases[z_sites[0]] = -coefficient
        elif len(z_sites) == 2:
            couplings.append(ZZCoupling(z_sites[0], z_sites[1], -coefficient))
    return tuple(biases), tuple(couplings)


def fit_power(rows: list[dict], field: str) -> dict:
    usable = [row for row in rows if abs(row[field]) > 1e-12]
    x = np.log(np.asarray([row["g"] for row in usable]))
    y = np.log(np.asarray([abs(row[field]) for row in usable]))
    power, intercept = np.polyfit(x, y, 1)
    fitted = power * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "power": float(power), "prefactor": float(np.exp(intercept)), "usable_points": len(usable),
        "r_squared_log": float(1.0 - np.sum((y - fitted) ** 2) / total) if total > 1e-30 else 1.0,
    }


def summary(audit: dict) -> dict:
    state = audit["state_vector"]
    alignment = state["target_alignment"]
    other_over_observed = float(np.sqrt((1.0 - alignment) / alignment)) if alignment > 0.0 else float("inf")
    return {
        "xxxx_coefficient": state["observed_target_coefficient"],
        "target_alignment": alignment,
        "operator_residual": state["fixed_scale_spectral_algebraic_residual"],
        "unwanted_pauli_norm_over_fixed_target_scale": state["unwanted_pauli_norm_over_target"],
        "other_pauli_norm_over_observed_xxxx": other_over_observed,
        "unwanted_pauli_norm_over_target": state["unwanted_pauli_norm_over_target"],
        "minimum_monomer_overlap": state["minimum_monomer_overlap_singular_value"],
        "low_high_gap": state["low_to_high_gap"],
        "reward": audit["reward"],
        "top_unwanted_paulis": audit["top_unwanted_paulis"],
    }


def main() -> None:
    rows = []
    reference_counterterms = None
    for g in G_VALUES:
        base = evaluate_local_candidate(base_candidate(g), "XXXX", include_pauli_coefficients=True)
        biases, zz = static_counterterms(base["full_traceless_pauli_coefficients"])
        compensated_candidate = MicroscopicCandidate2D(
            mott_u=U, maximum_sw_ratio=0.15, channels=crossed_channels(g), rail_biases=biases, zz_couplings=zz,
        )
        compensated = evaluate_local_candidate(compensated_candidate, "XXXX")
        rows.append({
            "g": g,
            "g_over_delta": g / DELTA,
            "base": summary(base),
            "counterterms": {
                "rail_biases": list(biases),
                "zz_couplings": [coupling.__dict__ for coupling in zz],
            },
            "compensated": summary(compensated),
        })
        if g == 0.5:
            reference_counterterms = rows[-1]
    base_rows = [{"g": row["g"], "coefficient": row["base"]["xxxx_coefficient"], "unwanted": row["base"]["unwanted_pauli_norm_over_target"]} for row in rows]
    compensated_rows = [{"g": row["g"], "coefficient": row["compensated"]["xxxx_coefficient"], "unwanted": row["compensated"]["unwanted_pauli_norm_over_target"]} for row in rows]
    connectivity = interaction_connectivity(base_candidate(0.5))
    if not connectivity["patch_connected"]:
        raise RuntimeError("crossed shared-mediator primitive must be treated as a connected plaquette hyperedge")
    out = {
        "schema": "antler.phase7c.crossed-charge2-cancellation-audit.v1",
        "primitive": {
            "description": "two charge-two mediators, each coherently coupled to two disjoint same-branch pairs on one four-edge plaquette",
            "parameters": {"U": U, "Delta": DELTA, "g_values": G_VALUES},
            "connectivity": connectivity,
        },
        "rows": rows,
        "power_fits": {
            "base_xxxx": fit_power(base_rows, "coefficient"),
            "base_unwanted": fit_power(base_rows, "unwanted"),
            "compensated_xxxx": fit_power(compensated_rows, "coefficient"),
            "compensated_unwanted": fit_power(compensated_rows, "unwanted"),
        },
        "g_0_5_reference": reference_counterterms,
        "decision": (
            "The crossed mediator is a genuine new local primitive: it yields a much larger four-body component and calibrated "
            "static Z/ZZ counterterms suppress its second-order diagonal content. It is not promoted because the remaining "
            "four-body operator is not pure XXXX and its absolute coefficient is still far below -0.5."
        ),
        "claim_boundary": (
            "Counterterms are fitted from the same finite block at each g. This is a numerical local calibration, not an "
            "independent analytic SW proof, a fixed microscopic implementation, a scalable parent, topological order or braiding."
        ),
    }
    path = ROOT / "results" / "phase7" / "crossed_charge2_cancellation_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"power_fits": out["power_fits"], "g_0_5_compensated": reference_counterterms["compensated"]}, indent=2))


if __name__ == "__main__":
    main()
