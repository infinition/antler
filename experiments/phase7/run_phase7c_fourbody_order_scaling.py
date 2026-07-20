"""Measure the coupling order of representative connected four-edge channels.

The calculation tests the order-counting claim independently of a language
model.  It uses one best fixed-scale bare topology from the exhaustive catalog
and the submitted alternating mixed ring (with its parity metadata corrected),
over the registered ``g/Delta <= 0.15`` window.
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
    MicroscopicCandidate2D,
    PairConversionChannel,
    evaluate_local_candidate,
    rail_mode,
)


G_VALUES = (0.30, 0.40, 0.50, 0.70, 0.90, 1.10, 1.30, 1.50)
U, DELTA = 20.0, 10.0


def channel(kind: str, left: int, right: int, g: float) -> PairConversionChannel:
    if kind == "E":
        pairs = ((rail_mode(left, 0), rail_mode(right, 0), 1.0), (rail_mode(left, 1), rail_mode(right, 1), 1.0))
    elif kind == "M":
        pairs = ((rail_mode(left, 0), rail_mode(right, 1), 1.0), (rail_mode(left, 1), rail_mode(right, 0), 1.0))
    else:
        raise ValueError("unknown channel kind")
    return PairConversionChannel(f"{kind}_{left}{right}", DELTA, g, 0.0, pairs)


TOPOLOGIES = {
    "best_bare_catalog": (("E", 0, 2), ("E", 1, 2), ("E", 1, 3), ("E", 2, 3)),
    "mixed_species_ring_counterfactual": (("E", 0, 1), ("M", 1, 2), ("E", 2, 3), ("M", 3, 0)),
}


def fit_power(rows: list[dict], field: str, floor: float = 1e-12) -> dict:
    usable = [row for row in rows if abs(row[field]) > floor]
    if len(usable) < 3:
        return {"status": "unresolved_below_numerical_floor", "usable_points": len(usable)}
    x = np.log(np.asarray([row["g"] for row in usable]))
    y = np.log(np.asarray([abs(row[field]) for row in usable]))
    power, intercept = np.polyfit(x, y, 1)
    fitted = power * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "status": "fit",
        "usable_points": len(usable),
        "power": float(power),
        "prefactor": float(np.exp(intercept)),
        "r_squared_log": float(1.0 - np.sum((y - fitted) ** 2) / total) if total > 1e-30 else 1.0,
    }


def main() -> None:
    output = {}
    for name, topology in TOPOLOGIES.items():
        rows = []
        for g in G_VALUES:
            candidate = MicroscopicCandidate2D(
                mott_u=U, maximum_sw_ratio=0.15,
                channels=tuple(channel(kind, left, right, g) for kind, left, right in topology),
            )
            audit = evaluate_local_candidate(candidate, target_label="XXXX")
            state = audit["state_vector"]
            rows.append({
                "g": g,
                "g_over_delta": g / DELTA,
                "xxxx_coefficient": state["observed_target_coefficient"],
                "unwanted_pauli_norm_over_target": state["unwanted_pauli_norm_over_target"],
                "operator_residual": state["fixed_scale_spectral_algebraic_residual"],
                "reward": audit["reward"],
            })
        output[name] = {
            "topology": topology,
            "rows": rows,
            "xxxx_power_fit": fit_power(rows, "xxxx_coefficient"),
            "unwanted_norm_power_fit": fit_power(rows, "unwanted_pauli_norm_over_target"),
        }
    out = {
        "schema": "antler.phase7c.fourbody-order-scaling.v1",
        "parameters": {"U": U, "Delta": DELTA, "g_values": G_VALUES, "maximum_g_over_delta": 0.15},
        "topologies": output,
        "decision": (
            "This is an independent finite-block order measurement. A power near eight supports the independent-mediator "
            "ring order count; it does not by itself establish a universal all-parameter no-go bound."
        ),
        "claim_boundary": (
            "Only two representative topologies are scaled. This does not prove a global native-ANTLER no-go, a controlled "
            "Schrieffer-Wolff theorem, a tiled parent, topological order or a braid result."
        ),
    }
    path = ROOT / "results" / "phase7" / "fourbody_order_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        name: {"xxxx": item["xxxx_power_fit"], "unwanted": item["unwanted_norm_power_fit"]}
        for name, item in output.items()
    }, indent=2))


if __name__ == "__main__":
    main()
