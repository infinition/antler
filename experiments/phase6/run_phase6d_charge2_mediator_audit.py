"""Phase 6D: exact branch-parity audit for a charge-two mediator candidate."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.native_charge2 import (
    branch_parity_operator,
    build_charge2_mediator_block,
    low_pair_masks,
)


PARAMETERS = {"Delta": 5.0, "V_mixed": 3.0}
G_VALUES = (0.05, 0.075, 0.10, 0.15, 0.20, 0.30)


def n1_cross_norm(g: float) -> float:
    H, states, _ = build_charge2_mediator_block(total_charge=1, g=g, **PARAMETERS)
    low_a = [i for i, state in enumerate(states) if int(state) in (1, 2)]
    low_b = [i for i, state in enumerate(states) if int(state) in (4, 8)]
    return float(np.linalg.norm(H[np.ix_(low_a, low_b)], ord="fro"))


def row(g: float) -> dict:
    H, states, index = build_charge2_mediator_block(total_charge=2, g=g, **PARAMETERS)
    values, vectors = eigh(H, subset_by_index=[0, 2], driver="evr")
    pair_a, pair_b = low_pair_masks()
    captures = [
        float(abs(vectors[index[pair_a], column]) ** 2 + abs(vectors[index[pair_b], column]) ** 2)
        for column in range(2)
    ]
    commutators = {
        branch: float(np.linalg.norm(H @ parity - parity @ H, ord="fro"))
        for branch, parity in (("a", branch_parity_operator(states, "a")), ("b", branch_parity_operator(states, "b")))
    }
    return {
        "g": g,
        "lowest_energies": values.tolist(),
        "pair_splitting": float(values[1] - values[0]),
        "pair_isolation_gap": float(values[2] - values[1]),
        "pair_subspace_captures": captures,
        "minimum_pair_subspace_capture": float(min(captures)),
        "branch_parity_commutator_norms": commutators,
        "one_particle_cross_norm": n1_cross_norm(g),
    }


def power_fit(rows: list[dict]) -> dict:
    x = np.log(np.asarray([item["g"] for item in rows]))
    y = np.log(np.asarray([item["pair_splitting"] for item in rows]))
    power, intercept = np.polyfit(x, y, 1)
    fit = power * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "power": float(power), "prefactor": float(np.exp(intercept)),
        "r_squared_log": float(1.0 - np.sum((y - fit) ** 2) / total),
    }


def main() -> None:
    rows = [row(g) for g in G_VALUES]
    fit = power_fit(rows)
    qualifies = (
        max(max(item["branch_parity_commutator_norms"].values()) for item in rows) < 1e-12
        and max(item["one_particle_cross_norm"] for item in rows) < 1e-12
        and min(item["pair_splitting"] for item in rows) > 1e-8
        and min(item["pair_isolation_gap"] for item in rows) > 0.1
        and min(item["minimum_pair_subspace_capture"] for item in rows) > 0.99
        and 1.5 < fit["power"] < 2.5
    )
    out = {
        "schema": "antler.phase6.charge2-mediator-local-audit.v1",
        "model": {
            "name": "charge-two molecular mediator with local pair conversion",
            "modes": ["a0", "a1", "b0", "b1", "d"],
            "mode_charges": [1, 1, 1, 1, 2],
            "parameters": PARAMETERS,
            "interaction": "-g[d^dagger(a0 a1 + b0 b1)+h.c.]",
        },
        "requirements": {
            "exact_branch_parity": "both branch-parity commutator norms below 1e-12",
            "no_single_particle_transfer": "charge-one a/b cross norm below 1e-12",
            "derived_pair_transfer": "nonzero isolated charge-two pair splitting with second-order g power",
            "qualification": bool(qualifies),
        },
        "rows": rows,
        "pair_splitting_power_fit": fit,
        "claim_boundary": (
            "This is a local candidate extension with a new charge-two degree of freedom. It is not yet an ANTLER lattice derivation, a protected code, "
            "a topological phase, a braid calculation, or an experimental realization of atom--molecule conversion."
        ),
        "decision": (
            "local exact-parity mechanism qualifies: next derive a tiled bond-mediator ladder and run the protection preflight"
            if qualifies else
            "local exact-parity mechanism fails; do not tile it into a ladder"
        ),
    }
    path = ROOT / "results" / "phase6" / "charge2_mediator_local_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "qualification": qualifies, "decision": out["decision"], "pair_power": fit,
        "max_parity_commutator": max(max(item["branch_parity_commutator_norms"].values()) for item in rows),
        "max_single_particle_cross": max(item["one_particle_cross_norm"] for item in rows),
    }, indent=2))


if __name__ == "__main__":
    main()
