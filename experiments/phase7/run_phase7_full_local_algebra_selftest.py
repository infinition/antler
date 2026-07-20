"""Self-test the exhaustive local algebra gate on the archived Ising benchmark."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_ising_parent import build_weighted_basis_fast, code_frame
from antler.phase7_local_algebra import (
    charge_conserving_local_probes,
    exhaustive_local_code_metric,
)


def main() -> None:
    L = 4
    states, index = build_weighted_basis_fast(L, L)
    G = code_frame(L, states, index)
    charges = (1,) * (2 * L) + (2,) * (L - 1)
    fermions = frozenset(range(2 * L))
    one_rung = charge_conserving_local_probes(
        states, index, charges, (0, 1), fermions, label_prefix="rung0",
    )
    rung_and_mediator = charge_conserving_local_probes(
        states, index, charges, (0, 1, 2 * L), fermions, label_prefix="cell0",
    )
    rung_metric = exhaustive_local_code_metric(G, one_rung)
    cell_metric = exhaustive_local_code_metric(G, rung_and_mediator)
    hermiticity = max(
        float(np.linalg.norm(operator - operator.conj().T))
        for operator in (*one_rung.values(), *rung_and_mediator.values())
    )
    passed = (
        hermiticity < 1e-12
        and rung_metric["worst_projected_non_scalar_frobenius"] > 1e-3
        and cell_metric["worst_projected_non_scalar_frobenius"] > 1e-3
    )
    if not passed:
        raise RuntimeError("full local algebra self-test did not expose the Ising logical operator")
    out = {
        "schema": "antler.phase7.full-local-algebra-selftest.v1",
        "reference": "archived fixed-point fixed-point Ising/cat benchmark at L=4,Q=4",
        "one_rung_complete_charge_conserving_algebra": rung_metric,
        "rung_plus_charge_two_mediator_complete_charge_conserving_algebra": cell_metric,
        "max_probe_hermiticity_error": hermiticity,
        "passes_negative_control": bool(passed),
        "decision": (
            "The exhaustive physical local algebra exposes a non-scalar logical probe already on one rung. "
            "Future candidates must pass this gate before any gap or edge-mode result is promoted."
        ),
        "claim_boundary": "This is an audit-infrastructure control, not a new Hamiltonian or a topological-code claim.",
    }
    path = ROOT / "results" / "phase7" / "full_local_algebra_selftest.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
