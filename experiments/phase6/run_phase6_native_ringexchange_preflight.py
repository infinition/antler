"""Phase 6B: scalar three-leg plaquette ring-exchange protection scan.

The added term is an oriented local exchange around native three-leg plaquettes.
It is deliberately not an Iemini-like two-wire pair-transfer parent term.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.native_threeleg import build_native_threeleg_hamiltonian, local_density_operator


def analyse(L: int, N: int, parameters: dict) -> dict:
    H, states, _ = build_native_threeleg_hamiltonian(L=L, N=N, **parameters)
    values, vectors = eigh(H, subset_by_index=[0, 2], driver="evr")
    code = vectors[:, :2]
    probes = []
    for site in range(3 * L):
        density = local_density_operator(states, site)
        projected = code.conj().T @ (density[:, None] * code)
        probes.append(float(np.linalg.norm(projected - 0.5 * np.trace(projected) * np.eye(2))))
    split = float(values[1] - values[0])
    gap = float(values[2] - values[1])
    return {
        "L": L, "N": N, "parameters": parameters,
        "lowest_energies": values.tolist(), "logical_split": split,
        "isolation_gap": gap,
        "split_over_gap": float(split / gap) if gap > 0 else None,
        "max_site_density_local_distinguishability": float(max(probes)),
        "site_density_local_distinguishability": probes,
    }


def score(row: dict) -> tuple[float, float, float]:
    return (row["split_over_gap"] if row["split_over_gap"] is not None else np.inf,
            row["max_site_density_local_distinguishability"], -row["isolation_gap"])


def main() -> None:
    family = [
        {"theta": 0.3, "J1": J1, "J2": 1.0, "Jperp": Jperp,
         "U_rung": U_rung, "V_leg": V_leg, "K_ring": K_ring,
         "phi_ring": phi_ring}
        for J1 in (0.2, 0.4)
        for Jperp in (0.05, 0.2)
        for U_rung in (0.0, 1.0)
        for V_leg in (-0.5, 0.5)
        for K_ring in (0.5, 1.0)
        for phi_ring in (0.3, np.pi / 2.0)
    ]
    rows = [analyse(L, N, parameters) for L, N in ((4, 2), (5, 3)) for parameters in family]
    candidates = [
        row for row in rows
        if row["split_over_gap"] is not None and row["split_over_gap"] < 1e-3
        and row["max_site_density_local_distinguishability"] < 0.1
        and row["isolation_gap"] > 1e-2
    ]
    best = sorted(rows, key=score)[:8]
    out = {
        "schema": "antler.phase6.native-ringexchange-preflight.v1",
        "model": {
            "name": "native scalar three-leg correlated hopping plus chiral plaquette exchange",
            "conserved_quantity": "total U(1) particle number",
            "not_iemini": "the plaquette exchange preserves the native three-leg geometry and is not the two-wire parent pair-transfer term",
        },
        "claim_boundary": (
            "This scans a new local scalar ansatz. It is not a microscopic derivation, a protected code, or a braid calculation. "
            "A failed local-density criterion rejects a candidate immediately."
        ),
        "criteria": {"split_over_gap_below": 1e-3, "max_site_density_local_distinguishability_below": 0.1, "isolation_gap_above": 1e-2},
        "rows": rows, "best_rows_by_preflight_score": best, "candidates": candidates,
        "decision": (
            "candidate(s) only pass the first local-density filter; require L scaling and a complete local-operator basis before braid work"
            if candidates else
            "no protected-doublet candidate in the scanned native ring-exchange family; do not launch braid or gate calculations"
        ),
    }
    path = Path("results/phase6/native_ringexchange_preflight.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "best_rows": best, "candidate_count": len(candidates)}, indent=2))


if __name__ == "__main__":
    main()
