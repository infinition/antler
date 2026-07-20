"""Phase 6A: initial protection scan for a native scalar three-leg ANTLER family.

The scan is intentionally performed before any Wilson line, braid or gate is
introduced.  A candidate must first exhibit an isolated ground doublet that
is not locally distinguishable by the physical site densities.
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
    split = float(values[1] - values[0])
    isolation_gap = float(values[2] - values[1])
    probes = []
    for site in range(3 * L):
        density = local_density_operator(states, site)
        projected = code.conj().T @ (density[:, None] * code)
        traceless = projected - 0.5 * np.trace(projected) * np.eye(2)
        probes.append(float(np.linalg.norm(traceless)))
    return {
        "L": L,
        "N": N,
        "parameters": parameters,
        "lowest_energies": values.tolist(),
        "logical_split": split,
        "isolation_gap": isolation_gap,
        "split_over_gap": float(split / isolation_gap) if isolation_gap > 0 else None,
        "site_density_local_distinguishability": probes,
        "max_site_density_local_distinguishability": float(max(probes)),
    }


def score(row: dict) -> tuple[float, float, float]:
    ratio = row["split_over_gap"] if row["split_over_gap"] is not None else np.inf
    return ratio, row["max_site_density_local_distinguishability"], -row["isolation_gap"]


def main() -> None:
    # No edge traps or imposed matrix links are used.  These are short-range,
    # scalar native terms selected before looking at spectra.
    family = [
        {"theta": 0.3, "J1": J1, "J2": 1.0, "Jperp": Jperp,
         "U_rung": U_rung, "V_leg": V_leg}
        for J1 in (0.2, 0.4)
        for Jperp in (0.05, 0.2)
        for U_rung in (0.0, 1.0)
        for V_leg in (-0.5, 0.5)
    ]
    rows = []
    for L, N in ((4, 2), (5, 3)):
        rows.extend(analyse(L, N, parameters) for parameters in family)
    candidates = [
        row for row in rows
        if row["split_over_gap"] is not None and row["split_over_gap"] < 1e-3
        and row["max_site_density_local_distinguishability"] < 0.1
        and row["isolation_gap"] > 1e-2
    ]
    best = sorted(rows, key=score)[:8]
    out = {
        "schema": "antler.phase6.native-threeleg-preflight.v1",
        "model": {
            "name": "native scalar three-leg correlated-hopping extension",
            "conserved_quantity": "total U(1) particle number",
            "ordering": "rung-major (i,leg=0,1,2)",
            "terms": ["scalar correlated hopping", "nearest rung hopping", "local rung density interaction", "nearest-leg density interaction"],
        },
        "claim_boundary": (
            "This is an initial native ANTLER candidate family, not a derivation from a microscopic experiment and not a braid implementation. "
            "A local-density preflight can falsify a protection claim but cannot by itself prove topological order."
        ),
        "criteria": {
            "split_over_gap_below": 1e-3,
            "max_site_density_local_distinguishability_below": 0.1,
            "isolation_gap_above": 1e-2,
        },
        "rows": rows,
        "best_rows_by_preflight_score": best,
        "candidates": candidates,
        "decision": (
            "native candidate(s) pass the first local-density filter only; next require L-scaling and a broader local-operator audit before any braid"
            if candidates else
            "no protected-doublet candidate in this pre-registered minimal native three-leg family; do not launch braid or gate calculations"
        ),
    }
    path = Path("results/phase6/native_threeleg_preflight.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "best_rows": best, "candidate_count": len(candidates)}, indent=2))


if __name__ == "__main__":
    main()
