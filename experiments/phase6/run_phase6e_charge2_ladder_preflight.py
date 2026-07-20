"""Phase 6E: pre-registered protection scan of the charge-two mediator ladder.

The code is selected from the two parity sectors allowed at even total charge.
No braid or dynamics is attempted unless this static protection gate passes.
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

from antler.native_charge2_ladder import (
    branch_parities,
    build_charge2_mediator_ladder,
    local_density,
)


def sector_data(H: np.ndarray, states: np.ndarray, L: int, sector: tuple[int, int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = np.asarray([i for i, state in enumerate(states) if branch_parities(int(state), L) == sector], dtype=int)
    values, vectors = eigh(H[np.ix_(selected, selected)], subset_by_index=[0, 1], driver="evr")
    return selected, values, vectors


def analyse(L: int, total_charge: int, parameters: dict) -> dict:
    H, states, _ = build_charge2_mediator_ladder(L=L, total_charge=total_charge, **parameters)
    sector_even, energy_even, vector_even = sector_data(H, states, L, (1, 1))
    sector_odd, energy_odd, vector_odd = sector_data(H, states, L, (-1, -1))
    full_even = np.zeros(len(states), dtype=complex)
    full_odd = np.zeros(len(states), dtype=complex)
    full_even[sector_even] = vector_even[:, 0]
    full_odd[sector_odd] = vector_odd[:, 0]
    probes = [
        float(abs(np.vdot(full_even, density * full_even) - np.vdot(full_odd, density * full_odd)))
        for density in (local_density(states, site) for site in range(3 * L - 1))
    ]
    code_energies = (float(energy_even[0]), float(energy_odd[0]))
    split = abs(code_energies[1] - code_energies[0])
    isolation = min(float(energy_even[1]), float(energy_odd[1])) - max(code_energies)
    parity_commutator_error = 0.0
    for row, state_row in enumerate(states):
        for column, state_column in enumerate(states):
            if abs(H[row, column]) > 1e-12 and branch_parities(int(state_row), L) != branch_parities(int(state_column), L):
                parity_commutator_error = max(parity_commutator_error, float(abs(H[row, column])))
    return {
        "L": L, "total_charge": total_charge, "parameters": parameters,
        "logical_sectors": {"even_even": [1, 1], "odd_odd": [-1, -1]},
        "logical_energies": list(code_energies), "logical_split": split,
        "isolation_gap": isolation,
        "split_over_gap": float(split / isolation) if isolation > 0 else None,
        "max_site_density_local_distinguishability": max(probes),
        "site_density_local_distinguishability": probes,
        "maximum_cross_parity_matrix_element": parity_commutator_error,
        "sector_dimensions": {"even_even": int(len(sector_even)), "odd_odd": int(len(sector_odd))},
    }


def score(item: dict) -> tuple[float, float, float]:
    return (
        item["split_over_gap"] if item["split_over_gap"] is not None else np.inf,
        item["max_site_density_local_distinguishability"],
        -item["isolation_gap"],
    )


def main() -> None:
    family = [
        {"t_leg": t_leg, "g": g, "Delta": Delta, "V_rung": V_rung, "V_leg": V_leg}
        for t_leg in (0.25, 0.5, 1.0)
        for g in (0.3, 0.6)
        for Delta in (2.0, 5.0)
        for V_rung in (0.0, 3.0)
        for V_leg in (-1.0, 0.0, 1.0)
    ]
    rows = [analyse(L, total_charge, parameters) for L, total_charge in ((4, 4), (5, 4)) for parameters in family]
    candidates = [
        item for item in rows
        if item["split_over_gap"] is not None and item["split_over_gap"] < 1e-3
        and item["max_site_density_local_distinguishability"] < 0.1
        and item["isolation_gap"] > 1e-2
        and item["maximum_cross_parity_matrix_element"] < 1e-12
    ]
    best = sorted(rows, key=score)[:10]
    out = {
        "schema": "antler.phase6.charge2-mediator-ladder-preflight.v1",
        "model": {
            "name": "two-rail ladder with explicit charge-two bond mediators",
            "conserved_quantity": "sum(n_a+n_b+2 n_d)",
            "exact_symmetries": ["(-1)^N_a", "(-1)^N_b"],
            "not_iemini": "the pair transfer is produced through explicit bond mediators rather than inserted as a two-wire parent-Hamiltonian term",
        },
        "claim_boundary": (
            "This is an initial finite-size protection preflight of a new candidate extension. A positive row would require larger-L scaling, "
            "a complete local-operator basis and a microscopic physical implementation before braid work. A negative result rejects only this registered scan."
        ),
        "criteria": {
            "split_over_gap_below": 1e-3,
            "max_site_density_local_distinguishability_below": 0.1,
            "isolation_gap_above": 1e-2,
            "maximum_cross_parity_matrix_element_below": 1e-12,
        },
        "rows": rows, "best_rows_by_preflight_score": best, "candidates": candidates,
        "decision": (
            "candidate(s) pass only the initial finite-size gate; require scaling and a complete local-operator audit before any braid calculation"
            if candidates else
            "no protected-doublet candidate in the registered charge-two mediator ladder scan; do not launch braid or gate calculations"
        ),
    }
    path = ROOT / "results" / "phase6" / "charge2_mediator_ladder_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "candidate_count": len(candidates), "best_rows": best}, indent=2))


if __name__ == "__main__":
    main()
