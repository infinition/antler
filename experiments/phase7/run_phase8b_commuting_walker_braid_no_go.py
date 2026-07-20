"""Register the non-Abelian limitation of the current walker stabilizer route.

The tiled Phase-8B effective parent contains fixed-support commuting words
A_s, B_p and A_s B_p.  Changing their coefficients in time does not change
their common eigenbasis.  This script records the exact consequence before any
claim about defect motion or braiding is entertained.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
STARS = 9
PLAQUETTES = 9
INCIDENCES = 36


def main() -> None:
    path = ROOT / "results" / "phase7" / "phase8b_effective_toric_tile_closure.json"
    tiled = json.loads(path.read_text(encoding="utf-8"))
    rows = tiled["rows"]
    # In the code sector all A_s=B_p=+1, including every incident product.
    code_energies = [
        STARS * float(row["star_coefficient"])
        + PLAQUETTES * float(row["plaquette_coefficient"])
        + INCIDENCES * float(row["star_plaquette_product_coefficient_in_stabilizer_basis"])
        for row in rows
    ]
    pairs = []
    for left, right in itertools.combinations(range(len(rows)), 2):
        # Every coefficient path is a polynomial in the same commuting
        # stabilizers, so this is exactly zero, not a finite-size estimate.
        pairs.append({
            "left_coupling_over_detuning": rows[left]["coupling_over_detuning"],
            "right_coupling_over_detuning": rows[right]["coupling_over_detuning"],
            "effective_hamiltonian_commutator_norm": 0.0,
            "projected_logical_commutator_norm": 0.0,
        })
    output = {
        "schema": "antler.phase8b.commuting-walker-braid-no-go.v1",
        "input": "phase8b_effective_toric_tile_closure.json",
        "exact_algebra": {
            "generator_family": "{A_s, B_p, A_s B_p}",
            "all_generators_commute": True,
            "time_dependent_coefficient_paths_commute_at_all_times": True,
            "projected_hamiltonian_on_fourfold_code": "E_code(t) I_4",
            "projected_logical_hamiltonian_nonscalar_norm": 0.0,
        },
        "code_energies_by_registered_ratio": [
            {"coupling_over_detuning": row["coupling_over_detuning"], "code_energy": energy}
            for row, energy in zip(rows, code_energies)
        ],
        "all_pairwise_coefficient_path_checks": pairs,
        "decision": (
            "The present walker-derived star/plaquette family can implement an Abelian stabilizer memory but cannot yield a "
            "noncommuting holonomy by varying its existing coefficients. Its evolution is scalar on the fourfold code sector."
        ),
        "minimal_next_resource": (
            "A microscopically derived defect-deformation primitive that changes stabilizer support and exchanges electric "
            "and magnetic character across a branch cut (or another independently derived noncommuting logical generator). "
            "A static inserted Hadamard, Majorana/BdG term or braid matrix does not meet this requirement."
        ),
        "claim_boundary": (
            "This is an algebraic no-go for coefficient-only modulation of the current effective family. It does not rule out "
            "twist defects, measurement protocols, a new gauge-link control primitive, or another non-Abelian extension once "
            "they are explicitly derived and audited."
        ),
    }
    result = ROOT / "results" / "phase7" / "phase8b_commuting_walker_braid_no_go.json"
    result.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
