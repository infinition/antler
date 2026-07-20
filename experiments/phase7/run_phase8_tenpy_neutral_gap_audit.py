"""Neutral-gap audit for the fixed-density Phase-8 Floquet candidate.

For each *explicitly conserved* branch-parity sector we first obtain its
ground state and then repeat finite-MPS DMRG with that state in
``orthogonal_to``.  Consequently this measures the first neutral excitation
in the same fixed-N, fixed-Z2 sector.  L=8 is registered against exact
diagonalization before the longer MPS values are retained.

It does not turn a finite-size candidate into a thermodynamic phase claim.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from phase8_nc_majorana_tenpy_model import DynamicNumberConservingLadder
from run_phase8_tenpy_dmrg_parity_scaling import product_state


U0, ALPHA, CHI_MAX = -2.0, 0.5, 384
CASES = ((8, 4), (12, 6), (16, 8))
DMRG_PARAMS = {
    "mixer": True,
    "max_E_err": 1e-11,
    "max_sweeps": 60,
    "trunc_params": {"chi_max": CHI_MAX, "svd_min": 1e-12},
}


def excited_seed(length: int, particle_number: int, parity: int) -> list[str]:
    """A same-charge, same-parity seed spatially distinct from the GS seed."""
    empty, a, b = "empty_A empty_B", "full_A empty_B", "empty_A full_B"
    state = [empty] * length
    # Spread the particles throughout the chain instead of filling its left end.
    positions = [round(index * (length - 1) / (particle_number - 1)) for index in range(particle_number)]
    for index in positions:
        state[index] = a
    if parity == 1:
        state[positions[0]] = b
    return state


def model_for(length: int) -> DynamicNumberConservingLadder:
    return DynamicNumberConservingLadder({
        "L": length,
        "bc_MPS": "finite",
        "u0": U0,
        "alpha": ALPHA,
        "t_leg": 1.0,
        "conserve_branch_parity": True,
    })


def sector_gap(length: int, particle_number: int, parity: int) -> dict:
    model = model_for(length)
    sites = model.lat.mps_sites()
    gs = MPS.from_product_state(
        sites, product_state(length, particle_number, parity), bc="finite", unit_cell_width=length
    )
    ground_info = dmrg.run(gs, model, DMRG_PARAMS)
    excited = MPS.from_product_state(
        sites, excited_seed(length, particle_number, parity), bc="finite", unit_cell_width=length
    )
    excited_info = dmrg.run(excited, model, DMRG_PARAMS, orthogonal_to=[gs])
    overlap = gs.overlap(excited)
    branch_parity = excited.expectation_value_term(
        [("PA", site) for site in range(length)], autoJW=False
    )
    return {
        "ground_energy": float(ground_info["E"]),
        "first_neutral_excitation": float(excited_info["E"]),
        "neutral_gap": float(excited_info["E"] - ground_info["E"]),
        "ground_excited_overlap_abs": float(abs(overlap)),
        "excited_branch_parity": float(np.real_if_close(branch_parity)),
        "ground_max_bond_dimension": int(max(gs.chi)),
        "excited_max_bond_dimension": int(max(excited.chi)),
    }


def l8_exact_gaps() -> dict[str, float]:
    payload = json.loads((ROOT / "results" / "phase7" / "nc_majorana_l8_sparse_audit.json").read_text(encoding="utf-8"))
    row = next(item for item in payload["constant_density_targets"] if (item["L"], item["N"]) == (8, 4))
    return {parity: float(data["neutral_gap"]) for parity, data in row["parity_sectors"].items()}


def main() -> None:
    exact = l8_exact_gaps()
    rows = []
    for length, particle_number in CASES:
        sectors = {str(parity): sector_gap(length, particle_number, parity) for parity in (0, 1)}
        row = {
            "L": length,
            "N": particle_number,
            "filling_N_over_2L": particle_number / (2.0 * length),
            "sectors": sectors,
            "smallest_neutral_gap": min(data["neutral_gap"] for data in sectors.values()),
        }
        if (length, particle_number) == (8, 4):
            residuals = {parity: abs(sectors[parity]["neutral_gap"] - exact[parity]) for parity in ("0", "1")}
            row["exact_ed_gap_residual_by_sector"] = residuals
            row["exact_ed_gap_threshold"] = 1e-7
            if max(residuals.values()) > row["exact_ed_gap_threshold"]:
                raise RuntimeError(f"L=8 neutral-gap validation failed: {row}")
        rows.append(row)
    out = {
        "schema": "antler.phase8.tenpy-neutral-gap-audit.v1",
        "parameters": {"u0_attractive_nn": U0, "alpha": ALPHA, "chi_max": CHI_MAX},
        "rows": rows,
        "decision": "Neutral excitations are variationally resolved in fixed U(1) x Z2 sectors; finite-size scaling remains required.",
        "claim_boundary": "No thermodynamic phase, ANTLER-native realization, braid, or non-Abelian conclusion follows.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_neutral_gap_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
