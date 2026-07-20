"""Bond-dimension convergence check for the registered Phase-8 candidate.

This script intentionally repeats the two explicitly conserved branch-parity
sectors.  It does *not* infer a topological phase from a small splitting; it
only establishes whether the reported energy difference is stable against the
MPS bond dimension.
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


U0, ALPHA = -2.0, 0.5
CASES = ((12, 6), (16, 8))
CHI_VALUES = (128, 256, 384)


def sector_energy(length: int, particle_number: int, parity: int, chi_max: int) -> dict:
    model = DynamicNumberConservingLadder({
        "L": length,
        "bc_MPS": "finite",
        "u0": U0,
        "alpha": ALPHA,
        "t_leg": 1.0,
        "conserve_branch_parity": True,
    })
    psi = MPS.from_product_state(
        model.lat.mps_sites(),
        product_state(length, particle_number, parity),
        bc="finite",
        unit_cell_width=length,
    )
    info = dmrg.run(psi, model, {
        "mixer": True,
        "max_E_err": 1e-11,
        "max_sweeps": 50,
        "trunc_params": {"chi_max": chi_max, "svd_min": 1e-12},
    })
    parity_value = psi.expectation_value_term(
        [("PA", site) for site in range(length)], autoJW=False
    )
    return {
        "energy": float(info["E"]),
        "measured_branch_parity": float(np.real_if_close(parity_value)),
        "max_bond_dimension": int(max(psi.chi)),
    }


def main() -> None:
    rows = []
    for length, particle_number in CASES:
        by_chi = {}
        for chi in CHI_VALUES:
            sectors = {str(parity): sector_energy(length, particle_number, parity, chi) for parity in (0, 1)}
            by_chi[str(chi)] = {
                "sectors": sectors,
                "parity_sector_split": abs(sectors["0"]["energy"] - sectors["1"]["energy"]),
            }
        reference = by_chi[str(max(CHI_VALUES))]["parity_sector_split"]
        rows.append({
            "L": length,
            "N": particle_number,
            "filling_N_over_2L": particle_number / (2.0 * length),
            "by_chi": by_chi,
            "split_change_from_chi256_to_chi384": abs(
                by_chi["256"]["parity_sector_split"] - reference
            ),
        })
    out = {
        "schema": "antler.phase8.tenpy-dmrg-bond-convergence.v1",
        "parameters": {"u0_attractive_nn": U0, "alpha": ALPHA, "chi_values": list(CHI_VALUES)},
        "rows": rows,
        "decision": "This is a numerical convergence audit of sector energies, not a thermodynamic phase certificate.",
        "claim_boundary": "No bulk neutral-gap, correlation-length, native-ANTLER, or non-Abelian claim is made.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_dmrg_bond_convergence.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
