"""Targeted bond-dimension convergence of the L=18 published-point split."""
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


LENGTH, PARTICLE_NUMBER, U0, ALPHA = 18, 12, -1.5, 0.5
CHI_VALUES = (256, 384, 512)


def run_sector(parity: int, chi_max: int) -> dict:
    model = DynamicNumberConservingLadder({
        "L": LENGTH, "bc_MPS": "finite", "u0": U0, "alpha": ALPHA, "t_leg": 1.0,
        "conserve_branch_parity": True,
    })
    psi = MPS.from_product_state(
        model.lat.mps_sites(), product_state(LENGTH, PARTICLE_NUMBER, parity),
        bc="finite", unit_cell_width=LENGTH,
    )
    info = dmrg.run(psi, model, {
        "mixer": True,
        "max_E_err": 1e-11,
        "max_sweeps": 80,
        "trunc_params": {"chi_max": chi_max, "svd_min": 1e-12},
    })
    pa = psi.expectation_value_term([("PA", site) for site in range(LENGTH)], autoJW=False)
    return {
        "energy": float(info["E"]),
        "measured_branch_parity": float(np.real_if_close(pa)),
        "max_bond_dimension": int(max(psi.chi)),
    }


def main() -> None:
    by_chi = {}
    for chi in CHI_VALUES:
        sectors = {str(parity): run_sector(parity, chi) for parity in (0, 1)}
        by_chi[str(chi)] = {
            "sectors": sectors,
            "parity_sector_split": abs(sectors["0"]["energy"] - sectors["1"]["energy"]),
        }
    split_384 = by_chi["384"]["parity_sector_split"]
    split_512 = by_chi["512"]["parity_sector_split"]
    out = {
        "schema": "antler.phase8.tenpy-published-l18-convergence.v1",
        "parameters": {"L": LENGTH, "N": PARTICLE_NUMBER, "u0_attractive_nn": U0, "alpha": ALPHA, "chi_values": list(CHI_VALUES)},
        "by_chi": by_chi,
        "absolute_split_change_chi384_to_chi512": abs(split_384 - split_512),
        "relative_split_change_chi384_to_chi512": abs(split_384 - split_512) / split_512,
        "decision": "Numerical convergence check only; it does not establish the thermodynamic asymptotics of the splitting.",
        "claim_boundary": "No native-ANTLER, braid, non-Abelian, or fault-tolerance claim follows.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_published_l18_convergence.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
