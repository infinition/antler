"""Canonical MPS parity-sector scaling at the published Floquet point.

This reproduces the *parameter regime* used by Defossez et al. rather than
promoting the earlier U0=-2 exploratory point (which is close to the reported
phase-separation boundary).  Total charge and branch parity are both exact
tensor charges.  The L=6 result is benchmarked against exact diagonalization.
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
from run_phase8_nc_majorana_l8_sparse_audit import analyze as exact_analyze
from run_phase8_tenpy_dmrg_parity_scaling import product_state


U0, ALPHA, CHI_MAX = -1.5, 0.5, 384
CASES = ((6, 4), (12, 8), (18, 12))  # nu=N/(2L)=1/3.
DMRG_PARAMS = {
    "mixer": True,
    "max_E_err": 1e-11,
    "max_sweeps": 60,
    "trunc_params": {"chi_max": CHI_MAX, "svd_min": 1e-12},
}


def sector(length: int, particle_number: int, parity: int) -> dict:
    model = DynamicNumberConservingLadder({
        "L": length,
        "bc_MPS": "finite",
        "u0": U0,
        "alpha": ALPHA,
        "t_leg": 1.0,
        "conserve_branch_parity": True,
    })
    psi = MPS.from_product_state(
        model.lat.mps_sites(), product_state(length, particle_number, parity),
        bc="finite", unit_cell_width=length,
    )
    info = dmrg.run(psi, model, DMRG_PARAMS)
    pa = psi.expectation_value_term([("PA", site) for site in range(length)], autoJW=False)
    return {
        "energy": float(info["E"]),
        "measured_branch_parity": float(np.real_if_close(pa)),
        "max_bond_dimension": int(max(psi.chi)),
    }


def main() -> None:
    exact = exact_analyze(6, 4, u0=U0, alpha=ALPHA)
    exact_energy = {parity: float(data["ground_energy"]) for parity, data in exact["parity_sectors"].items()}
    rows = []
    for length, particle_number in CASES:
        sectors = {str(parity): sector(length, particle_number, parity) for parity in (0, 1)}
        row = {
            "L": length,
            "N": particle_number,
            "filling_N_over_2L": particle_number / (2.0 * length),
            "sectors": sectors,
            "parity_sector_split": abs(sectors["0"]["energy"] - sectors["1"]["energy"]),
        }
        if (length, particle_number) == (6, 4):
            residuals = {key: abs(sectors[key]["energy"] - exact_energy[key]) for key in ("0", "1")}
            row["exact_ed_residual_by_sector"] = residuals
            row["exact_ed_threshold"] = 1e-8
            if max(residuals.values()) > row["exact_ed_threshold"]:
                raise RuntimeError(f"L=6 published-point DMRG validation failed: {row}")
        rows.append(row)
    out = {
        "schema": "antler.phase8.tenpy-published-point-parity-scaling.v1",
        "citation": "Defossez et al., arXiv:2412.14886v2 (2025)",
        "parameters": {"u0_attractive_nn": U0, "alpha": ALPHA, "chi_max": CHI_MAX},
        "rows": rows,
        "decision": "Fixed-density open-boundary parity-splitting scaling at the published effective-model point.",
        "claim_boundary": "The neutral sector is expected to be gapless. This is an external-model finite-MPS audit, not a thermodynamic, native-ANTLER, or non-Abelian proof.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_published_point_parity_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
