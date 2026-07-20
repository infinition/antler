"""U(1)-canonical DMRG scaling of the two branch-parity sectors.

The model itself conserves branch parity. We seed the two disconnected sectors
with product states of opposite PA=prod_j(-1)^n_A,j and verify that value from
the converged MPS. L=8 is checked against exact ED before L=12 is reported.
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


U0, ALPHA = -2.0, 0.5
# Keep the filling fixed and use even N: this is the only comparison for which
# the two branch-parity sectors are independent symmetry sectors rather than
# being related by the rail-exchange symmetry of an odd-N calculation.
CASES = ((8, 4), (12, 6), (16, 8))  # nu=N/(2L)=1/4.
DMRG_PARAMS = {
    "mixer": True,
    "max_E_err": 1e-11,
    "max_sweeps": 40,
    "trunc_params": {"chi_max": 256, "svd_min": 1e-12},
}


def product_state(length: int, particle_number: int, branch_parity: int) -> list[str]:
    if particle_number % 2:
        raise ValueError("This registered scaling uses even N only")
    empty, a, b = "empty_A empty_B", "full_A empty_B", "empty_A full_B"
    if branch_parity == 0:
        occupied = [a] * particle_number
    elif branch_parity == 1:
        occupied = [b] + [a] * (particle_number - 1)
    else:
        raise ValueError("branch_parity must be 0 or 1")
    return occupied + [empty] * (length - particle_number)


def exact_l8_sectors() -> dict[str, float]:
    payload = json.loads((ROOT / "results" / "phase7" / "nc_majorana_l8_sparse_audit.json").read_text(encoding="utf-8"))
    target = next(row for row in payload["constant_density_targets"] if (row["L"], row["N"]) == (8, 4))
    return {parity: float(data["ground_energy"]) for parity, data in target["parity_sectors"].items()}


def run_sector(length: int, particle_number: int, parity: int) -> dict:
    model = DynamicNumberConservingLadder({
        "L": length, "bc_MPS": "finite", "u0": U0, "alpha": ALPHA, "t_leg": 1.0,
        "conserve_branch_parity": True,
    })
    psi = MPS.from_product_state(
        model.lat.mps_sites(),
        product_state(length, particle_number, parity),
        bc="finite",
        unit_cell_width=length,
    )
    info = dmrg.run(psi, model, DMRG_PARAMS)
    pa = psi.expectation_value_term([("PA", site) for site in range(length)], autoJW=False)
    entropies = psi.entanglement_entropy()
    spectra = psi.entanglement_spectrum(by_charge=False)
    center = length // 2 - 1
    center_spectrum = np.asarray(spectra[center], dtype=float)
    return {
        "seeded_branch_parity": parity,
        "energy": float(info["E"]),
        "measured_branch_parity": float(np.real_if_close(pa)),
        "center_entanglement_entropy": float(entropies[center]),
        "center_entanglement_spectrum_first_12": [float(value) for value in center_spectrum[:12]],
    }


def main() -> None:
    rows = []
    exact = exact_l8_sectors()
    for length, particle_number in CASES:
        sectors = {str(parity): run_sector(length, particle_number, parity) for parity in (0, 1)}
        split = abs(sectors["0"]["energy"] - sectors["1"]["energy"])
        row = {
            "L": length,
            "N": particle_number,
            "filling_N_over_2L": particle_number / (2.0 * length),
            "sectors": sectors,
            "parity_sector_split": split,
        }
        if (length, particle_number) == (8, 4):
            residuals = {parity: abs(sectors[parity]["energy"] - exact[parity]) for parity in ("0", "1")}
            row["exact_ed_residual_by_sector"] = residuals
            row["exact_ed_threshold"] = 1e-8
            if max(residuals.values()) > row["exact_ed_threshold"]:
                raise RuntimeError(f"L=8 parity-sector DMRG validation failed: {row}")
        rows.append(row)
    out = {
        "schema": "antler.phase8.tenpy-dmrg-parity-scaling.v1",
        "parameters": {"u0_attractive_nn": U0, "alpha": ALPHA, "chi_max": 256},
        "rows": rows,
        "decision": (
            "U(1)-canonical parity-sector DMRG scaling. A decreasing split is finite-size evidence only; neutral-gap and "
            "correlation-length scaling remain required before a phase claim."
        ),
        "claim_boundary": "This remains an external effective model; no ANTLER-native or non-Abelian claim follows.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_dmrg_parity_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
