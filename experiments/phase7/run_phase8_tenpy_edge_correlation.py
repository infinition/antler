"""Finite-MPS edge-correlation diagnostic at the published Floquet point.

The single-particle correlator <a_0^dagger a_j> should decay away from the
left edge and revive at the far edge in the number-conserving Majorana ladder.
We keep total U(1) and branch Z2 explicitly; this is an external-model
signature check, not a proof of a native ANTLER edge qubit.
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


LENGTH, PARTICLE_NUMBER, U0, ALPHA, CHI_MAX = 12, 8, -1.5, 0.5, 384


def main() -> None:
    model = DynamicNumberConservingLadder({
        "L": LENGTH,
        "bc_MPS": "finite",
        "u0": U0,
        "alpha": ALPHA,
        "t_leg": 1.0,
        "conserve_branch_parity": True,
    })
    psi = MPS.from_product_state(
        model.lat.mps_sites(), product_state(LENGTH, PARTICLE_NUMBER, 0),
        bc="finite", unit_cell_width=LENGTH,
    )
    info = dmrg.run(psi, model, {
        "mixer": True,
        "max_E_err": 1e-11,
        "max_sweeps": 60,
        "trunc_params": {"chi_max": CHI_MAX, "svd_min": 1e-12},
    })
    sites = list(range(1, LENGTH))
    values = np.asarray(
        psi.correlation_function("CdA", "CA", sites1=[0], sites2=sites, autoJW=True), dtype=complex
    ).reshape(-1)
    absolute = np.abs(values)
    # Exclude the first and last two points; this is only a descriptive
    # revival ratio, never an invariant or a pass/fail topological test.
    interior = absolute[2:-2]
    endpoint = float(absolute[-1])
    out = {
        "schema": "antler.phase8.tenpy-edge-correlation.v1",
        "parameters": {
            "L": LENGTH,
            "N": PARTICLE_NUMBER,
            "filling_N_over_2L": PARTICLE_NUMBER / (2.0 * LENGTH),
            "u0_attractive_nn": U0,
            "alpha": ALPHA,
            "chi_max": CHI_MAX,
            "branch_parity_sector": 0,
        },
        "energy": float(info["E"]),
        "correlator_a0dagger_aj": [
            {"j": int(site), "real": float(value.real), "imag": float(value.imag), "abs": float(abs(value))}
            for site, value in zip(sites, values)
        ],
        "far_edge_abs": endpoint,
        "interior_abs_min": float(np.min(interior)),
        "interior_abs_max": float(np.max(interior)),
        "far_edge_over_interior_min": float(endpoint / np.min(interior)),
        "decision": "A finite-chain edge-correlation diagnostic at the published effective-model point.",
        "claim_boundary": "A correlation revival alone does not establish a thermodynamic topological phase, native ANTLER realization, braid, non-Abelian statistics, or fault tolerance.",
    }
    path = ROOT / "results" / "phase7" / "tenpy_edge_correlation.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
