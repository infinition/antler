"""Canonical U(1) DMRG reproduction of the exact L=8 Floquet candidate energy."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from phase8_nc_majorana_tenpy_model import DynamicNumberConservingLadder


L, N, U0, ALPHA = 8, 4, -2.0, 0.5
DMRG_PARAMS = {
    "mixer": True,
    "max_E_err": 1e-11,
    "max_sweeps": 30,
    "trunc_params": {"chi_max": 128, "svd_min": 1e-12},
    "verbose": 0,
}


def exact_reference_energy() -> float:
    payload = json.loads((ROOT / "results" / "phase7" / "nc_majorana_l8_sparse_audit.json").read_text(encoding="utf-8"))
    target = next(row for row in payload["constant_density_targets"] if (row["L"], row["N"]) == (L, N))
    return min(sector["ground_energy"] for sector in target["parity_sectors"].values())


def main() -> None:
    model = DynamicNumberConservingLadder({"L": L, "bc_MPS": "finite", "u0": U0, "alpha": ALPHA, "t_leg": 1.0})
    # Four A-rail particles: total N=4 and even branch parity. DMRG's U(1)
    # blocks preserve total charge exactly; the exact ground state is in this
    # parity block for the registered candidate.
    product_state = ["full_A empty_B"] * N + ["empty_A empty_B"] * (L - N)
    psi = MPS.from_product_state(model.lat.mps_sites(), product_state, bc="finite")
    info = dmrg.run(psi, model, DMRG_PARAMS)
    energy = float(info["E"])
    exact_energy = exact_reference_energy()
    residual = abs(energy - exact_energy)
    out = {
        "schema": "antler.phase8.tenpy-dmrg-l8-validation.v1",
        "parameters": {"L": L, "N": N, "u0_attractive_nn": U0, "alpha": ALPHA, "chi_max": 128},
        "dmrg_energy": energy,
        "exact_fixed_N_energy": exact_energy,
        "absolute_energy_residual": residual,
        "threshold": 1e-8,
        "sweeps": int(info.get("sweeps", -1)),
        "max_truncation_error": float(info.get("max_trunc_err", 0.0)),
        "decision": "Canonical DMRG is eligible for L>8 scaling only if it reproduces exact L=8 energy below threshold.",
    }
    if residual > out["threshold"]:
        raise RuntimeError(json.dumps(out, indent=2))
    path = ROOT / "results" / "phase7" / "tenpy_dmrg_l8_validation.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
