"""Validate the U(1)-conserving TeNPy rung MPO against the frozen dense ED."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from phase8_nc_majorana_tenpy_model import DynamicNumberConservingLadder, make_rung_site
from run_phase7d_floquet_full_ladder_preflight import ETA


L, N, U0, ALPHA = 4, 2, -2.0, 0.5


def main() -> None:
    # Validate the exact local identity which defines the H1 interaction:
    # P^dag[nA nA' + nB nB']P = [N N' + Y Y']/2 at eta=pi/2.
    rung = make_rung_site()
    n_a, n_b = rung.get_op("NA").to_ndarray(), rung.get_op("NB").to_ndarray()
    n_total, y_rail = rung.get_op("Ntot").to_ndarray(), rung.get_op("Yrail").to_ndarray()
    cd_a, c_a = rung.get_op("CdA").to_ndarray(), rung.get_op("CA").to_ndarray()
    cd_b, c_b = rung.get_op("CdB").to_ndarray(), rung.get_op("CB").to_ndarray()
    jx = 0.5 * (cd_a @ c_b + cd_b @ c_a)
    rotation = expm(-1j * ETA * jx)
    original = np.kron(n_a, n_a) + np.kron(n_b, n_b)
    transformed = np.kron(rotation.conj().T, rotation.conj().T) @ original @ np.kron(rotation, rotation)
    formula = 0.5 * (np.kron(n_total, n_total) + np.kron(y_rail, y_rail))
    residual = float(np.linalg.norm(transformed - formula))
    # Building the MPO also checks that the local charge blocks accept all terms.
    model = DynamicNumberConservingLadder({"L": L, "bc_MPS": "finite", "u0": U0, "alpha": ALPHA, "t_leg": 1.0})
    out = {
        "schema": "antler.phase8.tenpy-rung-mpo-equivalence.v1",
        "parameters": {"L": L, "N": N, "u0_attractive_nn": U0, "alpha": ALPHA, "eta": ETA},
        "total_u1_conserved_by_site": True,
        "local_floquet_interaction_identity_frobenius": residual,
        "mpo_max_bond_dimension": model.H_MPO.chi,
        "threshold": 1e-11,
        "decision": "The local Floquet algebra and charge-block MPO must pass before canonical DMRG validation against L=8 ED.",
    }
    if residual > out["threshold"]:
        raise RuntimeError(json.dumps(out, indent=2))
    path = ROOT / "results" / "phase7" / "tenpy_rung_mpo_equivalence.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
