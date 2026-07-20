"""Exact periodic-boundary charged-gap audit for the external Floquet ladder.

The number-conserving Majorana ladder has a gapless neutral (total-charge)
mode, so its relevant topological mass is the charged-sector curvature

    Delta_topo = [E(N + 1) + E(N - 1) - 2 E(N)] / 2

on a periodic ladder.  This finite-ED implementation is deliberately limited
to small systems, but avoids open-boundary Majorana edge states contaminating
the charged diagnostic.  It is a reproduction/calibration of the published
effective model, not a native-ANTLER result.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply, wire_a_parity
from run_phase8_nc_majorana_l8_sparse_audit import ETA, global_rail_rotation_sparse, sparse_from_entries


T_HOP, ALPHA = 1.0, 0.5
# The first two are the published density nu=1/3 and interaction scale U0=-1.5.
# The final two record the earlier U0=-2, nu=1/4 finite-size candidate separately.
TARGETS = (
    {"L": 6, "N": 4, "u0": -1.5, "label": "published_density_L6"},
    {"L": 9, "N": 6, "u0": -1.5, "label": "published_density_L9"},
    {"L": 4, "N": 2, "u0": -2.0, "label": "registered_L4"},
    {"L": 8, "N": 4, "u0": -2.0, "label": "registered_L8"},
)


def h0_periodic_sparse(length: int, particle_number: int, interaction: float) -> tuple[csr_matrix, np.ndarray, dict[int, int]]:
    states, index = build_basis(2 * length, particle_number)
    entries: dict[tuple[int, int], complex] = {}
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        diagonal = 0.0
        for rung in range(length):
            next_rung = (rung + 1) % length
            for rail in (0, 1):
                left, right = site_index(rung, rail), site_index(next_rung, rail)
                for operations in ((("ann", right), ("create", left)), (("ann", left), ("create", right))):
                    item = _apply(state, operations)
                    if item is not None:
                        new_state, amplitude = item
                        key = (index[new_state], column)
                        entries[key] = entries.get(key, 0.0) - T_HOP * amplitude
                diagonal += interaction * (((state >> left) & 1) * ((state >> right) & 1))
        entries[(column, column)] = entries.get((column, column), 0.0) + diagonal
    return sparse_from_entries(entries, len(states)), states, index


def effective_periodic_sparse(length: int, particle_number: int, interaction: float) -> tuple[csr_matrix, np.ndarray]:
    h0, states, index = h0_periodic_sparse(length, particle_number, interaction)
    rotation = global_rail_rotation_sparse(length, states, index, ETA)
    return (ALPHA * h0 + (1.0 - ALPHA) * (rotation.conj().T @ h0 @ rotation)).tocsr(), states


def ground_energy(length: int, particle_number: int, interaction: float) -> dict:
    hamiltonian, states = effective_periodic_sparse(length, particle_number, interaction)
    value, vector = eigsh(hamiltonian, k=1, which="SA", tol=1e-10, maxiter=200_000)
    state = vector[:, 0]
    probabilities = np.abs(state) ** 2
    parity_weight = {
        str(parity): float(sum(probabilities[index] for index, raw in enumerate(states) if wire_a_parity(int(raw), length) == parity))
        for parity in (0, 1)
    }
    return {
        "N": particle_number,
        "hilbert_dimension": int(len(states)),
        "ground_energy": float(value[0]),
        "branch_parity_weight": parity_weight,
    }


def analyze(target: dict) -> dict:
    length, number, interaction = target["L"], target["N"], target["u0"]
    energies = {str(n): ground_energy(length, n, interaction) for n in (number - 1, number, number + 1)}
    lower, center, upper = (energies[str(n)]["ground_energy"] for n in (number - 1, number, number + 1))
    return {
        **target,
        "filling_N_over_2L": number / (2.0 * length),
        "boundary_condition": "periodic",
        "energies": energies,
        "charged_gap_curvature": float((upper + lower - 2.0 * center) / 2.0),
        "diagnostic_note": "Finite periodic-system estimate. It is not the infinite-MPS charged excitation gap used for a phase certificate.",
    }


def main() -> None:
    rows = [analyze(target) for target in TARGETS]
    out = {
        "schema": "antler.phase8.pbc-charged-gap-ed.v1",
        "citation": "Defossez et al., arXiv:2412.14886v2 (2025)",
        "parameters": {"alpha": ALPHA, "eta": ETA, "t_leg": T_HOP},
        "rows": rows,
        "decision": "Small-system periodic charged-gap calibration. The neutral gap is intentionally not used as a rejection condition.",
        "claim_boundary": "No thermodynamic phase, native ANTLER realization, protected hardware qubit, braid, or non-Abelian conclusion follows.",
    }
    path = ROOT / "results" / "phase7" / "pbc_charged_gap_ed.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
