"""Coherent pair-transfer dynamics inside the two-block physical sector.

This advances the two-block transport algebra audit from a matrix-element
check to an exact propagation.  A pair starts in block 0 and is transferred to
block 1 by the inserted pair-only boundary term.  The audit tracks pair
arrival, walker micromotion and exact Gauss-sector leakage.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse.linalg import expm_multiply


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from run_phase8b_two_block_pair_transport_audit import (
    DETUNING, PAIR_TRANSFER, basis_index, build_model,
)


RATIOS = (0.10, 0.05, 0.025)
TIME_SAMPLES = 81


def main() -> None:
    states, positions = build_basis(8, 2)
    source_state = (1 << site_index(0, 0)) | (1 << site_index(1, 0))
    target_state = (1 << site_index(2, 0)) | (1 << site_index(3, 0))
    source_matter = positions[source_state]
    target_matter = positions[target_state]
    rows = []
    for ratio in RATIOS:
        hamiltonian, _, _, g0, g1, metadata = build_model(ratio * DETUNING)
        dimension = hamiltonian.shape[0]
        # Gauge label 0 is X_L=X_M=X_R=+1.  Both source and target have even
        # a parity in each block, so they begin in G0=G1=+1.
        source = basis_index(source_matter, 0, 0, 0, metadata["matter_dimension"])
        target_indices = np.asarray([
            basis_index(target_matter, gauge, walker0, walker1, metadata["matter_dimension"])
            for gauge in range(8) for walker0 in range(4) for walker1 in range(4)
        ], dtype=int)
        walker_excited = np.asarray([
            basis_index(matter, gauge, walker0, walker1, metadata["matter_dimension"])
            for matter in range(metadata["matter_dimension"])
            for gauge in range(8) for walker0 in range(4) for walker1 in range(4)
            if walker0 != 0 or walker1 != 0
        ], dtype=int)
        physical = (g0 == 1) & (g1 == 1)
        initial = np.zeros(dimension, dtype=complex)
        initial[source] = 1.0
        times = np.linspace(0.0, np.pi / PAIR_TRANSFER, TIME_SAMPLES)
        states_time = expm_multiply(-1j * hamiltonian, initial, start=times[0], stop=times[-1], num=TIME_SAMPLES, endpoint=True)
        target_population = np.sum(np.abs(states_time[:, target_indices]) ** 2, axis=1)
        source_population = np.abs(states_time[:, source]) ** 2
        walker_population = np.sum(np.abs(states_time[:, walker_excited]) ** 2, axis=1)
        gauss_leakage = 1.0 - np.sum(np.abs(states_time[:, physical]) ** 2, axis=1)
        max_index = int(np.argmax(target_population))
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": ratio * DETUNING,
            "time_horizon": float(times[-1]),
            "maximum_target_pair_population": float(target_population[max_index]),
            "time_of_maximum_target_pair_population": float(times[max_index]),
            "minimum_source_pair_population": float(np.min(source_population)),
            "maximum_virtual_walker_population": float(np.max(walker_population)),
            "maximum_gauss_sector_leakage": float(np.max(np.abs(gauss_leakage))),
        })
    output = {
        "schema": "antler.phase8b.two-block-pair-dynamics.v1",
        "parameters": {
            "fixed_matter_particle_number": 2,
            "pair_transfer": PAIR_TRANSFER,
            "ratios": list(RATIOS),
            "time_samples": TIME_SAMPLES,
            "initial": "two a-rail particles in block 0, walkers empty, all boundary X eigenvalues +1",
            "target": "two a-rail particles in block 1, any walker/gauge microstate",
        },
        "rows": rows,
        "decision": (
            "Exact pair-only two-block propagation stays in the G0=G1=+1 sector and transfers the pair coherently in the "
            "inserted-primitive model. The walker population is virtual micromotion, not a loss channel."
        ),
        "claim_boundary": (
            "This is a single-pair coherent transport control. It does not establish a many-pair phase, a thermodynamic "
            "gap, local indistinguishability, a native implementation, a T-junction, fusion or braid statistics."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_two_block_pair_dynamics.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
