"""Audit of a *common* microscopic rail-parity symmetry for direct H0/H1.

The direct Phase-8 implementation reuses two physical charge-two mediator
slots per link while rotating their pair channels between H0 and H1.  A claim
of exact stroboscopic Z2 protection requires one fixed parity operator on the
full microscopic Hilbert space, not a different mediator-parity assignment in
each half-period.  We exhaust every mediator Z2 assignment on the registered
L=3,N=2 block.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import site_index
from run_phase8_native_micro_floquet_l3 import LENGTH, PARTICLE_NUMBER, RATIO, TARGET_U0
from run_phase8_native_direct_h1_closure import build_micro, mediator_mode
from run_phase8_direct_channel_ramp_audit import channel_path


def parity_operator(states: np.ndarray, assignment: tuple[int, ...]) -> np.ndarray:
    """(-1)^(N_a + sum_m assignment_m n_d,m) on the weighted Fock basis."""
    diagonal = []
    for raw_state in states:
        state = int(raw_state)
        exponent = sum((state >> site_index(rung, 0)) & 1 for rung in range(LENGTH))
        for link in range(LENGTH - 1):
            for slot in (0, 1):
                position = 2 * link + slot
                exponent += assignment[position] * ((state >> mediator_mode(LENGTH, link, slot)) & 1)
        diagonal.append(-1.0 if exponent % 2 else 1.0)
    return np.diag(diagonal)


def normalized_commutator(hamiltonian: np.ndarray, parity: np.ndarray) -> float:
    denominator = np.linalg.norm(hamiltonian, ord="fro")
    return float(np.linalg.norm(hamiltonian @ parity - parity @ hamiltonian, ord="fro") / denominator)


def main() -> None:
    detuning = abs(TARGET_U0) / RATIO**2
    g = RATIO * detuning
    h0, states, _ = build_micro(LENGTH, PARTICLE_NUMBER, channel_path(0.0), g, detuning)
    h1, states_h1, _ = build_micro(LENGTH, PARTICLE_NUMBER, channel_path(np.pi / 2.0), g, detuning)
    if not np.array_equal(states, states_h1):
        raise RuntimeError("segment bases differ")
    rows = []
    for assignment in product((0, 1), repeat=2 * (LENGTH - 1)):
        parity = parity_operator(states, assignment)
        rows.append({
            "mediator_z2_charges_by_link_slot": list(assignment),
            "normalized_h0_commutator_frobenius": normalized_commutator(h0, parity),
            "normalized_h1_commutator_frobenius": normalized_commutator(h1, parity),
        })
    exact_h0 = [row["mediator_z2_charges_by_link_slot"] for row in rows if row["normalized_h0_commutator_frobenius"] < 1e-12]
    exact_h1 = [row["mediator_z2_charges_by_link_slot"] for row in rows if row["normalized_h1_commutator_frobenius"] < 1e-12]
    common = [row["mediator_z2_charges_by_link_slot"] for row in rows if max(
        row["normalized_h0_commutator_frobenius"], row["normalized_h1_commutator_frobenius"]
    ) < 1e-12]
    best = min(rows, key=lambda row: max(
        row["normalized_h0_commutator_frobenius"], row["normalized_h1_commutator_frobenius"]
    ))
    out = {
        "schema": "antler.phase8.direct-channel-common-parity-audit.v1",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0,
            "g_over_detuning": RATIO, "detuning": detuning, "g": g,
            "mediator_slots": "[link0_slot0, link0_slot1, link1_slot0, link1_slot1]",
            "segments": {
                "H0": "(aa-bb)/sqrt(2), (aa+bb)/sqrt(2)",
                "H1": "(aa-bb)/sqrt(2), (ab+ba)/sqrt(2)",
            },
        },
        "all_assignments": rows,
        "exact_assignments_for_h0": exact_h0,
        "exact_assignments_for_h1": exact_h1,
        "common_exact_assignments": common,
        "best_compromise_assignment": best,
        "decision": (
            "A common exact microscopic rail-parity symmetry exists only if the same fixed mediator Z2 assignment "
            "commutes with both H0 and H1. An empty intersection is a resource-accounting obstruction for the "
            "reused-two-mediator direct-channel implementation, not a no-go for an architecture with separately "
            "parity-charged mediator species or another symmetry mechanism."
        ),
        "claim_boundary": (
            "This audit addresses the explicit direct-channel microscopic block only. It does not decide the existence "
            "of a T-junction fusion space, a thermodynamic phase, a braid, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "direct_channel_common_parity_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "exact_assignments_for_h0": exact_h0,
        "exact_assignments_for_h1": exact_h1,
        "common_exact_assignments": common,
        "best_compromise_assignment": best,
    }, indent=2))


if __name__ == "__main__":
    main()
