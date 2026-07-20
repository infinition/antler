"""Algebraic Phase 7B preflight for a 2D charge-frozen surface-code control.

This is not a microscopic ANTLER simulation.  It checks the exact stabilizer
reference that a future two-dimensional ANTLER Hamiltonian would have to
realize, and records the missing low-body mediator derivation explicitly.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_stabilizer_algebra import toric_code_preflight


def main() -> None:
    audit = toric_code_preflight(lx=3, ly=3)
    passed = (
        audit["stabilizer_algebra"]["all_stabilizers_commute"]
        and audit["code"]["encoded_qubits"] == 2
        and audit["code"]["ground_space_degeneracy"] == 4
        and audit["code"]["minimum_logical_weight"] == 3
        and audit["complete_projected_local_pauli_gate"]["all_tested_probes_project_to_scalars"]
    )
    if not passed:
        raise RuntimeError("2D abstract stabilizer control failed its algebraic preflight")
    out = {
        "schema": "antler.phase7.2d-stabilizer-preflight.v1",
        "reference_model": (
            "abstract charge-frozen square-torus stabilizer parent; each edge is a one-particle rung qubit"
        ),
        "physical_embedding_contract": {
            "mott_constraint": "C_e=(n_a,e+n_b,e-1)^2, P_e=1-C_e",
            "rung_paulis_in_monomer_sector": {
                "X_e": "a_e^dagger b_e + b_e^dagger a_e",
                "Z_e": "n_a,e - n_b,e",
            },
            "commuting_parent_terms": {
                "star": "Q_s=(prod_(e in star s) P_e - prod_(e in star s) X_e)/2",
                "plaquette": "Q_p=(prod_(e in boundary p) P_e - prod_(e in boundary p) Z_e)/2",
                "hamiltonian": "H=U sum_e C_e + J_s sum_s Q_s + J_p sum_p Q_p",
            },
            "exact_symmetries": ["total U(1) charge", "global rail parities P_a and P_b"],
        },
        "algebraic_audit": audit,
        "passes_exact_2d_reference_gate": bool(passed),
        "decision": (
            "The exact 2D reference has the expected commuting stabilizer algebra, four torus ground states, "
            "and no non-scalar projected Pauli on fewer than three physical edges.  This is a reference-control "
            "result only; it is not a native ANTLER Hamiltonian."
        ),
        "mandatory_next_gate": (
            "Derive, from explicitly specified ANTLER local degrees of freedom and mediator couplings, a low-body "
            "Hamiltonian whose controlled low-energy theory produces these four-rung star and plaquette terms. "
            "Until that Schrieffer-Wolff/error analysis passes, the parent is imposed rather than emergent."
        ),
        "claim_boundary": (
            "This establishes neither a microscopic mediator realization, a dynamical exchange protocol, anyon braiding, "
            "non-Abelian statistics, universality, nor fault tolerance.  The reference stabilizer code is Abelian."
        ),
    }
    path = ROOT / "results" / "phase7" / "2d_surface_code_preflight.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
