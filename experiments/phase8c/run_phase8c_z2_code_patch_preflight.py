"""Phase 8C-T2: pure-gauge Z2 code-patch preflight on a 3x3 torus.

This is deliberately the pure-gauge branch of the separately declared neutral
link resource.  There is no mobile matter in this test: the local Gauss
generator is therefore the star check A_v=product_(e incident to v) tau^x_e.
The calculation is exact binary stabilizer algebra, not an imposed braid or a
microscopic derivation from the frozen ANTLER ladder.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_stabilizer_algebra import (  # noqa: E402
    BinaryPauli,
    _gf2_basis,
    _packed,
    commute,
    gf2_in_span,
    local_paulis,
    symplectic_parity,
    toric_code_geometry,
    toric_code_preflight,
)


LX = 3
LY = 3
STAR_COUPLING = 1.0
PLAQUETTE_COUPLING = 1.0


def pauli_labels(labels: tuple[str, ...], pauli: BinaryPauli) -> list[str]:
    """Return a phase-free human-readable Pauli word."""
    word = []
    for edge, label in enumerate(labels):
        x_bit = (pauli.x >> edge) & 1
        z_bit = (pauli.z >> edge) & 1
        if x_bit or z_bit:
            factor = "Y" if x_bit and z_bit else "X" if x_bit else "Z"
            word.append(f"{factor}[{label}]")
    return word


def find_anticommuting_logical_pair() -> tuple[BinaryPauli, BinaryPauli]:
    """Find two minimum-weight logical strings with nonzero mutual symplectic form."""
    geometry = toric_code_geometry(LX, LY)
    stabilizer_basis = _gf2_basis(
        tuple(_packed(stabilizer, geometry.n_edges) for stabilizer in geometry.stabilizers)
    )
    logicals = []
    for candidate in local_paulis(geometry.n_edges, min(LX, LY)):
        if all(commute(candidate, stabilizer) for stabilizer in geometry.stabilizers):
            if not gf2_in_span(_packed(candidate, geometry.n_edges), stabilizer_basis):
                logicals.append(candidate)
    for first in logicals:
        for second in logicals:
            if symplectic_parity(first, second):
                return first, second
    raise RuntimeError("expected an anticommuting logical-string pair on the torus")


def main() -> None:
    audit = toric_code_preflight(lx=LX, ly=LY)
    geometry = toric_code_geometry(LX, LY)
    logical_z, logical_x = find_anticommuting_logical_pair()

    independent_star_rank = len(_gf2_basis([star.x for star in geometry.stars]))
    independent_plaquette_rank = len(_gf2_basis([plaquette.z for plaquette in geometry.plaquettes]))
    total_stabilizer_rank = audit["stabilizer_algebra"]["independent_stabilizer_rank_over_GF2"]
    physical_gauge_dimension = 1 << (geometry.n_edges - independent_star_rank)
    ground_space_degeneracy = audit["code"]["ground_space_degeneracy"]

    # For H=-Js sum_v A_v-Jp sum_p B_p, the independent star and plaquette
    # syndromes occur in pairs.  A pair costs 4J, so this is the exact first
    # syndrome gap on a torus, not a numerical diagonalisation estimate.
    first_star_pair_cost = 4.0 * STAR_COUPLING
    first_plaquette_pair_cost = 4.0 * PLAQUETTE_COUPLING
    syndrome_gap = min(first_star_pair_cost, first_plaquette_pair_cost)
    ground_energy = -STAR_COUPLING * len(geometry.stars) - PLAQUETTE_COUPLING * len(geometry.plaquettes)

    local_gate = audit["complete_projected_local_pauli_gate"]
    passes = (
        audit["stabilizer_algebra"]["all_stabilizers_commute"]
        and independent_star_rank == len(geometry.stars) - 1
        and independent_plaquette_rank == len(geometry.plaquettes) - 1
        and total_stabilizer_rank == 2 * LX * LY - 2
        and audit["code"]["encoded_qubits"] == 2
        and ground_space_degeneracy == 4
        and audit["code"]["minimum_logical_weight"] == 3
        and local_gate["all_tested_probes_project_to_scalars"]
        and local_gate["nontrivial_logical"] == 0
        and symplectic_parity(logical_z, logical_x) == 1
    )
    if not passes:
        raise RuntimeError("Phase 8C-T2 pure-gauge code patch failed its pre-registered gate")

    output = {
        "schema": "antler.phase8c.z2-code-patch-preflight.v1",
        "parameters": {
            "new_declared_resource": "one neutral Z2 link qubit tau_e per edge; not an ANTLER charge-two mediator",
            "matter_sector": "pure gauge: matter absent/frozen, so G_v=A_v",
            "gauss_generator": "A_v=product_(e incident to v) tau^x_e",
            "magnetic_plaquette": "B_p=product_(e in boundary p) tau^z_e",
            "reference_hamiltonian": "H=-J_s sum_v A_v-J_p sum_p B_p",
            "star_coupling": STAR_COUPLING,
            "plaquette_coupling": PLAQUETTE_COUPLING,
        },
        "geometry": {
            "name": "3x3 square torus, pure Z2 gauge reference",
            "lx": LX,
            "ly": LY,
            "boundary": "periodic in both directions",
            "link_qubits": geometry.n_edges,
            "vertices": LX * LY,
            "plaquettes": LX * LY,
            "full_link_hilbert_dimension": 1 << geometry.n_edges,
            "physical_gauss_sector_dimension": physical_gauge_dimension,
        },
        "algebra": {
            "star_count": len(geometry.stars),
            "plaquette_count": len(geometry.plaquettes),
            "pair_anticommutations": audit["stabilizer_algebra"]["pair_anticommutations"],
            "all_star_plaquette_checks_commute": audit["stabilizer_algebra"]["all_stabilizers_commute"],
            "independent_gauss_rank_over_GF2": independent_star_rank,
            "independent_plaquette_rank_over_GF2": independent_plaquette_rank,
            "total_stabilizer_rank_over_GF2": total_stabilizer_rank,
            "global_star_relation": "product_v A_v=I",
            "global_plaquette_relation": "product_p B_p=I",
        },
        "code": {
            "encoded_qubits": audit["code"]["encoded_qubits"],
            "ground_space_degeneracy": ground_space_degeneracy,
            "minimum_logical_weight": audit["code"]["minimum_logical_weight"],
            "logical_z_representative": pauli_labels(geometry.edge_labels, logical_z),
            "logical_x_representative": pauli_labels(geometry.edge_labels, logical_x),
            "logical_pair_symplectic_parity": symplectic_parity(logical_z, logical_x),
        },
        "local_indistinguishability": {
            "operator_basis": "all nonidentity Pauli words X/Y/Z on link supports of weight <= distance-1",
            "maximum_tested_weight": local_gate["maximum_tested_weight"],
            "tested_nonidentity_paulis": local_gate["tested_nonidentity_paulis"],
            "anticommutes_with_a_check_and_projects_to_zero": local_gate["anticommutes_with_stabilizer"],
            "stabilizer_scalar": local_gate["stabilizer_scalar"],
            "nontrivial_logical_below_distance": local_gate["nontrivial_logical"],
            "all_tested_local_probes_project_to_scalars": local_gate["all_tested_probes_project_to_scalars"],
        },
        "spectrum": {
            "ground_energy": ground_energy,
            "ground_space_degeneracy": ground_space_degeneracy,
            "first_star_pair_cost": first_star_pair_cost,
            "first_plaquette_pair_cost": first_plaquette_pair_cost,
            "exact_first_syndrome_gap": syndrome_gap,
            "gap_interpretation": "commuting-projector syndrome gap of the imposed pure-gauge reference, not a native ANTLER gap",
        },
        "decision": (
            "PASS T2 for the pure-gauge reference only. The 3x3 torus has an exact commuting Z2 stabilizer algebra, "
            "a fourfold ground space, distance three, a nonzero syndrome gap, and no non-scalar projected Pauli probe "
            "on any support below distance. It is therefore a qualified Abelian code benchmark on which an e<->m "
            "domain-wall test can be posed."
        ),
        "next_gate": (
            "T3 must introduce and audit an actual support-changing e<->m domain wall: it must map electric and magnetic "
            "string operators into each other across the wall. A pi hopping phase, a background B_p=-1 flux, or a static "
            "sign line fails that definition."
        ),
        "claim_boundary": (
            "The neutral link qubits, star checks and plaquette checks are imposed reference resources. This result neither "
            "derives them from frozen ANTLER charge-two mediators nor introduces matter defects, a twist, non-Abelian fusion, "
            "a braid, universality, noise resilience or fault tolerance. The certified parent is the Abelian D(Z2) toric code."
        ),
    }
    path = ROOT / "results" / "phase8c" / "z2_code_patch_preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
