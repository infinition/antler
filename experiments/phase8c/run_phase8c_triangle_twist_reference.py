"""Phase 8C-T3b: minimal triangle twist-to-boundary reference.

This is a seven-neutral-qubit, distance-three non-CSS stabilizer patch.  It is
specified geometrically as the smallest triangular patch with three plaquettes
and three boundary loops.  Its central mixed plaquette supplies a finite branch
cut terminating at the boundary.  The audit verifies the stabilizer algebra,
local protection and a string-deformation witness in which a pure-Z boundary
representative is continued as an X/Y word across the mixed check.

It is deliberately a reference patch: it has a twist-to-boundary, not two
separated twist endpoints, a fusion space, or a braid.  The neutral qubits and
checks are inserted resources and are not derived from frozen ANTLER.
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
    gf2_rank,
    local_paulis,
    symplectic_parity,
)


QUBIT_COORDINATES = (
    (0, 0, 0),  # central twist corner
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
)

# The two CSS plaquettes live on the xz and yz faces. The xy face contains a
# Y at the central corner and is the non-CSS twist-to-boundary check. The last
# three words are boundary loops. All supports follow the r=s=t=2 triangular
# incidence pattern; no term is inferred from the ANTLER ladder.
CHECK_LABELS = (
    "XXIXIXI",  # xz plaquette, X type
    "ZIZZIIZ",  # yz plaquette, Z type
    "YZXIXII",  # xy mixed plaquette / twist-to-boundary junction
    "IZIIIZI",  # boundary loop
    "IIZIZII",  # boundary loop
    "IIIXIIX",  # boundary loop
)
TWIST_CHECK_INDEX = 2
Z_SIDE_LABEL = "ZZZIIII"
X_SIDE_LABEL = "IIXXXII"


def parse_pauli(label: str) -> BinaryPauli:
    x = sum(1 << qubit for qubit, factor in enumerate(label) if factor in "XY")
    z = sum(1 << qubit for qubit, factor in enumerate(label) if factor in "YZ")
    return BinaryPauli(x=x, z=z)


def label_of(pauli: BinaryPauli) -> str:
    out = []
    for qubit in range(len(QUBIT_COORDINATES)):
        x_bit = (pauli.x >> qubit) & 1
        z_bit = (pauli.z >> qubit) & 1
        out.append("Y" if x_bit and z_bit else "X" if x_bit else "Z" if z_bit else "I")
    return "".join(out)


def multiply_mod_phase(left: BinaryPauli, right: BinaryPauli) -> BinaryPauli:
    return BinaryPauli(x=left.x ^ right.x, z=left.z ^ right.z)


def is_logical(pauli: BinaryPauli, checks: tuple[BinaryPauli, ...], basis: dict[int, int]) -> bool:
    return all(commute(pauli, check) for check in checks) and not gf2_in_span(
        _packed(pauli, len(QUBIT_COORDINATES)), basis
    )


def minimum_distance(checks: tuple[BinaryPauli, ...], basis: dict[int, int]) -> int:
    for weight in range(1, len(QUBIT_COORDINATES) + 1):
        if any(is_logical(candidate, checks, basis) for candidate in local_paulis(len(QUBIT_COORDINATES), weight)):
            return weight
    raise RuntimeError("no logical Pauli found")


def main() -> None:
    checks = tuple(parse_pauli(label) for label in CHECK_LABELS)
    twist = checks[TWIST_CHECK_INDEX]
    packed = tuple(_packed(check, len(QUBIT_COORDINATES)) for check in checks)
    basis = _gf2_basis(packed)
    rank = gf2_rank(packed)
    pair_anticommutations = sum(
        not commute(first, second)
        for index, first in enumerate(checks)
        for second in checks[index + 1:]
    )
    distance = minimum_distance(checks, basis)

    z_side = parse_pauli(Z_SIDE_LABEL)
    x_side = parse_pauli(X_SIDE_LABEL)
    z_continuation = multiply_mod_phase(z_side, twist)
    x_continuation = multiply_mod_phase(x_side, twist)

    local_counts = {"projects_to_zero": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
    for candidate in local_paulis(len(QUBIT_COORDINATES), distance - 1):
        if not all(commute(candidate, check) for check in checks):
            local_counts["projects_to_zero"] += 1
        elif gf2_in_span(_packed(candidate, len(QUBIT_COORDINATES)), basis):
            local_counts["stabilizer_scalar"] += 1
        else:
            local_counts["nontrivial_logical"] += 1

    z_witness = (
        is_logical(z_side, checks, basis)
        and is_logical(z_continuation, checks, basis)
        and z_side.z != 0
        and z_side.x == 0
        and z_continuation.x != 0
        and multiply_mod_phase(z_side, z_continuation) == twist
    )
    x_witness = (
        is_logical(x_side, checks, basis)
        and is_logical(x_continuation, checks, basis)
        and x_side.x != 0
        and x_side.z == 0
        and x_continuation.z != 0
        and multiply_mod_phase(x_side, x_continuation) == twist
    )
    passes = (
        pair_anticommutations == 0
        and rank == 6
        and distance == 3
        and local_counts["nontrivial_logical"] == 0
        and symplectic_parity(z_side, x_side) == 1
        and z_witness
        and x_witness
    )
    if not passes:
        raise RuntimeError("minimal triangle twist reference failed its registered algebraic gate")

    output = {
        "schema": "antler.phase8c.triangle-twist-reference.v1",
        "parameters": {
            "new_declared_resource": "seven neutral qubits and local non-CSS stabilizer checks; not derived from ANTLER charge-two mediators",
            "geometry": "r=s=t=2 triangular stabilizer patch, with a twist-to-boundary reference junction",
            "qubit_coordinates": [list(coordinate) for coordinate in QUBIT_COORDINATES],
            "check_labels": list(CHECK_LABELS),
            "twist_check_index": TWIST_CHECK_INDEX,
            "twist_check": CHECK_LABELS[TWIST_CHECK_INDEX],
        },
        "stabilizer_algebra": {
            "check_count": len(checks),
            "pair_anticommutations": pair_anticommutations,
            "all_checks_commute": pair_anticommutations == 0,
            "independent_rank_over_GF2": rank,
            "encoded_qubits": len(QUBIT_COORDINATES) - rank,
            "ground_space_degeneracy": 1 << (len(QUBIT_COORDINATES) - rank),
            "exact_syndrome_gap_in_units_of_J": 2.0,
        },
        "local_protection": {
            "minimum_logical_weight": distance,
            "maximum_tested_weight": distance - 1,
            "tested_nonidentity_paulis": sum(local_counts.values()),
            **local_counts,
            "all_below_distance_probes_scalar_or_zero": local_counts["nontrivial_logical"] == 0,
        },
        "twist_to_boundary_string_witness": {
            "pure_z_side_representative": label_of(z_side),
            "z_side_continued_through_mixed_check": label_of(z_continuation),
            "pure_x_side_representative": label_of(x_side),
            "x_side_continued_through_mixed_check": label_of(x_continuation),
            "z_continuation_differs_by_twist_check": multiply_mod_phase(z_side, z_continuation) == twist,
            "x_continuation_differs_by_twist_check": multiply_mod_phase(x_side, x_continuation) == twist,
            "z_side_and_x_side_anticommute": symplectic_parity(z_side, x_side) == 1,
            "z_to_mixed_continuation_witness": z_witness,
            "x_to_mixed_continuation_witness": x_witness,
            "interpretation": "the central non-CSS check permits a stabilizer-equivalent string deformation that changes Pauli type across the twist-to-boundary junction",
        },
        "decision": (
            "PASS as a minimal local twist-to-boundary reference: the specified non-CSS triangle patch is a commuting [[7,1,3]] "
            "stabilizer code and supplies explicit X/Z-to-mixed string-deformation witnesses through its central junction. "
            "It is a calibration of the local e<->m branch-cut algebra, not a separated-defect fusion or braid result."
        ),
        "next_gate": (
            "T4 must embed two separated, geometrically explicit twist endpoints in one growing patch, prove their fusion-space "
            "dimension and local indistinguishability, and only then define defect transports."
        ),
        "claim_boundary": (
            "This finite reference inserts non-CSS stabilizer checks and a boundary. It does not derive a neutral link/check "
            "Hamiltonian from frozen ANTLER, demonstrate a thermodynamic wall, two separated twists, non-Abelian fusion, defect "
            "motion, a noncommutative braid, universality, noise resilience or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "triangle_twist_reference.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
