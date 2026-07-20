"""Phase 8C-T4: four corner-twist fusion-space calibration.

An odd-distance rotated planar surface-code patch has alternating X/Z boundary
types.  Their four corners are boundary twists: with fixed total topological
charge they calibrate a two-dimensional Ising-like fusion sector.  This script
constructs the CSS stabilizer family for d=3 and d=5, checks its full binary
stabilizer algebra, and exhaustively rejects all logical Pauli actions below
the geometric separation/distance.

It is an imposed neutral-qubit reference.  It calibrates fusion-space
protection and scaling before any microscopic ANTLER bridge or physical defect
motion is attempted; it does not manufacture a braid matrix.
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
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
    symplectic_parity,
)


DISTANCES = (3, 5)


def qubit(x: int, y: int, distance: int) -> int:
    return y * distance + x


def pauli(kind: str, support: tuple[int, ...] | list[int]) -> BinaryPauli:
    mask = sum(1 << position for position in support)
    return BinaryPauli(x=mask if kind == "X" else 0, z=mask if kind == "Z" else 0)


def label_of(word: BinaryPauli, qubits: int) -> str:
    out = []
    for position in range(qubits):
        x_bit = (word.x >> position) & 1
        z_bit = (word.z >> position) & 1
        out.append("Y" if x_bit and z_bit else "X" if x_bit else "Z" if z_bit else "I")
    return "".join(out)


def rotated_patch(distance: int) -> tuple[tuple[BinaryPauli, ...], list[dict[str, object]]]:
    """Return the conventional odd-distance rotated CSS patch.

    Interior checks are alternating X/Z four-body plaquettes.  Weight-two
    boundary checks alternate so each corner joins an X and a Z boundary.
    """
    checks: list[BinaryPauli] = []
    metadata: list[dict[str, object]] = []
    for y in range(distance - 1):
        for x in range(distance - 1):
            kind = "X" if (x + y) % 2 == 0 else "Z"
            support = (qubit(x, y, distance), qubit(x + 1, y, distance), qubit(x, y + 1, distance), qubit(x + 1, y + 1, distance))
            checks.append(pauli(kind, support))
            metadata.append({"kind": kind, "location": [x + 0.5, y + 0.5], "support": list(support), "boundary": False})
    # The four alternating boundary families close the planar patch.  For odd
    # distance they give exactly d^2-1 independent candidate checks.
    for x in range(0, distance - 1, 2):
        support = (qubit(x, 0, distance), qubit(x + 1, 0, distance))
        checks.append(pauli("Z", support))
        metadata.append({"kind": "Z", "location": [x + 0.5, 0.0], "support": list(support), "boundary": True, "side": "top"})
    for x in range(1, distance - 1, 2):
        support = (qubit(x, distance - 1, distance), qubit(x + 1, distance - 1, distance))
        checks.append(pauli("Z", support))
        metadata.append({"kind": "Z", "location": [x + 0.5, distance - 1.0], "support": list(support), "boundary": True, "side": "bottom"})
    for y in range(1, distance - 1, 2):
        support = (qubit(0, y, distance), qubit(0, y + 1, distance))
        checks.append(pauli("X", support))
        metadata.append({"kind": "X", "location": [0.0, y + 0.5], "support": list(support), "boundary": True, "side": "left"})
    for y in range(0, distance - 1, 2):
        support = (qubit(distance - 1, y, distance), qubit(distance - 1, y + 1, distance))
        checks.append(pauli("X", support))
        metadata.append({"kind": "X", "location": [distance - 1.0, y + 0.5], "support": list(support), "boundary": True, "side": "right"})
    return tuple(checks), metadata


def local_paulis(qubits: int, maximum_weight: int):
    for weight in range(1, maximum_weight + 1):
        for support in combinations(range(qubits), weight):
            for factors in product(((1, 0), (0, 1), (1, 1)), repeat=weight):
                x = sum(x_bit << site for site, (x_bit, _) in zip(support, factors))
                z = sum(z_bit << site for site, (_, z_bit) in zip(support, factors))
                yield BinaryPauli(x=x, z=z)


def local_gate(checks: tuple[BinaryPauli, ...], qubits: int, distance: int, basis: dict[int, int]) -> dict[str, int | bool]:
    counts: dict[str, int | bool] = {"projects_to_zero": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
    for candidate in local_paulis(qubits, distance - 1):
        if any(not commute(candidate, check) for check in checks):
            counts["projects_to_zero"] += 1
        elif gf2_in_span(_packed(candidate, qubits), basis):
            counts["stabilizer_scalar"] += 1
        else:
            counts["nontrivial_logical"] += 1
    counts["tested_nonidentity_paulis"] = sum(int(counts[key]) for key in ("projects_to_zero", "stabilizer_scalar", "nontrivial_logical"))
    counts["all_below_distance_probes_scalar_or_zero"] = counts["nontrivial_logical"] == 0
    return counts


def audit_distance(distance: int) -> dict[str, object]:
    checks, metadata = rotated_patch(distance)
    qubits = distance * distance
    packed = tuple(_packed(check, qubits) for check in checks)
    basis = _gf2_basis(packed)
    rank = gf2_rank(packed)
    pair_anticommutations = sum(
        not commute(first, second)
        for index, first in enumerate(checks)
        for second in checks[index + 1:]
    )
    logical_x = pauli("X", tuple(qubit(x, 0, distance) for x in range(distance)))
    logical_z = pauli("Z", tuple(qubit(0, y, distance) for y in range(distance)))
    local = local_gate(checks, qubits, distance, basis)
    corner_coordinates = [[0, 0], [distance - 1, 0], [distance - 1, distance - 1], [0, distance - 1]]
    return {
        "distance": distance,
        "qubits": qubits,
        "check_count": len(checks),
        "check_labels": [label_of(check, qubits) for check in checks],
        "check_metadata": metadata,
        "pair_anticommutations": pair_anticommutations,
        "independent_rank_over_GF2": rank,
        "encoded_qubits": qubits - rank,
        "ground_space_degeneracy": 1 << (qubits - rank),
        "exact_syndrome_gap_in_units_of_J": 2.0,
        "corner_twists": corner_coordinates,
        "nearest_corner_separation": distance - 1,
        "logical_x": label_of(logical_x, qubits),
        "logical_z": label_of(logical_z, qubits),
        "logical_pair_symplectic_parity": symplectic_parity(logical_x, logical_z),
        "logical_x_commutes_with_checks": all(commute(logical_x, check) for check in checks),
        "logical_z_commutes_with_checks": all(commute(logical_z, check) for check in checks),
        "logical_x_outside_stabilizers": not gf2_in_span(_packed(logical_x, qubits), basis),
        "logical_z_outside_stabilizers": not gf2_in_span(_packed(logical_z, qubits), basis),
        "local_protection": local,
    }


def main() -> None:
    rows = [audit_distance(distance) for distance in DISTANCES]
    for row in rows:
        distance = int(row["distance"])
        if not (
            int(row["pair_anticommutations"]) == 0
            and int(row["independent_rank_over_GF2"]) == int(row["qubits"]) - 1
            and int(row["encoded_qubits"]) == 1
            and int(row["ground_space_degeneracy"]) == 2
            and int(row["logical_pair_symplectic_parity"]) == 1
            and bool(row["logical_x_commutes_with_checks"])
            and bool(row["logical_z_commutes_with_checks"])
            and bool(row["logical_x_outside_stabilizers"])
            and bool(row["logical_z_outside_stabilizers"])
            and bool(row["local_protection"]["all_below_distance_probes_scalar_or_zero"])
        ):
            raise RuntimeError(f"corner-twist calibration failed at distance {distance}")

    output = {
        "schema": "antler.phase8c.corner-twist-fusion-calibration.v1",
        "parameters": {
            "new_declared_resource": "neutral stabilizer qubits with imposed local X/Z checks; not derived from ANTLER charge-two mediators",
            "family": "odd-distance rotated planar surface-code reference with alternating X/Z boundaries",
            "corner_interpretation": "each meeting of X and Z boundary is a boundary-twist calibration; the four-corner code doublet matches the fixed-sector dimension expected of a four-twist encoding",
            "distances": list(DISTANCES),
        },
        "rows": rows,
        "scaling": {
            "corner_separations": [row["nearest_corner_separation"] for row in rows],
            "code_distances": [row["distance"] for row in rows],
            "ground_space_dimensions": [row["ground_space_degeneracy"] for row in rows],
            "conclusion": "the protected doublet persists while the minimum logical support grows with the corner separation in this imposed reference family",
        },
        "decision": (
            "PASS T4 as an external four-boundary-twist code-space calibration. At d=3 and d=5 the locally indistinguishable "
            "two-dimensional code space has minimum logical support growing from 3 to 5. This calibrates fusion-space dimension "
            "and local protection before any defect-motion or braid protocol; it is not a physical fusion measurement or an "
            "ANTLER-derived anyon result."
        ),
        "next_gate": (
            "T5 must replace boundary-fixed corner twists by an explicitly deformable interior twist-pair geometry, derive local "
            "deformation Hamiltonians, and calculate two adjacent exchange holonomies with a nonzero commutator before reporting "
            "Yang-Baxter."
        ),
        "claim_boundary": (
            "This is a CSS stabilizer reference with boundary twists, not a microscopic ANTLER Hamiltonian. It does not demonstrate "
            "mobile interior defects, physical fusion measurement, non-Abelian braid matrices, universality, material noise or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "corner_twist_fusion_calibration.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
