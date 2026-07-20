"""Exact algebraic completion of a minimal mixed-check stabilizer patch.

The Phase-8B local audits derive YXZX and XZXZ walker words, including their
shared-link compatibility.  Here they seed a deliberately small completion on
seven gauge links.  All checks have four-link support and use only the
phase-programmable X/Y/Z walker vocabulary.  The completion is searched
algebraically subject to commuting, rank-six and distance-three gates, then
frozen below as an independently reproducible reference.

This is an abstract finite stabilizer-code control.  Only the first adjacent
mixed-check pair has a joint microscopic walker downfolding so far; the other
overlaps and a scalable geometry remain future gates.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_stabilizer_algebra import (
    BinaryPauli,
    _gf2_basis,
    _packed,
    commute,
    gf2_in_span,
    gf2_rank,
    local_paulis,
    symplectic_parity,
)


N_QUBITS = 7
CHECK_LABELS = (
    "YXZXIII",  # audited local YXZX loop
    "XZIIXZI",  # audited neighboring XZXZ loop
    "XZXYIII",
    "YXIIYXI",
    "XIXIXIX",
    "YIZIYIY",
)
SEEDED_OVERLAP = (CHECK_LABELS[0], CHECK_LABELS[1])


def parse_pauli(label: str) -> BinaryPauli:
    if len(label) != N_QUBITS:
        raise ValueError(label)
    x = sum((1 << qubit) for qubit, factor in enumerate(label) if factor in "XY")
    z = sum((1 << qubit) for qubit, factor in enumerate(label) if factor in "YZ")
    return BinaryPauli(x=x, z=z)


def label_of(pauli: BinaryPauli) -> str:
    symbols = []
    for qubit in range(N_QUBITS):
        x_bit, z_bit = (pauli.x >> qubit) & 1, (pauli.z >> qubit) & 1
        symbols.append("Y" if x_bit and z_bit else "X" if x_bit else "Z" if z_bit else "I")
    return "".join(symbols)


def minimum_logicals(checks: tuple[BinaryPauli, ...], basis: dict[int, int]) -> tuple[BinaryPauli, BinaryPauli]:
    centralizer = []
    for weight in range(1, N_QUBITS + 1):
        for candidate in local_paulis(N_QUBITS, weight):
            if all(commute(candidate, check) for check in checks):
                if not gf2_in_span(_packed(candidate, N_QUBITS), basis):
                    centralizer.append(candidate)
        if centralizer:
            break
    logical_x = centralizer[0]
    logical_z = next(candidate for candidate in centralizer if symplectic_parity(logical_x, candidate) == 1)
    return logical_x, logical_z


def main() -> None:
    checks = tuple(parse_pauli(label) for label in CHECK_LABELS)
    pair_anticommutations = [
        (left, right)
        for left, first in enumerate(checks)
        for right, second in enumerate(checks[left + 1:], start=left + 1)
        if not commute(first, second)
    ]
    packed = tuple(_packed(check, N_QUBITS) for check in checks)
    rank = gf2_rank(packed)
    basis = _gf2_basis(packed)
    logical_x, logical_z = minimum_logicals(checks, basis)
    distance = logical_x.weight()
    local_counts = {"anticommutes_with_check": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
    for candidate in local_paulis(N_QUBITS, distance - 1):
        if not all(commute(candidate, check) for check in checks):
            local_counts["anticommutes_with_check"] += 1
        elif gf2_in_span(_packed(candidate, N_QUBITS), basis):
            local_counts["stabilizer_scalar"] += 1
        else:
            local_counts["nontrivial_logical"] += 1
    output = {
        "schema": "antler.phase8b.mixed-check-patch-algebra.v1",
        "parameters": {
            "gauge_links": N_QUBITS,
            "stabilizer_labels": list(CHECK_LABELS),
            "all_check_weights": [check.weight() for check in checks],
            "seeded_jointly_audited_overlap": list(SEEDED_OVERLAP),
            "remaining_check_words": list(CHECK_LABELS[2:]),
            "candidate_rule": (
                "minimal seven-link, six-check completion seeded by the audited YXZX/XZXZ pair; "
                "all words are four-link products in the declared phase-programmable walker vocabulary"
            ),
        },
        "stabilizer_algebra": {
            "pair_anticommutations": pair_anticommutations,
            "all_checks_commute": not pair_anticommutations,
            "independent_rank_over_GF2": rank,
        },
        "code": {
            "encoded_qubits": N_QUBITS - rank,
            "ground_space_degeneracy": 1 << (N_QUBITS - rank),
            "minimum_logical_weight": distance,
            "one_logical_x_representative": label_of(logical_x),
            "one_anticommuting_logical_z_representative": label_of(logical_z),
            "logical_pair_symplectic_parity": symplectic_parity(logical_x, logical_z),
            "commuting_parent": "H_patch = -J sum_a S_a",
            "exact_syndrome_gap_in_units_of_J": 2.0,
        },
        "complete_projected_local_pauli_gate": {
            "maximum_tested_weight": distance - 1,
            "tested_nonidentity_paulis": sum(local_counts.values()),
            **local_counts,
            "all_tested_probes_project_to_scalars_or_zero": local_counts["nontrivial_logical"] == 0,
        },
        "decision": (
            "The seeded mixed-check vocabulary admits a finite seven-link commuting completion with one encoded qubit, "
            "distance three and a nonzero exact stabilizer-parent syndrome gap."
        ),
        "claim_boundary": (
            "This is a finite algebraic reference completion, not a complete dislocation/twist geometry. Only its first "
            "mixed-check overlap is microscopically walker-audited; all other overlaps, multi-walker crosstalk, physical "
            "U(1) embedding, scalable topological order, twist fusion, defect motion, non-Abelian braid, universality and "
            "fault tolerance remain unestablished."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_mixed_check_patch_algebra.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
