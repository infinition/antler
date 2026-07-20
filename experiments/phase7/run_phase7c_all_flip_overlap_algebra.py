"""Exact overlap algebra for the four-rung all-flip interaction.

The crossed mediator generates F_p = prod S+ + prod S- rather than XXXX.
This finite logical-space control determines whether two such identical terms
can be repeated on supports with 0..4 shared rungs.  It is an algebraic filter,
not a microscopic Schrieffer--Wolff calculation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


RUNG_COUNT = 8
FIRST_SUPPORT = (0, 1, 2, 3)
SECOND_SUPPORTS = {
    0: (4, 5, 6, 7),
    1: (3, 4, 5, 6),
    2: (2, 3, 4, 5),
    3: (1, 2, 3, 4),
    4: (0, 1, 2, 3),
}


def all_flip_operator(support: tuple[int, ...]) -> np.ndarray:
    raise_one = np.asarray([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    lower_one = raise_one.conj().T
    identity = np.eye(2, dtype=complex)
    raise_all = np.asarray([[1.0]], dtype=complex)
    lower_all = np.asarray([[1.0]], dtype=complex)
    for rung in range(RUNG_COUNT):
        raise_all = np.kron(raise_all, raise_one if rung in support else identity)
        lower_all = np.kron(lower_all, lower_one if rung in support else identity)
    return raise_all + lower_all


def main() -> None:
    first = all_flip_operator(FIRST_SUPPORT)
    rows = []
    for overlap, second_support in SECOND_SUPPORTS.items():
        second = all_flip_operator(second_support)
        commutator = first @ second - second @ first
        rows.append({
            "shared_rungs": overlap,
            "second_support": list(second_support),
            "commutator_frobenius": float(np.linalg.norm(commutator)),
            "commutator_spectral": float(np.linalg.norm(commutator, 2)),
        })
    out = {
        "schema": "antler.phase7c.all-flip-overlap-algebra.v1",
        "first_support": list(FIRST_SUPPORT),
        "rows": rows,
        "decision": (
            "Identical all-flip terms commute only when their supports are disjoint or identical. "
            "A nontrivial overlapping repetition therefore cannot serve as a commuting stabilizer family."
        ),
        "claim_boundary": (
            "This is an exact algebraic statement for the identified all-flip operator only. It does not exclude "
            "a different microscopic primitive, mixed stabilizer families, a noncommuting parent phase, or topological order."
        ),
    }
    path = ROOT / "results" / "phase7" / "all_flip_overlap_algebra.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
