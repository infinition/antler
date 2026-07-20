"""Phase 7 harness self-test on a deliberately unprotected classical code.

The toy Hamiltonian is a commuting-projector ferromagnetic chain.  Its two
ground states have an exactly conserved edge Z operator, yet a local density
reads the code directly.  Passing the algebra/edge checks while failing local
indistinguishability is the required negative control for the Phase 7 audit.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_parent_audit import (
    local_indistinguishability,
    local_projector_algebra,
    projected_edge_metrics,
    sector_spectrum,
    symmetry_audit,
)


def site_z(site: int, sites: int) -> np.ndarray:
    return np.diag([
        1.0 if ((state >> site) & 1) == 0 else -1.0
        for state in range(1 << sites)
    ])


def site_flip(site: int, sites: int) -> np.ndarray:
    dimension = 1 << sites
    out = np.zeros((dimension, dimension), dtype=complex)
    for state in range(dimension):
        out[state ^ (1 << site), state] = 1.0
    return out


def domain_wall_projector(left: int, right: int, sites: int) -> np.ndarray:
    return np.diag([
        float(((state >> left) & 1) != ((state >> right) & 1))
        for state in range(1 << sites)
    ])


def main() -> None:
    sites = 3
    terms = [
        domain_wall_projector(0, 1, sites),
        domain_wall_projector(1, 2, sites),
    ]
    H = sum(terms, start=np.zeros_like(terms[0]))
    frame = np.zeros((1 << sites, 2), dtype=complex)
    frame[0, 0] = 1.0       # |000>
    frame[(1 << sites) - 1, 1] = 1.0  # |111>
    total_charge = [state.bit_count() for state in range(1 << sites)]
    total_parity = [1.0 if charge % 2 == 0 else -1.0 for charge in total_charge]
    algebra = local_projector_algebra(terms, supports=((0, 1), (1, 2)), names=("Pi_01", "Pi_12"))
    good_edge = projected_edge_metrics(H, frame, site_z(0, sites))
    bad_edge = projected_edge_metrics(H, frame, site_flip(0, sites))
    local = local_indistinguishability(frame, {"left_Z": site_z(0, sites)})
    symmetries = symmetry_audit(H, {"total_charge": total_charge, "total_parity": total_parity})
    charge_one = sector_spectrum(H, total_charge, sector=1, code_multiplicity=1)
    passed = (
        algebra["all_terms_are_projectors_at_1e_minus_10"]
        and algebra["all_terms_commute_at_1e_minus_10"]
        and symmetries["all_commute_at_1e_minus_10"]
        and good_edge["code_commutator_action_normalized"] < 1e-12
        and bad_edge["code_commutator_action_normalized"] > 1e-3
        and local["worst_projected_non_scalar_frobenius"] > 1e-3
    )
    if not passed:
        raise RuntimeError("Phase 7 harness self-test failed to separate the controls")
    out = {
        "schema": "antler.phase7.harness-selftest.v1",
        "model": "three-site diagonal domain-wall commuting-projector toy",
        "purpose": "negative control: exact edge conservation without local indistinguishability",
        "local_projector_algebra": algebra,
        "symmetry_audit": symmetries,
        "exactly_conserved_but_unprotected_left_Z": good_edge,
        "nonconserved_edge_negative_control_left_X": bad_edge,
        "local_indistinguishability": local,
        "charge_one_sector": charge_one,
        "passes_harness_selftest": bool(passed),
        "decision": (
            "The harness accepts exact local projector algebra and edge conservation, but rejects the toy code because a local edge probe distinguishes its ground states. "
            "A Phase 7 candidate must pass both categories."
        ),
        "claim_boundary": (
            "This is a software control only. The toy has no U(1)-fixed topological code, no protected edge mode, and no braid interpretation."
        ),
    }
    path = ROOT / "results" / "phase7" / "harness_selftest.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
