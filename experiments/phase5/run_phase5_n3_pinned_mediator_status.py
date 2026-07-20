"""Consolidate the two independent no-go diagnostics for the pinned N=3 mediator.

The static screen reports a near-degenerate logical pair and a small but finite
subspace-isolation gap.  The holonomy run separately tests (i) whether the two
exchanges are noncommuting and (ii) whether its overlap-tracked logical
subspace remains continuous enough to support the proposed adiabatic control.
The result deliberately does not turn a raw Yang--Baxter residual into evidence
when the commutator is below its stated threshold.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--static", type=Path,
                        default=Path("results/phase5/n3_mediator_gap_optimization.json"))
    parser.add_argument("--holonomy", type=Path,
                        default=Path("results/phase5/n3_local_exchange_holonomy.json"))
    parser.add_argument("--continuity-threshold", type=float, default=0.95,
                        help="minimum link singular value required by this audit")
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/n3_pinned_mediator_no_go.json"))
    args = parser.parse_args()
    if not 0.0 < args.continuity_threshold <= 1.0:
        raise ValueError("continuity-threshold must lie in (0, 1]")
    static, holonomy = load(args.static), load(args.holonomy)
    best = static["best"]
    spectra = best["spectra"]
    max_split = max(row["logical_split"] for row in spectra)
    min_gap = min(row["isolation_gap_to_low_energy_complement"] for row in spectra)
    min_link = min(
        holonomy["audit_A"]["minimum_parallel_transport_link_singular_value"],
        holonomy["audit_B"]["minimum_parallel_transport_link_singular_value"],
    )
    algebra_nonabelian = bool(holonomy["braid_relation_interpretable"])
    continuity_pass = min_link >= args.continuity_threshold
    output = {
        "schema": "antler.phase5.n3-pinned-mediator-status.v1",
        "claim_boundary": (
            "This rejects the tested static pinned-mediator route.  It does not prove "
            "that every mobile-mediator or enlarged-graph N=3 model is impossible."
        ),
        "static_near_crossing": {
            "logical_split_max": max_split,
            "subspace_isolation_gap_min": min_gap,
            "gap_over_split": min_gap / max_split if max_split else None,
            "interpretation": (
                "The ~1e-6 internal splitting is a near-degeneracy, not the gap to "
                "the complement.  It therefore requires subspace transport rather than "
                "independent-eigenstate adiabatic following."
            ),
        },
        "algebra_gate": {
            "commutator_norm": holonomy["commutator_norm"],
            "noncommuting_threshold": holonomy["noncommuting_threshold"],
            "passes": algebra_nonabelian,
            "braid_relation_status": holonomy["braid_relation_status"],
            "raw_braid_relation_residual": holonomy["raw_braid_relation_residual"],
        },
        "adiabatic_continuity_gate": {
            "minimum_parallel_transport_link_singular_value": min_link,
            "required_minimum": args.continuity_threshold,
            "passes": continuity_pass,
            "interpretation": (
                "This is a quality gate for the sampled control path, not a universal "
                "no-go theorem.  The selected branch is not smooth enough at this "
                "resolution to promote a finite-time braid."
            ),
        },
        "decision": (
            "reject_pinned_mediator_route__commuting_controls_and_unqualified_adiabatic_continuity"
            if not algebra_nonabelian and not continuity_pass else
            "do_not_promote_without_resolving_failed_gate"
        ),
        "sources": {"static": str(args.static), "holonomy": str(args.holonomy)},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({
        "decision": output["decision"],
        "commutator_norm": output["algebra_gate"]["commutator_norm"],
        "minimum_link_singular_value": min_link,
    }, indent=2))


if __name__ == "__main__":
    main()
