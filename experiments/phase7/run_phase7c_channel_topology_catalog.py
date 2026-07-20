"""Enumerate the registered four-edge charge-two channel topologies.

This is the discrete half of the classical baseline.  At fixed perturbative
parameters it screens every connected three- or four-link graph whose channel
on a link is either ``E`` (aa/bb, bare rail parities) or ``M`` (ab/ba,
mediator-dressed parities).  It does not tune continuous parameters or make a
global-phase claim.
"""
from __future__ import annotations

from itertools import combinations, product
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    EDGE_COUNT,
    MicroscopicCandidate2D,
    PairConversionChannel,
    evaluate_local_candidate,
    rail_mode,
)


LINKS = tuple(combinations(range(EDGE_COUNT), 2))
FIXED = {"mott_u": 20.0, "detuning": 10.0, "coupling": 0.5, "maximum_sw_ratio": 0.15}


def connected(links: tuple[tuple[int, int], ...]) -> bool:
    reached, frontier = {0}, [0]
    adjacency = {edge: set() for edge in range(EDGE_COUNT)}
    for left, right in links:
        adjacency[left].add(right)
        adjacency[right].add(left)
    while frontier:
        edge = frontier.pop()
        unseen = adjacency[edge] - reached
        reached |= unseen
        frontier.extend(unseen)
    return len(reached) == EDGE_COUNT


def channel(left: int, right: int, kind: str) -> PairConversionChannel:
    if kind == "E":
        pairs = ((rail_mode(left, 0), rail_mode(right, 0), 1.0), (rail_mode(left, 1), rail_mode(right, 1), 1.0))
    elif kind == "M":
        pairs = ((rail_mode(left, 0), rail_mode(right, 1), 1.0), (rail_mode(left, 1), rail_mode(right, 0), 1.0))
    else:
        raise ValueError("channel type must be E or M")
    return PairConversionChannel(
        name=f"{kind}_{left}{right}", detuning=FIXED["detuning"], coupling=FIXED["coupling"], phase=0.0,
        pair_terms=pairs,
    )


def row(links: tuple[tuple[int, int], ...], kinds: tuple[str, ...]) -> dict:
    candidate = MicroscopicCandidate2D(
        mott_u=FIXED["mott_u"], maximum_sw_ratio=FIXED["maximum_sw_ratio"],
        channels=tuple(channel(left, right, kind) for (left, right), kind in zip(links, kinds)),
    )
    audit = evaluate_local_candidate(candidate, target_label="XXXX", target_strength=1.0)
    state = audit["state_vector"]
    return {
        "topology": ",".join(f"{kind}{left}{right}" for (left, right), kind in zip(links, kinds)),
        "link_count": len(links),
        "links": [list(link) for link in links],
        "channel_kinds": list(kinds),
        "parity_type": candidate.parity_type(),
        "observed_xxxx_coefficient": state["observed_target_coefficient"],
        "sign_correct_xxxx_strength": max(0.0, -state["observed_target_coefficient"]),
        "target_alignment": state["target_alignment"],
        "operator_residual": state["fixed_scale_spectral_algebraic_residual"],
        "unwanted_pauli_norm_over_target": state["unwanted_pauli_norm_over_target"],
        "minimum_monomer_overlap": state["minimum_monomer_overlap_singular_value"],
        "low_high_gap": state["low_to_high_gap"],
        "reward": audit["reward"],
        "hard_failures": audit["hard_failures"],
    }


def summary(rows: list[dict], parity_type: str) -> dict:
    subset = [item for item in rows if item["parity_type"] == parity_type]
    ordered = sorted(subset, key=lambda item: item["sign_correct_xxxx_strength"], reverse=True)
    return {
        "count": len(subset),
        "best_five_by_sign_correct_xxxx": ordered[:5],
        "maximum_sign_correct_xxxx_strength": ordered[0]["sign_correct_xxxx_strength"] if ordered else 0.0,
        "candidates_above_1e_minus_7": sum(item["sign_correct_xxxx_strength"] > 1e-7 for item in subset),
    }


def main() -> None:
    rows = []
    for link_count in (3, 4):
        for links in combinations(LINKS, link_count):
            if not connected(links):
                continue
            for kinds in product(("E", "M"), repeat=link_count):
                rows.append(row(links, kinds))
    expected = 368
    if len(rows) != expected:
        raise RuntimeError(f"topology enumeration changed unexpectedly: {len(rows)} != {expected}")
    out = {
        "schema": "antler.phase7c.fixed-scale-channel-topology-catalog.v1",
        "fixed_parameters": FIXED,
        "scope": "all connected three- and four-link E/M channel graphs on four edges, evaluated against XXXX",
        "topology_count": len(rows),
        "bare_rail_parity_class": summary(rows, "bare_rail_parities"),
        "mediator_dressed_parity_class": summary(rows, "mediator_dressed_parities"),
        "rows": rows,
        "decision": (
            "This is a discrete classical prefilter. Only a topology with a sign-correct coefficient above numerical noise "
            "may enter continuous optimization; a high local gap or graph connectivity alone does not qualify it."
        ),
        "claim_boundary": (
            "A fixed-scale four-edge catalog neither proves a controlled Schrieffer-Wolff expansion nor a tiled two-dimensional "
            "parent, topological order, braiding, non-Abelian statistics or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "channel_topology_catalog.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "topology_count": len(rows),
        "bare": out["bare_rail_parity_class"],
        "dressed": out["mediator_dressed_parity_class"],
    }, indent=2))


if __name__ == "__main__":
    main()
