"""Run the registered Phase 7C local microscopic reward negative control."""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    compact_optimizer_observation,
    computational_budget,
    evaluate_local_candidate,
    registered_search_space,
    seeded_perturbative_candidate,
)


def main() -> None:
    candidate = seeded_perturbative_candidate()
    star = evaluate_local_candidate(candidate, target_label="XXXX", target_strength=1.0)
    plaquette = evaluate_local_candidate(candidate, target_label="ZZZZ", target_strength=1.0)
    if star["hard_failures"] or plaquette["hard_failures"]:
        raise RuntimeError("registered microscopic seed should satisfy structural, not target, gates")
    if abs(star["state_vector"]["observed_target_coefficient"]) > 1e-10:
        raise RuntimeError("disconnected pair seed unexpectedly generated a four-edge XXXX term")
    if abs(plaquette["state_vector"]["observed_target_coefficient"]) > 1e-10:
        raise RuntimeError("disconnected pair seed unexpectedly generated a four-edge ZZZZ term")
    out = {
        "schema": "antler.phase7.microscopic-reward-baseline.v1",
        "computational_budget": computational_budget(),
        "star_xxxx_negative_control": star,
        "plaquette_zzzz_negative_control": plaquette,
        "optimizer_interface": {
            "registered_search_space": registered_search_space(),
            "star_compact_observation": compact_optimizer_observation(star),
            "plaquette_compact_observation": compact_optimizer_observation(plaquette),
        },
        "decision": (
            "The symmetry-constrained, perturbative two-pair seed has a clean isolated monomer block but no connected "
            "four-edge stabilizer coefficient. The reward correctly records this as a non-solution rather than rewarding "
            "its unrelated pair processes."
        ),
        "claim_boundary": (
            "This is an optimizer-environment negative control. It does not identify a microscopic ANTLER realization."
        ),
    }
    path = ROOT / "results" / "phase7" / "microscopic_reward_baseline.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "local_dimensions": out["computational_budget"]["local_four_edge_fixed_charge_dimensions"],
        "global_3x3_no_mediator_dimension": out["computational_budget"]["global_3x3_without_mediators_fixed_charge_dimension"],
        "star_target_coefficient": star["state_vector"]["observed_target_coefficient"],
        "plaquette_target_coefficient": plaquette["state_vector"]["observed_target_coefficient"],
        "star_reward": star["reward"],
        "plaquette_reward": plaquette["reward"],
    }, indent=2))


if __name__ == "__main__":
    main()
