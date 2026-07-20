"""Emit the complete the submitted candidate-ready Phase 7C initial observation S_0."""
from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    candidate_payload,
    compact_optimizer_observation,
    evaluate_local_candidate,
    registered_search_space,
    seeded_perturbative_candidate,
)
from antler.phase7_optimizer_monitor import OptimizationMonitor


def main() -> None:
    candidate = seeded_perturbative_candidate()
    audit = evaluate_local_candidate(candidate, target_label="XXXX", target_strength=1.0)
    monitor = OptimizationMonitor(ROOT / "results" / "phase7" / "optimizer_monitor", reset=True)
    monitor.record(iteration=0, audit=audit, tag="seeded_perturbative_negative_control")
    out = {
        "schema": "antler.phase7c.initial-optimizer-observation.v1",
        "target": {
            "local_stabilizer": "XXXX",
            "required_effective_term": "-0.5 * XXXX up to an additive identity for J_star=1",
            "promotion_condition": "connected virtual graph, sign-correct target coefficient, then selectivity on XXXX and ZZZZ",
        },
        "candidate_payload": candidate_payload(candidate),
        "action_space": registered_search_space(),
        "initial_loss": {
            "name": "fixed_scale_spectral_algebraic_residual",
            "value": audit["state_vector"]["fixed_scale_spectral_algebraic_residual"],
            "secondary_operator_loss": audit["state_vector"]["best_scale_operator_residual"],
            "reward": audit["reward"],
        },
        "s0": compact_optimizer_observation(audit),
        "monitoring": {
            "history_jsonl": "results/phase7/optimizer_monitor/history.jsonl",
            "history_csv": "results/phase7/optimizer_monitor/history.csv",
            "dashboard_png": "results/phase7/optimizer_monitor/dashboard.png",
        },
        "claim_boundary": (
            "S_0 is a local optimizer state and negative control. It does not establish a microscopic stabilizer, "
            "a controlled SW derivation, a tiled parent, topological order, braiding, or non-Abelian computation."
        ),
    }
    path = ROOT / "results" / "phase7" / "optimizer_s0_xxxx.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
