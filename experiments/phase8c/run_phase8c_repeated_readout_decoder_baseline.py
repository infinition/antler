"""Phase 8C-T5l: exact independent-readout repetition/majority baseline."""
from __future__ import annotations

import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGES = 24
READOUT_FLIP_PROBABILITIES = (0.001, 0.01, 0.05)
REPETITIONS = (1, 3, 5, 7)


def majority_failure(probability: float, repetitions: int) -> float:
    return sum(comb(repetitions, errors) * probability**errors * (1.0 - probability) ** (repetitions - errors) for errors in range(repetitions // 2 + 1, repetitions + 1))


def main() -> None:
    rows = []
    for probability in READOUT_FLIP_PROBABILITIES:
        for repetitions in REPETITIONS:
            per_check = majority_failure(probability, repetitions)
            rows.append({
                "independent_readout_flip_probability": probability,
                "repetitions_per_check": repetitions,
                "majority_decoder_per_check_failure": per_check,
                "probability_at_least_one_wrong_decoded_outcome_in_24_stage_protocol": 1.0 - (1.0 - per_check) ** STAGES,
            })
    output = {
        "schema": "antler.phase8c.repeated-readout-decoder-baseline.v1",
        "parameters": {"measurement_stages": STAGES, "noise_model": "independent classical outcome-bit flips only; data faults, correlated errors, measurement backaction and time correlations excluded", "decoder": "odd-repetition majority vote"},
        "rows": rows,
        "decision": "PASS as a readout-only baseline: odd repetition suppresses the independent outcome-flip channel polynomially, but this is not a syndrome decoder for the full circuit and cannot establish fault tolerance.",
        "next_gate": "Integrate repeated checks with propagated data faults and a time-resolved decoder, then test correlated circuit noise before attempting any microscopic ANTLER derivation.",
        "claim_boundary": "This is an exact binomial calculation for an inserted classical readout-noise model. It does not include physical measurement dynamics, code thresholds, topological protection, non-Abelian braiding or an ANTLER realization.",
    }
    result = ROOT / "results" / "phase8c" / "repeated_readout_decoder_baseline.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
