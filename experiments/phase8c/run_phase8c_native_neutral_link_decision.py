"""Phase 8C-D1: bounded native neutral-link derivation verdict."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    source = json.loads((ROOT / "results" / "phase7" / "phase8b_jw_string_conditional_rung_bridge.json").read_text())
    rows = source["convention_rows"]
    frozen = [row for row in rows if int(row["string_weight"]) in {0, 2}]
    odd = [row for row in rows if int(row["string_weight"]) == 1]
    if not (all(abs(float(row["effective_xz_coefficient"])) < 1e-15 for row in frozen) and all(abs(float(row["effective_xz_coefficient"]) + 0.2) < 1e-15 for row in odd)):
        raise RuntimeError("registered JW resource distinction is not reproduced")
    output = {
        "schema": "antler.phase8c.native-neutral-link-decision.v1",
        "scope": "one bounded microscopic derivation attempt: a local generalized JW crossing bridge using only the frozen rung-major convention and charge-two mediators",
        "frozen_allowed_cases": [{"string_weight": int(row["string_weight"]), "mott_q": row["mott_q"], "effective_xz_coefficient": row["effective_xz_coefficient"], "largest_unwanted": row["largest_unwanted_non_scalar_pauli"]} for row in frozen],
        "counterfactual_required_case": [{"string_weight": int(row["string_weight"]), "mott_q": row["mott_q"], "effective_xz_coefficient": row["effective_xz_coefficient"], "unwanted_over_xz": row.get("unwanted_over_xz")} for row in odd],
        "decision": "NO-GO for this bounded direct derivation: frozen adjacent-rung string weight 0 and physical charge-two weight 2 both give exactly zero X_rail Z_link coefficient at every registered Mott depth, whereas the clean coefficient -0.2 requires the artificial odd weight 1. The required neutral/odd Z2 link is therefore not supplied by this ANTLER grammar.",
        "minimal_missing_resource": "an explicitly derived neutral odd Z2 link/statistical degree of freedom, or a different microscopic substrate with physical odd string weight",
        "next_gate": "Do not optimize pulses inside this closed grammar. Either declare and derive a new neutral-link resource, or close the native non-Abelian branch with this scoped no-go and publish the resource boundary.",
        "claim_boundary": "This closes one direct JW-crossing construction only. It does not prove that every possible ANTLER extension is impossible, and establishes no physical code, defect, braid, universality or topological quantum computer.",
    }
    result = ROOT / "results" / "phase8c" / "native_neutral_link_decision.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
