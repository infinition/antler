"""Archive and audit the submitted first mixed-species ring proposal.

The raw proposal labels its mixed ``a-b`` conversion mediators as parity-even.
The environment must reject that metadata claim.  A second, explicitly marked
counterfactual replaces only the labels by the parity actually implied by the
same pair terms, so that the underlying Hamiltonian can be audited without
silently repairing the submitted payload.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import (
    candidate_from_payload,
    compact_optimizer_observation,
    evaluate_local_candidate,
)


RAW_SUBMISSION = {
    "mott_u": 20.0,
    "maximum_sw_ratio": 0.15,
    "channels": [
        {"name": "even_pair_01", "detuning": 10.0, "coupling": 0.5, "phase": 0.0,
         "mediator_parity_signature": [0, 0],
         "pair_terms": [["a0", "a1", [1.0, 0.0]], ["b0", "b1", [1.0, 0.0]]]},
        {"name": "even_pair_12", "detuning": 10.0, "coupling": 0.5, "phase": 0.0,
         "mediator_parity_signature": [0, 0],
         "pair_terms": [["a1", "b2", [1.0, 0.0]], ["b1", "a2", [1.0, 0.0]]]},
        {"name": "even_pair_23", "detuning": 10.0, "coupling": 0.5, "phase": 0.0,
         "mediator_parity_signature": [0, 0],
         "pair_terms": [["a2", "a3", [1.0, 0.0]], ["b2", "b3", [1.0, 0.0]]]},
        {"name": "even_pair_30", "detuning": 10.0, "coupling": 0.5, "phase": 0.0,
         "mediator_parity_signature": [0, 0],
         "pair_terms": [["a3", "b0", [1.0, 0.0]], ["b3", "a0", [1.0, 0.0]]]},
    ],
    "rail_hops": [], "zz_couplings": [], "rail_biases": [0.0, 0.0, 0.0, 0.0],
}


def main() -> None:
    raw_error = None
    try:
        candidate_from_payload(RAW_SUBMISSION)
    except ValueError as error:
        raw_error = str(error)
    if raw_error is None:
        raise RuntimeError("raw the submitted candidate payload should fail its false mediator-parity metadata")
    corrected = deepcopy(RAW_SUBMISSION)
    corrected["channels"][1]["mediator_parity_signature"] = [1, 1]
    corrected["channels"][3]["mediator_parity_signature"] = [1, 1]
    candidate = candidate_from_payload(corrected)
    rows = {}
    for target in ("XXXX", "ZZZZ"):
        audit = evaluate_local_candidate(candidate, target_label=target)
        observation = compact_optimizer_observation(audit)
        rows[target] = {
            "compact_observation": observation,
            "observed_target_coefficient": audit["state_vector"]["observed_target_coefficient"],
            "reward": audit["reward"],
        }
    if not rows["XXXX"]["compact_observation"]["interaction_connectivity"]["patch_connected"]:
        raise RuntimeError("counterfactual mixed ring should connect the four-edge patch")
    if abs(rows["XXXX"]["observed_target_coefficient"]) > 1e-7:
        raise RuntimeError("unexpected large stabilizer coefficient; inspect this candidate before classification")
    out = {
        "schema": "antler.phase7c.mixed-species-ring-audit.v1",
        "raw_submission_status": "rejected",
        "raw_submission_rejection": raw_error,
        "counterfactual_metadata_corrected_only": {
            "description": (
                "The pair terms are unchanged. Only the two mixed mediator signatures are corrected from [0,0] to [1,1], "
                "which changes the exact microscopic symmetry from bare rail parities to mediator-dressed parities."
            ),
            "candidate_parity_type": candidate.parity_type(),
            "audits": rows,
        },
        "decision": (
            "Reject as a target primitive. The submitted payload has false parity metadata; after metadata-only correction, "
            "the connected mixed ring still has XXXX and ZZZZ coefficients of order 1e-10 rather than -0.5, while two-body "
            "Pauli terms remain of order 2.5e-3."
        ),
        "claim_boundary": (
            "This is a local counterfactual audit, not a repair or promotion of the submitted proposal. It establishes no "
            "microscopic stabilizer, global parent, topological order, braid or non-Abelian result."
        ),
    }
    path = ROOT / "results" / "phase7" / "mixed_species_ring_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "raw_submission_status": out["raw_submission_status"],
        "raw_submission_rejection": raw_error,
        "counterfactual_parity_type": candidate.parity_type(),
        "xxxx_coefficient": rows["XXXX"]["observed_target_coefficient"],
        "xxxx_reward": rows["XXXX"]["reward"],
    }, indent=2))


if __name__ == "__main__":
    main()
