"""Verify the final public numerical closure without running long simulations."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
payload = json.loads((ROOT / "results" / "phase4_7" / "publication_closure.json").read_text())
assert payload["gates"] == {
    "all_points_dt_0p125": True,
    "leakage_below_1e-4": True,
    "sigma_min_above_0p9999": True,
    "positive_handoff_gap": True,
}
assert abs(payload["fit"]["p"] - 2.286237051902575) < 1e-12
assert payload["fit"]["r2"] > 0.9999
print("PASS: public Abelian Phase 4.7 closure")
