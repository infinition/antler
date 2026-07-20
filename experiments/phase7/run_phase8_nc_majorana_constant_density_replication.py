"""Constant-density small-ED replication gate for the dynamic NC-Majorana ladder.

This corrects a limitation of the earlier Phase 7D preflight: its L=4,N=2 to
L=6,N=4 continuation changed the filling from 1/4 to 1/3 and used a local
transfer filter that is not the bulk-gap/entanglement-spectrum diagnostic of
Defossez et al. (arXiv:2412.14886v2).

The present script is deliberately a small-system reproduction control at the
published illustrative point U0=-0.7, alpha=1/3, eta=pi/2.  It preserves
filling within two sequences and serializes the central-cut Schmidt spectrum.
The nu=1/4 sequence has odd N at L=6, so rail exchange forces equality of the
two branch-parity spectra; that split is explicitly non-diagnostic.  The
second, even-N nu=1/3 sequence is included so that parity-sector splitting is
not symmetry-forced. It is not a thermodynamic MPS calculation and cannot
establish a phase by itself.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_floquet_full_ladder_preflight import ETA, build_h0_and_rotation
from antler.number_conserving_pairwire import wire_a_parity


CASES = (
    ("published_illustrative_nu_one_quarter_odd_N", 4, 2, -0.7, 1.0 / 3.0),
    ("published_illustrative_nu_one_quarter_odd_N", 6, 3, -0.7, 1.0 / 3.0),
    ("published_illustrative_nu_one_third_even_N", 3, 2, -0.7, 1.0 / 3.0),
    ("published_illustrative_nu_one_third_even_N", 6, 4, -0.7, 1.0 / 3.0),
    # These three points were already registered in the Phase 7D small-ED
    # scan.  The alpha=0.5 point had the best L=6 split/gap there, but was
    # rejected only by the non-literature local-transfer proxy.
    ("phase7d_low_split_candidate_nu_one_third_even_N", 6, 4, -2.0, 0.25),
    ("phase7d_low_split_candidate_nu_one_third_even_N", 6, 4, -2.0, 0.50),
    ("phase7d_low_split_candidate_nu_one_third_even_N", 6, 4, -2.0, 0.75),
)
BUILD_CACHE: dict[tuple[int, int, float], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}


def cached_h0_and_rotation(L: int, N: int, u0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Avoid rebuilding/exponentiating the same microscopic block across alpha."""
    key = (L, N, u0)
    if key not in BUILD_CACHE:
        h0, rotation, states, _ = build_h0_and_rotation(L, N, u0)
        BUILD_CACHE[key] = h0, rotation, states
    return BUILD_CACHE[key]


def central_cut_spectrum(vector: np.ndarray, states: np.ndarray, L: int) -> dict:
    """Return the nonzero central-cut Schmidt weights without an ES criterion."""
    left_rungs = L // 2
    left_modes = 2 * left_rungs
    left_dimension = 1 << left_modes
    right_dimension = 1 << (2 * L - left_modes)
    amplitude = np.zeros((left_dimension, right_dimension), dtype=complex)
    for coefficient, raw_state in zip(vector, states):
        state = int(raw_state)
        left = state & ((1 << left_modes) - 1)
        right = state >> left_modes
        amplitude[left, right] = coefficient
    weights = np.linalg.eigvalsh(amplitude @ amplitude.conj().T)
    weights = np.sort(weights[weights > 1e-13])[::-1]
    levels = -np.log(weights)
    pair_differences = [float(abs(levels[i] - levels[i + 1])) for i in range(0, len(levels) - 1, 2)]
    entropy = float(-np.sum(weights * np.log(weights)))
    return {
        "cut_after_rungs": left_rungs,
        "schmidt_weights_descending": [float(value) for value in weights],
        "entanglement_levels_ascending": [float(value) for value in levels],
        "adjacent_pair_splittings": pair_differences,
        "von_neumann_entropy": entropy,
    }


def analyze(sequence: str, L: int, N: int, u0: float, alpha: float) -> dict:
    h0, rotation, states = cached_h0_and_rotation(L, N, u0)
    effective = alpha * h0 + (1.0 - alpha) * (rotation.conj().T @ h0 @ rotation)
    parity_matrix = np.diag([(-1.0 if wire_a_parity(int(state), L) else 1.0) for state in states])
    parity_residual = float(np.linalg.norm(effective @ parity_matrix - parity_matrix @ effective))
    sectors = {}
    for parity in (0, 1):
        rows = np.asarray(
            [i for i, state in enumerate(states) if wire_a_parity(int(state), L) == parity], dtype=int
        )
        values, vectors = eigh(effective[np.ix_(rows, rows)], subset_by_index=[0, 1], driver="evr")
        ground = np.zeros(len(states), dtype=complex)
        ground[rows] = vectors[:, 0]
        sectors[str(parity)] = {
            "ground_energy": float(values[0]),
            "first_neutral_excitation": float(values[1]),
            "neutral_gap": float(values[1] - values[0]),
            "central_cut": central_cut_spectrum(ground, states, L),
        }
    return {
        "constant_density_sequence": sequence,
        "L": L,
        "N": N,
        "filling_N_over_2L": N / (2.0 * L),
        "parameters": {"u0_attractive_nn": u0, "alpha": alpha, "eta": ETA, "t_leg": 1.0},
        "hilbert_dimension": len(states),
        "parity_commutator_frobenius": parity_residual,
        "parity_sector_ground_split": abs(sectors["0"]["ground_energy"] - sectors["1"]["ground_energy"]),
        "parity_split_interpretation": (
            "non-diagnostic: for odd total N, rail-exchange symmetry maps the two branch-parity sectors into one another"
            if N & 1 else "diagnostic candidate: even total N does not force the two branch-parity spectra to coincide"
        ),
        "smallest_neutral_gap": min(sectors["0"]["neutral_gap"], sectors["1"]["neutral_gap"]),
        "parity_sectors": sectors,
    }


def main() -> None:
    rows = [analyze(sequence, L, N, u0, alpha) for sequence, L, N, u0, alpha in CASES]
    out = {
        "schema": "antler.phase8.nc-majorana-constant-density-replication.v1",
        "citation": "Defossez et al., arXiv:2412.14886v2 (2025)",
        "purpose": "constant-density small-ED replication gate using bulk-gap and entanglement-spectrum observables",
        "model": "H_eff=alpha H0+(1-alpha) P^dag H0 P, P=exp(-i eta Jx)",
        "rows": rows,
        "decision": (
            "Small-ED evidence only. The next valid promotion is an independent constant-density large-L tensor-network "
            "replication with the paper's diagnostics; neither a positive nor a negative two-size result establishes a phase."
        ),
        "claim_boundary": (
            "This is an external-model replication control, not a derivation from ANTLER's frozen charge-two mediator hardware. "
            "It establishes no ANTLER phase, edge qubit, braid, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "nc_majorana_constant_density_replication.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
