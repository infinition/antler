"""Phase 5K: static-disorder audit of the four-zero-mode Kitaev T-junction.

This is a microscopic BdG target with a phase-biased ``bb`` junction link.
For every realization it independently perturbs onsite potentials, the
normal/pairing amplitude of every arm bond, and the junction-link amplitude.
The audit asks whether the lowest four BdG modes remain spectrally isolated;
it is not a claim of a material-specific or fault-tolerance threshold.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def junction_edges(arm_length: int) -> tuple[int, list[tuple[int, int]]]:
    n = 1 + 3 * arm_length
    edges = []
    for arm in range(3):
        previous = 0
        for distance in range(arm_length):
            current = 1 + arm * arm_length + distance
            edges.append((previous, current))
            previous = current
    return n, edges


def disordered_bdg(arm_length: int, sigma: float, rng: np.random.Generator,
                   t: float, delta: float, junction_chiral: float) -> np.ndarray:
    n, edges = junction_edges(arm_length)
    normal = np.diag(-sigma * abs(t) * rng.standard_normal(n)).astype(complex)
    pairing = np.zeros((n, n), complex)
    for source, target in edges:
        hop = t * (1.0 + sigma * rng.standard_normal())
        pair = delta * (1.0 + sigma * rng.standard_normal())
        normal[source, target] = normal[target, source] = -hop
        pairing[source, target] = pair
        pairing[target, source] = -pair
    # The bb channel is the one that removes the accidental two-mode
    # degeneracy of the real trivalent p-wave star.
    first_arm_1, first_arm_2 = 1, 1 + arm_length
    kappa = junction_chiral * (1.0 + sigma * rng.standard_normal())
    normal[first_arm_1, first_arm_2] += 1j * kappa
    normal[first_arm_2, first_arm_1] -= 1j * kappa
    pairing[first_arm_1, first_arm_2] -= 1j * kappa
    pairing[first_arm_2, first_arm_1] += 1j * kappa
    return np.block([[normal, pairing], [-pairing.conj(), -normal.T]])


def statistics(values: np.ndarray) -> tuple[float, float]:
    return float(values.mean()), float(values.var(ddof=1))


def ensemble(arm_length: int, sigma: float, samples: int, rng: np.random.Generator,
             t: float, delta: float, junction_chiral: float) -> dict:
    rows = []
    for sample in range(samples):
        energies = np.linalg.eigvalsh(disordered_bdg(
            arm_length, sigma, rng, t, delta, junction_chiral
        ))
        magnitude = np.sort(np.abs(energies))
        low_split = float(magnitude[3])
        isolation_gap = float(magnitude[4])
        ratio = low_split / isolation_gap if isolation_gap else np.inf
        rows.append({
            "sample": sample,
            "four_mode_splitting": low_split,
            "isolation_gap": isolation_gap,
            "splitting_to_gap_ratio": ratio,
            "well_separated": bool(ratio < 0.1),
        })
    split = np.asarray([row["four_mode_splitting"] for row in rows])
    gap = np.asarray([row["isolation_gap"] for row in rows])
    ratio = np.asarray([row["splitting_to_gap_ratio"] for row in rows])
    split_mean, split_var = statistics(split)
    gap_mean, gap_var = statistics(gap)
    ratio_mean, ratio_var = statistics(ratio)
    return {
        "sigma": sigma,
        "samples": samples,
        "four_mode_splitting_mean": split_mean,
        "four_mode_splitting_variance": split_var,
        "isolation_gap_mean": gap_mean,
        "isolation_gap_variance": gap_var,
        "isolation_gap_min": float(gap.min()),
        "splitting_to_gap_ratio_mean": ratio_mean,
        "splitting_to_gap_ratio_variance": ratio_var,
        "splitting_to_gap_ratio_max": float(ratio.max()),
        "well_separated_fraction": float(np.mean(ratio < 0.1)),
        "realizations": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm-length", type=int, default=6)
    parser.add_argument("--sigmas", type=float, nargs="+",
                        default=[0.01, 0.03, 0.05, 0.10, 0.20])
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--delta", type=float, default=1.0)
    parser.add_argument("--junction-chiral", type=float, default=1.0)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/tjunction_kitaev_static_noise.json"))
    args = parser.parse_args()
    if args.arm_length < 2 or args.samples < 2 or any(sigma < 0 for sigma in args.sigmas):
        raise ValueError("arm length >= 2, samples >= 2, and non-negative sigmas are required")
    rng = np.random.default_rng(args.seed)
    ensembles = [ensemble(args.arm_length, sigma, args.samples, rng, args.t,
                          args.delta, args.junction_chiral) for sigma in args.sigmas]
    output = {
        "schema": "antler.phase5.tjunction-kitaev-static-noise.v1",
        "claim_boundary": (
            "This is a static BdG spectral-isolation audit for a phase-biased Kitaev "
            "T-junction target.  Pairing and the junction term are not derived from "
            "the frozen number-conserving ANTLER ladder."
        ),
        "noise_model": (
            "independent Gaussian perturbations on onsite potentials, each arm normal "
            "hopping, each pairing amplitude, and the bb junction coupling"
        ),
        "arm_length": args.arm_length, "t": args.t, "delta": args.delta,
        "junction_chiral": args.junction_chiral, "seed": args.seed,
        "ensembles": ensembles,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps([{
        key: row[key] for key in ("sigma", "four_mode_splitting_mean", "isolation_gap_mean",
                                   "well_separated_fraction")
    } for row in ensembles], indent=2))


if __name__ == "__main__":
    main()
