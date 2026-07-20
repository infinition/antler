"""Phase 5H: microscopic Kitaev T-junction zero-mode preflight.

This is the minimal number-nonconserving defect Hamiltonian needed by the
T-junction/fusion target.  It scans a three-arm p-wave junction, counts
near-zero BdG modes, measures the excitation gap, and records endpoint
localisation.  It is intentionally separate from the frozen number-conserving
ANTLER ladder; a future microscopic bridge must derive these pairing terms.
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


def bdg_tjunction(arm_length: int, mu: float, t: float, delta: float,
                  junction_chiral: float = 0.0,
                  junction_channel: str = "aa") -> tuple[np.ndarray, list[int]]:
    n, edges = junction_edges(arm_length)
    normal = -mu * np.eye(n, dtype=complex)
    pairing = np.zeros((n, n), complex)
    for source, target in edges:
        normal[source, target] = normal[target, source] = -t
        pairing[source, target] = delta
        pairing[target, source] = -delta
    # A trivalent real p-wave star has two accidental junction Majoranas in
    # addition to the centre and arm-end defects.  A phase-biased local
    # junction link gaps that pair while preserving the particle-hole BdG
    # structure.  It is an *extension target*, not a term in frozen ANTLER.
    if junction_chiral:
        if junction_channel not in {"aa", "bb"}:
            raise ValueError("junction_channel must be aa or bb")
        first_arm_1, first_arm_2 = 1, 1 + arm_length
        normal[first_arm_1, first_arm_2] += 1j * junction_chiral
        normal[first_arm_2, first_arm_1] -= 1j * junction_chiral
        sign = 1.0 if junction_channel == "aa" else -1.0
        pairing[first_arm_1, first_arm_2] += sign * 1j * junction_chiral
        pairing[first_arm_2, first_arm_1] -= sign * 1j * junction_chiral
    bdg = np.block([[normal, pairing], [-pairing.conj(), -normal.T]])
    endpoints = [arm * arm_length + arm_length for arm in range(3)]
    return bdg, endpoints


def analyse(arm_length: int, mu: float, t: float, delta: float, tolerance: float,
            junction_chiral: float = 0.0, junction_channel: str = "aa") -> dict:
    H, endpoints = bdg_tjunction(arm_length, mu, t, delta, junction_chiral,
                                 junction_channel)
    energies, vectors = np.linalg.eigh(H)
    zero = np.where(abs(energies) < tolerance)[0]
    nonzero = np.where(abs(energies) >= tolerance)[0]
    gap = float(np.min(abs(energies[nonzero]))) if len(nonzero) else 0.0
    profiles = []
    n = H.shape[0] // 2
    for index in zero:
        vector = vectors[:, index]
        weight = abs(vector[:n]) ** 2 + abs(vector[n:]) ** 2
        profiles.append({
            "energy": float(energies[index]),
            "endpoint_weight": float(np.sum(weight[endpoints])),
            "junction_weight": float(weight[0]),
            "max_site": int(np.argmax(weight)),
        })
    return {
        "arm_length": arm_length, "mu": mu, "t": t, "delta": delta,
        "junction_chiral": junction_chiral, "junction_channel": junction_channel,
        "zero_mode_count": int(len(zero)), "bulk_gap": gap,
        "endpoints": endpoints, "zero_mode_profiles": profiles,
        "fixed_parity_zero_mode_dimension": int(2 ** (len(zero) // 2 - 1))
        if len(zero) >= 2 and len(zero) % 2 == 0 else 0,
        "suitable_for_fixed_parity_fusion_space": bool(len(zero) == 4 and gap > 0.1 * abs(t)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm-lengths", type=int, nargs="+", default=[3, 4, 5, 6])
    parser.add_argument("--mus", type=float, nargs="+", default=[0.0, 0.1, 0.2])
    parser.add_argument("--deltas", type=float, nargs="+", default=[0.8, 1.0])
    parser.add_argument("--junction-chirals", type=float, nargs="+", default=[0.0, 0.25, 0.5, 1.0])
    parser.add_argument("--junction-channels", nargs="+", default=["aa", "bb"])
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--zero-tolerance", type=float, default=1e-8)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/tjunction_kitaev_preflight.json"))
    args = parser.parse_args()
    rows = [
        analyse(length, mu, args.t, delta, args.zero_tolerance, chiral, channel)
        for length in args.arm_lengths for mu in args.mus for delta in args.deltas
        for chiral in args.junction_chirals for channel in args.junction_channels
    ]
    rows.sort(key=lambda row: (
        not row["suitable_for_fixed_parity_fusion_space"],
        -row["zero_mode_count"], -row["bulk_gap"],
    ))
    viable = [row for row in rows if row["suitable_for_fixed_parity_fusion_space"]]
    out = {
        "schema": "antler.phase5.tjunction-kitaev-preflight.v1",
        "claim_boundary": (
            "This is a microscopic superconducting T-junction target, not yet a "
            "derivation from ANTLER's number-conserving correlated-hopping ladder."
        ),
        "best": rows[0], "viable_rows": viable, "rows": rows,
        "decision": (
            "derive_or_engineer_phase_biased_pairing_extension_then_simulate_braids" if viable else
            "no_viable_zero_mode_tjunction_in_scanned_range"
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"decision": out["decision"], "best": out["best"],
                      "viable_count": len(viable)}, indent=2))


if __name__ == "__main__":
    main()
