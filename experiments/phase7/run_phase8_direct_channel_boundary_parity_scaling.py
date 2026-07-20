"""Audit the claimed rail-parity flip at a direct-channel segment boundary.

The Phase-8 direct-channel bridge reuses two physical charge-two mediator
slots per link while changing their pair-conversion channels between H0 and
H1.  No fixed microscopic rail parity is then available.  This script asks a
narrower, falsifiable question: after a virtual-Rabi-closed H0 -> H1 cycle,
does the projected low-space rail-parity flip scale as (g/Delta)^2?

It includes a control with four *separate* mediator species per link: H0 uses
one fixed pair and H1 the other.  That control has an exact full-Hilbert-space
parity operator, so every parity-flip metric must vanish to numerical
precision.  Neither result establishes a fusion space or a braid.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from antler.basis import site_index
from antler.number_conserving_pairwire import _apply
from run_phase8_direct_channel_ramp_audit import channel_path
from run_phase8_native_direct_h1_closure import build_micro
from run_phase8_native_micro_floquet_l3 import polar_unitary


LENGTH = 3
PARTICLE_NUMBER = 2
TARGET_U0 = -1.5
ALPHA = 0.5
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025, 0.01875, 0.0125)
CYCLE_COUNTS = (1, 2, 4, 8)


def species_mode(length: int, link: int, family: int, slot: int) -> int:
    """Two charge-two channels in each of the H0 and H1 mediator families."""
    return 2 * length + 4 * link + 2 * family + slot


def species_weighted_basis(length: int, particle_number: int) -> tuple[np.ndarray, dict[int, int]]:
    charges = [1] * (2 * length) + [2] * (4 * (length - 1))
    states = np.asarray([
        state for state in range(1 << len(charges))
        if sum(((state >> mode) & 1) * charge for mode, charge in enumerate(charges)) == particle_number
    ], dtype=np.int64)
    return states, {int(state): position for position, state in enumerate(states)}


def build_separate_species(
    length: int, particle_number: int, channels: np.ndarray, family: int, g: float, detuning: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One segment with its own mediator family; inactive mediators remain detuned."""
    states, index = species_weighted_basis(length, particle_number)
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        for link in range(length - 1):
            for rail in (0, 1):
                first, second = site_index(link, rail), site_index(link + 1, rail)
                for operations in ((("ann", second), ("create", first)), (("ann", first), ("create", second))):
                    item = _apply(state, operations)
                    if item is not None:
                        new_state, amplitude = item
                        hamiltonian[index[new_state], column] += -amplitude
            for mediator_family in (0, 1):
                for slot in (0, 1):
                    mediator = species_mode(length, link, mediator_family, slot)
                    hamiltonian[column, column] += detuning * ((state >> mediator) & 1)
            for slot in (0, 1):
                mediator = species_mode(length, link, family, slot)
                if (state >> mediator) & 1:
                    continue
                pair_modes = (
                    (site_index(link, 0), site_index(link + 1, 0)),
                    (site_index(link, 0), site_index(link + 1, 1)),
                    (site_index(link, 1), site_index(link + 1, 0)),
                    (site_index(link, 1), site_index(link + 1, 1)),
                )
                for pair_slot, (first, second) in enumerate(pair_modes):
                    if not ((state >> first) & 1 and (state >> second) & 1):
                        continue
                    item = _apply(state, (("ann", second), ("ann", first), ("create", mediator)))
                    if item is None:
                        raise RuntimeError("valid pair conversion vanished")
                    new_state, amplitude = item
                    value = -g * channels[slot, pair_slot] * amplitude
                    hamiltonian[index[new_state], column] += value
                    hamiltonian[column, index[new_state]] += np.conj(value)
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-11):
        raise RuntimeError("separate-species Hamiltonian is not Hermitian")
    low = np.asarray([
        position for position, state in enumerate(states)
        if all(not ((int(state) >> species_mode(length, link, family, slot)) & 1)
               for link in range(length - 1) for family in (0, 1) for slot in (0, 1))
    ], dtype=int)
    frame = np.zeros((len(states), len(low)), dtype=complex)
    frame[low, np.arange(len(low))] = 1.0
    return hamiltonian, states, frame


def bare_low_parity(states: np.ndarray, frame: np.ndarray) -> np.ndarray:
    low_states = [int(states[int(np.argmax(np.abs(frame[:, column])))]) for column in range(frame.shape[1])]
    return np.diag([
        -1.0 if sum((state >> site_index(rung, 0)) & 1 for rung in range(LENGTH)) % 2 else 1.0
        for state in low_states
    ])


def full_separate_parity(states: np.ndarray) -> np.ndarray:
    """Q with fixed per-channel charges: H0=(even, even), H1=(even, odd)."""
    diagonal = []
    for raw_state in states:
        state = int(raw_state)
        exponent = sum((state >> site_index(rung, 0)) & 1 for rung in range(LENGTH))
        # H1 slot 0 is (aa-bb)/sqrt(2), hence even under rail parity;
        # H1 slot 1 is (ab+ba)/sqrt(2), hence odd.  The H0 slots are both
        # even.  This is the fixed charge assignment unavailable when the
        # physical mediator slots are reused across H0 and H1.
        exponent += sum(
            (state >> species_mode(LENGTH, link, 1, 1)) & 1
            for link in range(LENGTH - 1)
        )
        diagonal.append(-1.0 if exponent % 2 else 1.0)
    return np.diag(diagonal)


def low_metrics(unitary: np.ndarray, frame: np.ndarray, parity: np.ndarray) -> dict:
    raw = frame.conj().T @ unitary @ frame
    logical = polar_unitary(raw)
    even = np.flatnonzero(np.diag(parity) > 0)
    odd = np.flatnonzero(np.diag(parity) < 0)
    return {
        "raw_even_to_odd_operator_norm": float(np.linalg.norm(raw[np.ix_(even, odd)], ord=2)),
        "raw_odd_to_even_operator_norm": float(np.linalg.norm(raw[np.ix_(odd, even)], ord=2)),
        "polar_parity_commutator_normalized": float(np.linalg.norm(logical @ parity - parity @ logical) / np.sqrt(parity.shape[0])),
        "raw_min_singular_value": float(np.min(np.linalg.svd(raw, compute_uv=False))),
    }


def loglog_slope(rows: list[dict], field: str) -> float | None:
    selected = [(row["g_over_detuning"], row[field]) for row in rows if row[field] > 1e-13]
    if len(selected) < 3:
        return None
    x = np.log([point[0] for point in selected])
    y = np.log([point[1] for point in selected])
    return float(np.polyfit(x, y, 1)[0])


def deep_sw_multi_cycle_slope(rows: list[dict], cycles: int) -> float:
    selected = []
    for row in rows:
        if row["g_over_detuning"] > 0.075:
            continue
        metrics = next(item["reused_mediator"] for item in row["multi_cycle_metrics"] if item["cycles"] == cycles)
        selected.append((row["g_over_detuning"], metrics["polar_parity_commutator_normalized"]))
    return float(np.polyfit(np.log([x for x, _ in selected]), np.log([y for _, y in selected]), 1)[0])


def main() -> None:
    h0_channels = channel_path(0.0)
    h1_channels = channel_path(np.pi / 2.0)
    rows = []
    for ratio in RATIOS:
        detuning = abs(TARGET_U0) / ratio**2
        g = ratio * detuning
        omega = float(np.sqrt(detuning**2 + 4.0 * g**2))
        period = 4.0 * np.pi / omega

        h0_reused, states_reused, frame_reused = build_micro(LENGTH, PARTICLE_NUMBER, h0_channels, g, detuning)
        h1_reused, states_h1, frame_h1 = build_micro(LENGTH, PARTICLE_NUMBER, h1_channels, g, detuning)
        if not np.array_equal(states_reused, states_h1) or not np.allclose(frame_reused, frame_h1):
            raise RuntimeError("reused segment frames disagree")
        reused_cycle = expm(-1j * (1.0 - ALPHA) * period * h1_reused) @ expm(-1j * ALPHA * period * h0_reused)

        h0_sep, states_sep, frame_sep = build_separate_species(LENGTH, PARTICLE_NUMBER, h0_channels, 0, g, detuning)
        h1_sep, states_h1_sep, frame_h1_sep = build_separate_species(LENGTH, PARTICLE_NUMBER, h1_channels, 1, g, detuning)
        if not np.array_equal(states_sep, states_h1_sep) or not np.allclose(frame_sep, frame_h1_sep):
            raise RuntimeError("separate-species segment frames disagree")
        q_sep = full_separate_parity(states_sep)
        denominator = max(np.linalg.norm(h0_sep, ord="fro"), np.linalg.norm(h1_sep, ord="fro"))
        exact_q_residual = max(
            float(np.linalg.norm(h0_sep @ q_sep - q_sep @ h0_sep, ord="fro") / denominator),
            float(np.linalg.norm(h1_sep @ q_sep - q_sep @ h1_sep, ord="fro") / denominator),
        )
        separate_cycle = expm(-1j * (1.0 - ALPHA) * period * h1_sep) @ expm(-1j * ALPHA * period * h0_sep)

        pa_reused = bare_low_parity(states_reused, frame_reused)
        pa_sep = bare_low_parity(states_sep, frame_sep)
        if not np.array_equal(pa_reused, pa_sep):
            raise RuntimeError("physical low-space parity differs between primary and control")
        multi_cycle_metrics = []
        for cycles in CYCLE_COUNTS:
            multi_cycle_metrics.append({
                "cycles": cycles,
                "reused_mediator": low_metrics(np.linalg.matrix_power(reused_cycle, cycles), frame_reused, pa_reused),
                "separate_species_control": low_metrics(np.linalg.matrix_power(separate_cycle, cycles), frame_sep, pa_sep),
            })
        rows.append({
            "g_over_detuning": ratio,
            "detuning": detuning,
            "g": g,
            "virtual_rabi_closure_period": period,
            "reused_mediator": low_metrics(reused_cycle, frame_reused, pa_reused),
            "separate_species_control": {
                "full_microscopic_parity_commutator_normalized": exact_q_residual,
                **low_metrics(separate_cycle, frame_sep, pa_sep),
            },
            "multi_cycle_metrics": multi_cycle_metrics,
        })

    reused_flat = [dict(g_over_detuning=row["g_over_detuning"], **row["reused_mediator"]) for row in rows]
    out = {
        "schema": "antler.phase8.direct-channel-boundary-parity-scaling.v2",
        "parameters": {
            "L": LENGTH, "N": PARTICLE_NUMBER, "target_u0": TARGET_U0, "alpha": ALPHA,
            "ratios": list(RATIOS),
            "cycle_counts": list(CYCLE_COUNTS),
            "h0_channels": "(aa-bb)/sqrt(2), (aa+bb)/sqrt(2)",
            "h1_channels": "(aa-bb)/sqrt(2), (ab+ba)/sqrt(2)",
            "cycle": "exp[-i(1-alpha) T H1] exp[-i alpha T H0], T=4pi/sqrt(Delta^2+4g^2)",
        },
        "rows": rows,
        "reused_mediator_loglog_slopes": {
            "raw_even_to_odd_operator_norm": loglog_slope(reused_flat, "raw_even_to_odd_operator_norm"),
            "raw_odd_to_even_operator_norm": loglog_slope(reused_flat, "raw_odd_to_even_operator_norm"),
            "polar_parity_commutator_normalized": loglog_slope(reused_flat, "polar_parity_commutator_normalized"),
        },
        "deep_sw_multi_cycle_polar_parity_slopes": {
            str(cycles): deep_sw_multi_cycle_slope(rows, cycles) for cycles in CYCLE_COUNTS
        },
        "decision_rule": (
            "The claimed boundary mechanism is supported only if the reused-mediator low-space parity-flip metric has "
            "a stable deep-SW slope compatible with 2 while the separate-species control has exact microscopic and "
            "low-space parity residuals at numerical precision. Otherwise its stated order or mechanism is falsified."
        ),
        "claim_boundary": (
            "This is an L=3,N=2 ideal-pulse microscopic audit of one switched-channel block. It does not establish a "
            "T-junction, a thermodynamic phase, a protected fusion space, a braid, non-Abelian statistics, universality "
            "or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "direct_channel_boundary_parity_scaling.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
