"""Exact CUDA audit of systematic errors in the direct Peierls-sign echo.

The prior ramp audit establishes convergence of the numerical representation of
a finite phase waveform.  This independent test asks a different, physical
control question in the same 472-state microscopic block: how far can the
pi-phase plateau and the signed-time balance be miscalibrated before the local
four-pulse primitive misses its registered control targets?

It intentionally tests deterministic one-parameter offsets only.  It does not
model bandwidth, stochastic phase noise, correlations, calibration protocols,
or a many-body protected phase.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import numpy as np
import torch


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import (
    FRAME, LOGICAL_PA, LOGICAL_PB, PULSE_TIME, pair_gate, projected,
    pulse_hamiltonian, rail_rotation, remove_global_phase,
)
from run_phase7d_peierls_phase_ramp_audit import LEG_ZERO, leg_with_phase


if not torch.cuda.is_available():
    raise RuntimeError("CUDA is required for this registered switching-error audit")


DEVICE, DTYPE, REAL_DTYPE = torch.device("cuda:0"), torch.complex128, torch.float64
EPSILON, SUBCYCLES = 1e-2, 16
PHASE_ERROR_DEGREES = (0.0, 0.5, 1.0, 2.0, 5.0, 7.5, 10.0, 15.0)
TIME_IMBALANCE = (0.0, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 1e-1, 1.5e-1)
LOCAL_TARGET = {
    "max_monomer_leakage": 1e-4,
    "max_logical_parity_residual": 1e-4,
    "max_logical_deviation": 1e-4,
    "min_logical_singular_value": 0.9999,
}


def to_gpu(matrix: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(matrix, dtype=DTYPE, device=DEVICE)


def eigensystem(matrix) -> tuple[torch.Tensor, torch.Tensor]:
    dense = to_gpu(matrix.toarray() if hasattr(matrix, "toarray") else matrix)
    values, vectors = torch.linalg.eigh(dense)
    return values.to(REAL_DTYPE), vectors


def evolve(vectors: torch.Tensor, eig: tuple[torch.Tensor, torch.Tensor], duration: float) -> torch.Tensor:
    values, basis = eig
    coefficients = basis.mH @ vectors
    return basis @ (torch.exp((-1j * duration) * values).unsqueeze(1) * coefficients)


ROTATION_EIGENS = {
    rungs: {axis: eigensystem(rail_rotation(axis, rungs)) for axis in ("x", "y")}
    for rungs in ((0, 1, 2, 3), (1, 2))
}


def pulse_eigens(phase_error_rad: float) -> dict[tuple[tuple[int, ...], str, str], tuple[torch.Tensor, torch.Tensor]]:
    output = {}
    for active_links in ((0, 2), (1,)):
        for kind in ("same", "opposite"):
            h_pair = pulse_hamiltonian(active_links, kind, EPSILON, 0.0)
            output[(active_links, kind, "plus")] = eigensystem(h_pair + LEG_ZERO)
            output[(active_links, kind, "minus")] = eigensystem(h_pair + leg_with_phase(np.pi + phase_error_rad))
    return output


def direct_echo_pulse(
    vectors: torch.Tensor,
    active_links: tuple[int, ...],
    kind: str,
    eigens: dict[tuple[tuple[int, ...], str, str], tuple[torch.Tensor, torch.Tensor]],
    imbalance: float,
) -> torch.Tensor:
    """Keep each subcycle duration fixed while biasing its + and - plateaux."""
    delta = PULSE_TIME / SUBCYCLES
    plus_duration = delta * 0.25 * (1.0 + imbalance)
    minus_duration = delta * 0.50 * (1.0 - imbalance)
    for _ in range(SUBCYCLES):
        vectors = evolve(vectors, eigens[(active_links, kind, "plus")], plus_duration)
        vectors = evolve(vectors, eigens[(active_links, kind, "minus")], minus_duration)
        vectors = evolve(vectors, eigens[(active_links, kind, "plus")], plus_duration)
    return vectors


def direct_pair_gate(
    vectors: torch.Tensor,
    active_links: tuple[int, ...],
    eigens: dict[tuple[tuple[int, ...], str, str], tuple[torch.Tensor, torch.Tensor]],
    imbalance: float,
) -> torch.Tensor:
    rungs = (0, 1, 2, 3) if active_links == (0, 2) else (1, 2)
    rotation = ROTATION_EIGENS[rungs]
    vectors = evolve(vectors, rotation["y"], np.pi / 4.0)
    vectors = direct_echo_pulse(vectors, active_links, "same", eigens, imbalance)
    vectors = evolve(vectors, rotation["y"], -np.pi / 4.0)
    vectors = evolve(vectors, rotation["x"], np.pi / 4.0)
    vectors = direct_echo_pulse(vectors, active_links, "opposite", eigens, imbalance)
    return evolve(vectors, rotation["x"], -np.pi / 4.0)


def metrics(vectors: np.ndarray, reference: np.ndarray) -> dict[str, float | bool]:
    logical = projected(vectors)
    leakage = float(np.linalg.norm(vectors - FRAME @ logical, ord=2) ** 2)
    parity_a = float(np.linalg.norm(logical @ LOGICAL_PA - LOGICAL_PA @ logical, ord=2))
    parity_b = float(np.linalg.norm(logical @ LOGICAL_PB - LOGICAL_PB @ logical, ord=2))
    deviation = float(np.linalg.norm(logical - remove_global_phase(reference, logical), ord=2))
    singular_min = float(np.linalg.svd(logical, compute_uv=False)[-1])
    passed = (
        leakage <= LOCAL_TARGET["max_monomer_leakage"]
        and max(parity_a, parity_b) <= LOCAL_TARGET["max_logical_parity_residual"]
        and deviation <= LOCAL_TARGET["max_logical_deviation"]
        and singular_min >= LOCAL_TARGET["min_logical_singular_value"]
    )
    return {
        "monomer_leakage": leakage,
        "logical_deviation_from_zero_leg_schedule": deviation,
        "logical_parity_a_residual": parity_a,
        "logical_parity_b_residual": parity_b,
        "logical_singular_value_min": singular_min,
        "passes_registered_local_target": passed,
    }


def evaluate(
    frame: torch.Tensor,
    reference: np.ndarray,
    phase_error_rad: float,
    imbalance: float,
) -> dict[str, float | bool]:
    eigens = pulse_eigens(phase_error_rad)
    even = direct_pair_gate(frame.clone(), (0, 2), eigens, imbalance)
    complete = direct_pair_gate(even, (1,), eigens, imbalance)
    torch.cuda.synchronize()
    return metrics(complete.detach().cpu().numpy(), reference)


def tolerance_bracket(rows: list[dict], key: str) -> dict[str, float | None]:
    passing = [float(row[key]) for row in rows if row["passes_registered_local_target"]]
    failing = [float(row[key]) for row in rows if not row["passes_registered_local_target"]]
    return {
        "largest_tested_passing": max(passing) if passing else None,
        "smallest_tested_failing": min(failing) if failing else None,
    }


def main() -> None:
    torch.cuda.synchronize()
    started = time.perf_counter()
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    frame = to_gpu(FRAME)

    phase_rows = []
    for degrees in PHASE_ERROR_DEGREES:
        phase_rows.append({
            "systematic_pi_plateau_error_degrees": degrees,
            "systematic_pi_plateau_error_radians": float(np.deg2rad(degrees)),
            **evaluate(frame, reference, float(np.deg2rad(degrees)), 0.0),
        })

    nominal_eigens = pulse_eigens(0.0)
    timing_rows = []
    for imbalance in TIME_IMBALANCE:
        even = direct_pair_gate(frame.clone(), (0, 2), nominal_eigens, imbalance)
        complete = direct_pair_gate(even, (1,), nominal_eigens, imbalance)
        torch.cuda.synchronize()
        timing_rows.append({
            "signed_plateau_time_imbalance_fraction": imbalance,
            "plus_plateau_relative_error_percent": 100.0 * imbalance,
            "minus_plateau_relative_error_percent": -100.0 * imbalance,
            **metrics(complete.detach().cpu().numpy(), reference),
        })

    phase_bracket = tolerance_bracket(phase_rows, "systematic_pi_plateau_error_degrees")
    time_bracket = tolerance_bracket(timing_rows, "plus_plateau_relative_error_percent")
    torch.cuda.synchronize()
    out = {
        "schema": "antler.phase7d.peierls-switch-error-gpu-audit.v1",
        "method": "exact segment spectral propagation on CUDA, complex128",
        "device": torch.cuda.get_device_name(0),
        "parameters": {"inactive_channel_coupling_over_g": EPSILON, "subcycles": SUBCYCLES},
        "registered_local_target": LOCAL_TARGET,
        "systematic_phase_plateau_scan": phase_rows,
        "systematic_signed_time_imbalance_scan": timing_rows,
        "registered_tolerance_brackets": {
            "pi_plateau_error_degrees": phase_bracket,
            "signed_time_imbalance_percent": time_bracket,
        },
        "wall_clock_seconds_including_eigensystems": time.perf_counter() - started,
        "decision": (
            "For the registered local target, the tested deterministic plateau-phase tolerance passes through "
            f"{phase_bracket['largest_tested_passing']} degrees and first fails at "
            f"{phase_bracket['smallest_tested_failing']} degrees; the signed-time imbalance passes through "
            f"{time_bracket['largest_tested_passing']} percent and first fails at "
            f"{time_bracket['smallest_tested_failing']} percent. These are scan brackets, not interpolated specifications."
        ),
        "claim_boundary": (
            "This does not test bandwidth, stochastic or correlated noise, calibration, full-ladder dynamics, a protected phase, "
            "edge modes, a 2D code, braid, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "peierls_switch_error_gpu_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
