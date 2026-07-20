"""CUDA exact-eigenpropagation completion of the finite Peierls-ramp audit.

The CPU Krylov implementation is retained as an incomplete stress-test.  This
script diagonalizes each repeated Hermitian segment once on CUDA and applies
its exact spectral propagator to all 16 monomer-frame columns.  It therefore
tests the same 472-state Hamiltonian without reducing the ramp physics.
"""
from __future__ import annotations

import json
import os
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
    FRAME, L, LINKS, LOGICAL_PA, LOGICAL_PB, PULSE_TIME, pair_gate,
    projected, pulse_hamiltonian, rail_rotation, remove_global_phase,
)
from run_phase7d_peierls_phase_ramp_audit import LEG_PI, LEG_ZERO, leg_with_phase


if not torch.cuda.is_available():
    raise RuntimeError("CUDA is required for this registered GPU completion audit")


DEVICE, DTYPE = torch.device("cuda:0"), torch.complex128
REAL_DTYPE = torch.float64
LEG_HOPPING, EPSILON, SUBCYCLES = 1.0, 1e-2, 16
RAMP_STEPS = int(os.environ.get("ANTLER_RAMP_STEPS", "2"))
RAMP_FRACTIONS = tuple(float(item) for item in os.environ.get("ANTLER_RAMP_FRACTIONS", "0,0.005,0.02").split(","))
RAMP_PHASES = tuple(np.pi * (step + 0.5) / RAMP_STEPS for step in range(RAMP_STEPS))
OUTPUT_NAME = "peierls_phase_ramp_gpu_audit.json" if RAMP_STEPS == 2 else f"peierls_phase_ramp_gpu_{RAMP_STEPS}step_refinement.json"


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


def pair_eigens(active_links: tuple[int, ...], kind: str) -> dict[str, tuple[torch.Tensor, torch.Tensor]]:
    h_pair = pulse_hamiltonian(active_links, kind, EPSILON, 0.0)
    output = {
        "plus": eigensystem(h_pair + LEG_ZERO),
        "minus": eigensystem(h_pair + LEG_PI),
    }
    for step, phase in enumerate(RAMP_PHASES):
        output[f"ramp_{step}"] = eigensystem(h_pair + leg_with_phase(phase))
    return output


PAIR_EIGENS = {(active, kind): pair_eigens(active, kind) for active in ((0, 2), (1,)) for kind in ("same", "opposite")}
ROTATION_EIGENS = {
    rungs: {axis: eigensystem(rail_rotation(axis, rungs)) for axis in ("x", "y")}
    for rungs in ((0, 1, 2, 3), (1, 2))
}


def ramped_pulse(vectors: torch.Tensor, active_links: tuple[int, ...], kind: str, fraction: float) -> torch.Tensor:
    eig = PAIR_EIGENS[(active_links, kind)]
    delta = PULSE_TIME / SUBCYCLES
    plus_duration = delta * (0.25 - fraction / 2.0)
    minus_duration = delta * (0.50 - fraction)
    ramp_step_duration = delta * fraction / RAMP_STEPS
    for _ in range(SUBCYCLES):
        vectors = evolve(vectors, eig["plus"], plus_duration)
        if ramp_step_duration:
            for step in range(RAMP_STEPS):
                vectors = evolve(vectors, eig[f"ramp_{step}"], ramp_step_duration)
        vectors = evolve(vectors, eig["minus"], minus_duration)
        if ramp_step_duration:
            for step in reversed(range(RAMP_STEPS)):
                vectors = evolve(vectors, eig[f"ramp_{step}"], ramp_step_duration)
        vectors = evolve(vectors, eig["plus"], plus_duration)
    return vectors


def ramped_pair_gate(vectors: torch.Tensor, active_links: tuple[int, ...], fraction: float) -> torch.Tensor:
    rungs = (0, 1, 2, 3) if active_links == (0, 2) else (1, 2)
    rotation = ROTATION_EIGENS[rungs]
    vectors = evolve(vectors, rotation["y"], np.pi / 4.0)
    vectors = ramped_pulse(vectors, active_links, "same", fraction)
    vectors = evolve(vectors, rotation["y"], -np.pi / 4.0)
    vectors = evolve(vectors, rotation["x"], np.pi / 4.0)
    vectors = ramped_pulse(vectors, active_links, "opposite", fraction)
    return evolve(vectors, rotation["x"], -np.pi / 4.0)


def metrics(vectors: np.ndarray, reference: np.ndarray) -> dict:
    logical = projected(vectors)
    leakage = float(np.linalg.norm(vectors - FRAME @ logical, ord=2) ** 2)
    return {
        "monomer_leakage": leakage,
        "logical_deviation_from_zero_leg_schedule": float(np.linalg.norm(logical - remove_global_phase(reference, logical), ord=2)),
        "logical_parity_a_residual": float(np.linalg.norm(logical @ LOGICAL_PA - LOGICAL_PA @ logical, ord=2)),
        "logical_parity_b_residual": float(np.linalg.norm(logical @ LOGICAL_PB - LOGICAL_PB @ logical, ord=2)),
        "logical_singular_value_min": float(np.linalg.svd(logical, compute_uv=False)[-1]),
    }


def main() -> None:
    torch.cuda.synchronize()
    started = time.perf_counter()
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    frame_gpu = to_gpu(FRAME)
    rows = []
    for fraction in RAMP_FRACTIONS:
        even = ramped_pair_gate(frame_gpu.clone(), (0, 2), fraction)
        complete = ramped_pair_gate(even, (1,), fraction)
        torch.cuda.synchronize()
        rows.append({
            "ramp_fraction_of_each_subcycle": fraction,
            "ramp_duration": fraction * PULSE_TIME / SUBCYCLES,
            "ramp_steps": RAMP_STEPS,
            "total_phase_ramp_time_per_mediator_pulse": 2.0 * fraction * PULSE_TIME,
            **metrics(complete.detach().cpu().numpy(), reference),
        })
    torch.cuda.synchronize()
    out = {
        "schema": "antler.phase7d.peierls-phase-ramp-gpu-audit.v1",
        "method": "exact segment spectral propagation on CUDA, complex128",
        "device": torch.cuda.get_device_name(0),
        "parameters": {"leg_hopping": LEG_HOPPING, "inactive_channel_coupling_over_g": EPSILON, "subcycles": SUBCYCLES, "ramp_fraction_scan": list(RAMP_FRACTIONS)},
        "rows": rows,
        "wall_clock_seconds_including_cached_eigensystems": time.perf_counter() - started,
        "decision": "CUDA completion of the finite Peierls-ramp audit; numerical temporal convergence must be checked before promotion.",
        "claim_boundary": "This is a finite-block pulse control only; no hardware bandwidth, protected phase, 2D code, braid, non-Abelianity, universality or fault tolerance is established.",
    }
    path = ROOT / "results" / "phase7" / OUTPUT_NAME
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
