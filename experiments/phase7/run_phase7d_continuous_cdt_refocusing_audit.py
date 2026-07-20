"""Continuous coherent-destruction-of-tunnelling (CDT) control for Phase 7D.

Unlike the parity-only echo generator, K_cont assigns a charge-two mediator on
link j the potential weight 2j+1.  Hence every nearest-neighbour pair
conversion exactly commutes with K_cont, while each intraleg hopping changes
K_cont by one.  A sinusoidal onsite drive at the first J0 zero suppresses the
hopping in the high-frequency average without suppressing pair conversion.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.special import jn_zeros, jv
from scipy.sparse import diags
from scipy.sparse.linalg import norm as sparse_norm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase7d_four_rung_microscopic_pulse_audit import (
    CHANNELS, FRAME, L, LINKS, LOGICAL_PA, LOGICAL_PB, LOW_MODES, PULSE_TIME,
    STATES, evolve, leakage, pair_gate, projected, pulse_hamiltonian,
    rail_rotation, remove_global_phase,
)


LEG_HOPPING, EPSILON = 1.0, 1e-2
CYCLES = (1, 2, 4)
STEPS_PER_CYCLE = (16, 32)
XI = float(jn_zeros(0, 1)[0])


def continuous_generator_values() -> np.ndarray:
    values = []
    for raw_state in STATES:
        state = int(raw_state)
        rail_weight = sum(
            rung * (((state >> (2 * rung)) & 1) + ((state >> (2 * rung + 1)) & 1))
            for rung in range(L)
        )
        mediator_weight = sum(
            (2 * link + 1) * ((state >> (LOW_MODES + mediator)) & 1)
            for mediator, (link, _kind) in enumerate(CHANNELS)
        )
        values.append(rail_weight + mediator_weight)
    return np.asarray(values, dtype=float)


K_VALUES = continuous_generator_values()


def continuous_pulse(vectors: np.ndarray, active_links: tuple[int, ...], kind_prefix: str,
                     cycles: int, steps_per_cycle: int) -> np.ndarray:
    h_total = pulse_hamiltonian(active_links, kind_prefix, EPSILON, LEG_HOPPING)
    h_pair = pulse_hamiltonian(active_links, kind_prefix, EPSILON, 0.0)
    h_leg = h_total - h_pair
    total_steps = cycles * steps_per_cycle
    dt = PULSE_TIME / total_steps
    u_pair_half = expm((-0.5j * dt) * h_pair.toarray())
    u_leg = expm((-1j * dt) * h_leg.toarray())
    omega = 2.0 * np.pi * cycles / PULSE_TIME
    for step in range(total_steps):
        midpoint = (step + 0.5) * dt
        phase = XI * np.sin(omega * midpoint)
        rotating = np.exp(-1j * phase * K_VALUES)
        vectors = u_pair_half @ vectors
        vectors = rotating[:, None] * vectors
        vectors = u_leg @ vectors
        vectors = rotating.conj()[:, None] * vectors
        vectors = u_pair_half @ vectors
    return vectors


def continuous_pair_gate(vectors: np.ndarray, active_links: tuple[int, ...], cycles: int,
                         steps_per_cycle: int) -> np.ndarray:
    rungs = tuple(sorted({rung for link in active_links for rung in LINKS[link]}))
    ry, rx = rail_rotation("y", rungs), rail_rotation("x", rungs)
    vectors = evolve(vectors, ry, np.pi / 4.0)
    vectors = continuous_pulse(vectors, active_links, "same", cycles, steps_per_cycle)
    vectors = evolve(vectors, ry, -np.pi / 4.0)
    vectors = evolve(vectors, rx, np.pi / 4.0)
    vectors = continuous_pulse(vectors, active_links, "opposite", cycles, steps_per_cycle)
    return evolve(vectors, rx, -np.pi / 4.0)


def metrics(vectors: np.ndarray, reference: np.ndarray) -> dict:
    logical = projected(vectors)
    return {
        "monomer_leakage": leakage(vectors),
        "logical_deviation_from_zero_leg_schedule": float(np.linalg.norm(
            logical - remove_global_phase(reference, logical), ord=2
        )),
        "logical_parity_a_residual": float(np.linalg.norm(logical @ LOGICAL_PA - LOGICAL_PA @ logical, ord=2)),
        "logical_parity_b_residual": float(np.linalg.norm(logical @ LOGICAL_PB - LOGICAL_PB @ logical, ord=2)),
        "logical_singular_value_min": float(np.linalg.svd(logical, compute_uv=False)[-1]),
    }


def main() -> None:
    reference = projected(pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0), (1,), 0.0))
    identity_checks = []
    k_matrix = diags(K_VALUES, format="csr")
    for kind in ("same", "opposite"):
        h_pair = pulse_hamiltonian((1,), kind, EPSILON, 0.0)
        h_leg = pulse_hamiltonian((1,), kind, EPSILON, LEG_HOPPING) - h_pair
        identity_checks.append({
            "channel_kind": kind,
            "pair_commutator_frobenius": float(sparse_norm(h_pair @ k_matrix - k_matrix @ h_pair)),
            "leg_commutator_frobenius": float(sparse_norm(h_leg @ k_matrix - k_matrix @ h_leg)),
        })
    bare = pair_gate(pair_gate(FRAME.copy(), (0, 2), 0.0, LEG_HOPPING), (1,), 0.0, LEG_HOPPING)
    rows = [{"method": "bare_no_cdt", **metrics(bare, reference)}]
    for cycles in CYCLES:
        for steps in STEPS_PER_CYCLE:
            even = continuous_pair_gate(FRAME.copy(), (0, 2), cycles, steps)
            complete = continuous_pair_gate(even, (1,), cycles, steps)
            omega = 2.0 * np.pi * cycles / PULSE_TIME
            rows.append({
                "method": "continuous_cdt",
                "cycles_per_rabi_pulse": cycles,
                "steps_per_cycle": steps,
                "drive_frequency": omega,
                "dimensionless_amplitude_xi": XI,
                "base_potential_amplitude": XI * omega,
                "max_rail_onsite_amplitude": (L - 1) * XI * omega,
                "J0_xi": float(jv(0, XI)),
                **metrics(complete, reference),
            })
    out = {
        "schema": "antler.phase7d.continuous-cdt-refocusing-audit.v1",
        "identity": (
            "K_cont=sum_j j(n_aj+n_bj)+sum_link(2j+1)N_mediator_link: nearest-link pair conversion commutes exactly, "
            "whereas intraleg hopping changes K_cont by one."
        ),
        "parameters": {
            "leg_hopping_during_pulses": LEG_HOPPING,
            "inactive_channel_coupling_over_g": EPSILON,
            "first_bessel_zero_xi": XI,
            "cycle_scan": list(CYCLES), "steps_per_cycle_scan": list(STEPS_PER_CYCLE),
        },
        "identity_checks": identity_checks,
        "rows": rows,
        "decision": (
            "The pair-preserving CDT identity is exact and the registered 1..4-cycle dynamics suppresses leg-hopping errors, but no "
            "tested point is promoted as a high-fidelity pulse. Hardware bandwidth, finite-frequency optimization, pulse-error and "
            "extended-ladder phase audits remain required."
        ),
        "claim_boundary": (
            "The link-dependent mediator potentials and sinusoidal modulation are added ideal controls. This is not a derivation of "
            "their experimental availability and establishes no protected phase, 2D code, braid, non-Abelian statistics, universality or fault tolerance."
        ),
    }
    path = ROOT / "results" / "phase7" / "continuous_cdt_refocusing_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
