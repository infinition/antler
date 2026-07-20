"""Convention audit for a proposed JW-string conditional-rung extension.

The frozen rung-major ladder has an unconditioned physical rung hop.  Here we
test a *declared graph extension*: a local rail bridge whose Jordan--Wigner
path crosses one charge-two walker mode.  At theta=pi its microscopic hopping
is

    -J (q_a^dag exp(i*pi*n_d1) q_b + h.c.),

which projects directly to `-J X_rail Z_walker` in the one-rail-one-walker
code.  Unlike the rejected shared-matter sign echo, the desired conditional
term is present before any Floquet averaging.  Residual ordinary ANTLER leg
hops are retained and a Mott-depth scan measures their leakage and Pauli
contamination.

The frozen model counts only unit-charge spatial rail occupations and has no
string on an adjacent rung hop.  A charge-two mediator is not part of that
string convention.  This script therefore contrasts three explicitly labelled
extensions: string weight zero (frozen rung rule), weight two (physical-charge
counting for a molecule), and the tempting but *new* odd weight one rule.  The
last one is retained only as an algebraic counterfactual: a positive gate
there cannot be promoted without deriving a new statistical-string resource.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[2]
PHASE7 = ROOT / "experiments" / "phase7"
if str(PHASE7) not in sys.path:
    sys.path.insert(0, str(PHASE7))

from run_phase8b_shared_matter_conditional_link_sw_audit import (
    CHARGES,
    DELTA_PAIR,
    D_0,
    D_1,
    LEG_HOPPING,
    PAULIS,
    Q_A,
    Q_B,
    R_A,
    R_B,
    RUNG_HOPPING,
    THETA,
    add_directed_hop,
    annihilate,
    code_indices,
    create,
    pauli_coefficients,
    schur_effective,
    weighted_basis,
)
from run_phase8b_shared_matter_pulse_closure_audit import polar_unitary, phase_aligned_distance


MOTT_DEPTHS = (15.0, 30.0, 60.0, 120.0)
TARGET_ANGLES = (0.1, 0.3, 0.6)
RELATIVE_ERROR_TARGET = 0.10
LEAKAGE_TARGET = 1e-4
STRING_WEIGHTS = (0, 1, 2)


def add_weighted_string_hop(
    hamiltonian: np.ndarray,
    states: np.ndarray,
    positions: dict[int, int],
    destination: int,
    source: int,
    amplitude: complex,
    string_weight: int,
) -> None:
    """A proposed generalized JW bridge; `string_weight=1` is non-frozen."""
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        item = annihilate(state, source)
        if item is None:
            continue
        intermediate, sign = item
        item = create(intermediate, destination)
        if item is None:
            continue
        final, final_sign = item
        phase = np.exp(1j * THETA * string_weight * ((state >> D_1) & 1))
        hamiltonian[positions[final], column] += amplitude * phase * sign * final_sign


def build_jw_crossing_bridge(
    mott_q: float,
    string_weight: int,
) -> tuple[np.ndarray, np.ndarray, dict[int, int]]:
    states, positions = weighted_basis()
    hamiltonian = np.zeros((len(states), len(states)), dtype=complex)
    for position, raw_state in enumerate(states):
        state = int(raw_state)
        n_q = ((state >> Q_A) & 1) + ((state >> Q_B) & 1)
        hamiltonian[position, position] += mott_q * (n_q - 1) ** 2
        hamiltonian[position, position] += DELTA_PAIR * ((state >> R_A) & 1) * ((state >> R_B) & 1)
    # Proposed bridge.  Weight 0 is the frozen adjacent-rung rule; weight 2
    # is the physical charge count of the molecular mediator; only weight 1
    # produces the apparent XZ gate and is a new odd-string convention.
    add_weighted_string_hop(hamiltonian, states, positions, Q_A, Q_B, -RUNG_HOPPING, string_weight)
    add_directed_hop(hamiltonian, states, positions, Q_A, R_A, -LEG_HOPPING, conditional_mode=Q_B)
    add_directed_hop(hamiltonian, states, positions, Q_B, R_B, -LEG_HOPPING, conditional_mode=R_A)
    hamiltonian = hamiltonian + hamiltonian.conj().T
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("JW-crossing bridge is not Hermitian")
    return hamiltonian, states, positions


def main() -> None:
    rows = []
    xz = np.kron(PAULIS["X"], PAULIS["Z"])
    static_rows = []
    convention_rows = []
    for mott_q in MOTT_DEPTHS:
        for string_weight in STRING_WEIGHTS:
            hamiltonian, states, positions = build_jw_crossing_bridge(mott_q, string_weight)
            low = code_indices(positions)
            effective, capture, gap = schur_effective(hamiltonian, low)
            coefficients = pauli_coefficients(effective)
            unwanted = {label: value for label, value in coefficients.items() if label not in {"II", "XZ"}}
            unwanted_label, unwanted_value = max(unwanted.items(), key=lambda item: abs(item[1]))
            target_coefficient = coefficients["XZ"]
            convention_rows.append({
                "mott_q": mott_q,
                "string_weight": string_weight,
                "effective_xz_coefficient": target_coefficient,
                "maximum_unwanted_non_scalar_coefficient": abs(unwanted_value),
                "largest_unwanted_non_scalar_pauli": unwanted_label,
                "low_frame_capture": capture,
                "low_high_gap": gap,
            })
            if string_weight != 1:
                continue
            static_rows.append({
                "mott_q": mott_q,
                "string_weight": string_weight,
                "effective_xz_coefficient": target_coefficient,
                "maximum_unwanted_non_scalar_coefficient": abs(unwanted_value),
                "largest_unwanted_non_scalar_pauli": unwanted_label,
                "unwanted_over_xz": float(abs(unwanted_value / target_coefficient)),
                "low_frame_capture": capture,
                "low_high_gap": gap,
            })
            counterfactual_hamiltonian = hamiltonian
            counterfactual_states = states
            counterfactual_low = low
            counterfactual_effective = effective
            counterfactual_target_coefficient = target_coefficient
        frame = np.zeros((len(counterfactual_states), len(counterfactual_low)), dtype=complex)
        frame[counterfactual_low, np.arange(len(counterfactual_low))] = 1.0
        projector = frame @ frame.conj().T
        for angle in TARGET_ANGLES:
            duration = float(angle / abs(counterfactual_target_coefficient))
            physical = expm(-1j * duration * counterfactual_hamiltonian)
            sw_target = expm(-1j * duration * counterfactual_effective)
            ideal = expm(-1j * duration * counterfactual_target_coefficient * xz)
            raw = frame.conj().T @ physical @ frame
            logical = polar_unitary(raw)
            sw_signal = phase_aligned_distance(sw_target, np.eye(4, dtype=complex))
            ideal_signal = phase_aligned_distance(ideal, np.eye(4, dtype=complex))
            physical_to_sw = phase_aligned_distance(logical, sw_target)
            physical_to_ideal = phase_aligned_distance(logical, ideal)
            sw_to_ideal = phase_aligned_distance(sw_target, ideal)
            leakage = float(np.linalg.norm((np.eye(len(counterfactual_states)) - projector) @ physical @ frame, ord=2) ** 2)
            rows.append({
                "mott_q": mott_q,
                "target_xz_angle": angle,
                "duration": duration,
                "time_ordered_sw_signal_distance_from_identity": sw_signal,
                "time_ordered_sw_to_ideal_xz_distance": sw_to_ideal,
                "physical_polar_signal_distance_from_identity": phase_aligned_distance(logical, np.eye(4, dtype=complex)),
                "relative_physical_to_sw_error": float(physical_to_sw / sw_signal) if sw_signal > 1e-15 else None,
                "relative_physical_to_ideal_xz_error": float(physical_to_ideal / ideal_signal) if ideal_signal > 1e-15 else None,
                "low_frame_leakage_worst": leakage,
                "passes_local_screen": bool(
                    physical_to_sw / sw_signal < RELATIVE_ERROR_TARGET
                    and physical_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                    and sw_to_ideal / ideal_signal < RELATIVE_ERROR_TARGET
                    and leakage < LEAKAGE_TARGET
                ),
            })
    output = {
        "schema": "antler.phase8b.jw-string-conditional-rung-bridge.v1",
        "parameters": {
            "mode_order": ["q_a", "q_b", "r_a", "r_b", "d0", "d1"],
            "weighted_charges": list(CHARGES),
            "theta": THETA,
            "proposed_graph_resource": "local q_a--q_b bridge whose generalized string crosses d1",
            "string_weight_conventions": {
                "0": "frozen adjacent-rung convention: no string",
                "2": "physical U(1)-charge counting for a charge-two mediator",
                "1": "algebraic counterfactual: odd hard-core mediator string, not derived",
            },
            "target": "X_rail Z_walker",
            "mott_depths": list(MOTT_DEPTHS),
            "target_angles": list(TARGET_ANGLES),
            "screen": {
                "relative_error": RELATIVE_ERROR_TARGET,
                "low_frame_leakage_worst": LEAKAGE_TARGET,
            },
        },
        "convention_rows": convention_rows,
        "counterfactual_weight_one_static_rows": static_rows,
        "rows": rows,
        "passing_rows": [row for row in rows if row["passes_local_screen"]],
        "decision": "Pending execution.",
        "claim_boundary": "The weight-one dynamical rows are an algebraic counterfactual, not a valid frozen-ANTLER result. They cannot derive a graph edge, a walker/code, protection, defects, fusion, non-Abelian braiding, universality or fault tolerance.",
    }
    output["decision"] = (
        "The odd string-weight-one counterfactual has local XZ gate rows, but the frozen (weight 0) and physical-charge (weight 2) conventions both have XZ=0. A new odd-string mediator resource is therefore necessary and not derived."
    )
    path = ROOT / "results" / "phase7" / "phase8b_jw_string_conditional_rung_bridge.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
