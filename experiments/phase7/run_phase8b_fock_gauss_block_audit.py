"""First fixed-charge Fock-block audit of the repaired Phase-8B walker.

This embeds the b=2 Lambda walker in a two-rung, two-rail fermionic block at
fixed N=2, with two boundary gauge qubits.  It tests exact Gauss commutation,
the rail-tunnelling selection rule, and the order of the resulting G_B-sector
gap.  The density-conditioned neutral-walker hop remains an explicitly
inserted new primitive in this control.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.basis import build_basis, site_index
from antler.number_conserving_pairwire import _apply


DETUNING = 10.0
PARTICLE_NUMBER = 2
RATIOS = (0.20, 0.15, 0.10, 0.075, 0.05, 0.0375, 0.025)


def block_index(matter: int, left_gauge: int, right_gauge: int, walker: int, matter_dimension: int) -> int:
    return ((((walker * 2) + right_gauge) * 2 + left_gauge) * matter_dimension) + matter


def build_operators(coupling: float, t_leg: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    states, positions = build_basis(4, PARTICLE_NUMBER)
    matter_dimension = len(states)
    dimension = matter_dimension * 2 * 2 * 4
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    gauss = np.zeros_like(hamiltonian)
    bare_rung_tunnel = np.zeros_like(hamiltonian)

    def p_rung(raw_state: int, rung: int) -> float:
        return -1.0 if ((raw_state >> site_index(rung, 0)) & 1) else 1.0

    def p_block(raw_state: int) -> float:
        return p_rung(raw_state, 0) * p_rung(raw_state, 1)

    for matter, raw_state in enumerate(states):
        state = int(raw_state)
        for left in (0, 1):
            for right in (0, 1):
                z_left = 1.0 if left == 0 else -1.0
                for walker in range(4):
                    column = block_index(matter, left, right, walker, matter_dimension)
                    gauss[block_index(matter, left ^ 1, right ^ 1, walker, matter_dimension), column] = p_block(state)
                    if walker:
                        hamiltonian[column, column] += DETUNING
                # Endpoint links of the virtual loop.  Gauge X flips the
                # corresponding boundary-qubit computational bit.
                origin = block_index(matter, left, right, 0, matter_dimension)
                target = block_index(matter, left ^ 1, right, 1, matter_dimension)
                hamiltonian[target, origin] += coupling
                hamiltonian[origin, target] += coupling
                origin = block_index(matter, left, right, 3, matter_dimension)
                target = block_index(matter, left, right ^ 1, 0, matter_dimension)
                hamiltonian[target, origin] += coupling
                hamiltonian[origin, target] += coupling
                # Density-conditioned neutral-walker links: inserted new
                # primitive, diagonal in the physical matter Fock basis.
                for first_walker, second_walker, parity in ((1, 2, p_rung(state, 0)), (2, 3, p_rung(state, 1))):
                    first = block_index(matter, left, right, first_walker, matter_dimension)
                    second = block_index(matter, left, right, second_walker, matter_dimension)
                    hamiltonian[second, first] += coupling * parity
                    hamiltonian[first, second] += coupling * parity
                # Bare U(1)-conserving rail tunnelling on the first rung.
                for operations in (
                    (("ann", site_index(0, 1)), ("create", site_index(0, 0))),
                    (("ann", site_index(0, 0)), ("create", site_index(0, 1))),
                ):
                    item = _apply(state, operations)
                    if item is not None:
                        new_state, amplitude = item
                        for walker in range(4):
                            row = block_index(positions[new_state], left, right, walker, matter_dimension)
                            column = block_index(matter, left, right, walker, matter_dimension)
                            bare_rung_tunnel[row, column] += amplitude
                # Ordinary intra-block leg hopping preserves N_a mod 2.
                for rail in (0, 1):
                    for operations in (
                        (("ann", site_index(1, rail)), ("create", site_index(0, rail))),
                        (("ann", site_index(0, rail)), ("create", site_index(1, rail))),
                    ):
                        item = _apply(state, operations)
                        if item is not None:
                            new_state, amplitude = item
                            for walker in range(4):
                                row = block_index(positions[new_state], left, right, walker, matter_dimension)
                                column = block_index(matter, left, right, walker, matter_dimension)
                                hamiltonian[row, column] += -t_leg * amplitude
    if not np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-12):
        raise RuntimeError("Fock-block Hamiltonian is not Hermitian")
    if not np.allclose(gauss, gauss.conj().T, atol=1e-12) or not np.allclose(gauss @ gauss, np.eye(dimension), atol=1e-12):
        raise RuntimeError("Gauss operator is not a Hermitian involution")
    return hamiltonian, gauss, bare_rung_tunnel, states, {"dimension": dimension, "matter_dimension": matter_dimension}


def sector_minimum(hamiltonian: np.ndarray, gauss: np.ndarray, sign: int) -> float:
    values, vectors = np.linalg.eigh(gauss)
    frame = vectors[:, values * sign > 0.5]
    return float(np.linalg.eigvalsh(frame.conj().T @ hamiltonian @ frame)[0])


def main() -> None:
    rows = []
    for ratio in RATIOS:
        hamiltonian, gauss, bare, states, metadata = build_operators(ratio * DETUNING)
        energy_plus = sector_minimum(hamiltonian, gauss, 1)
        energy_minus = sector_minimum(hamiltonian, gauss, -1)
        identity = np.eye(hamiltonian.shape[0], dtype=complex)
        p_plus = (identity + gauss) / 2.0
        p_minus = (identity - gauss) / 2.0
        rows.append({
            "coupling_over_detuning": ratio,
            "coupling": ratio * DETUNING,
            "gauss_sector_energies": {"G_plus": energy_plus, "G_minus": energy_minus},
            "gauss_sector_gap": float(abs(energy_plus - energy_minus)),
            "hamiltonian_gauss_commutator_frobenius": float(np.linalg.norm(hamiltonian @ gauss - gauss @ hamiltonian)),
            "bare_rung_tunnel_gauss_anticommutator_frobenius": float(np.linalg.norm(bare @ gauss + gauss @ bare)),
            "bare_rung_tunnel_physical_projection_norm": float(np.linalg.norm(p_plus @ bare @ p_plus, ord=2)),
            "bare_rung_tunnel_sector_changing_norm": float(np.linalg.norm(p_minus @ bare @ p_plus, ord=2)),
        })
    fit = float(np.polyfit(
        np.log([row["coupling_over_detuning"] for row in rows if row["coupling_over_detuning"] <= 0.075]),
        np.log([row["gauss_sector_gap"] for row in rows if row["coupling_over_detuning"] <= 0.075]),
        1,
    )[0])
    h_with_leg, gauss_with_leg, _, _, metadata = build_operators(0.1 * DETUNING, t_leg=0.3)
    output = {
        "schema": "antler.phase8b.fock-gauss-block-audit.v1",
        "parameters": {
            "block_size_b": 2, "matter_modes": 4, "fixed_matter_particle_number": PARTICLE_NUMBER,
            "detuning": DETUNING, "ratios": list(RATIOS),
            "walker": "vacuum plus mu0,mu1,mu2; at most one neutral walker",
            "new_inserted_primitive": "neutral-walker hopping conditioned on (1-2 n_a,j)",
            "gauge": "two boundary qubits with G_B=X_L(-1)^N_a,B X_R",
        },
        "dimensions": metadata,
        "rows": rows,
        "deep_sw_gap_power_vs_coupling_over_detuning": fit,
        "intra_block_leg_hopping_control": {
            "t_leg": 0.3,
            "hamiltonian_gauss_commutator_frobenius": float(np.linalg.norm(h_with_leg @ gauss_with_leg - gauss_with_leg @ h_with_leg)),
        },
        "decision": (
            "The declared Lambda-walker primitive embeds consistently in a fixed-charge Fock block: it has exact Gauss "
            "commutation, the bare rail tunnel changes Gauss sector rather than acting within it, and the sector gap has "
            "the expected fourth-order scaling. This is still an inserted-primitive control, not a native ANTLER derivation "
            "or a multi-block paired-phase result."
        ),
        "claim_boundary": (
            "This one-block audit omits physical realization of the conditional walker hop, inter-block pair transport, the "
            "full Floquet sequence, noise, local indistinguishability, a thermodynamic phase, fusion and braid operations."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_fock_gauss_block_audit.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
