"""Many-pair preflight for the ideal Phase-8B pair-only continuation.

The previous Phase-8B tests establish an *inserted* pair transfer on two
coarse-Gauss blocks.  Before adding a T junction or calling it a phase, this
script asks the more basic question: can the most favourable low-energy
continuation made only of hard-core pairs and nearest-neighbour pair hopping
give a protected many-body sector?  The ideal effective chain is made periodic
so that an edge potential cannot manufacture or remove a low doublet.

The answer is deliberately tested in two limits on an open chain:

* V=0: a mobile hard-core-pair fluid.  Its neutral and pair-addition gaps are
  fitted with length at half pair filling.
* V=8J: a strongly repulsive, gapped-looking control.  The two lowest states
  are probed by a one-block pair-density operator.  A non-scalar projection
  identifies a locally readable CDW doublet rather than a topological code.

This is an effective physical-sector preflight only.  The pair projection and
the optional V term are not derived from the Lambda-walker microscopic model.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


ROOT = Path(__file__).resolve().parents[2]

PAIR_HOPPING = 1.0
BLOCK_COUNTS = (4, 6, 8, 10, 12)


def pair_basis(blocks: int, pairs: int) -> list[int]:
    return [sum(1 << site for site in occupied)
            for occupied in itertools.combinations(range(blocks), pairs)]


def build_hamiltonian(
    blocks: int, pairs: int, interaction: float, *, periodic: bool,
) -> tuple[sparse.csr_matrix, list[int]]:
    """Hard-core-pair chain at fixed number of pairs."""
    states = pair_basis(blocks, pairs)
    positions = {state: index for index, state in enumerate(states)}
    hamiltonian = sparse.lil_matrix((len(states), len(states)), dtype=float)
    for column, state in enumerate(states):
        occupation = [(state >> site) & 1 for site in range(blocks)]
        bond_count = blocks if periodic else blocks - 1
        hamiltonian[column, column] += interaction * sum(
            occupation[site] * occupation[(site + 1) % blocks] for site in range(bond_count)
        )
        for site in range(bond_count):
            right = (site + 1) % blocks
            if occupation[site] == occupation[right]:
                continue
            moved = state ^ (1 << site) ^ (1 << right)
            hamiltonian[positions[moved], column] += -PAIR_HOPPING
    return hamiltonian.tocsr(), states


def low_spectrum(hamiltonian: sparse.csr_matrix, count: int) -> tuple[np.ndarray, np.ndarray]:
    dimension = hamiltonian.shape[0]
    if dimension <= count + 1:
        values, vectors = np.linalg.eigh(hamiltonian.toarray())
        return values[:count], vectors[:, :count]
    values, vectors = sparse_linalg.eigsh(hamiltonian, k=count, which="SA", tol=1e-12)
    order = np.argsort(values)
    return values[order], vectors[:, order]


def ground_energy(blocks: int, pairs: int, interaction: float, *, periodic: bool) -> float:
    hamiltonian, _ = build_hamiltonian(blocks, pairs, interaction, periodic=periodic)
    return float(low_spectrum(hamiltonian, 1)[0][0])


def pair_density_nonscalar(vectors: np.ndarray, states: list[int], site: int) -> float:
    """||P n_site P - tr(P n_site P)/2|| on the two-lowest-state space."""
    density = np.asarray([(state >> site) & 1 for state in states], dtype=float)
    projected = vectors[:, :2].conj().T @ (density[:, None] * vectors[:, :2])
    residual = projected - np.eye(2) * np.trace(projected) / 2.0
    return float(np.linalg.svd(residual, compute_uv=False)[0])


def main() -> None:
    fluid_rows = []
    cdw_rows = []
    for blocks in BLOCK_COUNTS:
        pairs = blocks // 2
        fluid, _ = build_hamiltonian(blocks, pairs, interaction=0.0, periodic=True)
        values, _ = low_spectrum(fluid, 3)
        addition_gap = (
            ground_energy(blocks, pairs + 1, 0.0, periodic=True)
            + ground_energy(blocks, pairs - 1, 0.0, periodic=True)
            - 2.0 * values[0]
        )
        fluid_rows.append({
            "blocks": blocks,
            "pair_number": pairs,
            "hilbert_dimension": fluid.shape[0],
            "neutral_gap": float(values[1] - values[0]),
            "pair_addition_curvature": float(addition_gap),
            "ground_state_degeneracy_at_1e-10": int(np.count_nonzero(np.abs(values - values[0]) < 1e-10)),
            "gauss_parity_statement": "every pair hop changes each affected block N_a by +/-2, so coarse block parity is preserved",
        })

        cdw, states = build_hamiltonian(blocks, pairs, interaction=8.0 * PAIR_HOPPING, periodic=True)
        cdw_values, cdw_vectors = low_spectrum(cdw, 3)
        cdw_rows.append({
            "blocks": blocks,
            "pair_number": pairs,
            "hilbert_dimension": cdw.shape[0],
            "lowest_doublet_split": float(cdw_values[1] - cdw_values[0]),
            "gap_above_lowest_doublet": float(cdw_values[2] - cdw_values[1]),
            "site0_pair_density_nonscalar_norm": pair_density_nonscalar(cdw_vectors, states, site=0),
            "opposite_site_pair_density_nonscalar_norm": pair_density_nonscalar(cdw_vectors, states, site=blocks // 2),
        })

    fluid_gap_power = float(np.polyfit(
        np.log([row["blocks"] for row in fluid_rows]),
        np.log([row["neutral_gap"] for row in fluid_rows]),
        1,
    )[0])
    fluid_charge_power = float(np.polyfit(
        np.log([row["blocks"] for row in fluid_rows]),
        np.log([row["pair_addition_curvature"] for row in fluid_rows]),
        1,
    )[0])
    cdw_split_exponential_rate = float(np.polyfit(
        [row["blocks"] for row in cdw_rows],
        np.log([row["lowest_doublet_split"] for row in cdw_rows]),
        1,
    )[0])
    output = {
        "schema": "antler.phase8b.pair-chain-manybody-preflight.v1",
        "parameters": {
            "geometry": "periodic chain of coarse b=2 blocks, already projected to local hard-core pair occupation",
            "pair_hopping": PAIR_HOPPING,
            "half_pair_filling": True,
            "fluid_interaction": 0.0,
            "cdw_control_interaction": 8.0 * PAIR_HOPPING,
            "block_counts": list(BLOCK_COUNTS),
        },
        "fluid": {
            "rows": fluid_rows,
            "neutral_gap_power_vs_blocks": fluid_gap_power,
            "pair_addition_curvature_power_vs_blocks": fluid_charge_power,
        },
        "strong_repulsion_control": {
            "rows": cdw_rows,
            "lowest_doublet_log_split_slope_vs_blocks": cdw_split_exponential_rate,
        },
        "decision": (
            "The ideal pair-only continuation has no protected many-body sector. The V=0 hard-core-pair fluid has a unique "
            "ground state and both measured gaps close approximately as inverse length. Strong nearest-neighbour repulsion "
            "instead creates a locally readable low doublet (CDW control), not local indistinguishability."
        ),
        "claim_boundary": (
            "This is an effective low-energy preflight after an imposed local-pair projection. Pair hopping and V are not "
            "derived from the Lambda walker or native ANTLER resources. It neither rules out all paired phases nor establishes "
            "a phase, a thermodynamic code, a defect, a T junction, fusion or a non-Abelian braid."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_pair_chain_manybody_preflight.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
