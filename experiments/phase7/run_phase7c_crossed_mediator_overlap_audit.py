"""Stress-test two overlapping crossed charge-two plaquettes.

The local crossed primitive is not useful for a two-dimensional parent unless
it survives repetition.  This script constructs two four-rung motifs on the
six-rung supports (0,1,2,3) and (2,3,4,5), adds their *separately calibrated*
Z/ZZ counterterms, and exactly downfolds the fixed-charge block.  The overlap
of two rungs is a deliberately compact compatibility stress-test; it is not a
claim that this support is already a complete surface-code tiling.
"""
from __future__ import annotations

from functools import lru_cache
import json
from itertools import product
from math import comb
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_microscopic_optimizer import PAULI_MATRICES
from experiments.phase7.run_phase7c_crossed_mediator_cancellation_audit import (
    U, DELTA, base_candidate, static_counterterms,
)


RUNG_COUNT = 6
LOW_MODE_COUNT = 2 * RUNG_COUNT
LOW_DIMENSION = 1 << RUNG_COUNT
TOTAL_CHARGE = RUNG_COUNT
G = 0.5
SUPPORTS = ((0, 1, 2, 3), (2, 3, 4, 5))


def rail_mode(rung: int, rail: int) -> int:
    return 2 * rung + rail


def fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if (state & ((1 << mode) - 1)).bit_count() & 1 else 1.0


def annihilate_pair(state: int, first: int, second: int) -> tuple[int, float] | None:
    first, second = sorted((first, second))
    if not ((state >> first) & 1):
        return None
    intermediate = state ^ (1 << first)
    sign = fermionic_sign(state, first)
    if not ((intermediate >> second) & 1):
        return None
    return intermediate ^ (1 << second), sign * fermionic_sign(intermediate, second)


def crossed_channels() -> tuple[tuple[str, tuple[tuple[int, int], ...]], ...]:
    channels = []
    for plaquette, support in enumerate(SUPPORTS):
        left, right, third, fourth = support
        channels.extend((
            (f"p{plaquette}_aa_left_bb_right", ((rail_mode(left, 0), rail_mode(right, 0)), (rail_mode(third, 1), rail_mode(fourth, 1)))),
            (f"p{plaquette}_bb_left_aa_right", ((rail_mode(left, 1), rail_mode(right, 1)), (rail_mode(third, 0), rail_mode(fourth, 0)))),
        ))
    return tuple(channels)


def fixed_charge_basis(channel_count: int) -> tuple[np.ndarray, dict[int, int]]:
    states = np.asarray([
        state for state in range(1 << (LOW_MODE_COUNT + channel_count))
        if (state & ((1 << LOW_MODE_COUNT) - 1)).bit_count()
        + 2 * (state >> LOW_MODE_COUNT).bit_count() == TOTAL_CHARGE
    ], dtype=np.int64)
    expected = sum(
        comb(channel_count, mediators) * comb(LOW_MODE_COUNT, TOTAL_CHARGE - 2 * mediators)
        for mediators in range(min(channel_count, TOTAL_CHARGE // 2) + 1)
    )
    if len(states) != expected:
        raise RuntimeError("fixed-charge basis count mismatch")
    return states, {int(state): index for index, state in enumerate(states)}


def monomer_embedding(states: np.ndarray, index: dict[int, int]) -> np.ndarray:
    frame = np.zeros((len(states), LOW_DIMENSION), dtype=complex)
    for column, rails in enumerate(product((0, 1), repeat=RUNG_COUNT)):
        state = sum(1 << rail_mode(rung, rail) for rung, rail in enumerate(rails))
        frame[index[state], column] = 1.0
    return frame


def local_counterterms() -> tuple[tuple[float, ...], dict[tuple[int, int], float]]:
    """Embed the independently fitted one-plaquette Z/ZZ cancellation twice."""
    audit = __import__("antler.phase7_microscopic_optimizer", fromlist=["evaluate_local_candidate"]).evaluate_local_candidate(
        base_candidate(G), "XXXX", include_pauli_coefficients=True
    )
    local_biases, local_zz = static_counterterms(audit["full_traceless_pauli_coefficients"])
    biases = [0.0] * RUNG_COUNT
    zz: dict[tuple[int, int], float] = {}
    for support in SUPPORTS:
        for local_rung, strength in enumerate(local_biases):
            biases[support[local_rung]] += strength
        for coupling in local_zz:
            key = tuple(sorted((support[coupling.left_edge], support[coupling.right_edge])))
            zz[key] = zz.get(key, 0.0) + coupling.strength
    return tuple(biases), zz


def build_hamiltonian(include_counterterms: bool) -> tuple[coo_matrix, np.ndarray, dict[int, int], np.ndarray, np.ndarray, dict]:
    channels = crossed_channels()
    states, index = fixed_charge_basis(len(channels))
    biases, zz = local_counterterms() if include_counterterms else ((0.0,) * RUNG_COUNT, {})
    rows, columns, data = [], [], []
    mott_cost = np.zeros(len(states), dtype=float)
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        diagonal = 0.0
        for rung in range(RUNG_COUNT):
            n_a, n_b = (state >> rail_mode(rung, 0)) & 1, (state >> rail_mode(rung, 1)) & 1
            cell_cost = float((n_a + n_b - 1) ** 2)
            mott_cost[column] += cell_cost
            diagonal += U * cell_cost + biases[rung] * (n_a - n_b)
        for (left, right), strength in zz.items():
            z_left = ((state >> rail_mode(left, 0)) & 1) - ((state >> rail_mode(left, 1)) & 1)
            z_right = ((state >> rail_mode(right, 0)) & 1) - ((state >> rail_mode(right, 1)) & 1)
            diagonal += strength * z_left * z_right
        for channel_index, (_, pairs) in enumerate(channels):
            mediator = LOW_MODE_COUNT + channel_index
            diagonal += DELTA * ((state >> mediator) & 1)
        rows.append(column)
        columns.append(column)
        data.append(diagonal)
        for channel_index, (_, pairs) in enumerate(channels):
            mediator = LOW_MODE_COUNT + channel_index
            if (state >> mediator) & 1:
                continue
            for first, second in pairs:
                item = annihilate_pair(state, first, second)
                if item is None:
                    continue
                low_state, sign = item
                new_state = low_state | (1 << mediator)
                amplitude = -G * sign
                rows.extend((index[new_state], column))
                columns.extend((column, index[new_state]))
                data.extend((amplitude, amplitude))
    matrix = coo_matrix((np.asarray(data, dtype=complex), (rows, columns)), shape=(len(states), len(states))).tocsr()
    if np.linalg.norm((matrix - matrix.getH()).data) > 1e-10:
        raise RuntimeError("overlap compiler produced a non-Hermitian Hamiltonian")
    metadata = {
        "block_dimension": len(states), "channel_count": len(channels), "channels": [name for name, _ in channels],
        "rail_biases": list(biases),
        "zz_couplings": [{"rungs": list(key), "strength": value} for key, value in sorted(zz.items())],
    }
    return matrix, states, index, monomer_embedding(states, index), mott_cost, metadata


@lru_cache(maxsize=None)
def pauli_word(label: str) -> np.ndarray:
    operator = np.asarray([[1.0]], dtype=complex)
    for letter in label:
        operator = np.kron(operator, PAULI_MATRICES[letter])
    return operator


def pauli_coefficients(matrix: np.ndarray) -> dict[str, float]:
    rows = {}
    for letters in product("IXYZ", repeat=RUNG_COUNT):
        label = "".join(letters)
        coefficient = np.trace(pauli_word(label) @ matrix) / LOW_DIMENSION
        if abs(coefficient.imag) > 1e-8:
            raise RuntimeError(f"complex Pauli coefficient for {label}")
        rows[label] = float(coefficient.real)
    return rows


def all_flip_operator(support: tuple[int, ...]) -> np.ndarray:
    """Return prod S+ + prod S- on one four-rung support.

    This is the complete four-body combination generated by the crossed
    primitive: it contains XXXX but is not equal to XXXX.
    """
    raise_one = np.asarray([[0.0, 1.0], [0.0, 0.0]], dtype=complex)
    lower_one = raise_one.conj().T
    identity = np.eye(2, dtype=complex)
    raise_all = np.asarray([[1.0]], dtype=complex)
    lower_all = np.asarray([[1.0]], dtype=complex)
    for rung in range(RUNG_COUNT):
        raise_all = np.kron(raise_all, raise_one if rung in support else identity)
        lower_all = np.kron(lower_all, lower_one if rung in support else identity)
    return raise_all + lower_all


def diagonalize_by_bare_parity(
    hamiltonian: coo_matrix, states: np.ndarray, monomers: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Resolve exact parity sectors before extracting a degenerate low band.

    A generic iterative eigensolver can omit vectors from an exactly degenerate
    zero-energy sector.  The bare parities are exact here, so diagonalizing the
    two finite blocks with a Hermitian subset solver is both cheaper and a
    reliable degeneracy control.
    """
    state_sector = []
    for state in states:
        n_a = sum((int(state) >> rail_mode(rung, 0)) & 1 for rung in range(RUNG_COUNT))
        n_b = sum((int(state) >> rail_mode(rung, 1)) & 1 for rung in range(RUNG_COUNT))
        state_sector.append((n_a & 1, n_b & 1))
    monomer_sector = []
    for rails in product((0, 1), repeat=RUNG_COUNT):
        n_a = sum(rail == 0 for rail in rails)
        monomer_sector.append((n_a & 1, (RUNG_COUNT - n_a) & 1))
    low_values, low_vectors, high_values = [], [], []
    for sector in sorted(set(monomer_sector)):
        state_indices = np.asarray([index for index, value in enumerate(state_sector) if value == sector], dtype=int)
        logical_count = monomer_sector.count(sector)
        if not len(state_indices) > logical_count:
            raise RuntimeError("parity sector is too small to isolate its monomer subspace")
        sector_values, sector_vectors = eigh(
            hamiltonian[state_indices][:, state_indices].toarray(),
            subset_by_index=[0, logical_count], driver="evr", check_finite=False,
        )
        embedded = np.zeros((hamiltonian.shape[0], logical_count), dtype=complex)
        embedded[state_indices, :] = sector_vectors[:, :logical_count]
        low_values.extend(sector_values[:logical_count])
        low_vectors.append(embedded)
        high_values.append(float(sector_values[logical_count]))
    low_values_array = np.asarray(low_values)
    low_vectors_array = np.hstack(low_vectors)
    order = np.argsort(low_values_array)
    return low_values_array[order], low_vectors_array[:, order], min(high_values)


def audit(include_counterterms: bool) -> dict:
    hamiltonian, states, _, monomers, mott_cost, metadata = build_hamiltonian(include_counterterms)
    low_values, low_vectors, first_high_value = diagonalize_by_bare_parity(hamiltonian, states, monomers)
    overlap = monomers.conj().T @ low_vectors
    left, singular_values, right_dagger = np.linalg.svd(overlap)
    aligned = left @ right_dagger
    effective = aligned @ np.diag(low_values) @ aligned.conj().T
    effective = 0.5 * (effective + effective.conj().T)
    traceless = effective - np.trace(effective).real / LOW_DIMENSION * np.eye(LOW_DIMENSION)
    coefficients = pauli_coefficients(traceless)
    target_labels = ("XXXXII", "IIXXXX")
    target_coefficients = {label: coefficients[label] for label in target_labels}
    target = sum((-0.5 * pauli_word(label) for label in target_labels), np.zeros_like(traceless))
    residual = float(np.linalg.norm(traceless - target) / np.linalg.norm(target))
    nonidentity = {label: value for label, value in coefficients.items() if label != "I" * RUNG_COUNT}
    target_square = sum(value * value for label, value in nonidentity.items() if label in target_labels)
    total_square = sum(value * value for value in nonidentity.values())
    parity_a = np.asarray([(-1.0 if sum((int(state) >> rail_mode(rung, 0)) & 1 for rung in range(RUNG_COUNT)) & 1 else 1.0) for state in states])
    parity_b = np.asarray([(-1.0 if sum((int(state) >> rail_mode(rung, 1)) & 1 for rung in range(RUNG_COUNT)) & 1 else 1.0) for state in states])
    low_mott = np.real(np.sum(np.abs(low_vectors) ** 2 * mott_cost[:, None], axis=0))
    target_alignment = float(target_square / total_square) if total_square > 1e-24 else 0.0
    full_flips = {f"P{index}": all_flip_operator(support) for index, support in enumerate(SUPPORTS)}
    full_flip_coefficients = {
        label: float(np.trace(operator @ traceless).real / np.trace(operator @ operator).real)
        for label, operator in full_flips.items()
    }
    full_flip_square = sum(
        coefficient * coefficient * np.trace(full_flips[label] @ full_flips[label]).real
        for label, coefficient in full_flip_coefficients.items()
    )
    return {
        "metadata": metadata,
        "minimum_monomer_overlap_singular_value": float(np.min(singular_values)),
        "low_to_high_gap": float(first_high_value - low_values[-1]),
        "maximum_low_state_mott_violation": float(np.max(low_mott)),
        "bare_parity_commutator_frobenius_max": max(
            float(np.linalg.norm(hamiltonian.multiply(parity_a[None, :] - parity_a[:, None]).data)),
            float(np.linalg.norm(hamiltonian.multiply(parity_b[None, :] - parity_b[:, None]).data)),
        ),
        "target_coefficients": target_coefficients,
        "two_plaquette_target_alignment": target_alignment,
        "all_flip_coefficients": full_flip_coefficients,
        "all_flip_alignment": float(full_flip_square / np.linalg.norm(traceless) ** 2),
        "fixed_scale_two_plaquette_residual": residual,
        "top_non_target_paulis": sorted(
            ({"label": label, "coefficient": value} for label, value in nonidentity.items() if label not in target_labels),
            key=lambda row: abs(row["coefficient"]), reverse=True,
        )[:12],
    }


def main() -> None:
    base = audit(include_counterterms=False)
    compensated = audit(include_counterterms=True)
    first_flip, second_flip = (all_flip_operator(support) for support in SUPPORTS)
    flip_commutator = first_flip @ second_flip - second_flip @ first_flip
    out = {
        "schema": "antler.phase7c.crossed-charge2-overlap-audit.v1",
        "geometry": {
            "rung_count": RUNG_COUNT, "logical_monomer_dimension": LOW_DIMENSION,
            "plaquette_supports": [list(support) for support in SUPPORTS], "shared_rungs": [2, 3],
            "all_flip_commutator_frobenius": float(np.linalg.norm(flip_commutator)),
            "all_flip_commutator_spectral": float(np.linalg.norm(flip_commutator, 2)),
        },
        "parameters": {"U": U, "Delta": DELTA, "g": G, "g_over_delta": G / DELTA},
        "base": base,
        "separately_calibrated_counterterms": compensated,
        "decision": (
            "The separately calibrated motifs preserve a well-isolated monomer manifold and retain their local all-flip "
            "component. However, the two all-flip operators have a nonzero exact commutator on this overlap, so they are "
            "not a commuting stabilizer parent in this geometry. A selective repeated stabilizer and an independent SW "
            "derivation remain required."
        ),
        "claim_boundary": (
            "Two overlapping supports are not a surface-code lattice, and counterterms are imported from a one-plaquette "
            "finite-block fit. This audit cannot establish topological order, a protected global code, braiding or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase7" / "crossed_charge2_overlap_audit.json"
    result.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "base": {key: base[key] for key in ("minimum_monomer_overlap_singular_value", "low_to_high_gap", "target_coefficients", "two_plaquette_target_alignment", "fixed_scale_two_plaquette_residual")},
        "compensated": {key: compensated[key] for key in ("minimum_monomer_overlap_singular_value", "low_to_high_gap", "target_coefficients", "two_plaquette_target_alignment", "fixed_scale_two_plaquette_residual")},
    }, indent=2))


if __name__ == "__main__":
    main()
