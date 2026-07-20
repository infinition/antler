"""Exhaustive charge-conserving local probes for Phase 7 code audits.

The Phase 7 gate cannot rely on a hand-picked set of densities: a local rail
transfer may be the very operator that exposes a false code.  This module
constructs a Hermitian basis of all local weighted-charge-conserving matrix
units on a chosen set of hard-core modes, embedded in a fixed-charge Fock
basis. Rail modes may be fermionic while charge-two mediators are declared
hard-core commuting modes.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np


def _local_charge(pattern: int, charges: tuple[int, ...]) -> int:
    return sum(((pattern >> site) & 1) * charges[site] for site in range(len(charges)))


def _fermionic_sign(state: int, mode: int) -> float:
    return -1.0 if ((state & ((1 << mode) - 1)).bit_count() & 1) else 1.0


def _transition_matrix(
    states: np.ndarray,
    index: dict[int, int],
    support: tuple[int, ...],
    source_pattern: int,
    target_pattern: int,
    fermionic_modes: frozenset[int],
) -> np.ndarray:
    """Embed one local occupation matrix unit with canonical Fock signs."""
    dimension = len(states)
    out = np.zeros((dimension, dimension), dtype=complex)
    source_mask = sum(((source_pattern >> position) & 1) << mode for position, mode in enumerate(support))
    target_mask = sum(((target_pattern >> position) & 1) << mode for position, mode in enumerate(support))
    support_mask = sum(1 << mode for mode in support)
    remove = tuple(mode for mode in support if (source_mask >> mode) & 1 and not ((target_mask >> mode) & 1))
    add = tuple(mode for mode in support if (target_mask >> mode) & 1 and not ((source_mask >> mode) & 1))
    for column, raw_state in enumerate(states):
        state = int(raw_state)
        if state & support_mask != source_mask:
            continue
        current, amplitude = state, 1.0
        # The operations implement a canonical local monomial.  The eventual
        # Hermitian probes are explicitly symmetrized with their adjoints.
        for mode in remove:
            if mode in fermionic_modes:
                amplitude *= _fermionic_sign(current, mode)
            current ^= 1 << mode
        for mode in add:
            if mode in fermionic_modes:
                amplitude *= _fermionic_sign(current, mode)
            current |= 1 << mode
        out[index[current], column] = amplitude
    return out


def charge_conserving_local_probes(
    states: np.ndarray,
    index: dict[int, int],
    mode_charges: tuple[int, ...],
    support_modes: tuple[int, ...],
    fermionic_modes: frozenset[int],
    label_prefix: str = "local",
) -> dict[str, np.ndarray]:
    """Return a Hermitian basis for the complete local charge-preserving algebra.

    This includes diagonal occupation projectors and both Hermitian quadratures
    of every same-weighted-charge transition.  All returned probes have
    operator norm at most one.
    """
    if not support_modes or len(set(support_modes)) != len(support_modes):
        raise ValueError("support modes must be nonempty and unique")
    if any(mode < 0 or mode >= len(mode_charges) for mode in support_modes):
        raise ValueError("support mode outside mode-charge table")
    local_charges = tuple(mode_charges[mode] for mode in support_modes)
    groups: dict[int, list[int]] = defaultdict(list)
    for pattern in range(1 << len(support_modes)):
        groups[_local_charge(pattern, local_charges)].append(pattern)
    probes: dict[str, np.ndarray] = {}
    for charge, patterns in sorted(groups.items()):
        for source in patterns:
            unit = _transition_matrix(
                states, index, support_modes, source, source, fermionic_modes,
            )
            probes[f"{label_prefix}:q{charge}:diag:{source:0{len(support_modes)}b}"] = unit
        for source_index, source in enumerate(patterns):
            for target in patterns[source_index + 1:]:
                forward = _transition_matrix(
                    states, index, support_modes, source, target, fermionic_modes,
                )
                # Explicit symmetrization makes the construction insensitive to
                # the chosen canonical monomial phase for the reverse channel.
                hermitian_x = forward + forward.conj().T
                hermitian_y = -1j * (forward - forward.conj().T)
                stem = f"{label_prefix}:q{charge}:{source:0{len(support_modes)}b}->{target:0{len(support_modes)}b}"
                probes[f"{stem}:X"] = hermitian_x
                probes[f"{stem}:Y"] = hermitian_y
    return probes


def exhaustive_local_code_metric(frame: np.ndarray, probes: dict[str, np.ndarray]) -> dict:
    """Return the worst non-scalar code projection over a complete probe basis."""
    G = np.asarray(frame, dtype=complex)
    if G.ndim != 2 or not np.allclose(G.conj().T @ G, np.eye(G.shape[1]), atol=1e-11):
        raise ValueError("frame must have orthonormal columns")
    dimension, code_dimension = G.shape
    rows = []
    for label, operator in probes.items():
        if operator.shape != (dimension, dimension):
            raise ValueError(f"probe {label} has incompatible shape")
        logical = G.conj().T @ operator @ G
        scalar = np.trace(logical) / code_dimension
        rows.append({
            "label": label,
            "projected_non_scalar_frobenius": float(
                np.linalg.norm(logical - scalar * np.eye(code_dimension))
            ),
            "logical_action_frobenius": float(np.linalg.norm(logical)),
        })
    worst = max(rows, key=lambda row: row["projected_non_scalar_frobenius"])
    return {
        "probe_count": len(rows),
        "worst_projected_non_scalar_frobenius": worst["projected_non_scalar_frobenius"],
        "worst_probe": worst,
        "all_probes_scalar_at_1e_minus_10": bool(worst["projected_non_scalar_frobenius"] < 1e-10),
    }
