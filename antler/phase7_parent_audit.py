"""Model-agnostic dense audits for Phase 7 parent-Hamiltonian proposals.

This module intentionally does not define a new physical Hamiltonian.  It is a
small audit contract for a theory-supplied finite-size parent: local terms,
their supports, exact symmetry labels, a proposed code frame, local probes,
and proposed edge operators.  The same diagnostics must be passed before a
candidate is promoted from algebra to a larger sparse or matrix-free study.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np


def _as_square(name: str, matrix: np.ndarray) -> np.ndarray:
    array = np.asarray(matrix, dtype=complex)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be square")
    return array


def _as_frame(frame: np.ndarray, dimension: int) -> np.ndarray:
    array = np.asarray(frame, dtype=complex)
    if array.ndim != 2 or array.shape[0] != dimension or array.shape[1] < 1:
        raise ValueError("code frame has incompatible shape")
    if not np.allclose(array.conj().T @ array, np.eye(array.shape[1]), atol=1e-11):
        raise ValueError("code frame must have orthonormal columns")
    return array


def frobenius_commutator(left: np.ndarray, right: np.ndarray) -> float:
    """Return ``||[left,right]||_F`` after shape validation."""
    A, B = _as_square("left", left), _as_square("right", right)
    if A.shape != B.shape:
        raise ValueError("commutator operands have incompatible shapes")
    return float(np.linalg.norm(A @ B - B @ A))


def local_projector_algebra(
    terms: Sequence[np.ndarray], supports: Sequence[Sequence[int]], names: Sequence[str] | None = None,
) -> dict:
    """Audit Hermiticity, idempotency, and pairwise local commutativity.

    A theory may deliberately provide noncommuting parent terms.  In that case
    this function records the failure rather than declaring a commuting-parent
    proof.  Supports are labels only; they let the report separate disjoint
    and overlapping pairs.
    """
    if not terms or len(terms) != len(supports):
        raise ValueError("terms and supports must be nonempty and have equal length")
    matrices = [_as_square(f"term[{index}]", term) for index, term in enumerate(terms)]
    dimension = matrices[0].shape[0]
    if any(term.shape != (dimension, dimension) for term in matrices):
        raise ValueError("all local terms must share a Hilbert-space dimension")
    labels = list(names) if names is not None else [f"Pi_{index}" for index in range(len(matrices))]
    if len(labels) != len(matrices):
        raise ValueError("names must have one entry per term")
    hermiticity = [float(np.linalg.norm(term - term.conj().T)) for term in matrices]
    idempotency = [float(np.linalg.norm(term @ term - term)) for term in matrices]
    pair_rows = []
    for left in range(len(matrices)):
        for right in range(left + 1, len(matrices)):
            shared = sorted(set(supports[left]) & set(supports[right]))
            pair_rows.append({
                "left": labels[left],
                "right": labels[right],
                "overlapping_support": bool(shared),
                "shared_sites": shared,
                "commutator_frobenius": frobenius_commutator(matrices[left], matrices[right]),
            })
    overlapping = [row["commutator_frobenius"] for row in pair_rows if row["overlapping_support"]]
    disjoint = [row["commutator_frobenius"] for row in pair_rows if not row["overlapping_support"]]
    return {
        "term_count": len(matrices),
        "dimension": dimension,
        "max_term_hermiticity_error": float(max(hermiticity)),
        "max_term_idempotency_error": float(max(idempotency)),
        "pairwise_commutators": pair_rows,
        "max_overlapping_commutator": float(max(overlapping, default=0.0)),
        "max_disjoint_commutator": float(max(disjoint, default=0.0)),
        "all_terms_are_projectors_at_1e_minus_10": bool(max(idempotency) < 1e-10),
        "all_terms_commute_at_1e_minus_10": bool(
            max((row["commutator_frobenius"] for row in pair_rows), default=0.0) < 1e-10
        ),
    }


def symmetry_audit(H: np.ndarray, diagonal_labels: Mapping[str, Sequence[float | int]]) -> dict:
    """Audit symmetries supplied as diagonal eigenvalue labels in this basis."""
    matrix = _as_square("H", H)
    rows = []
    for name, labels in diagonal_labels.items():
        values = np.asarray(labels, dtype=complex)
        if values.shape != (matrix.shape[0],):
            raise ValueError(f"symmetry labels for {name} have incompatible length")
        generator = np.diag(values)
        rows.append({"name": name, "commutator_frobenius": frobenius_commutator(matrix, generator)})
    return {
        "symmetries": rows,
        "max_commutator_frobenius": float(max((row["commutator_frobenius"] for row in rows), default=0.0)),
        "all_commute_at_1e_minus_10": bool(
            max((row["commutator_frobenius"] for row in rows), default=0.0) < 1e-10
        ),
    }


def projected_edge_metrics(H: np.ndarray, frame: np.ndarray, operator: np.ndarray) -> dict:
    """Compute logical action, leakage, and ``||(I-P)[H,O]P||/||OP||``."""
    matrix = _as_square("H", H)
    G = _as_frame(frame, matrix.shape[0])
    O = _as_square("operator", operator)
    if O.shape != matrix.shape:
        raise ValueError("operator has incompatible shape")
    operated = O @ G
    denominator = float(np.linalg.norm(operated))
    if denominator <= 1e-15:
        raise ValueError("operator annihilates the proposed code frame")
    logical = G.conj().T @ operated
    leakage = operated - G @ logical
    commutator_on_code = matrix @ operated - O @ (matrix @ G)
    projected_commutator = commutator_on_code - G @ (G.conj().T @ commutator_on_code)
    return {
        "operator_on_code_frobenius": denominator,
        "logical_action_frobenius": float(np.linalg.norm(logical)),
        "logical_leakage_amplitude_frobenius": float(np.linalg.norm(leakage)),
        "full_commutator_on_code_frobenius": float(np.linalg.norm(commutator_on_code)),
        "code_commutator_action_frobenius": float(np.linalg.norm(projected_commutator)),
        "code_commutator_action_normalized": float(np.linalg.norm(projected_commutator) / denominator),
    }


def local_indistinguishability(
    frame: np.ndarray, operators: Mapping[str, np.ndarray],
) -> dict:
    """Measure the non-scalar action of each proposed local probe on the code."""
    first = next(iter(operators.values()), None)
    if first is None:
        raise ValueError("at least one local operator is required")
    dimension = _as_square("local operator", first).shape[0]
    G = _as_frame(frame, dimension)
    code_dimension = G.shape[1]
    rows = []
    for name, operator in operators.items():
        O = _as_square(name, operator)
        if O.shape != (dimension, dimension):
            raise ValueError(f"local operator {name} has incompatible shape")
        projected = G.conj().T @ O @ G
        scalar = np.trace(projected) / code_dimension
        rows.append({
            "name": name,
            "projected_non_scalar_frobenius": float(
                np.linalg.norm(projected - scalar * np.eye(code_dimension))
            ),
        })
    return {
        "operators": rows,
        "worst_projected_non_scalar_frobenius": float(
            max(row["projected_non_scalar_frobenius"] for row in rows)
        ),
    }


def sector_spectrum(
    H: np.ndarray,
    labels: Sequence[float | int],
    sector: float | int,
    code_multiplicity: int = 1,
) -> dict:
    """Dense finite-size spectrum in one exact diagonal symmetry sector."""
    matrix = _as_square("H", H)
    values = np.asarray(labels)
    if values.shape != (matrix.shape[0],):
        raise ValueError("sector labels have incompatible length")
    if code_multiplicity < 1:
        raise ValueError("code_multiplicity must be positive")
    rows = np.flatnonzero(values == sector)
    if len(rows) <= code_multiplicity:
        raise ValueError("sector has no excitation above requested code multiplicity")
    block = matrix[np.ix_(rows, rows)]
    energies = np.linalg.eigvalsh(block)
    return {
        "sector": float(sector),
        "sector_dimension": int(len(rows)),
        "lowest_energies": [float(value) for value in energies[: code_multiplicity + 1]],
        "ground_energy": float(energies[0]),
        "gap_above_requested_multiplicity": float(energies[code_multiplicity] - energies[code_multiplicity - 1]),
    }
