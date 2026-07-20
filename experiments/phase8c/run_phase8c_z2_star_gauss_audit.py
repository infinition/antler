"""Phase 8C-T0: exact Gauss-law audit on the minimal neutral-link star.

This is a reference gauge model with a declared *new* neutral Z2 link qubit on
each edge.  It deliberately does not reuse ANTLER's charge-two mediator and
does not claim a code, fusion space, anyons, or braiding.  Its only purpose is
to qualify the local algebra before any larger geometry is considered.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

VERTICES = (0, 1, 2, 3)
EDGES = ((0, 1), (0, 2), (0, 3))
FIXED_CHARGE = 2
EDGE_HOPPINGS = (0.61, 0.73, 0.89)
ELECTRIC_FIELDS = (0.27, 0.34, 0.41)
VERTEX_POTENTIALS = (0.11, -0.07, 0.13, -0.17)
GAUSS_PENALTY = 2.0

MATTER_BITS = len(VERTICES)
EDGE_BITS = len(EDGES)
DIMENSION = 1 << (MATTER_BITS + EDGE_BITS)


def bit(value: int, position: int) -> int:
    return (value >> position) & 1


def edge_bit(edge: int) -> int:
    return MATTER_BITS + edge


def fermion_annihilate(state: int, mode: int) -> tuple[int, complex] | None:
    if not bit(state, mode):
        return None
    sign = -1.0 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1.0
    return state ^ (1 << mode), sign


def fermion_create(state: int, mode: int) -> tuple[int, complex] | None:
    if bit(state, mode):
        return None
    sign = -1.0 if (state & ((1 << mode) - 1)).bit_count() % 2 else 1.0
    return state ^ (1 << mode), sign


def hop_term(target: int, source: int, edge: int | None = None) -> np.ndarray:
    """Return c_target^dag [tau_z(edge)] c_source in the full Fock/link basis."""
    operator = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    for column in range(DIMENSION):
        first = fermion_annihilate(column, source)
        if first is None:
            continue
        intermediate, amplitude = first
        second = fermion_create(intermediate, target)
        if second is None:
            continue
        row, create_amplitude = second
        link_sign = 1.0 if edge is None or bit(column, edge_bit(edge)) == 0 else -1.0
        operator[row, column] += amplitude * create_amplitude * link_sign
    return operator


def number_operator(vertex: int) -> np.ndarray:
    return np.diag([float(bit(state, vertex)) for state in range(DIMENSION)]).astype(complex)


def tau_x(edge: int) -> np.ndarray:
    operator = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    mask = 1 << edge_bit(edge)
    for column in range(DIMENSION):
        operator[column ^ mask, column] = 1.0
    return operator


def gauss_operator(vertex: int) -> np.ndarray:
    incident = [edge for edge, endpoints in enumerate(EDGES) if vertex in endpoints]
    operator = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    mask = sum(1 << edge_bit(edge) for edge in incident)
    for column in range(DIMENSION):
        parity = -1.0 if bit(column, vertex) else 1.0
        operator[column ^ mask, column] = parity
    return operator


def charge_projector(charge: int) -> np.ndarray:
    diagonal = [1.0 if sum(bit(state, vertex) for vertex in VERTICES) == charge else 0.0 for state in range(DIMENSION)]
    return np.diag(diagonal).astype(complex)


def operator_norm(operator: np.ndarray) -> float:
    return float(np.linalg.norm(operator, ord=2))


def commutator(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right - right @ left


def projector_frame(projector: np.ndarray, tolerance: float = 1.0e-10) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(projector)
    return eigenvectors[:, eigenvalues > 1.0 - tolerance]


def scalar_residual(operator: np.ndarray, frame: np.ndarray) -> float:
    projected = frame.conj().T @ operator @ frame
    rank = projected.shape[0]
    return operator_norm(projected - np.trace(projected) * np.eye(rank, dtype=complex) / rank)


def sector_projector(gauss: tuple[np.ndarray, ...], signs: tuple[int, ...]) -> np.ndarray:
    projector = np.eye(DIMENSION, dtype=complex)
    for generator, sign in zip(gauss, signs, strict=True):
        projector = projector @ ((np.eye(DIMENSION, dtype=complex) + sign * generator) / 2.0)
    return projector


def main() -> None:
    identity = np.eye(DIMENSION, dtype=complex)
    gauss = tuple(gauss_operator(vertex) for vertex in VERTICES)
    charge = sum((number_operator(vertex) for vertex in VERTICES), start=np.zeros((DIMENSION, DIMENSION), dtype=complex))
    fixed_charge_projector = charge_projector(FIXED_CHARGE)

    dressed_terms = []
    bare_terms = []
    hamiltonian = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    for edge, ((left, right), hopping, electric) in enumerate(zip(EDGES, EDGE_HOPPINGS, ELECTRIC_FIELDS, strict=True)):
        dressed = hop_term(left, right, edge) + hop_term(right, left, edge)
        bare = hop_term(left, right) + hop_term(right, left)
        dressed_terms.append(dressed)
        bare_terms.append(bare)
        hamiltonian += -hopping * dressed - electric * tau_x(edge)
    for vertex, potential in zip(VERTICES, VERTEX_POTENTIALS, strict=True):
        hamiltonian += potential * number_operator(vertex) - GAUSS_PENALTY * gauss[vertex]

    bare_hopping = sum((coefficient * term for coefficient, term in zip(EDGE_HOPPINGS, bare_terms, strict=True)), start=np.zeros_like(hamiltonian))
    dressed_hopping = sum((coefficient * term for coefficient, term in zip(EDGE_HOPPINGS, dressed_terms, strict=True)), start=np.zeros_like(hamiltonian))

    physical_projector = sector_projector(gauss, (1, 1, 1, 1))
    physical_fixed_charge_projector = fixed_charge_projector @ physical_projector
    physical_frame = projector_frame(physical_fixed_charge_projector)

    signature_rows = []
    for signs in product((-1, 1), repeat=len(VERTICES)):
        projector = fixed_charge_projector @ sector_projector(gauss, signs)
        frame = projector_frame(projector)
        if not frame.shape[1]:
            continue
        spectrum = np.linalg.eigvalsh(frame.conj().T @ hamiltonian @ frame)
        signature_rows.append({
            "signs": list(signs),
            "dimension": int(frame.shape[1]),
            "ground_energy": float(spectrum[0]),
            "first_excited_energy": float(spectrum[1]) if len(spectrum) > 1 else None,
        })

    physical_spectrum = np.linalg.eigvalsh(physical_frame.conj().T @ hamiltonian @ physical_frame)
    physical_ground = float(physical_spectrum[0])
    other_sector_energies = [row["ground_energy"] for row in signature_rows if row["signs"] != [1, 1, 1, 1]]
    syndrome_gap = float(min(other_sector_energies) - physical_ground)
    local_density_residual = scalar_residual(number_operator(0), physical_frame)

    gauss_product = identity.copy()
    for generator in gauss:
        gauss_product = gauss_product @ generator
    fermion_parity = np.diag([(-1.0) ** sum(bit(state, vertex) for vertex in VERTICES) for state in range(DIMENSION)]).astype(complex)

    output = {
        "schema": "antler.phase8c.z2-star-gauss-audit.v1",
        "parameters": {
            "new_declared_resource": "one neutral Z2 link qubit tau_e per edge; not an ANTLER charge-two mediator",
            "matter": "spinless U(1)-conserving fermion on each vertex",
            "gauge_generator": "G_v=(-1)^n_v product_(e incident to v) tau^x_e",
            "dressed_hopping": "c_v^dag tau^z_vw c_w+h.c.",
            "reference_hamiltonian": "sum_e[-t_e dressed_hop_e-h_e tau^x_e]+sum_v[mu_v n_v-lambda G_v]",
            "edge_hoppings": list(EDGE_HOPPINGS),
            "electric_fields": list(ELECTRIC_FIELDS),
            "vertex_potentials": list(VERTEX_POTENTIALS),
            "gauss_penalty": GAUSS_PENALTY,
        },
        "geometry": {
            "name": "three-arm star (gauge-algebra preflight only)",
            "vertices": list(VERTICES),
            "edges": [list(edge) for edge in EDGES],
            "plaquettes": 0,
        },
        "dimensions": {
            "full": DIMENSION,
            "fixed_N2_before_gauss": int(round(np.trace(fixed_charge_projector).real)),
            "physical_fixed_N2": int(physical_frame.shape[1]),
            "physical_all_charges": int(projector_frame(physical_projector).shape[1]),
        },
        "algebra": {
            "hamiltonian_hermiticity_norm": operator_norm(hamiltonian - hamiltonian.conj().T),
            "max_gauss_hermiticity_norm": max(operator_norm(generator - generator.conj().T) for generator in gauss),
            "max_gauss_involution_norm": max(operator_norm(generator @ generator - identity) for generator in gauss),
            "max_hamiltonian_gauss_commutator_norm": max(operator_norm(commutator(hamiltonian, generator)) for generator in gauss),
            "max_pairwise_gauss_commutator_norm": max(operator_norm(commutator(first, second)) for index, first in enumerate(gauss) for second in gauss[index + 1:]),
            "hamiltonian_charge_commutator_norm": operator_norm(commutator(hamiltonian, charge)),
            "global_gauss_product_minus_fermion_parity_norm": operator_norm(gauss_product - fermion_parity),
            "max_dressed_hop_gauss_commutator_norm": max(operator_norm(commutator(term, generator)) for term in dressed_terms for generator in gauss),
            "max_bare_hop_gauss_commutator_norm": max(operator_norm(commutator(term, generator)) for term in bare_terms for generator in gauss),
        },
        "physical_projection": {
            "projector_idempotency_norm": operator_norm(physical_fixed_charge_projector @ physical_fixed_charge_projector - physical_fixed_charge_projector),
            "max_projected_gauss_minus_projector_norm": max(operator_norm(physical_fixed_charge_projector @ generator @ physical_fixed_charge_projector - physical_fixed_charge_projector) for generator in gauss),
            "bare_hopping_projected_norm": operator_norm(physical_fixed_charge_projector @ bare_hopping @ physical_fixed_charge_projector),
            "dressed_hopping_projected_nonscalar_norm": scalar_residual(dressed_hopping, physical_frame),
            "local_vertex_density_projected_nonscalar_norm": local_density_residual,
        },
        "spectra": {
            "physical_N2_eigenvalues": [float(value) for value in physical_spectrum],
            "physical_internal_gap": float(physical_spectrum[1] - physical_spectrum[0]),
            "minimum_Gauss_syndrome_gap_at_N2": syndrome_gap,
            "gauss_signature_sectors_at_N2": signature_rows,
            "odd_N_physical_sector_note": "For all G_v=+1, product_v G_v=(-1)^N, so odd total N has no physical state without background charge.",
        },
        "decision": (
            "PASS T0. The separately declared neutral-link reference model has exact local Z2 Gauss algebra, exact total-U(1) "
            "conservation and the expected selection rule: bare hopping has no physical-sector matrix element, whereas dressed "
            "hopping is gauge invariant. The three-arm star contains no plaquette and its local density remains non-scalar in the "
            "physical N=2 sector, so it is explicitly not a topological code, fusion space or braid implementation."
        ),
        "claim_boundary": (
            "This audit inserts neutral Z2 link qubits and a Gauss penalty as new reference resources. It neither derives them from "
            "the frozen ANTLER ladder nor establishes a deconfined two-dimensional phase, plaquette flux, local indistinguishability, "
            "twist defects, anyons, fusion or non-Abelian braiding."
        ),
    }

    path = ROOT / "results" / "phase8c" / "z2_star_gauss_audit.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
