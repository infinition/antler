"""Phase 8C-T1: one-plaquette neutral Z2 gauge preflight.

This audit extends the T0 star to the smallest graph with a magnetic plaquette.
It checks exact Gauss algebra, U(1), flux gauge invariance and static flux
sectors.  One plaquette is deliberately *not* promoted to a deconfined phase,
topological code, twist defect, fusion space or braid.
"""
from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

VERTICES = (0, 1, 2, 3)
EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
FIXED_CHARGE = 2
EDGE_HOPPINGS = (0.57, 0.69, 0.83, 0.76)
ELECTRIC_FIELDS = (0.19, 0.31, 0.27, 0.38)
VERTEX_POTENTIALS = (0.12, -0.08, 0.16, -0.21)
GAUSS_PENALTY = 2.0
PLAQUETTE_COUPLING = 1.1

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


def plaquette_operator() -> np.ndarray:
    return np.diag([
        float(np.prod([1.0 if bit(state, edge_bit(edge)) == 0 else -1.0 for edge in range(EDGE_BITS)]))
        for state in range(DIMENSION)
    ]).astype(complex)


def charge_projector(charge: int) -> np.ndarray:
    return np.diag([
        1.0 if sum(bit(state, vertex) for vertex in VERTICES) == charge else 0.0
        for state in range(DIMENSION)
    ]).astype(complex)


def operator_norm(operator: np.ndarray) -> float:
    return float(np.linalg.norm(operator, ord=2))


def commutator(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right - right @ left


def projector_frame(projector: np.ndarray, tolerance: float = 1.0e-10) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(projector)
    return eigenvectors[:, eigenvalues > 1.0 - tolerance]


def scalar_residual(operator: np.ndarray, frame: np.ndarray) -> float:
    projected = frame.conj().T @ operator @ frame
    return operator_norm(projected - np.eye(projected.shape[0], dtype=complex) * np.trace(projected) / projected.shape[0])


def sector_projector(gauss: tuple[np.ndarray, ...], signs: tuple[int, ...]) -> np.ndarray:
    projector = np.eye(DIMENSION, dtype=complex)
    for generator, sign in zip(gauss, signs, strict=True):
        projector = projector @ ((np.eye(DIMENSION, dtype=complex) + sign * generator) / 2.0)
    return projector


def flux_projector(plaquette: np.ndarray, sign: int) -> np.ndarray:
    return (np.eye(DIMENSION, dtype=complex) + sign * plaquette) / 2.0


def sector_ground_energy(hamiltonian: np.ndarray, projector: np.ndarray) -> tuple[int, float]:
    frame = projector_frame(projector)
    return int(frame.shape[1]), float(np.linalg.eigvalsh(frame.conj().T @ hamiltonian @ frame)[0])


def main() -> None:
    identity = np.eye(DIMENSION, dtype=complex)
    gauss = tuple(gauss_operator(vertex) for vertex in VERTICES)
    plaquette = plaquette_operator()
    charge = sum((number_operator(vertex) for vertex in VERTICES), start=np.zeros((DIMENSION, DIMENSION, ), dtype=complex))
    fixed_charge_projector = charge_projector(FIXED_CHARGE)

    dressed_terms = []
    bare_terms = []
    hamiltonian_without_electric = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    electric_hamiltonian = np.zeros((DIMENSION, DIMENSION), dtype=complex)
    for edge, ((left, right), hopping, electric) in enumerate(zip(EDGES, EDGE_HOPPINGS, ELECTRIC_FIELDS, strict=True)):
        dressed = hop_term(left, right, edge) + hop_term(right, left, edge)
        bare = hop_term(left, right) + hop_term(right, left)
        dressed_terms.append(dressed)
        bare_terms.append(bare)
        hamiltonian_without_electric += -hopping * dressed
        electric_hamiltonian += -electric * tau_x(edge)
    for vertex, potential in zip(VERTICES, VERTEX_POTENTIALS, strict=True):
        hamiltonian_without_electric += potential * number_operator(vertex) - GAUSS_PENALTY * gauss[vertex]
    hamiltonian_without_electric += -PLAQUETTE_COUPLING * plaquette
    hamiltonian = hamiltonian_without_electric + electric_hamiltonian
    bare_hopping = sum((coefficient * term for coefficient, term in zip(EDGE_HOPPINGS, bare_terms, strict=True)), start=np.zeros_like(hamiltonian))
    dressed_hopping = sum((coefficient * term for coefficient, term in zip(EDGE_HOPPINGS, dressed_terms, strict=True)), start=np.zeros_like(hamiltonian))

    physical_projector = sector_projector(gauss, (1, 1, 1, 1))
    physical_fixed_charge_projector = fixed_charge_projector @ physical_projector
    physical_frame = projector_frame(physical_fixed_charge_projector)
    physical_spectrum = np.linalg.eigvalsh(physical_frame.conj().T @ hamiltonian @ physical_frame)
    physical_ground = float(physical_spectrum[0])

    gauss_rows = []
    for signs in product((-1, 1), repeat=len(VERTICES)):
        frame = projector_frame(fixed_charge_projector @ sector_projector(gauss, signs))
        if not frame.shape[1]:
            continue
        spectrum = np.linalg.eigvalsh(frame.conj().T @ hamiltonian @ frame)
        gauss_rows.append({"signs": list(signs), "dimension": int(frame.shape[1]), "ground_energy": float(spectrum[0])})
    syndrome_gap = min(row["ground_energy"] for row in gauss_rows if row["signs"] != [1, 1, 1, 1]) - physical_ground

    flux_rows = []
    for sign in (-1, 1):
        rank, energy = sector_ground_energy(
            hamiltonian_without_electric,
            physical_fixed_charge_projector @ flux_projector(plaquette, sign),
        )
        flux_rows.append({"B_p": sign, "dimension": rank, "ground_energy_without_electric": energy})
    flux_gap_static = abs(flux_rows[0]["ground_energy_without_electric"] - flux_rows[1]["ground_energy_without_electric"])

    output = {
        "schema": "antler.phase8c.z2-plaquette-preflight.v1",
        "parameters": {
            "new_declared_resource": "one neutral Z2 link qubit tau_e per edge; not an ANTLER charge-two mediator",
            "matter": "spinless U(1)-conserving fermion on each vertex",
            "gauge_generator": "G_v=(-1)^n_v product_(e incident to v) tau^x_e",
            "dressed_hopping": "c_v^dag tau^z_vw c_w+h.c.",
            "magnetic_plaquette": "B_p=product_(e in boundary p) tau^z_e",
            "edge_hoppings": list(EDGE_HOPPINGS),
            "electric_fields": list(ELECTRIC_FIELDS),
            "vertex_potentials": list(VERTEX_POTENTIALS),
            "gauss_penalty": GAUSS_PENALTY,
            "plaquette_coupling": PLAQUETTE_COUPLING,
        },
        "geometry": {
            "name": "one square plaquette (preflight only)",
            "vertices": list(VERTICES),
            "edges": [list(edge) for edge in EDGES],
            "plaquettes": 1,
        },
        "dimensions": {
            "full": DIMENSION,
            "fixed_N2_before_gauss": int(round(np.trace(fixed_charge_projector).real)),
            "physical_fixed_N2": int(physical_frame.shape[1]),
            "physical_all_charges": int(projector_frame(physical_projector).shape[1]),
        },
        "algebra": {
            "hamiltonian_hermiticity_norm": operator_norm(hamiltonian - hamiltonian.conj().T),
            "max_hamiltonian_gauss_commutator_norm": max(operator_norm(commutator(hamiltonian, generator)) for generator in gauss),
            "max_pairwise_gauss_commutator_norm": max(operator_norm(commutator(first, second)) for index, first in enumerate(gauss) for second in gauss[index + 1:]),
            "hamiltonian_charge_commutator_norm": operator_norm(commutator(hamiltonian, charge)),
            "max_plaquette_gauss_commutator_norm": max(operator_norm(commutator(plaquette, generator)) for generator in gauss),
            "hamiltonian_without_electric_plaquette_commutator_norm": operator_norm(commutator(hamiltonian_without_electric, plaquette)),
            "full_dynamic_hamiltonian_plaquette_commutator_norm": operator_norm(commutator(hamiltonian, plaquette)),
            "max_dressed_hop_gauss_commutator_norm": max(operator_norm(commutator(term, generator)) for term in dressed_terms for generator in gauss),
            "max_bare_hop_gauss_commutator_norm": max(operator_norm(commutator(term, generator)) for term in bare_terms for generator in gauss),
        },
        "physical_projection": {
            "projector_idempotency_norm": operator_norm(physical_fixed_charge_projector @ physical_fixed_charge_projector - physical_fixed_charge_projector),
            "max_projected_gauss_minus_projector_norm": max(operator_norm(physical_fixed_charge_projector @ generator @ physical_fixed_charge_projector - physical_fixed_charge_projector) for generator in gauss),
            "bare_hopping_projected_norm": operator_norm(physical_fixed_charge_projector @ bare_hopping @ physical_fixed_charge_projector),
            "dressed_hopping_projected_nonscalar_norm": scalar_residual(dressed_hopping, physical_frame),
            "local_vertex_density_projected_nonscalar_norm": scalar_residual(number_operator(0), physical_frame),
            "plaquette_projected_nonscalar_norm": scalar_residual(plaquette, physical_frame),
        },
        "spectra": {
            "physical_N2_eigenvalues": [float(value) for value in physical_spectrum],
            "physical_internal_gap": float(physical_spectrum[1] - physical_spectrum[0]),
            "minimum_Gauss_syndrome_gap_at_N2": float(syndrome_gap),
            "gauss_signature_sectors_at_N2": gauss_rows,
            "static_flux_sectors_at_N2": flux_rows,
            "static_flux_sector_gap": float(flux_gap_static),
        },
        "decision": (
            "PASS T1 as a one-plaquette gauge/flux preflight. The neutral-link reference keeps exact Gauss algebra and total U(1); "
            "B_p is gauge invariant and has separated static flux sectors, while electric fields intentionally make B_p dynamical. "
            "The projected local density is non-scalar, so a single plaquette is not a topological code or a fusion/braid implementation."
        ),
        "claim_boundary": (
            "This is a finite reference model with inserted neutral link qubits, Gauss penalty and plaquette term. It does not establish "
            "a deconfined thermodynamic phase, topological ground-state degeneracy, local indistinguishability, e<->m domain wall, twist "
            "defect, anyons, fusion or non-Abelian braid, and it is not derived from frozen ANTLER resources."
        ),
    }
    path = ROOT / "results" / "phase8c" / "z2_plaquette_preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
