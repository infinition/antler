"""Exact binary-stabilizer checks for the Phase 7B two-dimensional control.

This module is deliberately algebraic.  It verifies the code properties of an
*abstract* square-lattice toric-code parent after each physical ANTLER edge has
been frozen to a one-particle rung qubit.  It neither supplies a low-body
mediator realization nor claims one follows from the original two-rail ladder.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product


@dataclass(frozen=True)
class BinaryPauli:
    """A Pauli modulo phase, encoded by its X and Z bit masks."""

    x: int
    z: int

    def weight(self) -> int:
        return (self.x | self.z).bit_count()


@dataclass(frozen=True)
class ToricCodeGeometry:
    """Square torus with qubits on oriented horizontal and vertical edges."""

    lx: int
    ly: int
    edge_labels: tuple[str, ...]
    stars: tuple[BinaryPauli, ...]
    plaquettes: tuple[BinaryPauli, ...]

    @property
    def n_edges(self) -> int:
        return len(self.edge_labels)

    @property
    def stabilizers(self) -> tuple[BinaryPauli, ...]:
        return self.stars + self.plaquettes


def symplectic_parity(left: BinaryPauli, right: BinaryPauli) -> int:
    """Return one exactly when two binary Paulis anticommute."""
    return (((left.x & right.z).bit_count() + (left.z & right.x).bit_count()) & 1)


def commute(left: BinaryPauli, right: BinaryPauli) -> bool:
    return symplectic_parity(left, right) == 0


def _packed(pauli: BinaryPauli, n_qubits: int) -> int:
    return pauli.x | (pauli.z << n_qubits)


def _gf2_basis(rows: tuple[int, ...] | list[int]) -> dict[int, int]:
    """Row-echelon GF(2) basis indexed by most-significant pivot bit."""
    basis: dict[int, int] = {}
    for row in rows:
        vector = int(row)
        while vector:
            pivot = vector.bit_length() - 1
            old = basis.get(pivot)
            if old is None:
                basis[pivot] = vector
                break
            vector ^= old
    return basis


def gf2_rank(rows: tuple[int, ...] | list[int]) -> int:
    return len(_gf2_basis(rows))


def gf2_in_span(vector: int, basis: dict[int, int]) -> bool:
    """Test membership in a row space represented by ``_gf2_basis``."""
    remainder = int(vector)
    while remainder:
        pivot = remainder.bit_length() - 1
        old = basis.get(pivot)
        if old is None:
            return False
        remainder ^= old
    return True


def toric_code_geometry(lx: int, ly: int) -> ToricCodeGeometry:
    """Return the conventional toric-code stabilizers for an ``lx`` by ``ly`` torus.

    Values below three are refused: a two-site periodic direction identifies
    nominally distinct edges and is not a faithful square-lattice control.
    """
    if lx < 3 or ly < 3:
        raise ValueError("toric-code preflight requires lx, ly >= 3")

    def horizontal(x: int, y: int) -> int:
        return (y % ly) * lx + (x % lx)

    def vertical(x: int, y: int) -> int:
        return lx * ly + (y % ly) * lx + (x % lx)

    labels = tuple(
        [f"h({x},{y})" for y in range(ly) for x in range(lx)]
        + [f"v({x},{y})" for y in range(ly) for x in range(lx)]
    )
    stars: list[BinaryPauli] = []
    plaquettes: list[BinaryPauli] = []
    for y in range(ly):
        for x in range(lx):
            star_edges = (
                horizontal(x, y), horizontal(x - 1, y),
                vertical(x, y), vertical(x, y - 1),
            )
            stars.append(BinaryPauli(sum(1 << edge for edge in star_edges), 0))
            plaquette_edges = (
                horizontal(x, y), vertical(x + 1, y),
                horizontal(x, y + 1), vertical(x, y),
            )
            plaquettes.append(BinaryPauli(0, sum(1 << edge for edge in plaquette_edges)))
    return ToricCodeGeometry(lx, ly, labels, tuple(stars), tuple(plaquettes))


def _pauli_name(geometry: ToricCodeGeometry, pauli: BinaryPauli) -> list[str]:
    out = []
    for edge, label in enumerate(geometry.edge_labels):
        x_bit, z_bit = (pauli.x >> edge) & 1, (pauli.z >> edge) & 1
        if x_bit or z_bit:
            out.append(f"{'Y' if x_bit and z_bit else 'X' if x_bit else 'Z'}[{label}]")
    return out


def local_paulis(n_qubits: int, maximum_weight: int):
    """Yield all nonidentity Paulis through a specified physical-edge weight."""
    if maximum_weight < 1:
        return
    for weight in range(1, maximum_weight + 1):
        for support in combinations(range(n_qubits), weight):
            for factors in product(((1, 0), (0, 1), (1, 1)), repeat=weight):
                x = sum(x_bit << edge for edge, (x_bit, _) in zip(support, factors))
                z = sum(z_bit << edge for edge, (_, z_bit) in zip(support, factors))
                yield BinaryPauli(x, z)


def _first_nontrivial_logical(
    geometry: ToricCodeGeometry,
    stabilizer_basis: dict[int, int],
) -> BinaryPauli:
    """Find a minimum-weight centralizer element outside the stabilizer span."""
    for weight in range(1, min(geometry.lx, geometry.ly) + 1):
        for candidate in local_paulis(geometry.n_edges, weight):
            if all(commute(candidate, stabilizer) for stabilizer in geometry.stabilizers):
                if not gf2_in_span(_packed(candidate, geometry.n_edges), stabilizer_basis):
                    return candidate
    raise RuntimeError("no logical found through the expected toric-code distance")


def toric_code_preflight(lx: int = 3, ly: int = 3) -> dict:
    """Certify the exact abstract 2D control and its local-Pauli code gate.

    In the charge-frozen ANTLER edge embedding, every projected physical local
    charge-conserving operator is a linear combination of these Pauli probes.
    Hence a probe that anticommutes with a stabilizer has zero code projection,
    and a stabilizer product has a scalar projection.  The exhaustive search
    below verifies that no non-scalar logical action occurs below the code
    distance.
    """
    geometry = toric_code_geometry(lx, ly)
    stabilizer_rows = tuple(_packed(stabilizer, geometry.n_edges) for stabilizer in geometry.stabilizers)
    stabilizer_basis = _gf2_basis(stabilizer_rows)
    pair_anticommutations = sum(
        not commute(left, right)
        for position, left in enumerate(geometry.stabilizers)
        for right in geometry.stabilizers[position + 1:]
    )
    rank = gf2_rank(stabilizer_rows)
    encoded_qubits = geometry.n_edges - rank
    logical = _first_nontrivial_logical(geometry, stabilizer_basis)
    distance = logical.weight()
    rows = {"anticommutes_with_stabilizer": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
    for candidate in local_paulis(geometry.n_edges, distance - 1):
        if not all(commute(candidate, stabilizer) for stabilizer in geometry.stabilizers):
            rows["anticommutes_with_stabilizer"] += 1
        elif gf2_in_span(_packed(candidate, geometry.n_edges), stabilizer_basis):
            rows["stabilizer_scalar"] += 1
        else:
            rows["nontrivial_logical"] += 1
    checked = sum(rows.values())
    return {
        "geometry": {"lx": lx, "ly": ly, "boundary": "torus", "physical_edge_qubits": geometry.n_edges},
        "stabilizer_algebra": {
            "star_count": len(geometry.stars),
            "plaquette_count": len(geometry.plaquettes),
            "independent_stabilizer_rank_over_GF2": rank,
            "pair_anticommutations": pair_anticommutations,
            "all_stabilizers_commute": pair_anticommutations == 0,
        },
        "code": {
            "encoded_qubits": encoded_qubits,
            "ground_space_degeneracy": 1 << encoded_qubits,
            "minimum_logical_weight": distance,
            "one_minimum_logical": _pauli_name(geometry, logical),
        },
        "complete_projected_local_pauli_gate": {
            "maximum_tested_weight": distance - 1,
            "tested_nonidentity_paulis": checked,
            **rows,
            "all_tested_probes_project_to_scalars": rows["nontrivial_logical"] == 0,
        },
    }
