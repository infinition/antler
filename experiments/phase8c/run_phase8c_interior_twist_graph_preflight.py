"""Phase 8C-T5a: interior graph-twist stabilizer preflight.

Sarkar and Yoder's graph construction associates qubits to vertices and face
checks to an embedded graph; odd-degree vertices are twist defects.  This
script constructs an explicit periodic square-graph reference with two
disjoint deleted edges.  The four endpoints have degree three, while every
other vertex retains degree four.  A cyclically anticommuting local Pauli list
(`X,Z,X,Z` at degree four and `X,Y,Z` at degree three) turns each oriented
face into a commuting stabilizer check.

The result is deliberately an *external vertex-qubit stabilizer reference*.
It is not the neutral-link gauge Hamiltonian of Phase 8C, is not derived from
the frozen ANTLER microscopic model, and contains neither a prescribed braid
matrix nor a physical defect-motion protocol.  Its purpose is narrower: make
the first interior-twist geometry and its code-space gates reproducible before
any attempted holonomy.
"""
from __future__ import annotations

import json
import sys
from collections import deque
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_stabilizer_algebra import (  # noqa: E402
    BinaryPauli,
    _gf2_basis,
    _packed,
    commute,
    gf2_in_span,
    gf2_rank,
    local_paulis,
)


DIRECTIONS = ("E", "N", "W", "S")
DELTA = {"E": (1, 0), "N": (0, 1), "W": (-1, 0), "S": (0, -1)}
OPPOSITE = {"E": "W", "W": "E", "N": "S", "S": "N"}
SIZES = (4, 6)


def neighbour(vertex: tuple[int, int], direction: str, size: int) -> tuple[int, int]:
    dx, dy = DELTA[direction]
    return ((vertex[0] + dx) % size, (vertex[1] + dy) % size)


def edge(left: tuple[int, int], right: tuple[int, int]) -> frozenset[tuple[int, int]]:
    return frozenset((left, right))


def edge_label(item: frozenset[tuple[int, int]]) -> list[list[int]]:
    return [list(vertex) for vertex in sorted(item)]


def all_grid_edges(size: int) -> set[frozenset[tuple[int, int]]]:
    vertices = [(x, y) for y in range(size) for x in range(size)]
    return {
        edge(vertex, neighbour(vertex, direction, size))
        for vertex in vertices
        for direction in ("E", "N")
    }


def direction_to(left: tuple[int, int], right: tuple[int, int], size: int) -> str:
    return next(direction for direction in DIRECTIONS if neighbour(left, direction, size) == right)


def pauli_label(word: BinaryPauli, qubits: int) -> str:
    symbols = []
    for qubit in range(qubits):
        x_bit = (word.x >> qubit) & 1
        z_bit = (word.z >> qubit) & 1
        symbols.append("Y" if x_bit and z_bit else "X" if x_bit else "Z" if z_bit else "I")
    return "".join(symbols)


def parse_pauli(label: str) -> BinaryPauli:
    return BinaryPauli(
        x=sum(1 << index for index, symbol in enumerate(label) if symbol in {"X", "Y"}),
        z=sum(1 << index for index, symbol in enumerate(label) if symbol in {"Z", "Y"}),
    )


def local_sector_pauli(degree: int, sector_index: int) -> str:
    """One cyclically anticommuting list for each allowed vertex degree."""
    if degree == 4:
        return "XZXZ"[sector_index]
    if degree == 3:
        return "XYZ"[sector_index]
    raise ValueError(f"unsupported degree {degree}; this preflight requires only degree three or four")


def graph_code(
    size: int,
    removed_edges: set[frozenset[tuple[int, int]]],
    *,
    full_local_audit: bool = True,
) -> dict[str, object]:
    vertices = [(x, y) for y in range(size) for x in range(size)]
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    active_edges = all_grid_edges(size) - removed_edges
    adjacency: dict[tuple[int, int], list[str]] = {vertex: [] for vertex in vertices}
    darts: dict[tuple[tuple[int, int], str], tuple[int, int]] = {}
    for item in active_edges:
        left, right = tuple(item)
        for origin, target in ((left, right), (right, left)):
            direction = direction_to(origin, target, size)
            adjacency[origin].append(direction)
            darts[(origin, direction)] = target
    for directions in adjacency.values():
        directions.sort(key=DIRECTIONS.index)
    if any(len(directions) not in {3, 4} for directions in adjacency.values()):
        raise RuntimeError("the graph contains a degree outside the registered interior-twist grammar")

    # Trace each oriented face.  At the arrival vertex, selecting the preceding
    # direction in cyclic order keeps the face on the left of the traversal.
    faces: list[list[tuple[tuple[int, int], int]]] = []
    visited: set[tuple[tuple[int, int], str]] = set()
    for start in darts:
        if start in visited:
            continue
        current = start
        face: list[tuple[tuple[int, int], int]] = []
        while current not in visited:
            visited.add(current)
            _, direction = current
            arrival = darts[current]
            incoming = OPPOSITE[direction]
            ordered = adjacency[arrival]
            sector_index = ordered.index(incoming)
            outgoing = ordered[(sector_index - 1) % len(ordered)]
            face.append((arrival, sector_index))
            current = (arrival, outgoing)
        faces.append(face)

    checks: list[BinaryPauli] = []
    face_metadata: list[dict[str, object]] = []
    for face in faces:
        x_mask = 0
        z_mask = 0
        sectors = []
        for vertex, sector_index in face:
            local_pauli = local_sector_pauli(len(adjacency[vertex]), sector_index)
            qubit = vertex_index[vertex]
            if local_pauli in {"X", "Y"}:
                x_mask |= 1 << qubit
            if local_pauli in {"Z", "Y"}:
                z_mask |= 1 << qubit
            sectors.append({"vertex": list(vertex), "sector_index": sector_index, "pauli": local_pauli})
        word = BinaryPauli(x=x_mask, z=z_mask)
        checks.append(word)
        face_metadata.append({"weight": word.weight(), "sectors": sectors, "check": pauli_label(word, size * size)})

    packed_checks = tuple(_packed(check, size * size) for check in checks)
    basis = _gf2_basis(packed_checks)
    rank = gf2_rank(packed_checks)
    pair_anticommutations = sum(
        not commute(left, right)
        for index, left in enumerate(checks)
        for right in checks[index + 1:]
    )
    minimum_logical: BinaryPauli | None = None
    minimum_weight: int | None = None
    local_counts: dict[str, int | bool] | None = None
    if full_local_audit:
        for weight in range(1, size * size + 1):
            for candidate in local_paulis(size * size, weight):
                if all(commute(candidate, check) for check in checks) and not gf2_in_span(_packed(candidate, size * size), basis):
                    minimum_logical = candidate
                    minimum_weight = weight
                    break
            if minimum_logical is not None:
                break
        if minimum_logical is None or minimum_weight is None:
            raise RuntimeError("no logical Pauli found")

        local_counts = {"projects_to_zero": 0, "stabilizer_scalar": 0, "nontrivial_logical": 0}
        for candidate in local_paulis(size * size, minimum_weight - 1):
            if any(not commute(candidate, check) for check in checks):
                local_counts["projects_to_zero"] += 1
            elif gf2_in_span(_packed(candidate, size * size), basis):
                local_counts["stabilizer_scalar"] += 1
            else:
                local_counts["nontrivial_logical"] += 1
        local_counts["tested_nonidentity_paulis"] = sum(int(value) for value in local_counts.values())
        local_counts["all_below_distance_probes_scalar_or_zero"] = local_counts["nontrivial_logical"] == 0

    degree_three = sorted((vertex for vertex, directions in adjacency.items() if len(directions) == 3), key=lambda item: (item[1], item[0]))
    return {
        "size": size,
        "qubits": size * size,
        "active_edges": len(active_edges),
        "face_count": len(faces),
        "face_weights": [metadata["weight"] for metadata in face_metadata],
        "face_checks": face_metadata,
        "degree_three_interior_twists": [list(vertex) for vertex in degree_three],
        "degree_four_vertices": sum(len(directions) == 4 for directions in adjacency.values()),
        "pair_anticommutations": pair_anticommutations,
        "independent_rank_over_GF2": rank,
        "encoded_qubits": size * size - rank,
        "ground_space_degeneracy": 1 << (size * size - rank),
        "exact_syndrome_gap_in_units_of_J": 2.0,
        "minimum_logical_weight": minimum_weight,
        "one_minimum_logical_representative": pauli_label(minimum_logical, size * size) if minimum_logical is not None else None,
        "local_protection": local_counts,
        "check_labels": [pauli_label(check, size * size) for check in checks],
    }


def graph_distance(size: int, origin: tuple[int, int], target: tuple[int, int]) -> int:
    frontier: deque[tuple[int, int]] = deque((origin,))
    distance = {origin: 0}
    while frontier:
        vertex = frontier.popleft()
        if vertex == target:
            return distance[vertex]
        for direction in DIRECTIONS:
            next_vertex = neighbour(vertex, direction, size)
            if next_vertex not in distance:
                distance[next_vertex] = distance[vertex] + 1
                frontier.append(next_vertex)
    raise RuntimeError("periodic grid is disconnected")


def deformation_comparison(size: int) -> dict[str, object]:
    """Compare two static graph codes differing by a one-edge local mutation."""
    fixed_pair = edge((size // 2, size // 2), (size // 2 + 1, size // 2))
    initial_edge = edge((0, 0), (1, 0))
    mutated_edge = edge((0, 0), (0, 1))
    initial = graph_code(size, {initial_edge, fixed_pair})
    final = graph_code(size, {mutated_edge, fixed_pair})
    initial_checks = set(initial["check_labels"])
    final_checks = set(final["check_labels"])
    changed_labels = (initial_checks - final_checks) | (final_checks - initial_checks)
    changed_qubits = sorted({index for label in changed_labels for index, symbol in enumerate(label) if symbol != "I"})
    changed_vertices = [(index % size, index // size) for index in changed_qubits]
    initial_paulis = tuple(parse_pauli(label) for label in initial["check_labels"])
    final_paulis = tuple(parse_pauli(label) for label in final["check_labels"])
    cross_anticommutations = sum(not commute(left, right) for left in initial_paulis for right in final_paulis)
    return {
        "initial_removed_edges": [edge_label(initial_edge), edge_label(fixed_pair)],
        "final_removed_edges": [edge_label(mutated_edge), edge_label(fixed_pair)],
        "initial_twists": initial["degree_three_interior_twists"],
        "final_twists": final["degree_three_interior_twists"],
        "shared_checks": len(initial_checks & final_checks),
        "initial_only_checks": len(initial_checks - final_checks),
        "final_only_checks": len(final_checks - initial_checks),
        "cross_configuration_anticommutations": cross_anticommutations,
        "changed_support_qubits": changed_qubits,
        "changed_support_vertices": [list(vertex) for vertex in changed_vertices],
        "maximum_periodic_graph_distance_from_mutation_pivot": max(graph_distance(size, (0, 0), vertex) for vertex in changed_vertices),
        "static_initial_code": {key: initial[key] for key in ("encoded_qubits", "ground_space_degeneracy", "minimum_logical_weight")},
        "static_final_code": {key: final[key] for key in ("encoded_qubits", "ground_space_degeneracy", "minimum_logical_weight")},
    }


def main() -> None:
    rows = []
    base_controls = []
    for size in SIZES:
        deleted = {
            edge((0, 0), (1, 0)),
            edge((size // 2, size // 2), (size // 2 + 1, size // 2)),
        }
        row = graph_code(size, deleted)
        base = graph_code(size, set(), full_local_audit=False)
        if not (
            row["pair_anticommutations"] == 0
            and len(row["degree_three_interior_twists"]) == 4
            and row["degree_four_vertices"] == size * size - 4
            and row["encoded_qubits"] == 3
            and row["ground_space_degeneracy"] == 8
            and row["minimum_logical_weight"] >= 3
            and row["local_protection"]["all_below_distance_probes_scalar_or_zero"]
            and base["encoded_qubits"] == 2
        ):
            raise RuntimeError(f"interior-twist graph preflight failed at size {size}")
        row["removed_edges"] = [edge_label(item) for item in sorted(deleted, key=edge_label)]
        rows.append(row)
        base_controls.append({
            "size": size,
            "undislocated_encoded_qubits": base["encoded_qubits"],
            "dislocated_encoded_qubits": row["encoded_qubits"],
            "extra_encoded_qubits": row["encoded_qubits"] - base["encoded_qubits"],
        })

    mutation = deformation_comparison(6)
    if not (
        mutation["initial_twists"] != mutation["final_twists"]
        and mutation["initial_only_checks"] == mutation["final_only_checks"] == 5
        and mutation["cross_configuration_anticommutations"] > 0
        and mutation["static_initial_code"] == mutation["static_final_code"]
    ):
        raise RuntimeError("the registered local graph mutation is not a nontrivial equal-code-parameter deformation candidate")

    output = {
        "schema": "antler.phase8c.interior-twist-graph-preflight.v1",
        "parameters": {
            "reference": "Sarkar-Yoder graph surface-code construction: vertex qubits, face stabilizers, odd-degree vertices as twists",
            "new_declared_resource": "external neutral vertex-qubit stabilizer graph; separate from Phase 8C edge links and not derived from ANTLER",
            "local_sector_lists": {"degree_4": ["X", "Z", "X", "Z"], "degree_3": ["X", "Y", "Z"]},
            "periodic_sizes": list(SIZES),
            "construction": "two disjoint deleted square-grid edges merge adjacent faces and leave four degree-three interior vertices",
        },
        "base_torus_controls": base_controls,
        "rows": rows,
        "one_edge_mutation_candidate": mutation,
        "decision": (
            "PASS T5a as an external interior-twist graph-code preflight. Each registered periodic graph has four degree-three "
            "interior vertices, a commuting face-check parent, an eight-dimensional ground space (one extra encoded qubit relative "
            "to the undislocated torus), a nonzero syndrome gap and no logical Pauli below its measured distance. The changed static "
            "graphs have equal code parameters but noncommuting changed checks, so a defect-motion unitary or measurement protocol "
            "has not been inferred from their mere existence."
        ),
        "next_gate": (
            "T5b must derive a local gapped or measurement-resolved code-deformation schedule for the registered one-edge mutation, "
            "track the three logical Pauli pairs through it, and independently audit leakage and the instantaneous gap. Only after two "
            "adjacent deformations yield a nonzero commutator may a Yang-Baxter residual be interpreted."
        ),
        "claim_boundary": (
            "This is an imposed graph-stabilizer reference, not a neutral-link ANTLER Hamiltonian and not a microscopic resource "
            "derivation. It does not prove physical e<->m transport, a fusion measurement, mobile-defect adiabaticity, a braid matrix, "
            "non-Abelian statistics, universality, experimental feasibility or fault tolerance."
        ),
    }
    result = ROOT / "results" / "phase8c" / "interior_twist_graph_preflight.json"
    result.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
