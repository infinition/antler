"""Phase 8C-T3a: distinguish an e<->m duality from false branch cuts.

The pure-gauge 3x3 torus from T2 has an exact global self-duality: a lattice
duality permutation followed by Hadamard on every neutral link maps stars to
plaquettes and Z strings to X strings.  That is a *calibration* of the desired
anyon permutation, not a local domain wall.

This audit then rejects two tempting but insufficient local candidates:
  * a static pi/sign cocycle, which changes no Pauli type; and
  * a bounded partial-Hadamard change of basis, which creates mixed checks but
    maps no star support to a plaquette support.

No twist, fusion space, or braid is inserted or claimed.  A future T3b must
give a local dislocation/check geometry whose interface realizes the calibrated
map and has physical endpoints.
"""
from __future__ import annotations

import json
import sys
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
    local_paulis,
    toric_code_geometry,
)


LX = 3
LY = 3


def horizontal(x: int, y: int) -> int:
    return (y % LY) * LX + (x % LX)


def vertical(x: int, y: int) -> int:
    return LX * LY + (y % LY) * LX + (x % LX)


def dual_edge_permutation() -> tuple[int, ...]:
    """Map each primal edge to the intersecting dual edge on the same torus.

    h(x,y) -> v(x,y-1), and v(x,y) -> h(x-1,y).  Together with H on every
    link, this maps an X star to a Z plaquette and conversely.
    """
    permutation = [0] * (2 * LX * LY)
    for y in range(LY):
        for x in range(LX):
            permutation[horizontal(x, y)] = vertical(x, y - 1)
            permutation[vertical(x, y)] = horizontal(x - 1, y)
    if sorted(permutation) != list(range(2 * LX * LY)):
        raise RuntimeError("duality edge map is not a permutation")
    return tuple(permutation)


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for source, target in enumerate(permutation):
        if (mask >> source) & 1:
            result |= 1 << target
    return result


def global_duality(pauli: BinaryPauli, permutation: tuple[int, ...]) -> BinaryPauli:
    """Dual lattice permutation followed by H on every link (X <-> Z)."""
    return BinaryPauli(
        x=permute_mask(pauli.z, permutation),
        z=permute_mask(pauli.x, permutation),
    )


def partial_hadamard(pauli: BinaryPauli, region_mask: int) -> BinaryPauli:
    """Conjugate a word by H only on a bounded collection of links."""
    x_inside, z_inside = pauli.x & region_mask, pauli.z & region_mask
    return BinaryPauli(
        x=(pauli.x & ~region_mask) | z_inside,
        z=(pauli.z & ~region_mask) | x_inside,
    )


def label(labels: tuple[str, ...], pauli: BinaryPauli) -> list[str]:
    out = []
    for edge, name in enumerate(labels):
        x_bit = (pauli.x >> edge) & 1
        z_bit = (pauli.z >> edge) & 1
        if x_bit or z_bit:
            out.append(f"{'Y' if x_bit and z_bit else 'X' if x_bit else 'Z'}[{name}]")
    return out


def is_logical(pauli: BinaryPauli, stabilizers: tuple[BinaryPauli, ...], basis: dict[int, int], n_qubits: int) -> bool:
    return all(commute(pauli, stabilizer) for stabilizer in stabilizers) and not gf2_in_span(
        _packed(pauli, n_qubits), basis
    )


def first_z_logical(geometry, basis: dict[int, int]) -> BinaryPauli:
    for candidate in local_paulis(geometry.n_edges, min(LX, LY)):
        if candidate.x == 0 and is_logical(candidate, geometry.stabilizers, basis, geometry.n_edges):
            return candidate
    raise RuntimeError("expected a minimum-weight Z logical")


def main() -> None:
    geometry = toric_code_geometry(LX, LY)
    stabilizers = geometry.stabilizers
    basis = _gf2_basis(tuple(_packed(stabilizer, geometry.n_edges) for stabilizer in stabilizers))
    permutation = dual_edge_permutation()

    star_words = set(geometry.stars)
    plaquette_words = set(geometry.plaquettes)
    mapped_stars = tuple(global_duality(star, permutation) for star in geometry.stars)
    mapped_plaquettes = tuple(global_duality(plaquette, permutation) for plaquette in geometry.plaquettes)
    star_to_plaquette = sum(word in plaquette_words for word in mapped_stars)
    plaquette_to_star = sum(word in star_words for word in mapped_plaquettes)

    z_logical = first_z_logical(geometry, basis)
    dual_z_logical = global_duality(z_logical, permutation)

    # A static pi cocycle only decorates coefficients. It has no action on a
    # binary Pauli word, so neither X nor Z support can become the other type.
    static_pi_type_change_count = 0

    # This is deliberately bounded: one column of horizontal links, not a
    # noncontractible duality transformation. Local H keeps commutation but
    # cannot transform the support of a vertex star into a plaquette support.
    bounded_region = sum(1 << horizontal(0, y) for y in range(LY))
    bounded_stars = tuple(partial_hadamard(star, bounded_region) for star in geometry.stars)
    bounded_plaquettes = tuple(partial_hadamard(plaquette, bounded_region) for plaquette in geometry.plaquettes)
    bounded_all = bounded_stars + bounded_plaquettes
    bounded_anticommutations = sum(
        not commute(first, second)
        for index, first in enumerate(bounded_all)
        for second in bounded_all[index + 1:]
    )
    bounded_star_to_pure_plaquette = sum(word in plaquette_words for word in bounded_stars)
    bounded_plaquette_to_pure_star = sum(word in star_words for word in bounded_plaquettes)
    bounded_mixed_checks = sum(bool(word.x and word.z) for word in bounded_all)

    global_pass = (
        star_to_plaquette == len(geometry.stars)
        and plaquette_to_star == len(geometry.plaquettes)
        and is_logical(dual_z_logical, stabilizers, basis, geometry.n_edges)
        and dual_z_logical.x != 0
        and dual_z_logical.z == 0
    )
    if not global_pass:
        raise RuntimeError("global e<->m duality calibration failed")
    if bounded_anticommutations != 0:
        raise RuntimeError("partial Clifford transformation must preserve stabilizer commutation")

    output = {
        "schema": "antler.phase8c.em-duality-wall-preflight.v1",
        "geometry": {
            "name": "3x3 pure-gauge square torus used for T2",
            "link_qubits": geometry.n_edges,
            "stars": len(geometry.stars),
            "plaquettes": len(geometry.plaquettes),
            "boundary": "periodic",
        },
        "global_duality_calibration": {
            "map": "D = H_on_all_links after h(x,y)->v(x,y-1), v(x,y)->h(x-1,y)",
            "edge_permutation": list(permutation),
            "stars_mapped_exactly_to_plaquettes": star_to_plaquette,
            "plaquettes_mapped_exactly_to_stars": plaquette_to_star,
            "z_string_before": label(geometry.edge_labels, z_logical),
            "image_of_z_string": label(geometry.edge_labels, dual_z_logical),
            "image_is_nontrivial_x_logical": bool(
                is_logical(dual_z_logical, stabilizers, basis, geometry.n_edges)
                and dual_z_logical.x != 0
                and dual_z_logical.z == 0
            ),
            "passes_exact_em_permutation_calibration": global_pass,
            "scope": "global lattice duality; not a spatially bounded wall and has no endpoints",
        },
        "false_wall_controls": {
            "static_pi_cocycle": {
                "pauli_type_change_count": static_pi_type_change_count,
                "verdict": "FAIL: signs/phases alter coefficients but do not map Z-string support to an X-string",
            },
            "bounded_partial_hadamard": {
                "region": "three horizontal links h(0,y), y=0,1,2",
                "checks_still_commute": bounded_anticommutations == 0,
                "mixed_check_count": bounded_mixed_checks,
                "stars_mapped_to_pure_plaquettes": bounded_star_to_pure_plaquette,
                "plaquettes_mapped_to_pure_stars": bounded_plaquette_to_pure_star,
                "verdict": "FAIL: a bounded basis change creates mixed checks but no star/plaquette support exchange",
            },
        },
        "decision": (
            "T3a PASS as a discriminating calibration only: the 3x3 neutral-link code has an exact global e<->m symmetry, "
            "and the registered static-pi and bounded-mixed-check controls fail the support-changing wall criterion. "
            "T3 itself remains open because no local dislocation/domain-wall Hamiltonian with endpoints has been constructed."
        ),
        "next_gate": (
            "T3b must give an explicit local cellulation/check set with a finite e<->m branch cut, endpoint checks and string "
            "operators whose syndrome labels are exchanged on crossing. It must then re-run commutation, rank, gap and local "
            "indistinguishability audits before any fusion or braid calculation."
        ),
        "claim_boundary": (
            "This is an algebraic symmetry calibration on an imposed Abelian pure-gauge reference. It neither constructs a local "
            "domain wall nor establishes twist defects, a fusion space, defect transport, non-Abelian braiding, universality, a "
            "microscopic ANTLER derivation or experimental feasibility."
        ),
    }
    path = ROOT / "results" / "phase8c" / "em_duality_wall_preflight.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
