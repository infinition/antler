"""Tile the audited Phase-8B walker coefficients on the 3x3 toric reference.

The local overlap audit produces, after downfolding, only I, A_s, B_p and the
commuting product A_s B_p.  This script asks a narrowly scoped question:
assuming independent copies of that already-audited local gadget can be tiled,
does the *resulting effective commuting Hamiltonian* retain the fourfold torus
ground space and a nonzero syndrome gap?

It enumerates all 2^16 independent stabilizer-syndrome assignments rather
than diagonalising 2^18 edge states.  It is not a microscopic multi-walker
tiling simulation; inter-gadget crosstalk is deliberately left open.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.phase7_stabilizer_algebra import toric_code_geometry, toric_code_preflight


def parity_constrained_signs(count: int) -> list[np.ndarray]:
    """All sign arrays with product +1, enforcing the torus stabilizer relation."""
    arrays = []
    for bits in itertools.product((-1.0, 1.0), repeat=count - 1):
        last = float(np.prod(bits))
        arrays.append(np.asarray((*bits, last), dtype=float))
    return arrays


def main() -> None:
    local_path = ROOT / "results" / "phase7" / "phase8b_star_plaquette_walker_overlap_audit.json"
    local = json.loads(local_path.read_text(encoding="utf-8"))
    geometry = toric_code_geometry(lx=3, ly=3)
    stars, plaquettes = geometry.stars, geometry.plaquettes
    incidences = [
        (star_index, plaquette_index)
        for star_index, star in enumerate(stars)
        for plaquette_index, plaquette in enumerate(plaquettes)
        if (star.x & plaquette.z).bit_count() == 2
    ]
    if len(incidences) != 36:
        raise RuntimeError(f"expected 36 toric star-plaquette incidences, got {len(incidences)}")
    sign_arrays = parity_constrained_signs(9)
    rows = []
    for local_row in local["rows"]:
        c_star = float(local_row["star_coefficient"])
        c_plaquette = float(local_row["plaquette_coefficient"])
        # A_s B_p = -YYXXZZ in the local Z-basis convention.
        c_product = -float(local_row["star_plaquette_product_coefficient"])
        energies: list[tuple[float, np.ndarray, np.ndarray]] = []
        for star_signs in sign_arrays:
            for plaquette_signs in sign_arrays:
                energy = (
                    c_star * float(np.sum(star_signs))
                    + c_plaquette * float(np.sum(plaquette_signs))
                    + c_product * sum(star_signs[star] * plaquette_signs[plaquette]
                                      for star, plaquette in incidences)
                )
                energies.append((float(energy), star_signs, plaquette_signs))
        energies.sort(key=lambda item: item[0])
        ground_energy = energies[0][0]
        tolerance = 1e-12 * max(1.0, abs(ground_energy))
        ground = [item for item in energies if abs(item[0] - ground_energy) <= tolerance]
        first_excited = next(item[0] for item in energies if item[0] > ground_energy + tolerance)
        all_plus_ground = all(
            np.all(item[1] == 1.0) and np.all(item[2] == 1.0) for item in ground
        )
        rows.append({
            "coupling_over_detuning": local_row["coupling_over_detuning"],
            "star_coefficient": c_star,
            "plaquette_coefficient": c_plaquette,
            "star_plaquette_product_coefficient_in_stabilizer_basis": c_product,
            "product_to_star_ratio": float(abs(c_product / c_star)),
            "ground_syndrome_degeneracy": len(ground),
            "effective_toric_ground_space_degeneracy": 4 * len(ground),
            "all_plus_syndrome_is_unique_ground": bool(all_plus_ground and len(ground) == 1),
            "syndrome_gap": float(first_excited - ground_energy),
        })
    deep = [row for row in rows if float(row["coupling_over_detuning"]) <= 0.075]
    gap_power = float(np.polyfit(
        np.log([float(row["coupling_over_detuning"]) for row in deep]),
        np.log([float(row["syndrome_gap"]) for row in deep]),
        1,
    )[0])
    reference = toric_code_preflight(lx=3, ly=3)
    output = {
        "schema": "antler.phase8b.effective-toric-tile-closure.v1",
        "parameters": {
            "geometry": "3x3 square torus",
            "edges": geometry.n_edges,
            "stars": len(stars),
            "plaquettes": len(plaquettes),
            "star_plaquette_incidences": len(incidences),
            "coefficients_source": "phase8b_star_plaquette_walker_overlap_audit.json",
            "effective_parent": "sum_s c_A A_s + sum_p c_B B_p + sum_(s,p incident) c_AB A_s B_p",
        },
        "rows": rows,
        "deep_sw_syndrome_gap_power": gap_power,
        "inherited_exact_stabilizer_reference": {
            "independent_stabilizer_rank": reference["stabilizer_algebra"]["independent_stabilizer_rank_over_GF2"],
            "encoded_qubits": reference["code"]["encoded_qubits"],
            "minimum_logical_weight": reference["code"]["minimum_logical_weight"],
            "all_weight_below_distance_probes_scalar_or_null": reference["complete_projected_local_pauli_gate"]["all_tested_probes_project_to_scalars"],
        },
        "decision": (
            "Conditioned on independent local-gadget tiling, the measured fourth-order star/plaquette coefficients and "
            "eighth-order commuting product retain a unique all-plus syndrome, a fourfold torus ground space and a nonzero "
            "effective syndrome gap on the 3x3 reference."
        ),
        "claim_boundary": (
            "This is an exact effective stabilizer-spectrum closure, not a microscopic 2D walker simulation. It assumes "
            "independent gadget tiling and therefore does not bound inter-gadget crosstalk, establish a thermodynamic gap "
            "or native ANTLER realization, or demonstrate anyon motion, defects, fusion, non-Abelian braid, universality "
            "or fault tolerance. The toric reference itself is Abelian."
        ),
    }
    path = ROOT / "results" / "phase7" / "phase8b_effective_toric_tile_closure.json"
    path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
