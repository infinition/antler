"""Phase 6I: distinguish fixed-charge code isolation from charge-sector gaps.

The external U(1) parent is a calibration reference.  A protected doublet in a
fixed charge sector and incompressibility under N -> N +/- 1 are separate
questions; neither may be substituted for the other.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.number_conserving_pairwire import build_iemini_hamiltonian, wire_a_parity


def ground_energy(L: int, N: int) -> float:
    sparse = L >= 8
    H, _, _ = build_iemini_hamiltonian(L, N, lam=1.0, sparse=sparse)
    if sparse:
        return float(eigsh(H, k=1, which="SA", return_eigenvectors=False, tol=1e-11)[0])
    return float(eigh(H, eigvals_only=True, subset_by_index=[0, 0], driver="evr")[0])


def neutral_code_gap(L: int, N: int) -> float:
    sparse = L >= 8
    H, states, _ = build_iemini_hamiltonian(L, N, lam=1.0, sparse=sparse)
    gaps = []
    for parity in (0, 1):
        rows = np.asarray([row for row, state in enumerate(states) if wire_a_parity(int(state), L) == parity], dtype=int)
        block = H[rows, :][:, rows] if sparse else H[np.ix_(rows, rows)]
        values = (
            np.sort(eigsh(block, k=2, which="SA", return_eigenvectors=False, tol=1e-11))
            if sparse else
            eigh(block, eigvals_only=True, subset_by_index=[0, 1], driver="evr")
        )
        gaps.append(float(values[1] - values[0]))
    return min(gaps)


def analyse(L: int) -> dict:
    N = L
    e_minus, e_zero, e_plus = (ground_energy(L, particles) for particles in (N - 1, N, N + 1))
    return {
        "L": L,
        "reference_charge_N": N,
        "pair_filling_N_over_2L": N / (2.0 * L),
        "ground_energies": {"N_minus_1": e_minus, "N": e_zero, "N_plus_1": e_plus},
        "addition_energy": e_plus - e_zero,
        "removal_energy": e_minus - e_zero,
        "fixed_N_neutral_code_gap": neutral_code_gap(L, N),
    }


def main() -> None:
    rows = [analyse(L) for L in (4, 6, 8)]
    out = {
        "schema": "antler.phase6.charge-sector-audit.v1",
        "reference": "external Iemini lambda=1 parent at N=L, corresponding to pair filling N/(2L)=1/2",
        "rows": rows,
        "interpretation": (
            "The parent has an exactly isolated fixed-N parity doublet, but its global charge-sector addition/removal energies must be inspected independently. "
            "A charge gap can be supplied by a conserved-sector preparation or an external charging energy; it cannot be silently inferred from fixed-N ED."
        ),
        "claim_boundary": (
            "This is a charge-versus-neutral-gap audit of an external parent only. It neither creates a native ANTLER qubit nor proves experimental protection against particle exchange with a reservoir."
        ),
    }
    path = ROOT / "results" / "phase6" / "charge_sector_audit.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
