"""Phase 5A: N=3 mobile-mediator logical-doublet preflight.

Before attempting two time-dependent exchanges, this sparse ED scan asks the
necessary question: can the unmodified two-leg correlated-hopping ladder host
a *two-dimensional, localized and spectrally isolated* code with a pinned
third particle?  A positive result only authorises dynamic braid searches; a
negative result is useful evidence that N=3 alone is insufficient and that a
synthetic dimension or T-junction must be added.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh


ROOT = Path(__file__).resolve().parents[2]
PHASE41 = ROOT / "experiments" / "phase4_1"
for path in (ROOT, PHASE41):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from antler.basis import build_basis
from antler.phase1 import hop_table
from run_phase4_1_logical_gate import build_occ


@dataclass(frozen=True)
class PreflightConfig:
    L: int = 14
    N: int = 3
    J1: float = 0.4
    J2: float = 1.0
    JPERP: float = 0.1
    code_depth: float = -4.0
    mediator_depth_scale: float = 1.0
    mediator_rung: int = 7
    mediator_leg: int = 0


def bare_mask(*sites: int) -> int:
    mask = 0
    for site in sites:
        mask |= 1 << site
    return mask


class N3Preflight:
    def __init__(self, cfg: PreflightConfig):
        self.cfg = cfg
        self.M = 2 * cfg.L
        self.states, self.index = build_basis(self.M, cfg.N)
        self.table = hop_table(cfg.L, cfg.J1, cfg.J2, cfg.JPERP,
                               self.states, self.index)
        self.occ = build_occ(self.states, self.M)
        self.mediator = 2 * cfg.mediator_rung + cfg.mediator_leg
        if self.mediator in {0, 1, self.M - 2, self.M - 1}:
            raise ValueError("mediator must not occupy one of the four code-edge sites")
        self.i_left = self.index[bare_mask(0, 1, self.mediator)]
        self.i_right = self.index[bare_mask(self.M - 2, self.M - 1, self.mediator)]

    def mu(self) -> np.ndarray:
        mu = np.zeros(self.M)
        mu[[0, 1, self.M - 2, self.M - 1]] = self.cfg.code_depth
        mu[self.mediator] = self.cfg.code_depth * self.cfg.mediator_depth_scale
        return mu

    def hamiltonian(self, theta: float) -> csr_matrix:
        rows, cols, mJ, nmid = self.table
        one = csr_matrix((mJ * np.exp(1j * theta * nmid), (rows, cols)),
                         shape=(len(self.states), len(self.states)))
        return one + one.conj().T + diags(self.occ @ self.mu())

    def analyse(self, theta: float, n_eigs: int) -> dict:
        H = self.hamiltonian(theta)
        # `which=SA` tests the physically relevant pinned-trap manifold.
        energies, vectors = eigsh(H, k=n_eigs, which="SA", tol=1e-10)
        order = np.argsort(energies)
        energies, vectors = energies[order], vectors[:, order]
        target_weights = abs(vectors[self.i_left, :]) ** 2 + abs(vectors[self.i_right, :]) ** 2
        pair = np.sort(np.argsort(-target_weights)[:2])
        outside = np.setdiff1d(np.arange(n_eigs), pair)
        isolation = float(np.min(abs(energies[pair, None] - energies[outside])))
        capture = float(target_weights[pair].sum())
        logical_split = float(abs(energies[pair[1]] - energies[pair[0]]))
        mediator_occ = np.real(abs(vectors[self._state_indices_with_mediator(), :]) ** 2).sum(axis=0)
        return {
            "theta": theta,
            "energies_lowest": energies.tolist(),
            "selected_eigenvalue_indices": pair.tolist(),
            "selected_energies": energies[pair].tolist(),
            "target_weights_selected": target_weights[pair].tolist(),
            "capture": capture,
            "logical_split": logical_split,
            "isolation_gap_to_low_energy_complement": isolation,
            "mediator_occupation_selected": mediator_occ[pair].tolist(),
        }

    def _state_indices_with_mediator(self) -> np.ndarray:
        return np.asarray(
            [i for i, state in enumerate(self.states) if (int(state) >> self.mediator) & 1],
            dtype=int,
        )


def assess(rows: list[dict]) -> dict:
    """Conservative gatekeeper for whether dynamics are worth compiling."""

    acceptable = []
    for row in rows:
        scans = row["spectra"]
        capture = min(scan["capture"] for scan in scans)
        occupation = min(min(scan["mediator_occupation_selected"]) for scan in scans)
        gap = min(scan["isolation_gap_to_low_energy_complement"] for scan in scans)
        split = max(scan["logical_split"] for scan in scans)
        candidate = capture >= 0.80 and occupation >= 0.80 and gap >= 5.0 * split
        row["gatekeeper"] = {
            "capture_min_over_theta": capture,
            "mediator_occupation_min_over_theta": occupation,
            "isolation_gap_min_over_theta": gap,
            "logical_split_max_over_theta": split,
            "passes_static_preflight": candidate,
        }
        if candidate:
            acceptable.append(row["label"])
    return {
        "pass_labels": acceptable,
        "decision": (
            "compile_two_exchange_dynamics" if acceptable else
            "do_not_interpret_N3_as_a_nonabelian_code; broaden_to_synthetic_dimension_or_T_junction"
        ),
        "note": (
            "A static pass is necessary but not sufficient.  It does not demonstrate "
            "non-commutativity; that requires two full logical propagations on the same "
            "tracked doublet."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--code-depth", type=float, default=-4.0,
                        help="negative depth of the four code-edge wells")
    parser.add_argument("--rungs", type=int, nargs="+", default=[3, 5, 7, 9])
    parser.add_argument("--depth-scales", type=float, nargs="+", default=[1.0, 1.25, 1.5])
    parser.add_argument("--n-eigs", type=int, default=12)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/n3_mediator_preflight.json"))
    args = parser.parse_args()
    if args.n_eigs < 4:
        raise ValueError("n-eigs must be at least four")
    if any(r < 1 or r > 12 for r in args.rungs):
        raise ValueError("mediator rungs must lie inside the 14-rung ladder")
    if any(scale <= 0 for scale in args.depth_scales):
        raise ValueError("mediator depth scales must be positive")
    if args.code_depth >= 0:
        raise ValueError("code-depth must be negative")

    started = time.time()
    rows = []
    for rung in args.rungs:
        for leg in (0, 1):
            for scale in args.depth_scales:
                cfg = PreflightConfig(mediator_rung=rung, mediator_leg=leg,
                                      mediator_depth_scale=scale,
                                      code_depth=args.code_depth)
                label = f"r{rung}_leg{leg}_scale{scale:g}"
                print(f"scan {label}", flush=True)
                system = N3Preflight(cfg)
                spectra = [system.analyse(theta, args.n_eigs) for theta in (0.0, args.theta)]
                rows.append({"label": label, "config": asdict(cfg), "spectra": spectra})
    gate = assess(rows)
    rows.sort(key=lambda row: (
        not row["gatekeeper"]["passes_static_preflight"],
        -row["gatekeeper"]["capture_min_over_theta"],
        -row["gatekeeper"]["isolation_gap_min_over_theta"],
    ))
    out = {
        "schema": "antler.phase5.n3-preflight.v1",
        "claim_boundary": (
            "This tests only a pinned-mediator static code candidate.  A positive "
            "result does not establish a non-Abelian braid; a negative result does "
            "not rule out a synthetic-dimension or T-junction extension."
        ),
        "dimension": int(len(N3Preflight(PreflightConfig()).states)),
        "theta_tested": [0.0, args.theta],
        "rows": rows,
        "gatekeeper": gate,
        "runtime_s": time.time() - started,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(gate, indent=2), flush=True)
    print(f"saved {args.out}", flush=True)


if __name__ == "__main__":
    main()
