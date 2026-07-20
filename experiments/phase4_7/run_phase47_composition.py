"""Phase 4.7C -- direct logical-gate composition for n=1,2,4,8.

For each n this runs the physically repeated exchange U_Z(theta)^n and an
independent one-cycle reference at the angle n*theta.  It therefore separates
the desired phase additivity from leakage growth and coherent-axis drift.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from phase47_common import (DigitalConfig, average_gate_fidelity, gate_axis,
                            matrix_from_json, remove_global, require_json,
                            run_gate)


def wrapped(angle: float) -> float:
    return float(np.angle(np.exp(1j * angle)))


def unwrap_near(angle: float, target: float) -> float:
    return float(angle + 2.0 * np.pi * round((target - angle) / (2.0 * np.pi)))


def compatible(cached: dict, cfg: DigitalConfig, theta: float, dt: float,
               cycles: int) -> bool:
    return (
        cached.get("schema") == "antler.phase47.gate.v1"
        and cached.get("config") == asdict(cfg)
        and float(cached.get("theta", np.nan)) == theta
        and float(cached.get("dt_requested", np.nan)) == dt
        and int(cached.get("cycles", 0)) == cycles
    )


def run_or_resume(path: Path, cfg: DigitalConfig, theta: float, dt: float,
                  cycles: int, gap_samples: int, resume: bool) -> dict:
    cached = require_json(path) if resume else None
    if cached is not None and compatible(cached, cfg, theta, dt, cycles):
        print(f"resume {path.name}", flush=True)
        return cached
    print(f"run {path.stem}: theta={theta:g}, cycles={cycles}", flush=True)
    started = time.time()
    result = run_gate(cfg, theta=theta, dt=dt, cycles=cycles,
                      include_gap=True, gap_samples=gap_samples)
    result["runtime_s"] = time.time() - started
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["metrics"], indent=2), flush=True)
    return result


def target_gate(angle: float) -> np.ndarray:
    return np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])


def coherent_distance(U: np.ndarray, V: np.ndarray) -> float:
    """Projective distance, minimizing over the physically irrelevant U(1) phase.

    ``remove_global`` fixes an SU(2) representative but retains its central
    sign.  The trace alignment here correctly identifies ``U`` and ``-U`` as
    the same logical operation, which is essential after branch selection.
    """
    W = V.conj().T @ U
    phase = float(np.angle(np.trace(W)))
    return float(np.linalg.norm(np.exp(-1j * phase) * W - np.eye(2)))


def relative_phase(U: np.ndarray) -> float:
    """Principal phase of U[0,0]/U[1,1] after global-phase removal."""
    U = remove_global(U)
    return float(np.angle(np.exp(1j * (np.angle(U[0, 0]) - np.angle(U[1, 1])))))


def balanced_square(result: dict) -> np.ndarray:
    """Rebuild Q = D(+theta) D(-theta)^dagger from cached raw transports.

    ``odd_gate_from_runs`` returns a principal square root of Q.  That root
    aliases when the composed relative phase crosses pi/2, so composition
    needs this raw object to select the continuous root independently.
    """
    runs = result["runs"]
    theta = float(result["theta"])
    def U(name: str, signed: float) -> np.ndarray:
        return matrix_from_json(runs[f"{name}_{signed:+.8f}"]["U"])
    plus = remove_global(U("rt", theta).conj().T @ U("ex", theta))
    minus = remove_global(U("rt", -theta).conj().T @ U("ex", -theta))
    return remove_global(plus @ minus.conj().T)


def branch_corrected_gate(result: dict, target_relative_phase: float) -> tuple[np.ndarray, dict]:
    """Select the square-root branch continuous with the requested cycle count.

    The correction is valid only for the near-diagonal logical gates already
    established in this campaign.  If Q is not near diagonal, a general
    matrix branch must be transported continuously instead of using this
    shortcut.
    """
    principal = matrix_from_json(result["Uodd"])
    Q = balanced_square(result)
    q_principal = relative_phase(Q)
    q_unwrapped = unwrap_near(q_principal, 2.0 * target_relative_phase)
    recovered = 0.5 * q_unwrapped
    principal_phase = relative_phase(principal)
    shift = wrapped(recovered - principal_phase)
    q_offdiag = float(np.linalg.norm(Q - np.diag(np.diag(Q))))
    if abs(shift) < 1e-7:
        gate, branch = principal, "principal"
    elif abs(abs(shift) - np.pi) < 1e-6 and q_offdiag < 1e-5:
        # The alternate SU(2) root of an effectively diagonal Q.  Its square
        # differs only by a global sign and hence represents the same Q.
        gate, branch = np.diag([1j, -1j]) @ principal, "alternate_diagonal_SU2_root"
    else:
        raise RuntimeError(
            "cannot select a safe composition branch: Q is not sufficiently diagonal "
            f"or the requested branch is not a square-root branch (shift={shift}, offdiag={q_offdiag})"
        )
    diagnostics = {
        "method": "raw_balanced_Q_then_target_anchored_square_root_branch",
        "Q_relative_phase_principal": q_principal,
        "Q_relative_phase_unwrapped": q_unwrapped,
        "recovered_relative_phase": recovered,
        "principal_root_relative_phase": principal_phase,
        "sqrt_branch": branch,
        "Q_offdiag_norm": q_offdiag,
    }
    return remove_global(gate), diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta", type=float, default=0.3)
    parser.add_argument("--cycles", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--depth", type=float, default=6.0,
                        help="positive trap-depth magnitude")
    parser.add_argument("--T", type=float, default=45_000.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--R", type=int, default=4)
    parser.add_argument("--gap-samples", type=int, default=61)
    parser.add_argument("--outdir", type=Path,
                        default=Path("results/phase4_7/composition"))
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.depth <= 0 or args.T <= 0 or args.dt <= 0 or args.theta <= 0:
        raise ValueError("theta, depth, T, and dt must be positive")
    if (sorted(set(args.cycles)) != sorted(args.cycles) or any(n < 1 for n in args.cycles)
            or args.cycles[0] != 1):
        raise ValueError("cycles must be a strictly increasing list of positive integers starting at 1")
    cfg = DigitalConfig(DEPTH=-args.depth, R_LOOP=args.R, T_TOTAL=args.T)
    started = time.time()
    results: dict[int, dict[str, dict]] = {}
    for n in args.cycles:
        composed_path = args.outdir / f"composed_n{n}.json"
        composed = run_or_resume(composed_path, cfg, args.theta, args.dt, n,
                                 args.gap_samples, resume=not args.no_resume)
        if n == 1:
            direct = composed
        else:
            direct_path = args.outdir / f"direct_theta_{n * args.theta:.6g}.json"
            direct = run_or_resume(direct_path, cfg, n * args.theta, args.dt, 1,
                                   args.gap_samples, resume=not args.no_resume)
        results[n] = {"composed": composed, "direct": direct}

    Uone, one_branch = branch_corrected_gate(
        results[min(args.cycles)]["composed"], -args.theta
    )
    leak_one = results[min(args.cycles)]["composed"]["metrics"]["leak_worst"]
    rows = []
    for n in args.cycles:
        composed, direct = results[n]["composed"], results[n]["direct"]
        Ucomp, comp_branch = branch_corrected_gate(composed, -n * args.theta)
        Udirect, direct_branch = branch_corrected_gate(direct, -n * args.theta)
        expected = target_gate(n * args.theta)
        Upower = remove_global(np.linalg.matrix_power(Uone, n))
        mcomp, mdirect = composed["metrics"], direct["metrics"]
        phase_comp = comp_branch["recovered_relative_phase"]
        phase_direct = direct_branch["recovered_relative_phase"]
        axis = gate_axis(Ucomp)
        direct_axis = gate_axis(Udirect)
        rows.append({
            "n": n,
            "target_phase": -n * args.theta,
            "composed_phase_unwrapped": phase_comp,
            "direct_phase_unwrapped": phase_direct,
            "phase_additivity_error": phase_comp - phase_direct,
            "phase_error_to_ideal": phase_comp + n * args.theta,
            "leak_worst": mcomp["leak_worst"],
            "leak_per_cycle": mcomp["leak_worst"] / n,
            "leak_relative_to_n_times_single": (
                mcomp["leak_worst"] / (n * leak_one) if leak_one > 0 else None
            ),
            "coherent_distance_to_ideal": coherent_distance(Ucomp, expected),
            "coherent_distance_to_direct": coherent_distance(Ucomp, Udirect),
            "coherent_distance_to_power_of_one_cycle": coherent_distance(Ucomp, Upower),
            "favg_to_ideal": average_gate_fidelity(Ucomp, expected),
            "odd_offdiag_norm": mcomp["odd_offdiag_norm"],
            "axis": axis,
            "direct_axis": direct_axis,
            "axis_drift_from_z": axis["z_axis_drift"],
            "axis_drift_relative_to_direct": float(
                np.linalg.norm(np.asarray(axis["axis"]) - np.asarray(direct_axis["axis"]))
            ),
            "minimum_handoff_isolation_gap": (composed.get("gap") or {}).get(
                "minimum_handoff_isolation_gap"
            ),
            "phase_extraction": {"composed": comp_branch, "direct": direct_branch},
            "composed_file": f"composed_n{n}.json",
            "direct_file": (
                f"composed_n{n}.json" if n == 1 else f"direct_theta_{n * args.theta:.6g}.json"
            ),
        })

    summary = {
        "schema": "antler.phase47.composition.v2",
        "claim_boundary": (
            "Composition is tested only for the present Abelian Z-family in the "
            "frozen correlated-hopping ladder.  Passing it does not make the family "
            "universal or non-commutative."
        ),
        "config": asdict(cfg),
        "theta": args.theta,
        "rows": rows,
        "leakage_model_diagnostic": (
            "leak_relative_to_n_times_single near one indicates approximately linear "
            "worst-case leakage; a systematic growth above one flags coherent or "
            "superlinear accumulation."
        ),
        "phase_branch_note": (
            "The phase is reconstructed from the cached raw +/-theta differentials. "
            "This removes the principal-square-root alias that otherwise appears for n=8; "
            "the alternate root is used only after verifying that Q is near diagonal."
        ),
        "one_cycle_phase_extraction": one_branch,
        "runtime_s": time.time() - started,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    path = args.outdir / "summary.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"saved {path}", flush=True)


if __name__ == "__main__":
    main()
