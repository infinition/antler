"""Shared, additive Phase 4.7 digital-quantization machinery.

Nothing in this module changes the frozen ED model, logical encoding, or
Phase 4.3 protocol implementation.  It reuses the same rung-major hopping
table and logical-frame construction, and makes the Phase 4.7 controls and
diagnostics explicit and serialisable.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
from scipy.linalg import polar, sqrtm
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix, diags


ROOT = Path(__file__).resolve().parents[2]
PHASE43 = ROOT / "experiments" / "phase4_1"
for _path in (ROOT, PHASE43):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from antler.basis import build_basis
from antler.phase1 import hop_table
from run_phase4_1_logical_gate import build_occ
from run_phase4_1_logical_gate_strang import exact_logical_frame


@dataclass(frozen=True)
class DigitalConfig:
    """Frozen physical model plus explicit Phase 4.7 control parameters."""

    L: int = 14
    N: int = 2
    J1: float = 0.4
    J2: float = 1.0
    JPERP: float = 0.1
    DEPTH: float = -4.0
    R_LOOP: int = 4
    T_TOTAL: float = 20_000.0
    ramp: str = "sin2"
    handoff_fraction: float = 0.10
    pause_fraction: float = 0.0
    parking_depth_scale: float = 1.0
    spectator_depth_scale: float = 1.0
    handoff_order: str = "left_then_rung"


def ramp(s: float, kind: str) -> float:
    """Monotone C0/C1 digital handoff profiles, normalised to [0, 1]."""

    x = float(np.clip(s, 0.0, 1.0))
    if kind == "linear":
        return x
    if kind == "sin2":
        return float(np.sin(0.5 * np.pi * x) ** 2)
    if kind == "smoothstep":
        return x * x * (3.0 - 2.0 * x)
    raise ValueError(f"unknown ramp: {kind}")


def add_discrete_trap(mu: np.ndarray, leg: int, x: float, depth: float, L: int,
                      profile: str) -> None:
    """Place one compact trap on a lattice site or cross-fade it across a bond."""

    x = float(np.clip(x, 0.0, L - 1.0))
    i = int(np.floor(x))
    if i >= L - 1:
        mu[2 * (L - 1) + leg] += depth
        return
    q = ramp(x - i, profile)
    mu[2 * i + leg] += depth * (1.0 - q)
    mu[2 * (i + 1) + leg] += depth * q


def _durations(cfg: DigitalConfig) -> tuple[float, float, float, float]:
    """Return outbound, handoff, handoff, return durations before pauses."""

    h = float(cfg.handoff_fraction)
    if not 0.0 < h < 0.5:
        raise ValueError("handoff_fraction must be in (0, 0.5)")
    travel = 1.0 - 2.0 * h
    # Preserve the original 0.35:0.45 outbound:return ratio.
    return travel * 0.35 / 0.80, h, h, travel * 0.45 / 0.80


def _protocol_stage(u: float, cfg: DigitalConfig) -> tuple[str, float]:
    """Map normalised time to a stage and its local coordinate.

    There are three optional constant-potential pauses: after outbound motion,
    between the two handoffs, and before return.  They deliberately preserve
    the worldline topology while changing the time parametrisation.
    """

    if not 0.0 <= cfg.pause_fraction < 1.0 / 3.0:
        raise ValueError("pause_fraction must be in [0, 1/3)")
    out, first, second, ret = _durations(cfg)
    active = 1.0 - 3.0 * cfg.pause_fraction
    entries = [
        ("outbound", active * out),
        ("pause_outbound", cfg.pause_fraction),
        ("first_handoff", active * first),
        ("pause_middle", cfg.pause_fraction),
        ("second_handoff", active * second),
        ("pause_return", cfg.pause_fraction),
        ("return", active * ret),
    ]
    x = float(np.clip(u, 0.0, 1.0))
    acc = 0.0
    for name, duration in entries:
        if x <= acc + duration + 1e-15:
            local = 1.0 if duration == 0.0 else (x - acc) / duration
            return name, float(np.clip(local, 0.0, 1.0))
        acc += duration
    return "return", 1.0


def protocol_description(cfg: DigitalConfig) -> dict:
    out, first, second, ret = _durations(cfg)
    return {
        "ramp": cfg.ramp,
        "handoff_order": cfg.handoff_order,
        "base_stage_fractions": {
            "outbound": out,
            "first_handoff": first,
            "second_handoff": second,
            "return": ret,
        },
        "pause_fraction_each": cfg.pause_fraction,
        "parking_depth_scale": cfg.parking_depth_scale,
        "spectator_depth_scale": cfg.spectator_depth_scale,
    }


def mu_digital(u: float, exchange: bool, cfg: DigitalConfig) -> np.ndarray:
    """Sequential digital protocol with controllable but topology-fixed timing.

    The default arguments reproduce the Phase 4.3b sequential protocol.  The
    only new degrees of freedom are deliberately exposed path deformations.
    In particular, the two swap handoffs always remain non-simultaneous.
    """

    if cfg.handoff_order not in {"left_then_rung", "rung_then_left"}:
        raise ValueError("handoff_order must be left_then_rung or rung_then_left")
    mu = np.zeros(2 * cfg.L)
    D, R = cfg.DEPTH, cfg.R_LOOP
    # The untouched right cat branch remains a spectator, as in Phase 4.3b.
    mu[-2] = D * cfg.spectator_depth_scale
    mu[-1] = D * cfg.spectator_depth_scale
    park = D * cfg.parking_depth_scale
    # A parking-well deformation must not change the endpoint Hamiltonian:
    # it is turned on after departure and turned off before recapture.
    def load_parking(x: float) -> float:
        q = ramp(x, cfg.ramp)
        return D * (1.0 - q) + park * q

    def unload_parking(x: float) -> float:
        q = ramp(x, cfg.ramp)
        return park * (1.0 - q) + D * q

    stage, s = _protocol_stage(u, cfg)
    first_is_left = cfg.handoff_order == "left_then_rung"
    if stage in {"pause_outbound", "pause_middle", "pause_return"}:
        s = 1.0 if stage != "pause_outbound" else 0.0
        if stage == "pause_outbound":
            stage = "first_handoff"
        elif stage == "pause_middle":
            stage = "first_handoff"
        else:
            stage = "second_handoff"

    if stage == "outbound":
        # The centre moves linearly through rung space; `add_discrete_trap`
        # supplies the selected local cross-fade.  For ramp='sin2' this is
        # exactly the frozen Phase 4.3b schedule.
        add_discrete_trap(mu, 0, R * s, D, cfg.L, cfg.ramp)
        mu[1] += load_parking(s)
        return mu

    # Round-trip has no rung transfers, but retains matched dwell intervals.
    if not exchange:
        if stage == "return":
            add_discrete_trap(mu, 0, R * (1.0 - s), D, cfg.L,
                              cfg.ramp)
        else:
            mu[2 * R] += D
        mu[1] += unload_parking(s) if stage == "return" else park
        return mu

    if stage == "first_handoff":
        doing_left = first_is_left
        q = ramp(s, cfg.ramp)
        if doing_left:
            mu[2 * R] += D
            mu[1] += park * (1.0 - q)
            mu[0] += park * q
        else:
            mu[1] += park
            mu[2 * R] += D * (1.0 - q)
            mu[2 * R + 1] += D * q
        return mu

    if stage == "second_handoff":
        doing_left = not first_is_left
        q = ramp(s, cfg.ramp)
        if doing_left:
            mu[2 * R + 1] += D
            mu[1] += park * (1.0 - q)
            mu[0] += park * q
        else:
            mu[0] += park
            mu[2 * R] += D * (1.0 - q)
            mu[2 * R + 1] += D * q
        return mu

    if stage == "return":
        add_discrete_trap(mu, 1, R * (1.0 - s), D, cfg.L,
                          cfg.ramp)
        mu[0] += unload_parking(s)
        return mu
    raise RuntimeError(f"unhandled protocol stage: {stage}")


def remove_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def matrix_json(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def matrix_from_json(data: dict) -> np.ndarray:
    return np.asarray(data["real"], dtype=float) + 1j * np.asarray(data["imag"], dtype=float)


def average_gate_fidelity(U: np.ndarray, V: np.ndarray) -> float:
    return float((abs(np.trace(V.conj().T @ U)) ** 2 + 2.0) / 6.0)


class DigitalSystem:
    """Cached fixed Hilbert-space data for one Phase 4.7 campaign."""

    def __init__(self, cfg: DigitalConfig):
        self.cfg = cfg
        self.M = 2 * cfg.L
        self.states, self.index = build_basis(self.M, cfg.N)
        self.table = hop_table(cfg.L, cfg.J1, cfg.J2, cfg.JPERP,
                               self.states, self.index)
        self.occ = build_occ(self.states, self.M)
        self._hop_cache: dict[float, tuple[csr_matrix, np.ndarray, np.ndarray]] = {}

    def hop(self, theta: float) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
        key = float(theta)
        if key not in self._hop_cache:
            rows, cols, mJ, nmid = self.table
            one = csr_matrix((mJ * np.exp(1j * theta * nmid), (rows, cols)),
                             shape=(len(self.states), len(self.states)))
            H = one + one.conj().T
            E, V = np.linalg.eigh(H.toarray())
            self._hop_cache[key] = H, E, V
        return self._hop_cache[key]

    def hamiltonian(self, theta: float, u: float, exchange: bool) -> csr_matrix:
        Hhop, _, _ = self.hop(theta)
        return Hhop + diags(self.occ @ mu_digital(u, exchange, self.cfg))

def _propagate(system: DigitalSystem, theta: float, exchange: bool, dt: float,
               cycles: int) -> dict:
    cfg = system.cfg
    Hhop, E, V = system.hop(theta)
    H0 = Hhop + diags(system.occ @ mu_digital(0.0, exchange, cfg))
    Hf = Hhop + diags(system.occ @ mu_digital(1.0, exchange, cfg))
    if np.linalg.norm((H0 - Hf).toarray()) > 1e-11:
        raise AssertionError("digital protocol endpoint Hamiltonians differ")
    U0, frame = exact_logical_frame(H0, system.index, system.M)
    nseg = int(round(cfg.T_TOTAL / dt))
    dt_eff = cfg.T_TOTAL / nseg
    Uhop = (V * np.exp(-1j * E * dt_eff)) @ V.conj().T
    psi = U0.copy()
    for _cycle in range(cycles):
        for a in range(nseg):
            u = (a + 0.5) / nseg
            potential = system.occ @ mu_digital(u, exchange, cfg)
            phase = np.exp(-0.5j * dt_eff * potential)[:, None]
            psi = phase * (Uhop @ (phase * psi))
    S = U0.conj().T @ psi
    return {
        "S": S,
        "U": polar(S)[0],
        "leak": 1.0 - np.sum(abs(S) ** 2, axis=0),
        "frame": frame,
        "nseg_per_cycle": nseg,
        "dt": dt_eff,
    }


def odd_gate_from_runs(runs: dict[tuple[str, float], dict], theta: float) -> tuple[np.ndarray, dict]:
    """Cancel even-in-theta backgrounds exactly as in validated Phase 4.3."""

    differential: dict[float, np.ndarray] = {}
    singular_values: list[float] = []
    unitarity_errors: list[float] = []
    for signed_theta in (theta, -theta):
        for name in ("rt", "ex"):
            S = runs[(name, signed_theta)]["S"]
            singular_values.extend(np.linalg.svd(S, compute_uv=False).tolist())
            unitarity_errors.append(float(np.linalg.norm(S.conj().T @ S - np.eye(2))))
        differential[signed_theta] = remove_global(
            runs[("rt", signed_theta)]["U"].conj().T @ runs[("ex", signed_theta)]["U"]
        )
    Q = remove_global(differential[theta] @ differential[-theta].conj().T)
    if np.linalg.norm(-Q - np.eye(2)) < np.linalg.norm(Q - np.eye(2)):
        Q = -Q
    Uodd = remove_global(polar(sqrtm(Q).astype(complex))[0])
    phase = float(np.angle(np.exp(1j * (
        np.angle(Uodd[0, 0]) - np.angle(Uodd[1, 1])
    ))))
    metrics = {
        "sigma_min": float(min(singular_values)),
        "sigma_max": float(max(singular_values)),
        "leak_worst": float(1.0 - min(singular_values) ** 2),
        "unitarity_frob_max": float(max(unitarity_errors)),
        "odd_phase": phase,
        "odd_slope": float(phase / theta),
        "odd_offdiag_norm": float(np.linalg.norm(Uodd - np.diag(np.diag(Uodd)))),
        "favg_target": average_gate_fidelity(
            Uodd, np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)])
        ),
    }
    return Uodd, metrics


def run_gate(cfg: DigitalConfig, theta: float, dt: float, cycles: int = 1,
             include_gap: bool = True, gap_samples: int = 61) -> dict:
    """Run matched exchange/round-trip controls at +/-theta and analyse them."""

    system = DigitalSystem(cfg)
    runs: dict[tuple[str, float], dict] = {}
    for signed_theta in (theta, -theta):
        for name, exchange in (("rt", False), ("ex", True)):
            runs[(name, signed_theta)] = _propagate(system, signed_theta, exchange,
                                                    dt, cycles)
    Uodd, metrics = odd_gate_from_runs(runs, theta)
    if cycles != 1:
        phase = metrics["odd_phase"]
        metrics["odd_phase_unwrapped_near_target"] = float(
            phase + 2.0 * np.pi * round((-cycles * theta - phase) / (2.0 * np.pi))
        )
        metrics["odd_slope_per_cycle"] = float(
            metrics["odd_phase_unwrapped_near_target"] / (cycles * theta)
        )
        target = np.diag([np.exp(-0.5j * cycles * theta),
                          np.exp(0.5j * cycles * theta)])
        metrics["favg_target"] = average_gate_fidelity(Uodd, target)
    gap = None
    if include_gap:
        gap = scan_handoff_gap(system, theta, gap_samples)
    serial_runs = {}
    for (name, signed_theta), run in runs.items():
        serial_runs[f"{name}_{signed_theta:+.8f}"] = {
            "S": matrix_json(run["S"]), "U": matrix_json(run["U"]),
            "leak": run["leak"].tolist(), "frame": run["frame"],
            "nseg_per_cycle": run["nseg_per_cycle"], "dt": run["dt"],
        }
    return {
        "schema": "antler.phase47.gate.v1",
        "config": asdict(cfg),
        "protocol": protocol_description(cfg),
        "theta": theta,
        "cycles": cycles,
        "dt_requested": dt,
        "metrics": metrics,
        "Uodd": matrix_json(Uodd),
        "gap": gap,
        "runs": serial_runs,
    }


def scan_handoff_gap(system: DigitalSystem, theta: float, samples: int = 61) -> dict:
    """Track the two-dimensional logical branch and its isolation gap.

    At each sample all ED eigenstates are available.  The two branch states
    are selected by maximum-overlap continuation from the preceding sample;
    the reported gap is their minimum distance to every other eigenvalue.
    This is a *subspace-isolation* gap, not the possibly much smaller internal
    logical splitting.
    """

    if samples < 3:
        raise ValueError("gap scan needs at least three samples")
    previous: np.ndarray | None = None
    rows: list[dict] = []
    for u in np.linspace(0.0, 1.0, samples):
        H = system.hamiltonian(theta, float(u), True).toarray()
        energies, vectors = np.linalg.eigh(H)
        if previous is None:
            frame, _ = exact_logical_frame(H, system.index, system.M)
            score = np.abs(frame.conj().T @ vectors) ** 2
        else:
            score = np.abs(previous.conj().T @ vectors) ** 2
        # Maximum-weight injective assignment of the two continuation vectors.
        selected_rows, selected = linear_sum_assignment(-score)
        if len(selected) != 2:
            raise RuntimeError("could not continue the two-state logical branch")
        selected = np.sort(selected)
        outside = np.setdiff1d(np.arange(len(energies)), selected)
        separation = float(np.min(np.abs(energies[selected, None] - energies[outside])))
        internal = float(abs(energies[selected[1]] - energies[selected[0]]))
        previous = vectors[:, selected]
        stage, local = _protocol_stage(float(u), system.cfg)
        rows.append({
            "u": float(u), "stage": stage, "stage_coordinate": local,
            "branch_energies": energies[selected].tolist(),
            "internal_split": internal, "isolation_gap": separation,
        })
    worst = min(rows, key=lambda row: row["isolation_gap"])
    handoff_rows = [row for row in rows if row["stage"] in {"first_handoff", "second_handoff"}]
    handoff_worst = min(handoff_rows, key=lambda row: row["isolation_gap"])
    return {
        "definition": "minimum spectral separation between the overlap-tracked "
                      "two-state logical branch and its complement",
        "theta": theta,
        "samples": samples,
        "minimum_protocol_isolation_gap": worst["isolation_gap"],
        "protocol_at_u": worst["u"], "protocol_stage": worst["stage"],
        "minimum_handoff_isolation_gap": handoff_worst["isolation_gap"],
        "handoff_at_u": handoff_worst["u"], "handoff_stage": handoff_worst["stage"],
        "internal_split_at_handoff_minimum": handoff_worst["internal_split"],
        "trace": rows,
    }


def gate_axis(U: np.ndarray) -> dict:
    """Return SU(2) angle and axis after stripping a global phase."""

    U = remove_global(U)
    pauli = (
        np.array([[0.0, 1.0], [1.0, 0.0]], complex),
        np.array([[0.0, -1j], [1j, 0.0]], complex),
        np.array([[1.0, 0.0], [0.0, -1.0]], complex),
    )
    coefficients = np.array([(0.5j * np.trace(P @ U)).real for P in pauli])
    sine = float(np.linalg.norm(coefficients))
    angle = float(2.0 * np.arctan2(sine, np.clip(0.5 * np.trace(U).real, -1.0, 1.0)))
    axis = (coefficients / sine).tolist() if sine > 1e-14 else [0.0, 0.0, 1.0]
    return {"angle": angle, "axis": axis, "z_axis_drift": float(np.hypot(axis[0], axis[1]))}


def require_json(path: Path) -> dict | None:
    """Load a completed cache item, or return None when a run must be launched."""

    if not path.exists():
        return None
    import json
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
