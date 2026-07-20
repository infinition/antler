"""Phase 5I: finite-time braid of a six-Majorana T-junction coupling network.

The central Majorana gamma_0 is sequentially coupled to three junction arms.
With the third arm parked initially and finally, the three-stage adiabatic
cycle exchanges the two zero modes on arms 1 and 2.  Two remote zero modes
complete a fixed-even-parity logical qubit.  The script measures convergence
to the braid (or its inverse), leakage, and logical fidelity versus duration.

This remains an effective defect-network model.  A microscopic ANTLER bridge
must still derive the tunable Majorana couplings from pairing physics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, polar


I2 = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Y = np.array([[0.0, -1j], [1j, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def kron3(a, b, c):
    return np.kron(np.kron(a, b), c)


def gammas() -> list[np.ndarray]:
    return [kron3(X, I2, I2), kron3(Y, I2, I2),
            kron3(Z, X, I2), kron3(Z, Y, I2),
            kron3(Z, Z, X), kron3(Z, Z, Y)]


GAMMA = gammas()
PARITY = kron3(Z, Z, Z)
EVEN = np.eye(8, dtype=complex)[:, [0, 3, 5, 6]]
COUPLING_GENERATORS = np.asarray([1j * GAMMA[0] @ GAMMA[arm]
                                  for arm in range(1, 4)])
EVEN_COUPLING_GENERATORS = np.asarray([
    EVEN.conj().T @ generator @ EVEN for generator in COUPLING_GENERATORS
])


def no_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def global_distance(U: np.ndarray, V: np.ndarray) -> float:
    relative = V.conj().T @ U
    overlap = np.trace(relative)
    phase = 1.0 if abs(overlap) < 1e-15 else np.exp(-1j * np.angle(overlap))
    return float(np.linalg.norm(phase * relative - np.eye(2)))


def fidelity(U: np.ndarray, V: np.ndarray) -> float:
    relative = V.conj().T @ U
    return float((abs(np.trace(relative)) ** 2 + 2.0) / 6.0)


def hamiltonian(lambdas: np.ndarray) -> np.ndarray:
    return 0.5 * np.tensordot(lambdas, COUPLING_GENERATORS, axes=(0, 0))


def hamiltonian_even(lambdas: np.ndarray) -> np.ndarray:
    """Hamiltonian in the exactly conserved even-parity sector."""

    return 0.5 * np.tensordot(lambdas, EVEN_COUPLING_GENERATORS, axes=(0, 0))


def exact_step(lambdas: np.ndarray, dt: float, *, even_parity: bool = False) -> np.ndarray:
    """Exact short-time propagator for a three-arm Majorana coupling.

    The three Hermitian generators ``i gamma_0 gamma_a`` mutually
    anticommute and square to one.  Therefore ``H(lambdas)**2`` is the
    scalar ``||lambdas||**2 / 4`` times the identity, so this expression is
    algebraically identical to ``expm(-1j * dt * H)``.  It removes a generic
    matrix-exponential approximation from the finite-time convergence and
    makes ensemble robustness scans inexpensive.
    """

    norm = float(np.linalg.norm(lambdas))
    if norm < 1e-15:
        return np.eye(4 if even_parity else 8, dtype=complex)
    H = hamiltonian_even(lambdas) if even_parity else hamiltonian(lambdas)
    half_angle = 0.5 * dt * norm
    return (np.cos(half_angle) * np.eye(H.shape[0], dtype=complex)
            - 2j * np.sin(half_angle) * H / norm)


def schedule(u: float) -> np.ndarray:
    """lambda_3 -> lambda_1 -> lambda_2 -> lambda_3, with constant gap."""

    x = float(np.clip(u, 0.0, 1.0))
    segment = min(int(3.0 * x), 2)
    s = 3.0 * x - segment
    q = np.sin(0.5 * np.pi * s) ** 2
    lam = np.zeros(3)
    if segment == 0:
        lam[2], lam[0] = 1.0 - q, q
    elif segment == 1:
        lam[0], lam[1] = 1.0 - q, q
    else:
        lam[1], lam[2] = 1.0 - q, q
    return lam


def ground_code() -> np.ndarray:
    H_even = hamiltonian_even(np.array([0.0, 0.0, 1.0]))
    energy, vectors = np.linalg.eigh(H_even)
    ground = EVEN @ vectors[:, :2]
    assert np.linalg.norm(ground.conj().T @ PARITY @ ground - np.eye(2)) < 1e-12
    return ground


def braid(T: float, dt: float, arm_scales: np.ndarray | None = None) -> dict:
    """Evolve one exchange, optionally with static independent arm scales."""

    scales = np.ones(3) if arm_scales is None else np.asarray(arm_scales, float)
    if scales.shape != (3,) or np.any(scales <= 0.0):
        raise ValueError("arm_scales must contain three positive values")
    frame = ground_code()
    frame_even = EVEN.conj().T @ frame
    nseg = int(round(T / dt)); dt = T / nseg
    psi = frame_even.copy()
    minimum_gap = float("inf")
    for step in range(nseg):
        u = (step + .5) / nseg
        lambdas = schedule(u) * scales
        # In this Clifford network, the excited-pair separation is exactly
        # ||lambda||; no numerical eigensolver is involved in the gap audit.
        minimum_gap = min(minimum_gap, float(np.linalg.norm(lambdas)))
        psi = exact_step(lambdas, dt, even_parity=True) @ psi
    S = frame_even.conj().T @ psi
    U = no_global(polar(S)[0])
    singular = np.linalg.svd(S, compute_uv=False)
    B12_full = expm(np.pi * GAMMA[2] @ GAMMA[1] / 4.0)
    B12 = no_global(frame.conj().T @ B12_full @ frame)
    candidates = {"B12": B12, "B12_inverse": B12.conj().T}
    label, target = min(candidates.items(), key=lambda item: global_distance(U, item[1]))
    return {
        "T": T, "dt": dt, "nseg": nseg, "arm_scales": scales.tolist(),
        "U": {"real": U.real.tolist(), "imag": U.imag.tolist()},
        "target_orientation": label,
        "matrix_error": global_distance(U, target),
        "favg_target": fidelity(U, target),
        "leak_worst": float(1.0 - min(singular) ** 2),
        "sigma_min": float(min(singular)),
        "minimum_even_parity_gap": minimum_gap,
    }


def braid_ensemble(T: float, dt: float, arm_scales: np.ndarray) -> list[dict]:
    """Vectorised fixed-parity evolution for many static coupling samples.

    All trajectories use the same temporal mesh, while each follows its own
    three-dimensional coupling path.  The only loop over samples is the final
    two-by-two polar decomposition used for logical diagnostics.
    """

    scales = np.asarray(arm_scales, float)
    if scales.ndim != 2 or scales.shape[1] != 3 or np.any(scales <= 0.0):
        raise ValueError("arm_scales must have shape (samples, 3) and be positive")
    frame = ground_code()
    frame_even = EVEN.conj().T @ frame
    nseg = int(round(T / dt)); dt = T / nseg
    samples = scales.shape[0]
    psi = np.broadcast_to(frame_even, (samples,) + frame_even.shape).copy()
    minimum_gap = np.full(samples, np.inf)
    eye = np.eye(4, dtype=complex)[None, :, :]
    for step in range(nseg):
        lambdas = schedule((step + .5) / nseg)[None, :] * scales
        norm = np.linalg.norm(lambdas, axis=1)
        H = 0.5 * np.einsum("sa,aij->sij", lambdas, EVEN_COUPLING_GENERATORS)
        half_angle = 0.5 * dt * norm
        propagator = (np.cos(half_angle)[:, None, None] * eye
                      - 2j * np.sin(half_angle)[:, None, None] * H / norm[:, None, None])
        psi = propagator @ psi
        minimum_gap = np.minimum(minimum_gap, norm)
    overlaps = np.einsum("ai,sib->sab", frame_even.conj().T, psi)
    B12_full = expm(np.pi * GAMMA[2] @ GAMMA[1] / 4.0)
    B12 = no_global(frame.conj().T @ B12_full @ frame)
    candidates = {"B12": B12, "B12_inverse": B12.conj().T}
    output = []
    for scale, S, gap in zip(scales, overlaps, minimum_gap):
        U = no_global(polar(S)[0])
        singular = np.linalg.svd(S, compute_uv=False)
        label, target = min(candidates.items(), key=lambda item: global_distance(U, item[1]))
        output.append({
            "T": T, "dt": dt, "nseg": nseg, "arm_scales": scale.tolist(),
            "target_orientation": label,
            "matrix_error": global_distance(U, target),
            "favg_target": fidelity(U, target),
            "leak_worst": float(1.0 - min(singular) ** 2),
            "sigma_min": float(min(singular)),
            "minimum_even_parity_gap": float(gap),
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--times", type=float, nargs="+", default=[10.0, 20.0, 40.0, 80.0])
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--out", type=Path,
                        default=Path("results/phase5/tjunction_braid_dynamics.json"))
    args = parser.parse_args()
    rows = [braid(T, args.dt) for T in args.times]
    out = {
        "schema": "antler.phase5.tjunction-braid-dynamics.v1",
        "claim_boundary": (
            "This is a finite-time effective Majorana-network braid.  The microscopic "
            "pairing/T-junction derivation from ANTLER remains a separate requirement."
        ),
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
