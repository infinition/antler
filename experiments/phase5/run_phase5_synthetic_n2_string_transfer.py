"""Phase 5E: exact N=2 correlated-string SU(2) transfer on a local ladder cell.

One hard-core spectator occupies the intermediate rung-major site while an
active pseudo-spin is digitally transferred over it.  The hopping therefore
contains the same density-dependent Jordan--Wigner phase as ANTLER, while the
synthetic link applies a noncommuting SU(2) matrix to the active spinor.

This is a local exact-ED validation of the new Hamiltonian ingredient; it is
not yet the protected L=14 logical braid protocol.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import polar

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from antler.synthetic import build_synthetic_hamiltonian, mode, su2


I = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def no_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def global_phase_distance(U: np.ndarray, V: np.ndarray) -> float:
    relative = V.conj().T @ U
    trace = np.trace(relative)
    phase = 1.0 if abs(trace) < 1e-15 else np.exp(-1j * np.angle(trace))
    return float(np.linalg.norm(phase * relative - I))


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def embeddings(states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    index = {int(state): i for i, state in enumerate(states)}
    source = np.zeros((len(states), 2), complex)
    target = np.zeros((len(states), 2), complex)
    spectator = 1 << mode(1, 0)
    for spin in range(2):
        source[index[spectator | (1 << mode(0, spin))], spin] = 1.0
        target[index[spectator | (1 << mode(2, spin))], spin] = 1.0
    return source, target


def align(frame: np.ndarray, embedding: np.ndarray) -> np.ndarray:
    return frame @ polar(embedding.conj().T @ frame)[0].conj().T


def frame(H: np.ndarray, embedding: np.ndarray) -> np.ndarray:
    energies, vectors = np.linalg.eigh(H)
    # The pinned spectator has an intentionally free internal component, so
    # the physical ground manifold is fourfold: spectator spin x active-spin
    # doublet.  Project the desired spectator-spin sector into that manifold
    # instead of selecting arbitrary vectors inside an exact degeneracy.
    manifold = vectors[:, :4]
    projected = manifold @ (manifold.conj().T @ embedding)
    Q, _ = np.linalg.qr(projected)
    return align(Q[:, :2], embedding)


def transfer(U_link: np.ndarray, theta: float = 0.3, T: float = 2400.0,
             dt: float = 0.05) -> dict:
    L, N, D_active, D_spectator = 2, 2, 10.0, 20.0
    zero = np.zeros(2 * L)
    Hhop, states, _ = build_synthetic_hamiltonian(
        L, N, theta=theta, J1=1.0, J2=0.0, Jperp=0.0, mu=zero,
        link_matrices={(0, 2): U_link},
    )
    Hhop = Hhop.toarray()
    source, target = embeddings(states)
    mu0 = zero.copy(); mu0[0] = -D_active; mu0[1] = -D_spectator
    muf = zero.copy(); muf[2] = -D_active; muf[1] = -D_spectator
    # State order in this compact basis is not mode-contiguous, so build the
    # potential diagonal from each hard-core basis mask directly.
    def potential(mu: np.ndarray) -> np.ndarray:
        out = np.zeros(len(states))
        for i, raw in enumerate(states):
            state = int(raw)
            for site in range(2 * L):
                if ((state >> mode(site, 0)) & 1) or ((state >> mode(site, 1)) & 1):
                    out[i] += mu[site]
        return out
    frame0 = frame(Hhop + np.diag(potential(mu0)), source)
    framef = frame(Hhop + np.diag(potential(muf)), target)
    energy, vectors = np.linalg.eigh(Hhop)
    nseg = int(round(T / dt)); dt = T / nseg
    Uhop = (vectors * np.exp(-1j * energy * dt)) @ vectors.conj().T
    psi = frame0.copy()
    for step in range(nseg):
        u = (step + 0.5) / nseg
        q = np.sin(0.5 * np.pi * u) ** 2
        mu = zero.copy(); mu[0] = -D_active * (1.0 - q); mu[2] = -D_active * q
        mu[1] = -D_spectator
        phase = np.exp(-0.5j * dt * potential(mu))[:, None]
        psi = phase * (Uhop @ (phase * psi))
    S = framef.conj().T @ psi
    Ueff = no_global(polar(S)[0])
    expected = no_global(U_link.conj().T)
    singular = np.linalg.svd(S, compute_uv=False)
    return {
        "Ueff": Ueff, "expected": expected,
        "matrix_error": global_phase_distance(Ueff, expected),
        "leak_worst": float(1.0 - min(singular) ** 2),
        "sigma_min": float(min(singular)),
        "T": T, "dt": dt,
    }


def main() -> None:
    B1 = su2(Z, np.pi / 2.0)
    B2 = su2(X, np.pi / 2.0)
    first, second = transfer(B1), transfer(B2)
    U1, U2 = first["Ueff"], second["Ueff"]
    out = {
        "schema": "antler.phase5.synthetic-n2-string-transfer.v1",
        "claim_boundary": (
            "This confirms a local N=2 density-dependent-string transfer with a "
            "matrix link.  A full L=14 protected logical braid, its gap, leakage, "
            "noise, and path audits remain required."
        ),
        "theta": 0.3,
        "string_occupation": 1,
        "B1": {**{k: v for k, v in first.items() if k not in {"Ueff", "expected"}},
               "Ueff": arr(U1), "expected": arr(first["expected"])},
        "B2": {**{k: v for k, v in second.items() if k not in {"Ueff", "expected"}},
               "Ueff": arr(U2), "expected": arr(second["expected"])},
        "commutator_norm": float(np.linalg.norm(U1 @ U2 - U2 @ U1)),
        "braid_relation_residual": float(np.linalg.norm(U1 @ U2 @ U1 - U2 @ U1 @ U2)),
    }
    path = Path("results/phase5/synthetic_n2_string_transfer.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({
        "B1_matrix_error": first["matrix_error"], "B1_leak": first["leak_worst"],
        "B2_matrix_error": second["matrix_error"], "B2_leak": second["leak_worst"],
        "commutator_norm": out["commutator_norm"],
        "braid_relation_residual": out["braid_relation_residual"],
    }, indent=2))


if __name__ == "__main__":
    main()
