"""Phase 5D: verify noncommuting matrix links inside an ANTLER ED Hamiltonian.

The synthetic SU(2) design is embedded in the new spinful ladder Hamiltonian.
Two independent digital handoffs are propagated with exact hopping and Strang
potential steps.  This is a one-particle integration test of the extension;
the next stage is the N=2 correlated-hopping logical-code audit.
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

from antler.synthetic import build_synthetic_hamiltonian, su2


I = np.eye(2, dtype=complex)
X = np.array([[0.0, 1.0], [1.0, 0.0]], complex)
Z = np.array([[1.0, 0.0], [0.0, -1.0]], complex)


def no_global(U: np.ndarray) -> np.ndarray:
    return U * np.exp(-0.5j * np.angle(np.linalg.det(U)))


def arr(A: np.ndarray) -> dict:
    return {"real": A.real.tolist(), "imag": A.imag.tolist()}


def global_phase_distance(U: np.ndarray, V: np.ndarray) -> float:
    """Distance after removing every U(1) phase, including the SU(2) centre -I."""

    relative = V.conj().T @ U
    overlap = np.trace(relative)
    phase = 1.0 if abs(overlap) < 1e-15 else np.exp(-1j * np.angle(overlap))
    return float(np.linalg.norm(phase * relative - I))


def align_to_embedding(frame: np.ndarray, embedding: np.ndarray) -> np.ndarray:
    """Fix the degenerate-frame gauge against a physical spin basis."""

    return frame @ polar(embedding.conj().T @ frame)[0].conj().T


def local_frame(H: np.ndarray, embedding: np.ndarray) -> np.ndarray:
    energies, vectors = np.linalg.eigh(H)
    scores = np.sum(abs(embedding.conj().T @ vectors) ** 2, axis=0)
    chosen = np.argsort(-scores)[:2]
    return align_to_embedding(vectors[:, chosen], embedding)


def transfer(U_link: np.ndarray, T: float = 2400.0, dt: float = 0.05) -> dict:
    """Run a closed-frame two-site digital handoff on the full L=14 ED basis."""

    L, N, D = 14, 1, 10.0
    mu_zero = np.zeros(2 * L)
    # Link (0,2) acts on the decreasing hop 2 -> 0; source 0 -> target 2
    # therefore implements its Hermitian conjugate.
    Hhop, states, _ = build_synthetic_hamiltonian(
        # Isolate the controlled (0,2) bond within the full L=14 spinful
        # Hilbert space.  The next N=2 audit restores the complete SSH graph.
        L, N, theta=0.0, J1=1.0, J2=0.0, Jperp=0.0, mu=mu_zero,
        link_matrices={(0, 2): U_link},
    )
    Hhop = Hhop.toarray()
    dim = len(states)
    source = np.zeros((dim, 2), complex)
    target = np.zeros((dim, 2), complex)
    for spin in range(2):
        source_index = next(i for i, state in enumerate(states) if int(state) == (1 << spin))
        target_index = next(i for i, state in enumerate(states) if int(state) == (1 << (4 + spin)))
        source[source_index, spin] = 1.0
        target[target_index, spin] = 1.0
    mu_initial = mu_zero.copy(); mu_initial[0] = -D
    mu_final = mu_zero.copy(); mu_final[2] = -D
    H0 = Hhop + np.diag(np.repeat(mu_initial, 2))
    Hf = Hhop + np.diag(np.repeat(mu_final, 2))
    frame0 = local_frame(H0, source)
    framef = local_frame(Hf, target)
    energies, vectors = np.linalg.eigh(Hhop)
    nseg = int(round(T / dt)); dt = T / nseg
    Uhop = (vectors * np.exp(-1j * energies * dt)) @ vectors.conj().T
    psi = frame0.copy()
    for step in range(nseg):
        u = (step + 0.5) / nseg
        q = np.sin(0.5 * np.pi * u) ** 2
        mu = mu_zero.copy(); mu[0] = -D * (1.0 - q); mu[2] = -D * q
        phase = np.exp(-0.5j * dt * np.repeat(mu, 2))[:, None]
        psi = phase * (Uhop @ (phase * psi))
    S = framef.conj().T @ psi
    Ueff = no_global(polar(S)[0])
    expected = no_global(U_link.conj().T)
    singular = np.linalg.svd(S, compute_uv=False)
    return {
        "Ueff": Ueff,
        "expected": expected,
        "matrix_error": global_phase_distance(Ueff, expected),
        "leak_worst": float(1.0 - min(singular) ** 2),
        "sigma_min": float(min(singular)),
        "offdiag_norm": float(np.linalg.norm(Ueff - np.diag(np.diag(Ueff)))),
        "T": T, "dt": dt,
    }


def main() -> None:
    B1 = su2(Z, np.pi / 2.0)
    B2 = su2(X, np.pi / 2.0)
    first, second = transfer(B1), transfer(B2)
    U1, U2 = first["Ueff"], second["Ueff"]
    out = {
        "schema": "antler.phase5.synthetic-ed-link-transfer.v1",
        "claim_boundary": (
            "This validates the new matrix-link Hamiltonian and its single-particle "
            "digital transfer.  It is not yet a two-particle correlated-hopping, "
            "protected logical braid demonstration."
        ),
        "B1": {**{k: v for k, v in first.items() if k not in {"Ueff", "expected"}},
               "Ueff": arr(U1), "expected": arr(first["expected"])},
        "B2": {**{k: v for k, v in second.items() if k not in {"Ueff", "expected"}},
               "Ueff": arr(U2), "expected": arr(second["expected"])},
        "commutator_norm": float(np.linalg.norm(U1 @ U2 - U2 @ U1)),
        "braid_relation_residual": float(np.linalg.norm(U1 @ U2 @ U1 - U2 @ U1 @ U2)),
    }
    path = Path("results/phase5/synthetic_ed_link_transfer.json")
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
