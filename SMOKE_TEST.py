from math import comb
import numpy as np

from antler.basis import build_basis
from antler.model import build_hamiltonian

L, N = 4, 2
states, index = build_basis(2 * L, N)
assert len(states) == comb(2 * L, N)
H, _, _ = build_hamiltonian(
    L=L, N=N, theta=0.37, J1=0.4, J2=1.0, Jperp=0.1,
    mu=np.linspace(-0.2, 0.2, 2 * L), basis=(states, index)
)
err = float(np.linalg.norm(H - H.conj().T))
assert err < 1e-12, err
print(f"PASS: basis dimension={len(states)}, hermiticity_error={err:.3e}")
