"""TeNPy U(1)-conserving rung-MPO for the external Floquet ladder candidate.

The two physical rails are grouped into one four-state fermionic rung.  This
makes the Floquet-conjugated density interaction a nearest-rung coupling while
retaining the total particle-number U(1) charge exactly.  The module is not a
native-ANTLER derivation; it is the canonical-MPS representation of the
already validated external H_eff benchmark.
"""
from __future__ import annotations

import numpy as np

from tenpy.linalg import np_conserved as npc
from tenpy.models.model import CouplingMPOModel
from tenpy.networks.site import FermionSite, GroupedSite


def make_rung_site(conserve_branch_parity: bool = False):
    """Four-state (empty, b, a, ab) rung with charge-preserving rail current."""
    rung = GroupedSite([FermionSite(conserve="N"), FermionSite(conserve="N")], labels=["A", "B"])
    cd_a, c_a = rung.get_op("CdA").to_ndarray(), rung.get_op("CA").to_ndarray()
    cd_b, c_b = rung.get_op("CdB").to_ndarray(), rung.get_op("CB").to_ndarray()
    n_a, n_b = rung.get_op("NA").to_ndarray(), rung.get_op("NB").to_ndarray()
    # P^dag (n_A-n_B) P = +/- Y at eta=pi/2. The sign cancels in Y_j Y_{j+1}.
    y_rail = -1j * (cd_a @ c_b) + 1j * (cd_b @ c_a)
    rung.add_op("Ntot", n_a + n_b, need_JW=False, hc="Ntot")
    rung.add_op("Yrail", y_rail, need_JW=False, hc="Yrail")
    rung.add_op("PA", np.eye(rung.dim) - 2.0 * n_a, need_JW=False, hc="PA")
    if conserve_branch_parity:
        # State order is (empty, B, A, AB). The charge is (N_total, N_A mod 2).
        chinfo = npc.ChargeInfo([1, 2], ["N", "PA"])
        qflat = np.asarray([[0, 0], [1, 0], [1, 1], [2, 1]], dtype=np.int64)
        rung.change_charge(npc.LegCharge.from_qflat(chinfo, qflat))
        # The Jordan-Wigner sign remains (-1)^N_total, i.e. the first charge.
        rung.charge_to_JW_parity = np.asarray([1, 0], dtype=np.int64)
    return rung


class DynamicNumberConservingLadder(CouplingMPOModel):
    """H_eff on a chain of grouped fermionic rungs with exact total U(1)."""

    def init_sites(self, model_params):
        conserve_branch_parity = model_params.get("conserve_branch_parity", False, bool)
        return make_rung_site(conserve_branch_parity)

    def init_terms(self, model_params):
        hopping = model_params.get("t_leg", 1.0, "real")
        interaction = model_params.get("u0", -2.0, "real")
        alpha = model_params.get("alpha", 0.5, "real")
        for u1, u2, dx in self.lat.pairs["nearest_neighbors"]:
            # Intraleg hopping is invariant under the global rail rotation P.
            self.add_coupling(-hopping, u1, "CdA", u2, "CA", dx, plus_hc=True)
            self.add_coupling(-hopping, u1, "CdB", u2, "CB", dx, plus_hc=True)
            # alpha H0: density interactions on the two physical legs.
            self.add_coupling(alpha * interaction, u1, "NA", u2, "NA", dx)
            self.add_coupling(alpha * interaction, u1, "NB", u2, "NB", dx)
            # (1-alpha) P^dag H0 P at eta=pi/2:
            # (U0/2) [N_j N_{j+1} + Y_j Y_{j+1}].
            self.add_coupling((1.0 - alpha) * interaction / 2.0, u1, "Ntot", u2, "Ntot", dx)
            self.add_coupling((1.0 - alpha) * interaction / 2.0, u1, "Yrail", u2, "Yrail", dx)
