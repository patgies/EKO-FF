"""D0-specific DGLAP evolution: BCFY and Kniehl-Kramer fragmentation
functions, built on top of ../core/evolve.py's generic evolve_flavor().
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "core"))
from evolve import DEFAULT_ZGRID, evolve_flavor  # noqa: E402

import physics as ph  # noqa: E402


def evolve_bcfy(q_values, zgrid=None, n_integration_cores=4, n_xgrid=100):
    """DGLAP-evolved BCFY D0 fragmentation function.

    D_D0(z, Q) = 0.168 * D_P_evolved(z, Q)
               + 0.39 * theta(mD/mDstar - z) * D_V_evolved(z * mDstar/mD, Q)

    P and V are evolved as two independent time-like DGLAP sets (both from
    mu0=mc) and combined after evolution.

    Returns (zgrid, {Q: D(z, Q) array}).
    """
    if zgrid is None:
        zgrid = DEFAULT_ZGRID
    z_scaled = ph.MASS_RATIO * zgrid
    step = (ph.MD / ph.MDSTAR - zgrid) >= 0.0

    _, dp = evolve_flavor(4, ph.bcfy_dp, ph.MC, q_values, zgrid, n_integration_cores, n_xgrid)
    _, dv = evolve_flavor(4, ph.bcfy_dv, ph.MC, q_values, z_scaled, n_integration_cores, n_xgrid)

    out = {q: 0.168 * dp[q] + np.where(step, 0.39 * dv[q], 0.0) for q in dp}
    return zgrid, out


def evolve_kk(q_values, zgrid=None, n_integration_cores=4, n_xgrid=100):
    """DGLAP-evolved Kniehl-Kramer D0 fragmentation function: charm channel
    (from mu0=mc) plus bottom channel (from mu0=MB_FF0, only defined/evolved
    for Q >= MB_FF0), summed linearly.

    Returns (zgrid, total, parts) where parts = {"charm": {Q: D(z,Q)},
    "bottom": {Q: D(z,Q)}} (bottom only has entries for Q >= MB_FF0).
    """
    if zgrid is None:
        zgrid = DEFAULT_ZGRID

    _, charm = evolve_flavor(4, ph.kk_charm_ic, ph.MC, q_values, zgrid, n_integration_cores, n_xgrid)

    bottom_qs = [q for q in q_values if q >= ph.MB_FF0]
    bottom = {}
    if bottom_qs:
        _, bottom = evolve_flavor(
            5, ph.kk_bottom_ic, ph.MB_FF0, bottom_qs, zgrid, n_integration_cores, n_xgrid
        )

    total = {q: charm[q] + bottom.get(q, np.zeros_like(zgrid)) for q in charm}
    return zgrid, total, {"charm": charm, "bottom": bottom}
