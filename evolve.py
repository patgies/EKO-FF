"""eko-based time-like DGLAP evolution of the BCFY/Kniehl-Kramer -> D0
fragmentation functions.
"""

import pathlib
import tempfile
from math import nan

import numpy as np
from eko import interpolation
from eko.io import EKO, runcards
from eko.io.types import ReferenceRunning
from eko.runner import managed
from ekobox import apply

import physics as ph

ALPHAS_MZ = 0.118
MZ = 91.1876
MASSES = (ph.MC, ph.MB_FF0, 173.0)
"""mc, mb, mt (GeV) -- standard VFNS thresholds. mb = MB_FF0 = 5.0 GeV,
matching hep-ph/0607306's own convention (mu0=mb=5 GeV is stated there as
where modern PDF sets place the flavour threshold), and the native-QCDNUM
port (KK_D0.cc), which uses the same value for both purposes."""

DEFAULT_ZGRID = np.linspace(0.05, 1.0 - 1e-6, 200)


class SeededFlavor:
    """lhapdf-like object seeding a single quark/antiquark pair (pid, -pid)
    with the same z-shape D(z) at the operator's initial scale."""

    def __init__(self, pid, func):
        self.pid = pid
        self.func = func

    def hasFlavor(self, pid):
        return pid in (self.pid, -self.pid)

    def xfxQ2(self, _pid, x, _q2):
        return x * float(self.func(np.array([x]))[0])


def nf_for_q(q, masses=MASSES):
    """Number of active flavours at scale q (GeV), given u/d/s always
    active and thresholds at `masses`."""
    nf = 3
    for m in masses:
        if q >= m:
            nf += 1
    return nf


def _build_cards(init_scale, mus, n_integration_cores, n_xgrid):
    mc, mb, mt = MASSES
    theory = runcards.TheoryCard.from_dict(
        dict(
            order=[1, 0],
            couplings=dict(alphas=ALPHAS_MZ, alphaem=0.007496, ref=(MZ, 5)),
            heavy=dict(
                masses=[
                    ReferenceRunning([mc, nan]),
                    ReferenceRunning([mb, nan]),
                    ReferenceRunning([mt, nan]),
                ],
                masses_scheme="POLE",
                matching_ratios=[1.0, 1.0, 1.0],
            ),
            xif=1.0,
            n3lo_ad_variation=(0, 0, 0, 0, 0, 0, 0),
            matching_order=[0, 0],
            use_fhmruvv=True,
        )
    )
    operator = runcards.OperatorCard.from_dict(
        dict(
            init=(init_scale, nf_for_q(init_scale)),
            mugrid=[(q, nf_for_q(q)) for q in mus],
            xgrid=interpolation.lambertgrid(n_xgrid, 1e-3, 1.0).tolist(),
            configs=dict(
                evolution_method="iterate-exact",
                ev_op_max_order=[10, 0],
                ev_op_iterations=10,
                interpolation_polynomial_degree=4,
                interpolation_is_log=True,
                scvar_method=None,
                inversion_method=None,
                n_integration_cores=n_integration_cores,
                polarized=False,
                time_like=True,
            ),
            debug=dict(skip_singlet=False, skip_non_singlet=False),
        )
    )
    return theory, operator


def evolve_flavor(
    pid, ic_func, init_scale, q_values, zgrid=None, n_integration_cores=4, n_xgrid=100
):
    """DGLAP-evolve a single quark-seeded fragmentation function `ic_func`
    (seeded on flavour `pid`, at scale `init_scale`) up to each of
    `q_values`.

    Returns (zgrid, {Q: D(z, Q) array}). Q values below `init_scale` are
    clamped to it (the fragmentation function isn't defined below its own
    starting scale).
    """
    if zgrid is None:
        zgrid = DEFAULT_ZGRID
    q_values = sorted({round(max(q, init_scale), 8) for q in q_values})
    theory, operator = _build_cards(init_scale, q_values, n_integration_cores, n_xgrid)

    with tempfile.TemporaryDirectory() as tmp:
        eko_path = pathlib.Path(tmp) / "eko.tar"
        managed.solve(theory, operator, path=eko_path)
        with EKO.read(eko_path) as eko_output:
            pdfs, _ = apply.apply_pdf(
                eko_output, SeededFlavor(pid, ic_func), targetgrid=zgrid
            )

    return zgrid, {q: pdfs[(q**2, nf_for_q(q))][pid] for q in q_values}


def evolve_bcfy(q_values, zgrid=None, n_integration_cores=4):
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

    _, dp = evolve_flavor(4, ph.bcfy_dp, ph.MC, q_values, zgrid, n_integration_cores)
    _, dv = evolve_flavor(4, ph.bcfy_dv, ph.MC, q_values, z_scaled, n_integration_cores)

    out = {q: 0.168 * dp[q] + np.where(step, 0.39 * dv[q], 0.0) for q in dp}
    return zgrid, out


def evolve_kk(q_values, zgrid=None, n_integration_cores=4):
    """DGLAP-evolved Kniehl-Kramer D0 fragmentation function: charm channel
    (from mu0=mc) plus bottom channel (from mu0=MB_FF0, only defined/evolved
    for Q >= MB_FF0), summed linearly.

    Returns (zgrid, total, parts) where parts = {"charm": {Q: D(z,Q)},
    "bottom": {Q: D(z,Q)}} (bottom only has entries for Q >= MB_FF0).
    """
    if zgrid is None:
        zgrid = DEFAULT_ZGRID

    _, charm = evolve_flavor(4, ph.kk_charm_ic, ph.MC, q_values, zgrid, n_integration_cores)

    bottom_qs = [q for q in q_values if q >= ph.MB_FF0]
    bottom = {}
    if bottom_qs:
        _, bottom = evolve_flavor(
            5, ph.kk_bottom_ic, ph.MB_FF0, bottom_qs, zgrid, n_integration_cores
        )

    total = {q: charm[q] + bottom.get(q, np.zeros_like(zgrid)) for q in charm}
    return zgrid, total, {"charm": charm, "bottom": bottom}
