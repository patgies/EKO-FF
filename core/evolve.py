import pathlib
import tempfile
from math import nan

import numpy as np
from eko import interpolation
from eko.io import EKO, runcards
from eko.io.types import ReferenceRunning
from eko.runner import managed
from ekobox import apply

ALPHAS_MZ = 0.118
MZ = 91.1876
MASSES = (1.5, 5.0, 173.0)
"""mc, mb, mt (GeV)"""

DEFAULT_ZGRID = np.linspace(0.05, 1.0 - 1e-6, 200)


class SeededFlavor:

    def __init__(self, pid, func):
        self.pid = pid
        self.func = func

    def hasFlavor(self, pid):
        return pid in (self.pid, -self.pid)

    def xfxQ2(self, _pid, x, _q2):
        return x * float(self.func(np.array([x]))[0])


def nf_for_q(q, masses=MASSES):
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
    """
    Returns (zgrid, {Q: D(z, Q) array}). Q values below `init_scale` are
    clamped to it.
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
