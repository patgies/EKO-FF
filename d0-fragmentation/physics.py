"""Analytic charm/bottom -> D0 fragmentation-function parametrizations.

BCFY: E. Braaten, K. Cheung, S. Fleming, T.C. Yuan, "Perturbative QCD
Fragmentation Functions as a Model for Heavy-Quark Fragmentation",
Phys. Rev. D51 (1995) 4819, arXiv:hep-ph/9409316. Gives the pseudoscalar
(c -> D0 direct) and vector (c -> D*0, feeding down to D0) channels
separately; they must be DGLAP-evolved independently and combined after
evolution (see evolve.evolve_bcfy).

Kniehl-Kramer: B.A. Kniehl, G. Kramer, hep-ph/0607306. LO charm fit
(Peterson-like) plus a separate bottom-quark fit; see evolve.evolve_kk for
how the two channels are combined.
"""

import numpy as np

R = 0.1
"""Light-quark/meson mass ratio in the BCFY formulas."""

MD, MDSTAR = 1.8648, 2.0067
"""D0, D*0 masses (GeV)."""

MASS_RATIO = MDSTAR / MD

MC = 1.5

MB_FF0 = 5.0
"""Kniehl-Kramer bottom-fragmentation starting scale mu0 (GeV, hep-ph/0607306)."""


def _bcfy_poly_p(z, r):
    return (
        6.0
        - 18.0 * (1.0 - 2.0 * r) * z
        + (21.0 - 74.0 * r + 68.0 * r**2) * z**2
        - 2.0 * (1.0 - r) * (6.0 - 19.0 * r + 18.0 * r**2) * z**3
        + 3.0 * (1.0 - r) ** 2 * (1.0 - 2.0 * r + 2.0 * r**2) * z**4
    )


def _bcfy_poly_v(z, r):
    return (
        2.0
        - 2.0 * (3.0 - 2.0 * r) * z
        + 3.0 * (3.0 - 2.0 * r + 4.0 * r**2) * z**2
        - 2.0 * (1.0 - r) * (4.0 - r + 2.0 * r**2) * z**3
        + (1.0 - r) ** 2 * (3.0 - 2.0 * r + 2.0 * r**2) * z**4
    )


def bcfy_dp(z, r=R):
    """Pseudoscalar (c -> D0 direct) BCFY fragmentation function, at mu0=mc."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    mask = (z > 0.0) & (z < 1.0)
    zz = z[mask]
    denom = 1.0 - (1.0 - r) * zz
    factor = r * zz * (1.0 - zz) ** 2 / denom**6
    out[mask] = factor * _bcfy_poly_p(zz, r)
    return out


def bcfy_dv(z, r=R):
    """Vector (c -> D*0) BCFY fragmentation function, at mu0=mc."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    mask = (z > 0.0) & (z < 1.0)
    zz = z[mask]
    denom = 1.0 - (1.0 - r) * zz
    factor = 3.0 * (r * zz * (1.0 - zz) ** 2) / denom**6
    out[mask] = factor * _bcfy_poly_v(zz, r)
    return out


def bcfy_d0(z, r=R):
    """Un-evolved c -> D0 fragmentation function (BCFY model), at mu=mc.

    D_D0(z) = 0.168 * D_P(z) + 0.39 * theta(mD/mDstar - z) * D_V(z * mDstar/mD)
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    dp_part = 0.168 * bcfy_dp(z, r)
    z_scaled = MASS_RATIO * z
    step = (MD / MDSTAR - z) >= 0.0
    dv_part = np.where(step, 0.39 * bcfy_dv(z_scaled, r) * MASS_RATIO, 0.0)
    return dp_part + dv_part


KK_CHARM_N, KK_CHARM_EPS = 0.694, 0.101


def kk_charm_ic(z):
    """Kniehl-Kramer charm -> D0 fragmentation function, at mu0=mc."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    mask = (z > 0.0) & (z < 1.0)
    zz = z[mask]
    denom = (1.0 - zz) ** 2 + KK_CHARM_EPS * zz
    out[mask] = KK_CHARM_N * zz * (1.0 - zz) ** 2 / denom**2
    return out


KK_BOTTOM_N, KK_BOTTOM_ALPHA, KK_BOTTOM_BETA = 81.7, 1.81, 4.95


def kk_bottom_ic(z):
    """Kniehl-Kramer bottom -> D0 fragmentation function, at mu0=MB_FF0."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    mask = (z > 0.0) & (z < 1.0)
    zz = z[mask]
    out[mask] = KK_BOTTOM_N * zz**KK_BOTTOM_ALPHA * (1.0 - zz) ** KK_BOTTOM_BETA
    return out
