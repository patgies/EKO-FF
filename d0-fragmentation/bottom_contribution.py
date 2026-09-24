"""Quantify the bottom-quark channel's contribution to the Kniehl-Kramer
D0 fragmentation function, at pT = 3 and 9 GeV.

Usage:
    python bottom_contribution.py
"""

import numpy as np

import physics as ph
from evolve_d0 import evolve_kk

PT_VALUES = [3.0, 9.0]
SCALE_FACTORS = [0.25, 0.5, 1.0, 2.0, 4.0]


def q_values_for_pt(pt):
    mt = np.sqrt(ph.MC**2 + pt**2)
    return [max(sf * mt, ph.MC) for sf in SCALE_FACTORS], mt


def main(cores=4):
    for pt in PT_VALUES:
        q_list, mt = q_values_for_pt(pt)
        zgrid, total, parts = evolve_kk(q_list, n_integration_cores=cores)
        mask = (zgrid > 0.1) & (zgrid < 0.9)

        print(f"\npT = {pt:.0f} GeV  (m_T = {mt:.3f} GeV, bottom channel active for Q >= {ph.MB_FF0} GeV)")
        header = f"{'Q (GeV)':>10s} {'charm D(0.5)':>13s} {'bottom D(0.5)':>14s} {'bottom/total @z=0.5':>21s} {'bottom/total, z-avg (0.1,0.9)':>30s}"
        print(header)
        print("-" * len(header))
        for q in sorted(total):
            c = parts["charm"][q]
            b = parts["bottom"].get(q, np.zeros_like(zgrid))
            mid = len(zgrid) // 2
            frac_mid = b[mid] / total[q][mid] if total[q][mid] != 0 else 0.0
            with np.errstate(invalid="ignore", divide="ignore"):
                frac = np.where(total[q][mask] != 0, b[mask] / total[q][mask], np.nan)
            frac_avg = np.nanmean(frac)
            has_bottom = "  (below bottom threshold)" if q not in parts["bottom"] else ""
            print(
                f"{q:10.3f} {c[mid]:13.5f} {b[mid]:14.5f} {frac_mid * 100:20.4f}% {frac_avg * 100:29.4f}%{has_bottom}"
            )


if __name__ == "__main__":
    main()
