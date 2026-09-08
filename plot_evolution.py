"""Plot the BCFY and Kniehl-Kramer D0 fragmentation functions, DGLAP-evolved
to a few reference scales at pT = 3 and 9 GeV (matching the pT points used
for the D0 spectrum in inclusive-D0-UPC).

Usage:
    python plot_evolution.py [--out ff_evolution.png]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from dglap_ff import physics as ph
from dglap_ff.evolve import evolve_bcfy, evolve_kk

PT_VALUES = [3.0, 9.0]
SCALE_FACTORS = [0.25, 0.5, 1.0, 2.0, 4.0]


def q_values_for_pt(pt):
    mt = np.sqrt(ph.MC**2 + pt**2)
    return [max(sf * mt, ph.MC) for sf in SCALE_FACTORS], mt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="ff_evolution.png")
    parser.add_argument("--cores", type=int, default=4)
    args = parser.parse_args()

    all_q = sorted({q for pt in PT_VALUES for q in q_values_for_pt(pt)[0]})

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    colors = plt.cm.Blues(np.linspace(0.35, 0.95, len(SCALE_FACTORS)))

    bcfy_z, bcfy_curves = evolve_bcfy(all_q, n_integration_cores=args.cores)
    kk_z, kk_curves, _ = evolve_kk(all_q, n_integration_cores=args.cores)

    rows = [
        ("BCFY", bcfy_z, bcfy_curves, ph.bcfy_d0),
        ("Kniehl–Kramer", kk_z, kk_curves, ph.kk_charm_ic),
    ]

    for row, (name, zgrid, curves, ic_fn) in enumerate(rows):
        for col, pt in enumerate(PT_VALUES):
            ax = axes[row, col]
            q_list, mt = q_values_for_pt(pt)
            ax.plot(zgrid, ic_fn(zgrid), color="0.55", ls=":", lw=1.3, label="IC (μ=m$_c$)")
            for i, (sf, q) in enumerate(zip(SCALE_FACTORS, q_list)):
                key = min(curves, key=lambda k: abs(k - q))
                ax.plot(zgrid, curves[key], color=colors[i], lw=1.8, label=f"{sf}×m$_T$={key:.3g} GeV")
            ax.set_title(f"{name}, p$_T$={pt:.0f} GeV", fontsize=11)
            ax.set_xlim(0.05, 1.0)
            if row == 1:
                ax.set_xlabel("z")
            if col == 0:
                ax.set_ylabel("D(z, Q)")
            if row == 0 and col == 1:
                ax.legend(fontsize=7.5, frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(args.out, dpi=160)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
