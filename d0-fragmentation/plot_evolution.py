"""Plot the BCFY and Kniehl-Kramer D0 fragmentation functions, DGLAP-evolved
to a few reference scales at pT = 3 and 9 GeV

Usage:
    python plot_evolution.py [--out ff_evolution.pdf] [--cache ff_evolution_data.npz]
"""

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np

import parametrizations as ph
from evolve_d0 import evolve_bcfy, evolve_kk

PT_VALUES = [3.0, 9.0]
SCALE_FACTORS = [0.25, 0.5, 1.0, 2.0, 4.0]

# House style (matches inclusive-D0-UPC/diffractive-D0-UPC python/*.py plots)
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 11,
        "axes.labelsize": 20,
        "axes.titlesize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 8,
        "ytick.major.size": 8,
        "xtick.minor.size": 4,
        "ytick.minor.size": 4,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
    }
)

# ColorBrewer Blues (5-class, dropping the palest step for visibility on
# white): wider perceptual spacing than an evenly-sampled colormap slice.
COLORS = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]


def q_values_for_pt(pt):
    mt = np.sqrt(ph.MC**2 + pt**2)
    return [max(sf * mt, ph.MC) for sf in SCALE_FACTORS], mt


def compute(cores):
    all_q = sorted({q for pt in PT_VALUES for q in q_values_for_pt(pt)[0]})
    bcfy_z, bcfy_curves = evolve_bcfy(all_q, n_integration_cores=cores)
    kk_z, kk_curves, _ = evolve_kk(all_q, n_integration_cores=cores)
    return bcfy_z, bcfy_curves, kk_z, kk_curves


def save_cache(path, bcfy_z, bcfy_curves, kk_z, kk_curves):
    payload = {"bcfy_z": bcfy_z, "kk_z": kk_z}
    for q, d in bcfy_curves.items():
        payload[f"bcfy_{q}"] = d
    for q, d in kk_curves.items():
        payload[f"kk_{q}"] = d
    np.savez(path, **payload)


def load_cache(path):
    data = np.load(path)
    bcfy_z = data["bcfy_z"]
    kk_z = data["kk_z"]
    bcfy_curves = {
        float(k.split("_", 1)[1]): data[k]
        for k in data.files
        if k.startswith("bcfy_") and k != "bcfy_z"
    }
    kk_curves = {
        float(k.split("_", 1)[1]): data[k]
        for k in data.files
        if k.startswith("kk_") and k != "kk_z"
    }
    return bcfy_z, bcfy_curves, kk_z, kk_curves


def plot(out, bcfy_z, bcfy_curves, kk_z, kk_curves):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True)

    rows = [
        ("BCFY", bcfy_z, bcfy_curves, ph.bcfy_d0),
        (r"Kniehl--Kramer", kk_z, kk_curves, ph.kk_charm_ic),
    ]

    for row, (name, zgrid, curves, ic_fn) in enumerate(rows):
        for col, pt in enumerate(PT_VALUES):
            ax = axes[row, col]
            q_list, mt = q_values_for_pt(pt)
            ax.plot(
                zgrid, ic_fn(zgrid), color="dimgray", ls=":", lw=1.5,
                label=r"IC ($\mu=m_c$)", zorder=10,
            )
            for i, (sf, q) in enumerate(zip(SCALE_FACTORS, q_list)):
                key = min(curves, key=lambda k: abs(k - q))
                ax.plot(
                    zgrid, curves[key], color=COLORS[i], lw=2.2,
                    label=rf"${sf:g}\times m_T={key:.3g}$ GeV",
                )
            ax.set_title(rf"{name}, $p_T={pt:.0f}$ GeV", pad=10)
            ax.set_xlim(0.05, 1.0)
            ax.set_ylim(bottom=0)
            if row == 1:
                ax.set_xlabel(r"$z$", labelpad=10)
            if col == 0:
                ax.set_ylabel(r"$D(z, Q^2)$", labelpad=10)
            if row == 0 and col == 0:
                ax.legend(fontsize=14, frameon=False, loc="upper left")
                # Separate caption, not a legend entry -- avoids the blank
                # handle indent a fake legend row would need.
                ax.text(
                    0.03, 0.50, r"$m_T=\sqrt{m_c^2+p_T^2}$",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=13, color="black",
                )

    fig.suptitle(
        r"DGLAP-evolved charm $\to D^0$ fragmentation functions (eko)",
        fontsize=18,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="output/ff_evolution.pdf")
    parser.add_argument("--cache", default="ff_evolution_data.npz")
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()

    cache_path = pathlib.Path(args.cache)
    if args.recompute or not cache_path.exists():
        data = compute(args.cores)
        save_cache(cache_path, *data)
    else:
        data = load_cache(cache_path)

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plot(out_path, *data)


if __name__ == "__main__":
    main()
