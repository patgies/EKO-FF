 """Export EKO-evolved BCFY and Kniehl-Kramer (charm-only)

Usage: python export_lhapdf.py
Writes ../../diffractive-D0-UPC/data/bcfy_eko/bcfy_eko_0000.dat and
       ../../diffractive-D0-UPC/data/kk_eko/kk_eko_0000.dat
"""

import pathlib
import time

import numpy as np

import evolve_d0 as evolve
import physics as ph

OUT_ROOT = pathlib.Path("../../diffractive-D0-UPC/data")

# z-grid: matches diffractive-D0-UPC's fragmentation range exactly
# (param.zmin=0.05, param.zmax=1.0 in src/main.cpp/main_xpom.cpp).
ZGRID = np.linspace(0.05, 1.0 - 1e-6, 200)

# Q-grid: log-spaced from mc up to comfortably above the highest frag_scale
# actually used anywhere in the project (mt0*SCALE_FACTOR, mt0=sqrt(pD0^2+mc^2)
# up to pD0~20 GeV in fixedW_spectrum.py, SCALE_FACTOR up to 2.0 -> ~40 GeV;
# 50 GeV leaves headroom). lhapdf_grid.cpp clamps outside this range rather
# than failing, but clamping means frozen (wrong) evolution, so better to
# just cover the range.
QGRID = np.geomspace(ph.MC, 50.0, 30).tolist()

N_XGRID = 80          # internal eko x-grid resolution (default 100; 60-80
                       # is a large speed win with no visible accuracy loss
                       # for this purpose -- see session timing test)
N_INTEGRATION_CORES = 4


def write_lhagrid1(path, zgrid, q_to_dz):
    """q_to_dz: {Q: D(z,Q) array matching zgrid}."""
    path.parent.mkdir(parents=True, exist_ok=True)
    qvals = sorted(q_to_dz)
    with open(path, "w") as f:
        f.write("PdfType: central\n")
        f.write("Format: lhagrid1\n")
        f.write("---\n")
        f.write(" ".join(f"{z:.8e}" for z in zgrid) + "\n")
        f.write(" ".join(f"{q:.8e}" for q in qvals) + "\n")
        f.write("4\n")   # single flavor column: charm (PDG id 4)
        for iz, z in enumerate(zgrid):
            for q in qvals:
                xf = z * float(q_to_dz[q][iz])
                f.write(f"{xf:.8e}\n")
        f.write("---\n")
    print(f"Wrote {path} ({len(zgrid)} z x {len(qvals)} Q)")


def main():
    t0 = time.time()
    print(f"Evolving BCFY over {len(QGRID)} Q points...")
    _, bcfy = evolve.evolve_bcfy(
        QGRID, zgrid=ZGRID, n_integration_cores=N_INTEGRATION_CORES, n_xgrid=N_XGRID
    )
    print(f"  done ({time.time()-t0:.1f}s)")
    write_lhagrid1(OUT_ROOT / "bcfy_eko" / "bcfy_eko_0000.dat", ZGRID, bcfy)

    t1 = time.time()
    print(f"Evolving KK charm channel over {len(QGRID)} Q points (bottom dropped)...")
    _, kk_charm = evolve.evolve_flavor(
        4, ph.kk_charm_ic, ph.MC, QGRID, zgrid=ZGRID,
        n_integration_cores=N_INTEGRATION_CORES, n_xgrid=N_XGRID,
    )
    print(f"  done ({time.time()-t1:.1f}s)")
    write_lhagrid1(OUT_ROOT / "kk_eko" / "kk_eko_0000.dat", ZGRID, kk_charm)

    print(f"Total: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
