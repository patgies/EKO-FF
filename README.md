# eko-based DGLAP evolution of fragmentation functions

Generic time-like DGLAP evolution (via [`eko`](https://github.com/NNPDF/eko))
of analytically-parametrized fragmentation functions, plus one project's
worth of physics built on top: charm/bottom → D0 fragmentation models used
for D0-meson production.

- **BCFY** — E. Braaten, K. Cheung, S. Fleming, T.C. Yuan,
  [Phys. Rev. D51 (1995) 4819](https://arxiv.org/abs/hep-ph/9409316):
  separate pseudoscalar (c → D0) and vector (c → D\*0 → D0) channels.
- **Kniehl–Kramer** — B.A. Kniehl, G. Kramer,
  [hep-ph/0607306](https://arxiv.org/abs/hep-ph/0607306): a charm channel
  plus a separate bottom-quark channel.

Both models are defined analytically at a starting scale μ₀ (m_c for the
charm channels, 5 GeV for the Kniehl–Kramer bottom channel) and evolved
with the LO time-like DGLAP equations.


## Layout

```
core/
  evolve.py    generic eko driver, no physics-specific content: seeds a
               single flavour at its starting scale (evolve_flavor()) and
               solves the time-like DGLAP evolution. Reusable as-is by any
               future project -- import this unchanged, don't fork it.

d0-fragmentation/   this project: charm/bottom -> D0 fragmentation
  physics.py         analytic BCFY/KK parametrizations (D0/D* masses, etc.)
  evolve_d0.py        D0-specific evolution: evolve_bcfy()/evolve_kk(),
                       built on core/evolve.py's evolve_flavor()
  export_lhapdf.py    writes eko-evolved BCFY/KK grids in LHAPDF lhagrid1
                       format, for diffractive-D0-UPC's src/bcfy_grid.cpp /
                       kk_grid.cpp to read (see that repo's README)
  plot_evolution.py   pT=3/9 GeV evolution plots for BCFY and Kniehl-Kramer
  bottom_contribution.py  quantifies the KK bottom channel vs. charm channel

A future project (different physics, e.g. a different meson/parametrization)
gets its own sibling folder next to d0-fragmentation/, importing the same
unchanged core/evolve.py -- no risk of two projects' physics colliding.
```

## Usage

```bash
pip install -r requirements.txt
cd d0-fragmentation
python plot_evolution.py             # writes output/ff_evolution.pdf
python bottom_contribution.py        # prints a table to stdout
python export_lhapdf.py              # writes the LHAPDF grids diffractive-D0-UPC reads
```



## Notes

- Both BCFY channels (P and V) start at μ₀ = m_c and are evolved as two
  *independent* DGLAP sets — the D\*0 → D0 decay is a hadron-level branching
  ratio (0.39) plus a momentum rescaling (z → z·m_D\*/m_D), not part of the
  QCD evolution, so it's applied only after both pieces are evolved.
- The Kniehl–Kramer bottom channel starts at μ₀ = 5 GeV (its own fitted
  scale) and is only defined/evolved for Q ≥ 5 GeV; below that it
  contributes zero by construction.
- All evolution is LO, time-like, with u/d/s always active, charm active
  from m_c=1.5 GeV, bottom from m_b=5.0 GeV, top from m_t=173 GeV. (m_b
  doubles as both the VFNS coupling/kernel threshold and the Kniehl–Kramer
  bottom channel's own starting scale — matching hep-ph/0607306's own
  convention and the native-QCDNUM port, KK_D0.cc.)

## Bottom channel contribution

`bottom_contribution.py` shows the Kniehl–Kramer bottom channel is
**not** negligible once Q crosses m_b: ~27–37% of the total D0
fragmentation function, both at pT=3 and pT=9 GeV. That's sizeable, and
consistent with hep-ph/0607306's own branching-fraction table, where a
bottom quark's branching to D0 (~53–58% at M_Z) is the same order as
charm's own (~66–68%).

