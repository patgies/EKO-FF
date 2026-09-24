# EDO DGLAP evolution of fragmentation functions

Generic time-like DGLAP evolution (via [`eko`](https://github.com/NNPDF/eko))

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
  evolve.py    generic eko driver

d0-fragmentation/   this project: charm/bottom -> D0 fragmentation
  parametrizations.py analytic BCFY/KK parametrizations (D0/D* masses, etc.)
  evolve_d0.py        D0-specific evolution: evolve_bcfy()/evolve_kk(),
                       built on core/evolve.py's evolve_flavor()
  export_grids.py    writes eko-evolved BCFY/KK grids in LHAPDF lhagrid1 format
  plot_evolution.py   pT=3/9 GeV evolution plots for BCFY and Kniehl-Kramer
  bottom_contribution.py  quantifies the KK bottom channel vs. charm channel

```

## Usage

```bash
pip install -r requirements.txt
cd d0-fragmentation
python plot_evolution.py             # writes output/ff_evolution.pdf
python bottom_contribution.py        # prints a table to stdout
python export_grids.py              # writes the LHAPDF grids 
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
  from m_c=1.5 GeV, bottom from m_b=5.0 GeV, top from m_t=173 GeV.


