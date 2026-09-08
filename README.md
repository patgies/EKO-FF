# DGLAP evolution of charm/bottom → D0 fragmentation functions

Time-like DGLAP evolution (via [`eko`](https://github.com/NNPDF/eko)) of two
charm-fragmentation models used for D0-meson production:

- **BCFY** — E. Braaten, K. Cheung, S. Fleming, T.C. Yuan,
  [Phys. Rev. D51 (1995) 4819](https://arxiv.org/abs/hep-ph/9409316):
  separate pseudoscalar (c → D0) and vector (c → D\*0 → D0) channels.
- **Kniehl–Kramer** — B.A. Kniehl, G. Kramer,
  [hep-ph/0607306](https://arxiv.org/abs/hep-ph/0607306): a charm channel
  plus a separate bottom-quark channel.

Both models are defined analytically at a starting scale μ₀ (m_c for the
charm channels, 5 GeV for the Kniehl–Kramer bottom channel) and evolved
upward with the LO time-like DGLAP equations.

## Why eko, not QCDNUM

This is a from-scratch reimplementation of an evolution originally done with
QCDNUM (see the `inclusive-D0-UPC` project). Two things came out of
cross-checking eko against that QCDNUM setup that are worth recording here:

1. **Coupling normalization.** QCDNUM's setup there pins n_f=4 as active
   forever above m_c (bottom/top thresholds pushed out of range, so the
   evolution never needs to cross them), which means its
   `setalf(0.118, M_Z²)` call ends up defining α_s(M_Z)=0.118 **at n_f=4**,
   not the usual n_f=5 world-average reading — an unintended side effect of
   freezing the flavour scheme for the DGLAP kernel, not a deliberate choice.
   This repo uses the standard **n_f=5 VFNS** reading instead (real b/t
   quark-mass thresholds), which is ~1–10% different in α_s over the
   1.5–35 GeV range relevant here.
2. **A residual few-percent gap** between eko and QCDNUM remains even after
   matching the coupling (median 3–9% relative, depending on channel and
   scale), most likely in the quark–gluon mixing sector of the timelike
   splitting functions. Not fully isolated — flagged here for anyone
   revisiting this.

## Layout

```
dglap_ff/
  physics.py   analytic fragmentation-function parametrizations
  evolve.py    eko driver: seeds a single flavour at its starting scale,
               solves the time-like DGLAP operator, applies it
plot_evolution.py      pT=3/9 GeV evolution plots for BCFY and Kniehl-Kramer
bottom_contribution.py quantifies the KK bottom channel's size vs. the charm channel
```

## Usage

```bash
pip install -r requirements.txt
python plot_evolution.py             # writes ff_evolution.png
python bottom_contribution.py        # prints a table to stdout
```

Each evolution point (one DGLAP operator solve per set of target scales) is
the expensive step — expect the full pT=3+9 GeV grid (5 scales each) to take
on the order of tens of minutes on a laptop; `--cores` controls how many
processes eko uses for the Mellin-space integration.

## Physics notes

- Both BCFY channels (P and V) start at μ₀ = m_c and are evolved as two
  *independent* DGLAP sets — the D\*0 → D0 decay is a hadron-level branching
  ratio (0.39) plus a momentum rescaling (z → z·m_D\*/m_D), not part of the
  QCD evolution, so it's applied only after both pieces are evolved.
- The Kniehl–Kramer bottom channel starts at μ₀ = 5 GeV (its own fitted
  scale) and is only defined/evolved for Q ≥ 5 GeV; below that it
  contributes zero by construction.
- All evolution is LO, time-like, with u/d/s always active, charm active
  from m_c=1.5 GeV, bottom from m_b=4.5 GeV, top from m_t=173 GeV.
