# Changelog

## Unreleased

- Optional protocol for enamel fragments measured at several angles
  (`eprdating.spectra.fragments`): orientation and angle from file names,
  orientation-balanced merged spectrum, angular profile, intensity of the
  merged spectrum, X-ray calibration (De fitted in seconds, converted with
  the calibration error). The isotropic correction is not included.
- Alpha escape at the enamel surfaces (`ToothSample(alpha_escape=True)`,
  optional, also in the interface): own alphas leaving the unstripped enamel
  and alphas entering from dentine, sediment or cementum, per emitter
  (`eprdating.alpha.escape_fractions`). Brings unstripped ROSY reference
  ages from up to 5-8 % to within 0.5-2.5 %.
- Interactive interface for Jupyter and Google Colab (`eprdating.gui.app`,
  extra `gui`): equivalent dose, sediment gamma spectrometry and age in
  three tabs, with downloads of results and settings.
- README: quick start from spectra to De with the functions; the interface
  is presented as optional (extra `gui`).
- `eprdating.spectra.empirical_template`: line-shape template from the
  strongest spectra of a series.
- Interface: the window centre can be given as a g value or as a field in
  G or mT; *Compare methods* gives De with the four intensity methods.
- Reader for Freiberg Instruments MS5000 spectra (ESRStudio `.xml`, its
  `.csv` export, or a folder of runs as scans): `read_ms5000`, also through
  `read_epr` and the interface. Repeats whose grids differ by a few points
  are combined. Checked on the public Calio dose series (Hakim et al. 2025):
  the published De is reproduced once the 8000 and 15 000 Gy labels of the
  archive are exchanged (`docs/validation.md`). The Calio spectra are
  included (`tests/data/calio_p4`, CC-BY 4.0) and tested (`tests/test_calio.py`).
- Double integral (`intensity(..., "double_integral")`): derivative baseline
  fitted to the sweep within one window width on each side of the window,
  absorption baseline through the outer 20 % of the window, noise injected
  over the window and both baseline regions. About 40 % less scatter than
  the previous end-point baselines. `double_integral(..., baseline="outside")`
  gives the same on its own.

## 0.1.1 — 2026-10-06

- Published on PyPI: `pip install eprdating`. Releases on GitHub publish to
  PyPI automatically (trusted publishing).
- Tests run on every push (Python 3.10 and 3.13).
- README renders on PyPI (absolute image and link URLs); Zenodo DOI and
  how to cite.
- Quickstart notebook: intensity window, repeated spectra and burial
  history; opens in Google Colab.

## 0.1.0 — 2026-10-06

First public release.

**EPR spectra and equivalent dose**
- Readers: Bruker BES3T (`.DSC`/`.DTA`) and ESP/WinEPR (`.par`/`.spc`),
  `.dat`/`.par` pairs, column text files.
- Baseline, power/gain/mass normalisation, field alignment; pseudo-modulation
  and time-constant broadening of simulated line shapes; EPRAYA backend
  (optional).
- Intensity window (default 100 G around g = 2.0023, configurable);
  intensities by template fit, peak-to-peak, T1–B2 or double integral, with
  noise-injection errors and a detection test (`p_noise`).
- Weighted averaging of repeated spectra of an aliquot with a repeatability
  check.
- Dose-response fits (LIN, SSE, EXPLIN, DSE), weighting, maximum-dose test,
  Birge-scaled or bootstrap De.

**Sediment gamma spectrometry**
- ORTEC `.Spe`/`.Chn`, N42 and text files; Canberra `.cnf`, ORTEC `.spc` and
  IEC 61455 through becquerel (optional).
- Automatic energy and resolution calibration, comparative method against
  IAEA RGU-1/RGTh-1/RGK-1, 226Ra/238U equilibrium check.

**Dose rates and ages**
- Conversion factors of Adamiec & Aitken (1998), Guérin et al. (2011) and
  Liritzis et al. (2013); water correction; cosmic dose rate (Prescott &
  Hutton 1994).
- U-series ingrowth with radon loss and measured or initial 234U/238U.
- One-group beta attenuation per emitter and U-series segment in
  sediment/cementum/enamel/dentine layers; configurable material
  compositions; DATA's single chain factor as an option.
- Energy-dependent alpha efficiency (ROSY's option).
- Piecewise histories of sediment water, burial depth (cosmic) and gamma.
- EU, LU and US uptake; US-ESR and CSUS-ESR; Monte Carlo over every input,
  geometry included.

**Validation**
- ROSY 2.0: EU, LU and CU ages within −0.6 to +0.7 %.
- DATA: EU/LU ages within −2.6 to +1.6 % (82 runs); US-ESR within −1.4 to
  +0.7 % (29 cases); CS-US within DATA's rounding.
- 58 published ages from nine studies within −9.5 to +3.4 %, 56 within the
  published 1σ.
- Conversion factors and cosmic dose rate against the R package
  Luminescence; file readers against EasySpin and becquerel test files.

**Tools**
- Drivers for ROSY (Wine) and DATA (DOSBox-X), without the programs.
- Synthetic additive-dose series to test analysis choices and measurement
  designs (`tools/synthetic_dose_series.py`).
