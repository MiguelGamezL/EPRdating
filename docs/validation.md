# Validation

What has been checked, against what, and where the details are. The tests
under `tests/validation` lock the agreement in; those that need external
files are skipped unless the files are available (see the end of this page).

## Dose rates and ages: ROSY 2.0

The original ROSY 2.0 program (Brennan et al. 1999) was run as a black box on
reference cases. Details: [`tests/validation/ROSY_FINDINGS.md`](../tests/validation/ROSY_FINDINGS.md).

| Piece | Agreement |
|---|---|
| Conversion factors (Adamiec & Aitken 1998), water correction | < 0.1 % |
| Cosmic dose rate | within 1 % |
| U-series ingrowth, radon loss | effective alpha efficiency predicted within −0.2 to +2.6 % |
| One-group beta attenuation (O'Brien et al. 1964, no fitted parameter) | dentine factors to three decimals; self-dose 1–2.5 %; sediment 2–5 % |
| Full ages, EU/LU/CU, 35 ka – 2 Ma, six real teeth of Brennan et al. (1997) | **−0.5 to +0.7 %** |
| Monte Carlo uncertainties vs ROSY's errors | within 25 % (geometry sampled) |

## US-ESR: published ages

Shao et al. (2015) and De Nadale et al. (2026). Details:
[`tests/validation/USESR_FINDINGS.md`](../tests/validation/USESR_FINDINGS.md).

- The uptake parameter *p* at the published age is reproduced within 0.015
  (Shao) and 0.03 (De Nadale, ages rounded to 1 ka).
- Ages agree within the published 1σ; the remaining differences come from
  the beta attenuation of the USESR program, which differs from the
  one-group method.
- A sample at the closed-system U-series bound (CN-55) is flagged as
  marginal by the Monte Carlo instead of being given a spurious age.
- The alternative CSUS-ESR model (Grün 2000) passes a synthetic round trip
  and gives ages 3–4 % older than US-ESR for the De Nadale teeth.

## Gamma spectrometry

- **Independent analysis of the same spectra** (UNAL, HPGe 40 %, sediment
  and IAEA RGU-1/RGTh-1/RGK-1): K, U and Th within 1 %.
- **Spectra from other laboratories** (becquerel project: ORTEC PopTop in a
  lead cave and with a Marinelli beaker; Canberra Falcon 5000 and ORTEC
  Trans-SPEC in the field): the automatic calibration, without any first
  guess, places the 352, 1461 and 2615 keV peaks within 0.5 keV.
  Scintillator spectra (NaI, CsI) are rejected with a clear message.

## File readers

- **Bruker** files from EasySpin's test collection: the same measurement
  stored in two formats (BES3T and ESP; ESP and WinEPR; big and little
  endian) gives the same spectrum.
- **Gamma**: ORTEC `.Spe`, Canberra `.cnf`, ORTEC `.spc` and IEC 61455 files
  read identically to becquerel; N42 files from Sandia's SpecUtils are read
  with their counts, times and calibration.
- ORTEC `.Chn` is tested on files written to the ORTEC specification only.

## EPR intensities

- Template-fit errors by noise injection reproduce the scatter of repeated
  simulated measurements, including time-constant correlated noise and the
  field-shift search on weak signals.
- Pseudo-modulation reproduces the optimum modulation amplitude of a
  Lorentzian line (3.5 ΔBpp).

## Running the validation tests

```bash
pytest -m validation -rs
```

| Variable | Data |
|---|---|
| `EPRDATING_NORM_DATA` | UNAL NORM spectra (not public) |
| `EPRDATING_BQ_SAMPLES` | `tests/samples` of [becquerel](https://github.com/lbl-anp/becquerel) |
| `EPRDATING_EASYSPIN_FILES` | `tests/eprfiles` of [EasySpin](https://github.com/StollLab/EasySpin) |

The ROSY and US-ESR reference data are in the repository.
