<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/MiguelGamezL/EPRdating/main/docs/img/logo/eprdating_logo_dark.png">
    <img src="https://raw.githubusercontent.com/MiguelGamezL/EPRdating/main/docs/img/logo/eprdating_logo.png" alt="EPRdating" width="560">
  </picture>
</p>

# EPRdating

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23192194.svg)](https://doi.org/10.5281/zenodo.23192194)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MiguelGamezL/EPRdating/blob/main/examples/quickstart.ipynb)

Open ESR (EPR) dating of tooth enamel in Python, from raw measurements to an
age with its uncertainty:

```
EPR spectra ──► intensities ──► dose response ──► De ──────┐
                                                           ├──► age: EU, LU, US, US-ESR, CSUS-ESR (Monte Carlo)
HPGe spectra ─► sediment U, Th, K ─► dose rate ────────────┘
                                         ▲
tooth U-series data ─► uptake model ─────┘
```

- **EPR spectra:** Bruker (BES3T, ESP/WinEPR), Freiberg MS5000 and text files; baseline,
  power/gain/mass normalisation, template fits (empirical or simulated with
  EPRAYA) with noise-injection errors; SSE, EXP+LIN, DSE and linear dose
  response with bootstrap De.
- **Sediment:** HPGe gamma spectra (ORTEC, Canberra, N42, text) by the
  comparative method against IAEA RGU-1/RGTh-1/RGK-1, with automatic energy
  calibration and a 226Ra/238U equilibrium check.
- **Dose rate:** conversion factors (Adamiec & Aitken 1998, Guérin et al.
  2011, Liritzis et al. 2013), water, cosmic, U-series ingrowth with radon
  loss, one-group beta attenuation per U-series segment in
  sediment/cementum/enamel/dentine layers, energy-dependent alpha efficiency.
- **Age:** EU, LU and US uptake, combined U-series/ESR (US-ESR) and CSUS-ESR;
  Monte Carlo over every input, geometry included.
- **Reproducibility:** every validation case and every recomputed published
  study is stored as JSON with the provenance of each value.

Inspired by ROSY (Brennan et al. 1997, 1999) and DATA (Grün 2009), written as
an independent, open implementation from the published equations, and
validated against both programs run as black boxes.

> **Version 0.1.1** (see `CHANGELOG.md`). Validation (details in
> `tests/validation/`):
>
> | Against | Agreement |
> |---|---|
> | ROSY 2.0: EU, LU and CU ages, incl. the six teeth of Brennan et al. (1997) | −0.6 to +0.7 % (`ROSY_FINDINGS.md`) |
> | DATA: EU and LU ages, one beta factor for the U chain as DATA | −2.6 to +1.6 % (`DATA_FINDINGS.md`) |
> | DATA: US-ESR and CS-US ages | −1.4 to +0.7 %; CS-US within DATA's 1 ka rounding |
> | 58 published ages from nine studies (DATA, USESR, ROSY) | −9.5 to +3.4 %, 56 within the published 1σ (`PUBLISHED_FINDINGS.md`) |
>
> An article describing the library is in preparation. See *Known
> limitations*.

## Install

```bash
pip install eprdating              # core: numpy + scipy
pip install "eprdating[plot]"      # + matplotlib, for the figures
pip install "eprdating[gui]"       # + ipywidgets, for the optional interface
```

The development version: `pip install "eprdating @ git+https://github.com/MiguelGamezL/EPRdating"`.

From a clone of the repository:

```bash
pip install -e .                 # core: numpy + scipy only
pip install -e '.[plot]'         # + matplotlib, for the figures
pip install -e '.[spectra]'      # + EPRAYA backend for simulated component shapes
pip install -e '.[gamma-formats]'  # + becquerel: Canberra .cnf, ORTEC .spc, IEC 61455
pip install -e '.[dev]'          # + pytest, ruff
```

EPRAYA (Física Aplicada group, Universidad Nacional de Colombia) is optional
because it pulls JAX, numba, ipywidgets and tkinter.

## Quick start

EPRdating can be used in two ways that run the same code: calling its
functions from any Python script or notebook, or through an optional
interface for Jupyter and Google Colab. The interface is a convenience
layer; everything it does is available as functions.

### With code

From spectra to De: read the spectra of each aliquot (repeated
measurements are averaged), measure the intensity in a field window and fit
the dose response.

```python
from eprdating import fit_dose_response
from eprdating.spectra import IntensityWindow, combined_intensity, read_epr

window = IntensityWindow(100, center_g=2.0023)   # the default; or center_mT=336.5, unit="mT"
aliquots = {                                     # added dose (Gy): files of that aliquot
    0: ["nat_rot1.DSC", "nat_rot2.DSC"],
    100: ["100Gy_rot1.DSC", "100Gy_rot2.DSC"],
    300: ["300Gy_rot1.DSC", "300Gy_rot2.DSC"],
    1000: ["1000Gy_rot1.DSC", "1000Gy_rot2.DSC"],
}
doses, I, S = [], [], []
for dose, files in aliquots.items():
    r = combined_intensity([read_epr(f) for f in files], "peak_to_peak", window=window)
    doses.append(dose); I.append(r.value); S.append(r.sigma)

drc = fit_dose_response(doses, I, model="SSE", sigma=S)
print(drc.De, drc.De_sigma)
```

`read_epr` reads Bruker, Freiberg MS5000, `.dat`/`.par` and text files. The
intensity can also be `"t1_b2"`, `"double_integral"` or `"template"` (a
line-shape fit; `empirical_template` builds the template from the strongest
spectra). `tests/validation/calio_ms5000.py` is a complete example on a
public dose series.

From De to age:

```python
from eprdating import ToothSample, ToothLayers, Sediment, LinearUptake, cosmic_dose_rate

sample = ToothSample(
    De=(drc.De, drc.De_sigma),
    enamel_U=(0.8, 0.08), dentine_U=(15.0, 1.5),
    sediment=Sediment(U=(2.1, 0.1), Th=(7.5, 0.4), K=(1.1, 0.05), water=(0.15, 0.05)),
    beta=ToothLayers(enamel_um=1100, dentine_um=2000, strip_outer_um=50, strip_inner_um=50),
    cosmic=cosmic_dose_rate(depth_m=1.5, density=1.9, lat_deg=4.6, lon_deg=-74.1, altitude_m=2600),
    uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake(),
    u234_u238_dentine=(1.25, 0.02), radon_loss_dentine=(0.3, 0.1),
)
print(sample.age().summary())
print(sample.age_mc(n=2000, seed=42).summary())
```

Every input accepts a number, a `(value, sigma)` tuple or a `Value`.

### Interface (optional)

For Jupyter and Google Colab, an interactive interface runs the whole chain
without code: load spectra, set the doses and the intensity window, tick
the points of the dose response, choose the model, compare the intensity
methods, analyse the HPGe spectra of the sediment and compute the age. The
De and the sediment flow into the age tab; every tab exports its results and
the settings (JSON) that reproduce them with the functions above. It needs
the `gui` extra; the rest of the package does not.

```python
%pip install "eprdating[gui]"
from eprdating.gui import app
app()
```

<img src="https://raw.githubusercontent.com/MiguelGamezL/EPRdating/main/docs/img/gui.png" alt="EPRdating interface: equivalent-dose tab" width="640">

The tabs are also objects that can be driven from code (`a = app()`, then
`a.de.load(files)`, `a.de.compute()`, `a.de.fit()`, `a.de.compare_methods()`).

Units: Gy, Gy/ka, ka, ppm (U, Th), % (K), mT, GHz.

**Documentation:** [overview](https://github.com/MiguelGamezL/EPRdating/blob/main/docs/index.md) · [user guide](https://github.com/MiguelGamezL/EPRdating/blob/main/docs/guide.md) ·
[figures](https://github.com/MiguelGamezL/EPRdating/blob/main/docs/plotting.md) · [API reference](https://github.com/MiguelGamezL/EPRdating/blob/main/docs/api/index.md) ·
[validation](https://github.com/MiguelGamezL/EPRdating/blob/main/docs/validation.md). Build the site locally with
`pip install -e '.[docs]' && mkdocs serve`.

See also
`examples/quickstart.py` and `examples/deconvolution_epraya.py`.

`examples/dose_series_dat.py` runs a real additive-dose series from raw
spectra to De: it reads `.dat`/`.par` files, builds the CO2- template (the
field-aligned average of the strong spectra, or the orthorhombic CO2- radical
simulated with EPRAYA and broadened by the time constant), averages the
repeated spectra of each aliquot, fits the amplitude inside the intensity
window (100 G around g = 2.0023 by default, `--window-G`, `--center-g`) with a
field-shift search and noise-injection errors, and fits the dose-response
line.

`examples/norm_gamma.py` goes from HPGe spectra of a sediment and the IAEA
reference materials to U, Th, K and the infinite-matrix dose rates.

## Modules

| module | content |
|---|---|
| `doseresponse` | additive-dose fits (SSE, EXPLIN, DSE, LIN), weighting, Dmax cut, covariance and bootstrap De |
| `dose_rate` | conversion factors (Adamiec & Aitken 1998, Guérin et al. 2011, Liritzis et al. 2013), water correction, alpha efficiency, cosmic dose rate (Prescott & Hutton 1994) |
| `uptake` | US model `U(t) = U_m (t/T)^(p+1)`; EU (p = −1), LU (p = 0) |
| `series` | U-series ingrowth after uptake (234U, 230Th, 231Pa), radon loss, measured or initial 234U/238U |
| `onegroup` | one-group (double-P0) beta transport in planar layers (O'Brien et al. 1964; Brennan et al. 1997), per emitter and per U-series segment |
| `beta` | fixed beta geometry factors, as an alternative to `onegroup` |
| `alpha` | energy-dependent alpha efficiency, k ∝ R(E)/E (ROSY's "varies with energy" option) |
| `usesr` | combined U-series/ESR (US-ESR): solves the age and the uptake parameter *p* of each tissue from its 230Th/234U and 234U/238U (Grün et al. 1988; Shao et al. 2015); CSUS-ESR (Grün 2000) with `age(model="CSUS")`; Monte Carlo with the fraction of draws that solve |
| `gamma` | HPGe gamma spectrometry of sediments: ORTEC `.Spe`/`.Chn`, N42, ASCII and column files (`.cnf`/`.spc`/IEC through becquerel); automatic energy/resolution calibration per spectrum (no first guess, absorbs gain drift), peak areas, comparative method against IAEA RGU-1/RGTh-1/RGK-1, 226Ra/238U equilibrium check, output as `Sediment` |
| `plot` | figures: stacked spectra with fits, dose-response with De, dose-rate budget, age distribution, gamma spectra and per-line contents |
| `age` | generic solver `∫₀ᵀ Ḋ(t) dt = De` and the `ToothSample` model with Monte Carlo; `beta_by_segment=False` applies one beta factor to the whole U chain, as DATA and USESR do |
| `spectra` | reading Bruker BES3T (`.DSC`/`.DTA`) and ESP/WinEPR (`.par`/`.spc`), Freiberg MS5000 (`.xml`), `.dat`/`.par` and column files (`read_epr`); baseline, power/gain/mass normalisation, field alignment, pseudo-modulation and time-constant broadening of simulated shapes; intensity window (default 100 G around g = 2.0023, configurable); intensities by template fit, peak-to-peak, T1–B2 or double integral, all with noise-injection errors; weighted average of repeated spectra of an aliquot with a repeatability check; template/component fits with field shift; EPRAYA backend |

## Known limitations

1. **U-series conventions differ between programs.** Daughter ingrowth
   (234U, 230Th, 231Pa), radon loss and a measured 234U/238U per tissue are
   modelled, with segment fractions from Adamiec & Aitken (1998)
   (`tools/derive_u_series_partition.py`); *p* comes from U-series data with
   `USESRSample` (US-ESR, CSUS-ESR). The p–T relation reproduces DATA and
   USESR. Ages differ from theirs only through the beta doses: EPRdating
   attenuates each U-series segment with the one-group method (as ROSY), so
   its ages come out 1–7 % younger than USESR's; DATA's single chain factor
   is available with `beta_by_segment=False` (see `DATA_FINDINGS.md`,
   `USESR_FINDINGS.md`, `PUBLISHED_FINDINGS.md`). Samples at the
   closed-system U-series bound may have no nominal solution; the Monte
   Carlo then reports the fraction of draws that solve and flags the age as
   *marginal* below 80 %.
2. **Material compositions are inputs, not uncertainties.** Enamel is
   hydroxyapatite, sediment silica and dentine 70 % mineral, 20 % collagen,
   10 % water by default; any composition can be given (`compound`,
   `mixture`, `sediment_material`, `dentine_material`). Thicknesses,
   stripping, densities and water contents are sampled in the Monte Carlo,
   compositions are not; their effect on the beta factors is ≤ 2–3 %
   (calcite instead of quartz: −2 % sediment factor; dentine mineral 60–80 %:
   < 1 %).
3. **Environmental histories are piecewise constant.** Sediment water,
   cosmic dose rate (from a burial-depth history, `cosmic_history`) and
   gamma can change during burial (`History`); each segment value is
   sampled in the Monte Carlo, the break times are fixed. The internal
   components (enamel, dentine, cementum water and the beta factors of the
   tooth) keep their present-day values.
4. **Reference data.** The three conversion-factor sets match the papers and
   an independent transcription (R package Luminescence); Guérin et al.
   (2011) give no uncertainties, so the relative ones of Adamiec & Aitken
   (1998) are carried, as in DRAC. The cosmic dose rate (Prescott & Hutton
   1994, DRAC's tables) agrees with Luminescence's within 2.5 %; under less
   than ~1.5 hg/cm² (about 0.8 m of sediment) the two fits of the soft
   component differ by up to 7 % (`tests/test_reference_data.py`).
5. **Alpha escape at layer surfaces** is optional (`alpha_escape=True`,
   needs a `ToothLayers` geometry; off by default, as in DATA). It matters
   only for enamel not stripped beyond the alpha range (12–40 µm): ~1 % of
   the enamel alpha dose for 300 µm, plus the alphas entering from dentine
   and sediment. Tracks are straight with constant stopping per emitter;
   the incoming alpha dose from dentine is ~12 % above ROSY's (a
   component below 1 % of the total in realistic cases).
6. **Radon loss from the sediment** is entered as a 226Ra deficit
   (`Sediment(U_ra226=...)`) and lowers only the 226Ra-onward part of the U
   chain (about 59 % of its beta dose). Published USESR values behave as if
   the loss removed the whole U-chain beta (Rising Star, `PUBLISHED_FINDINGS.md`).
7. **Dentine thickness.** When it is unknown the default is 2000 µm, within
   ~1 % of an infinitely thick dentine for the beta dose to the enamel.
8. **No JEOL readers** yet (Bruker BES3T, ESP/WinEPR, Freiberg MS5000 and text files are read).

## Roadmap

- JEOL spectrometer files.
- Releases on PyPI; article in *Quaternary Geochronology* (in preparation).

What each release contains: `CHANGELOG.md`.

## Validation against ROSY

`tools/rosy_harness` batch-runs the original ROSY 2.0 under Wine (bring your
own copy). Reference inputs and outputs live in
`tests/validation/rosy_reference/`; `tests/validation/test_rosy_reference.py`
locks in the agreement for gamma, cosmic, U-series ingrowth and end-to-end
ages, and `ROSY_FINDINGS.md` summarises the campaign, including the beta
geometry factors that the one-group solver has to reproduce.

## Validation against DATA

`tools/data_harness` batch-runs the original DATA (DOS) in DOSBox-X (bring
your own copy). Reference inputs, parsed outputs and raw printouts live in
`tests/validation/data_reference/`; `test_data_reference.py` checks gamma,
cosmic, internal and beta dose rates, DATA's printed beta factors and all 82
EU/LU ages. `DATA_FINDINGS.md` documents DATA's input conventions and the
one real difference: DATA attenuates the beta dose with ingrowth by one
chain factor, whereas EPRdating by default uses one factor per U-series
segment (`ToothSample(beta_by_segment=False)` reproduces DATA).
`run_useries.py` drives DATA's U-series/ESR screen (<F6>);
`test_data_useries.py` checks US-ESR ages and p values, CS-US ages and three
DATA limitations: the dentine p does not always converge, there is no result
near the closed-system bound, and CS-US ages can be younger than the uptake.

## Published ages

`tests/validation/published/` holds the inputs and published results of 18
studies, transcribed from the papers with the provenance of every value and
checked independently. `published_cases.py` recomputes the nine complete
ones with the conventions of the program each study used.
`test_published_ages.py` locks the agreement in, and
`tools/plot_published_ages.py` draws the figure.

## Testing

```bash
pytest                      # unit tests
pytest -m validation -rs    # published cases; skips list what data is missing
pytest -m epraya            # EPRAYA integration (needs the extra)
```

## References

- Brennan B.J. et al. (1997) Beta doses in tooth enamel by "one-group" theory and the ROSY ESR dating software. *Radiation Measurements* 27, 307–314.
- Grün R. (2009) The DATA program for the calculation of ESR age estimates on tooth enamel. *Quaternary Geochronology* 4, 231–232.
- O'Brien K., Samson S., Sanna R., McLaughlin J.E. (1964) The application of "one-group" transport theory to β-ray dosimetry. *Nuclear Science and Engineering* 18, 90–96.
- Prestwich W.V., Chan G.H. (2000) Beta dose scaling and the one group theory. *Radiation Physics and Chemistry* 59, 221–227.
- Grün R., Schwarcz H.P., Chadam J. (1988) ESR dating of tooth enamel: coupled correction for U-uptake and U-series disequilibrium. *Nuclear Tracks and Radiation Measurements* 14, 237–241.
- Durcan J.A., King G.E., Duller G.A.T. (2015) DRAC: Dose Rate and Age Calculator for trapped charge dating. *Quaternary Geochronology* 28, 54–61.
- Prescott J.R., Hutton J.T. (1994) Cosmic ray contributions to dose rates for luminescence and ESR dating. *Radiation Measurements* 23, 497–500.
- Adamiec G., Aitken M. (1998) Dose-rate conversion factors: update. *Ancient TL* 16, 37–50. doi:10.26034/la.atl.1998.292
- Guérin G., Mercier N., Adamiec G. (2011) Dose-rate conversion factors: update. *Ancient TL* 29, 5–8. doi:10.26034/la.atl.2011.443
- Brennan B.J., Rink W.J., Rule E.M., Schwarcz H.P., Prestwich W.V. (1999) The ROSY ESR dating program. *Ancient TL* 17, 45–53. doi:10.26034/la.atl.1999.307
- Grün R. (2000) An alternative model for open system U-series/ESR age calculations: (closed system U-series)-ESR, CSUS-ESR. *Ancient TL* 18, 1–4. doi:10.26034/la.atl.2000.313
- Shao Q., Bahain J.-J., Dolo J.-M., Falguères C. (2014) Monte Carlo approach to calculate US-ESR age and age uncertainty for tooth enamel. *Quaternary Geochronology* 22, 99–106.
- Shao Q., Chadam J., Grün R., Falguères C., Dolo J.-M., Bahain J.-J. (2015) The mathematical basis for the US-ESR dating method. *Quaternary Geochronology* 30, 1–8.
- Carvajal E., Montes L., Almanza O.A. (2011) Quaternary dating by electron spin resonance (ESR) applied to human tooth enamel. *Earth Sciences Research Journal* 15(2), 115–120. https://revistas.unal.edu.co/index.php/esrj/article/view/27715

## Logo

The logo is the CO2- spectrum of tooth enamel, the largest of an
additive-dose series, drawn without smoothing; its double minimum recalls a
molar. The mark of the Grupo de Física Aplicada (Universidad Nacional de
Colombia) sits in the corner. `tools/logo/make_logo.py` redraws it (Barlow
font, SIL Open Font License, `tools/logo/OFL.txt`).

## Authors

Miguel Enrique Gámez López ([ORCID](https://orcid.org/0000-0001-7831-3291)), Carol Jiseth Ospina Umaña,
Juan Sebastián Castro Millán and Ovidio Amado Almanza Montero —
Grupo de Física Aplicada, Departamento de Física, Universidad Nacional de
Colombia, Sede Bogotá. Contact: megamezl@unal.edu.co

## How to cite

An article describing the library is in preparation. Until then, cite the
software (GitHub's "Cite this repository" gives the same from `CITATION.cff`):

> Gámez López M.E., Ospina Umaña C.J., Castro Millán J.S., Almanza Montero O.A.
> (2026) EPRdating: open ESR dating of tooth enamel in Python, version 0.1.0.
> Zenodo. https://doi.org/10.5281/zenodo.23192195

https://doi.org/10.5281/zenodo.23192194 always points to the latest version;
cite the DOI of the version you used.

## License

MIT (see `LICENSE`).
