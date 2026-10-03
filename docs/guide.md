# EPRdating user guide

From raw measurements to an ESR age in three steps: the equivalent dose from
the EPR spectra, the sediment dose rate from gamma spectrometry, and the age
with its Monte Carlo uncertainty. Each step can be used on its own.

```bash
pip install eprdating                # core (numpy, scipy)
pip install "eprdating[spectra]"     # + EPRAYA, for simulated line shapes
pip install "eprdating[gamma-formats]"  # + Canberra .cnf, ORTEC .spc, IEC 61455 gamma files
```

Units everywhere: field in mT, frequency in GHz, dose in Gy, dose rate in
Gy/ka, age in ka, U and Th in µg/g (ppm), K in %. Any input can be a number,
a `(value, sigma)` tuple or a `Value`.

---

## 1. Equivalent dose from EPR spectra

### Reading spectra

```python
from eprdating.spectra import read_epr
s = read_epr("M18_3_19mW_4SCAN.dat")     # or .DSC/.DTA, .par/.spc, .csv
s.B, s.scans, s.y                        # field (mT), scans, scan average
s.freq_GHz, s.power_mW, s.gain, s.mod_amp_mT, s.time_constant_ms
```

| Files | Instrument / software |
|---|---|
| `.DSC` + `.DTA` | Bruker BES3T (Xepr: E500, E580, EMXplus, EMXmicro) |
| `.par` + `.spc` | Bruker ESP300 / WinEPR (EMX, ECS106) |
| `.dat` + `.par` | `KEY : value` parameters + five-column ASCII (UNAL X band) |
| `.txt` / `.csv` | a field column and one or more signal columns |

Bruker readers follow EasySpin's `eprload` and are checked on its test files.
Anything else can be loaded by hand: `Spectrum(B=..., scans=y[None, :])`.

### Intensities

Weak dating signals are best measured by fitting a line-shape template rather
than by peak-to-peak heights or double integrals of noisy spectra.

```python
from eprdating.spectra import ComponentBasis, normalise, subtract_baseline, aligned_average

y = normalise(s.y, power_mW=s.power_mW, ref_power_mW=19.0)      # also gain=, mass_mg=
win = (s.B > 333.1) & (s.B < 341.1)
basis = ComponentBasis(s.B[win], {"CO2-": template}, baseline_order=1)
noise = subtract_baseline(s.B, y, exclude=(331.4, 342.9))[s.B < 331.4]   # signal-free stretch
r = basis.fit(y[win], nonnegative=False, max_shift=0.6, noise=noise)
r.amplitudes["CO2-"], r.errors["CO2-"], r.shift
```

* `template`: an empirical average of the strongest spectra
  (`aligned_average`) or a simulation (`eprdating.spectra.epraya_backend`,
  broadened with `pseudo_modulation` and `time_constant_filter`).
* `max_shift` lets each spectrum find its own field position (tuning or
  frequency drift between measurements).
* `noise=` gives errors by noise injection, which account for the noise
  correlation of the time constant and for the shift search.
* `nonnegative=False` for dose-response work: the constraint biases weak
  amplitudes upwards.

### Dose-response and De

```python
from eprdating import fit_dose_response
drc = fit_dose_response(doses, amplitudes, "LIN", sigma=errors, De_min=-np.inf)
drc.De, drc.De_sigma, drc.chi2_red
```

Models `"LIN"`, `"SSE"`, `"EXPLIN"`, `"DSE"`. With `sigma`, errors are
inflated by the Birge ratio when `chi2_red > 1` (scatter between aliquots).
`De_min=-np.inf` shows where the data really extrapolate.

`examples/dose_series_dat.py` runs a complete series.

---

## 2. Sediment U, Th, K from gamma spectrometry (HPGe)

The comparative method: sample and reference materials measured in the same
container and geometry.

```python
from eprdating.gamma import analyse_files, IAEA_RGU_1, IAEA_RGTH_1, IAEA_RGK_1

r = analyse_files(
    "soil.Spe", mass_g=600.0,
    references={"U": ("RGU1.Spe", IAEA_RGU_1),
                "Th": ("RGTh1.Spe", IAEA_RGTH_1),
                "K": ("RGK1.Spe", IAEA_RGK_1)},
    background="empty_container.Spe",
)
print(r.summary())
sediment = r.sediment(water=0.10)   # -> ToothSample(sediment=...)
```

* **Files**: ORTEC `.Spe`/`.Chn`, ANSI N42.42 (`.n42`, `.xml`), ASCII
  spectra with a `# key: value` header, plain columns
  (`read_gamma(path, live_time=...)`), and through becquerel Canberra
  `.cnf`, ORTEC `.spc` and IEC 61455.
* **Calibration**: automatic, on each spectrum's own natural lines; no first
  guess is needed and gain drift between measurements is absorbed. A
  quadratic term is kept only when the data need it.
* **Lines**: K from 40K; U from 214Pb/214Bi (226Ra, equilibrium assumed);
  Th from 228Ac, 212Pb, 208Tl. 234Th and 234mPa give 238U directly, and
  `r.activity_ratio_ra226_u238` tests the equilibrium (seal the containers
  for 3–4 weeks before measuring so that radon grows back). If it departs
  from 1, `r.sediment(water, equilibrium=False)` uses 238U for the top of
  the chain and 226Ra for the rest (`Sediment(U=..., U_ra226=...)`); the
  dose rates then follow each part of the chain.
* **Masses**: `mass_g` is the dry mass of the sample; references default to
  500 g (`Reference(..., mass_g=...)`). Only mass is corrected: keep the
  fill height and density close to those of the references.
* **Water**: `water_content(wet_mass, dry_mass)` gives the water/dry mass
  ratio used by `Sediment`.
* **Not supported**: scintillators (NaI, CsI, LaBr3) — the line-by-line
  method needs HPGe resolution, and `auto_calibrate` refuses such spectra.

Your own reference materials: `Reference("my-standard", {"U": (12.3, 0.4)}, mass_g=450)`.

`examples/norm_gamma.py` prints contents and infinite-matrix dose rates.

---

## 3. Age

```python
from eprdating import ToothSample, ToothLayers, LinearUptake, cosmic_dose_rate

sample = ToothSample(
    De=(drc.De, drc.De_sigma),
    enamel_U=(0.8, 0.08), dentine_U=(15.0, 1.5),
    sediment=sediment,
    beta=ToothLayers(enamel_um=1100, strip_outer_um=50, strip_inner_um=50),
    cosmic=cosmic_dose_rate(depth_m=1.5, density=1.9, lat_deg=10.25, lon_deg=-73.4, altitude_m=150),
    uptake_enamel=LinearUptake(), uptake_dentine=LinearUptake(),
)
print(sample.age().summary())
print(sample.age_mc(n=2000, seed=1).summary())
```

With U-series data of the dental tissues use `USESRSample` (combined
U-series/ESR). See the README for the physics and its validation against
ROSY 2.0 and published US-ESR ages.

---

## Validation data

The tests under `tests/validation` that need external files are skipped
unless these variables point to them:

| Variable | Data |
|---|---|
| `EPRDATING_NORM_DATA` | UNAL NORM spectra (corte 0, IAEA references) |
| `EPRDATING_BQ_SAMPLES` | `tests/samples` of [becquerel](https://github.com/lbl-anp/becquerel) (HPGe from other laboratories) |
| `EPRDATING_EASYSPIN_FILES` | `tests/eprfiles` of [EasySpin](https://github.com/StollLab/EasySpin) (Bruker files) |
