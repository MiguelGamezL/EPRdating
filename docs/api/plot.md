# Plotting

`eprdating.plot`

Figures for every step: spectra, dose-response, dose rates, ages, gamma.

Needs matplotlib (`pip install "eprdating[plot]"`). Every function draws
on the axes passed as `ax` (or a new figure) and returns the axes, so the
figures can be combined and restyled with ordinary matplotlib calls.


```text
plot_spectra            a series of cw-EPR spectra, stacked
plot_spectrum           one spectrum with its individual scans
plot_dose_response      intensity vs added dose, fit and De
plot_dose_rate          contribution of each dose-rate component
plot_age_distribution   Monte Carlo age distribution
plot_gamma_spectrum     HPGe spectrum with the analysed lines
plot_gamma_lines        content from each gamma line, per group
```


- `Colours`: series follow a fixed categorical order, ordered series (doses) a
single-hue ramp, text stays in neutral inks, and grids are hairlines.

**Contents:** [`plot_age_distribution`](#plot_age_distribution), [`plot_dose_rate`](#plot_dose_rate), [`plot_dose_response`](#plot_dose_response), [`plot_gamma_lines`](#plot_gamma_lines), [`plot_gamma_spectrum`](#plot_gamma_spectrum), [`plot_spectra`](#plot_spectra), [`plot_spectrum`](#plot_spectrum), [`ramp`](#ramp)

### `plot_age_distribution`

`plot_age_distribution(mc, *, ax = None, bins: int = 50, unit: str = 'ka', reference: tuple | None = None, title: str | None = None)`

Histogram of Monte Carlo ages with the 68 % interval.

- `mc`: `AgeMC` or a US-ESR Monte Carlo result.
- `reference`: optional `(low, high)` range to compare with (e.g. the archaeological context), drawn as a band.

### `plot_dose_rate`

`plot_dose_rate(result, *, ax = None, unit: str = 'Gy/ka', title: str | None = None)`

Time-averaged dose rate of each component at the age found.

- `result`: `AgeResult` (`sample.age()`), or a `{name: rate}` mapping.

### `plot_dose_response`

`plot_dose_response(drc, *, excluded: Sequence[tuple] | None = None, ax = None, dose_unit: str = 'Gy', intensity_label: str = 'EPR intensity (a.u.)', show_De: bool = True, title: str | None = None)`

Additive-dose points, fitted curve extrapolated to zero and De.

- `drc`: `DoseResponseResult`.
- `excluded`: points left out of the fit, as `(dose, intensity[, sigma])`; drawn hollow.

### `plot_gamma_lines`

`plot_gamma_lines(result, *, ax = None, title: str | None = None)`

Content obtained from each line, grouped, with the group mean ± 1σ.

Lines that disagree with their group point to interferences or
self-absorption; the 238U group against the 226Ra group shows
disequilibrium.

### `plot_gamma_spectrum`

`plot_gamma_spectrum(spec, *, ax = None, lines = None, energy_range: tuple[float, float] | None = (30, 2700), log: bool = True, label_lines: bool = True, title: str | None = None)`

Calibrated HPGe spectrum with the analysed lines marked by group.

### `plot_spectra`

`plot_spectra(spectra: Sequence, labels: Sequence[str] | None = None, *, window: tuple[float, float] | None = None, offset: float | None = None, fits: Sequence | None = None, baseline: bool = True, ax = None, title: str | None = None)`

Stacked spectra of a dose series, light to dark with the dose.

- `spectra`: `Spectrum` objects or `(B, y)` pairs.
- `labels`: text at the right of each trace (e.g. `"40 Gy"`).
- `window`: field range (mT) to show.
- `offset`: vertical spacing; default 1.2 × the largest peak-to-peak.
- `fits`: optional fitted curves, `(B, y)` pairs or fit results with a `fitted` array on the same window (drawn over the data).
- `baseline`: subtract each trace's mean inside the window.

### `plot_spectrum`

`plot_spectrum(spectrum, *, window: tuple[float, float] | None = None, show_scans: bool = True, fit = None, ax = None, title: str | None = None)`

One spectrum: individual scans in grey, their average, and an
optional fitted curve (`(B, y)` or a fit result on the same window).

### `ramp`

`ramp(n: int) -> list[str]`

`n` colours from the single-hue ramp, light to dark (for doses).
