# EPR spectra: intensity of a spectrum

`eprdating.spectra.measure`

EPR intensity of one spectrum inside an intensity window, with its error.

`intensity` measures a `Spectrum` with one
of four methods, always inside an `IntensityWindow`
(default: 100 G around g = 2.0023):

* `"template"` (recommended): amplitude of a line-shape template fitted
  with a linear baseline and a small common field shift
  (`ComponentBasis`);
* `"peak_to_peak"`, `"t1_b2"` and `"double_integral"`, for comparison.
  For the double integral the derivative baseline is a line fitted to the
  signal-free sweep within one window width on each side of the window, and
  the absorption baseline a line through the outer 20 % of the window at
  each end.

The part of the sweep outside the window is taken as signal-free: after a
polynomial baseline it gives the noise, which is injected into the measured
(or fitted) spectrum to obtain the error of every method.

Detection. A template fit that searches the field position finds, in pure
noise, the place where the noise looks most like the signal, and returns a
positive amplitude. For weak signals this inflates the intensity. The
template method therefore also fits the template, with the same shift
search, to signal-free noise blocks of the same spectrum, and reports the
false-alarm probability `p_noise`: how often noise alone gives at least
the measured amplitude. The noise realisations are phase-randomised
surrogates of the signal-free sweep: they keep its spectrum, i.e. the
correlation from the time constant and the slow baseline wander, and give
independent draws even when the signal-free stretch is short. An aliquot
whose `p_noise` is not small (e.g.
`Intensity.detected(0.01)` is False) has no detectable signal, and its
intensity mostly measures the noise. For this the
longer signal-free side must hold at least as many points as the window;
otherwise the template error falls back to the residual autocovariance and
the other methods get no error (`nan`).

`combine_intensities` averages intensities of the same aliquot that
cannot be averaged as spectra (e.g. different sweep widths).

**Contents:** [`METHODS`](#methods), [`Intensity`](#intensity), [`combine_intensities`](#combine_intensities), [`empirical_template`](#empirical_template), [`intensity`](#intensity)

### `METHODS`

```python
METHODS = ('template', 'peak_to_peak', 't1_b2', 'double_integral')
```

### `Intensity`

*dataclass* `Intensity(value: float, sigma: float, method: str, window_mT: tuple[float, float], shift_mT: float = 0.0, fit: DeconvolutionResult | None = None, n_repeats: int = 1, chi2_red: float | None = None, repeats: list[Intensity] = <factory>, p_noise: float | None = None)`

An EPR intensity with its 1-sigma error and how it was obtained.

- `value`, `sigma`: intensity and error, in the units of the spectrum.
- `method`: one of `"template"`, `"peak_to_peak"`, `"t1_b2"`, `"double_integral"`.
- `window_mT`: field limits actually used.
- `shift_mT`: field shift found by the template fit.
- `fit`: the template fit (`method="template"`).
- `n_repeats`, `chi2_red`, `repeats`: for an intensity combined from repeated measurements, their number, the reduced chi-square of their scatter, and the individual intensities.
- `p_noise`: (`method="template"`) false-alarm probability: the fraction of signal-free noise blocks of the same spectrum that give at least this amplitude with the same fit and shift search. Small values mean the signal is detected; see `detected`.

**Members**

- `detected(self, alpha: float = 0.01) -> bool | None` — Whether the signal stands out of the noise at false-alarm level
`alpha` (None when `p_noise` is not available).

### `combine_intensities`

`combine_intensities(intensities: Sequence[Intensity]) -> Intensity`

Weighted mean of repeated intensities of the same aliquot.

The error is inflated by the Birge ratio, `sqrt(chi2_red)`, when the
repeats scatter more than their errors (repositioning in the cavity,
drift). Use it for repeats that cannot be averaged as spectra; otherwise
`combined_intensity` is better for weak
signals.

### `empirical_template`

`empirical_template(spectra: Sequence[Spectrum], window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), n_strongest: int = 3, max_shift: float = 0.6) -> tuple[np.ndarray, np.ndarray]`

Line-shape template from the strongest spectra of a series.

Each spectrum is normalised to the power and gain of the first, a cubic
baseline fitted outside the window is removed, and the `n_strongest`
by peak-to-peak inside the window are averaged after aligning them in
field (`aligned_average`). Only
spectra on the same field grid as the strongest one are used. Returns
`(B, shape)`, ready for `intensity` (`template=(B, shape)`).

### `intensity`

`intensity(spectrum: Spectrum, method: str = 'template', template: np.ndarray | tuple[np.ndarray, np.ndarray] | Callable | None = None, window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), *, ref_power_mW: float | None = None, mass_mg: float | None = None, max_shift: float = 0.6, baseline_order: int = 1, n_noise: int = 300, n_null: int = 1000, seed: int | None = 0) -> Intensity`

Intensity of `spectrum` (its scan average) inside `window`.

- `method`: `"template"` (default), `"peak_to_peak"`, `"t1_b2"` or `"double_integral"`.
- `template`: line shape for `"template"`: a callable `B -> shape`, a `(B, shape)` pair, or an array on the window's grid (e.g. from `aligned_average` or `epraya_backend`).
- `window`: `IntensityWindow` (default 100 G around g = 2.0023).
- `ref_power_mW`: if given, normalise to this microwave power (square-root law) and to unit receiver gain, from the spectrum's own values; `mass_mg` divides by the aliquot mass.
- `max_shift`: common field shift searched by the template fit (mT).
- `baseline_order`: polynomial baseline fitted with the template.
- `n_noise`, `seed`: noise-injection draws for the error.
- `n_null`: noise realisations for the detection test (`p_noise`).
