# EPR spectra: pre-processing

`eprdating.spectra.preprocess`

Pre-processing of cw-EPR spectra before intensity estimation.

* `subtract_baseline`  polynomial baseline fitted outside the signal.
* `aligned_average`    field-aligned average (empirical template).
* `normalise`          receiver gain, microwave power, mass.
* `time_constant_filter`  lock-in time constant (RC) along the sweep.
* `pseudo_modulation`  field-modulation broadening of a simulated
  derivative spectrum (Hyde, Pasenkiewicz-Gierula, Jesmanowicz & Antholine
  1990, *Appl. Magn. Reson.* 1, 483), so simulated shapes can be compared
  with spectra recorded with a large modulation amplitude.

**Contents:** [`aligned_average`](#aligned_average), [`noise_sigma`](#noise_sigma), [`normalise`](#normalise), [`pseudo_modulation`](#pseudo_modulation), [`subtract_baseline`](#subtract_baseline), [`time_constant_filter`](#time_constant_filter)

### `aligned_average`

`aligned_average(B, spectra, window: tuple[float, float], max_shift: float = 0.6, reference: int | None = None) -> tuple[np.ndarray, np.ndarray]`

Average of several spectra after aligning them in field.

Each spectrum is shifted (by whole field steps, up to `max_shift` mT)
to maximise its cross-correlation with the `reference` spectrum (the
one with the largest peak-to-peak inside `window` by default). Useful
to build an empirical line-shape template from the strongest spectra of
a dose series. Returns `(average, shifts_mT)`.

### `noise_sigma`

`noise_sigma(B, y, signal: tuple[float, float]) -> float`

Standard deviation of a baseline-corrected spectrum outside `signal`.

### `normalise`

`normalise(y, power_mW: float | None = None, gain: float | None = None, mass_mg: float | None = None, ref_power_mW: float = 1.0) -> np.ndarray`

Normalise a signal to unit gain, `ref_power_mW` and unit mass.

The signal of a non-saturated line grows as the square root of the
microwave power; the CO2- signal of enamel is not saturated below a few
tens of mW at room temperature. Intensities are divided by `gain` and
`mass_mg` when given.

### `pseudo_modulation`

`pseudo_modulation(B, y, mod_pp_mT: float) -> np.ndarray`

First-harmonic spectrum recorded with modulation amplitude `mod_pp_mT`.

`y` is the first-derivative spectrum for vanishing modulation on the
uniform grid `B` (mT). In Fourier space the modulated first harmonic is
the derivative multiplied by `2 J1(k a)/(k a)`, with `a` half the
peak-to-peak modulation amplitude (Hyde et al. 1990). The derivative
normalisation is kept, so a small modulation returns `y` unchanged;
the signal an instrument records is this times the modulation amplitude.

### `subtract_baseline`

`subtract_baseline(B, y, exclude: tuple[float, float] | None = None, order: int = 3) -> np.ndarray`

Subtract a polynomial fitted to the points outside `exclude` (mT).

`y` may be 1-D or 2-D (one row per scan).

### `time_constant_filter`

`time_constant_filter(y, tau_points: float) -> np.ndarray`

Single-pole (RC) low-pass filter along the sweep, as applied by the
lock-in time constant; `tau_points` is the time constant divided by
the time per point. The filter delays and broadens the line.
