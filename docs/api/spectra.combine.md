# EPR spectra: repeated spectra

`eprdating.spectra.combine`

Repeated spectra of the same aliquot: averaging before the intensity.

A weak signal measured several times (several files, several scans each) is
better measured once on the average of all its scans than separately on
- `each`: the template fit and its field-shift search need a visible signal.
`combine_spectra` builds that average:

1. every scan is normalised to a common microwave power (square-root law)
   and to unit receiver gain;
2. spectra recorded at different microwave frequencies are put on a common
   g-scale (the field axis is scaled by the frequency ratio);
3. optionally (`align=True`), each file (not each scan) is aligned in
   field to the others by cross-correlation inside the intensity window, to
   remove tuning shifts between measurements whose frequency was not
   recorded. Off by default: on weak signals the cross-correlation follows
   the noise, and the misaligned average loses 6-15 % of its amplitude in
   synthetic tests, more than the shifts it would correct; the template fit
   already searches a common shift. Use it for clear signals shifted by a
   sizeable fraction of the line width (a 0.3 mT shift of a 0.6 mT wide
   line costs 18 % of the amplitude without it, about 1 % with it);
4. every scan is weighted by `1/sigma**2`, with `sigma` its noise outside
   the intensity window, so noisier scans count less;
5. optionally, scans that deviate from the others by more than their noise
   allows (spikes, jumps) are rejected.

The result is a weighted **mean**, not a sum: a sum would raise the
intensity of aliquots measured more often and distort the dose-response
curve, while the mean keeps every aliquot on the same scale and still gains
`sqrt(N)` in signal-to-noise. No smoothing is applied: smoothing distorts
the line shape and amplitude and correlates the noise, and the template fit
already acts as the matched filter for a known line shape.

Only spectra recorded with the same sweep (number of points, field step,
modulation and time constant, when known) are averaged; the lock-in time
constant deforms the line differently for a different field step. Repeats
with different sweeps are measured separately and combined with
`combine_intensities`.

`combined_intensity` measures the average and checks the repeats
against each other: each file is also measured on its own, and if they
scatter more than their errors (repositioning of the tube in the cavity,
which matters for anisotropic enamel, or spectrometer drift) the error of
the combined intensity is inflated by the Birge ratio.

**Contents:** [`Combination`](#combination), [`combine_spectra`](#combine_spectra), [`combined_intensity`](#combined_intensity)

### `Combination`

*dataclass* `Combination(spectrum: Spectrum, members: list[str], weights: np.ndarray, noise: np.ndarray, shifts_mT: np.ndarray, rejected: list[str] = <factory>)`

Weighted average of repeated spectra, with what went into it.

- `spectrum`: the average, normalised to `ref_power_mW` and unit gain, on the field grid of the first spectrum.
- `members`: one label per scan, `"name#k"`.
- `weights`: normalised weight of each scan (0 for rejected ones).
- `noise`: noise of each scan outside the intensity window, after normalisation.
- `shifts_mT`: field shift applied to each file (alignment).
- `rejected`: labels of the rejected scans.

**Members**

- `summary(self) -> str`

### `combine_spectra`

`combine_spectra(spectra: Sequence[Spectrum], window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), *, ref_power_mW: float | None = None, align: bool = False, max_shift: float = 0.6, weighting: str = 'noise', reject_chi2: float | None = None) -> Combination`

Weighted average of repeated spectra of one aliquot (see the module notes).

- `window`: intensity window; the noise of each scan is measured outside it and the alignment is done inside it.
- `ref_power_mW`: common microwave power (default: that of the first spectrum; no power normalisation if it is unknown).
- `align`: align each file to the others by cross-correlation within `max_shift` mT (off by default; see the module notes).
- `weighting`: `"noise"` (`1/sigma**2`, default) or `"equal"`.
- `reject_chi2`: reject a scan whose mean squared deviation from the weighted mean of the other scans, inside the window and in units of its own noise, exceeds this value (e.g. 2; about 1 is expected). Scans are rejected one at a time, the worst first, and at least two are kept. `None` (default) keeps every scan.

### `combined_intensity`

`combined_intensity(spectra: Sequence[Spectrum], method: str = 'template', template: np.ndarray | tuple[np.ndarray, np.ndarray] | Callable | None = None, window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), *, ref_power_mW: float | None = None, align: bool = False, reject_chi2: float | None = None, weighting: str = 'noise', **intensity_kw) -> Intensity`

Intensity of repeated spectra of one aliquot, measured on their average.

With two or more spectra, each one is also measured on its own (with the
same normalisation); `chi2_red` of their scatter is reported and, when
above 1, the error of the combined intensity is multiplied by
`sqrt(chi2_red)`. `intensity_kw` go to
`intensity` (`max_shift`,
`baseline_order`, `n_noise`, `seed`, `mass_mg`).
