# EPR spectra: intensity window and scalar intensities

`eprdating.spectra.intensity`

Scalar ESR intensities from cw (first-derivative) spectra.

The magnetic field is in mT and the microwave frequency in GHz; spectra are
first-derivative cw-EPR traces sampled on increasing field.

The field region used for an intensity is an `IntensityWindow`: by
default 100 G (10 mT) centred on g = 2.0023, i.e. on the field of that g at
each spectrum's own microwave frequency. The window should hold the whole
dating signal with some signal-free margin on both sides for the baseline;
outside it, the spectrum is taken as signal-free and used to estimate the
noise. A wider window may take in other radicals (native signal, CO3-, SO2-,
methyl) that the intensity method does not describe; a narrower one leaves
few points for the baseline.

**Contents:** [`IntensityWindow`](#intensitywindow), [`double_integral`](#double_integral), [`field_for_g`](#field_for_g), [`g_for_field`](#g_for_field), [`outside_baseline`](#outside_baseline), [`peak_to_peak`](#peak_to_peak), [`t1_b2_amplitude`](#t1_b2_amplitude)

### `IntensityWindow`

*dataclass* `IntensityWindow(width: float = 100.0, unit: str = 'G', center_g: float = 2.0023, center_mT: float | None = None)`

Field window in which an EPR intensity is computed.

- `width`: full width of the window, in `unit` (default 100 G).
- `unit`: `"G"` or `"mT"`.
- `center_g`: g-value at the centre (default 2.0023); the centre field is computed from each spectrum's microwave frequency.
- `center_mT`: fixed centre field instead (overrides `center_g`).

For instance `IntensityWindow()` is 100 G around g = 2.0023,
`IntensityWindow(60)` 60 G around it, and
`IntensityWindow(8, unit="mT", center_mT=337.1)` 8 mT around 337.1 mT.

**Members**

- `bounds(self, freq_GHz: float | None) -> tuple[float, float]` — `(low, high)` field limits in mT.
- `center(self, freq_GHz: float | None) -> float` — Centre field (mT) for a spectrum recorded at `freq_GHz`.
- `describe(self, freq_GHz: float | None = None) -> str`
- `half_width_mT` *(property)*
- `mask(self, B, freq_GHz: float | None) -> np.ndarray` — Boolean mask of the points of `B` (mT) inside the window.

### `double_integral`

`double_integral(B, spectrum, baseline_points: int = 20, window: IntensityWindow | tuple[float, float] | None = None, freq_GHz: float | None = None, *, baseline: str = 'ends') -> float`

Double integral of a derivative spectrum (proportional to spin number).

- `baseline`: `"ends"` (default): a line through `baseline_points` at each end (of the `window`, if given) is subtracted before each integration. `"outside"` (needs a `window` with at least 5 points beyond it on each side): the derivative baseline is a line fitted to the sweep within one window width on each side of the window, which is far less noisy than a few end points and does not cut into the tails of the line (an error in it grows quadratically in the double integral); the absorption baseline is still a line through `baseline_points` at each end of the window.

### `field_for_g`

`field_for_g(g: float, freq_GHz: float) -> float`

Resonance field (mT) for a given g at microwave frequency `freq_GHz`.

### `g_for_field`

`g_for_field(B_mT, freq_GHz: float)`

g-value at field `B_mT` (mT).

### `outside_baseline`

`outside_baseline(B, y, m, reach: float = 1.0) -> np.ndarray | None`

`y` minus a line fitted to the points outside the window mask `m`
but within `reach` window widths of it, on both sides; None if either
side has fewer than 5 such points.

### `peak_to_peak`

`peak_to_peak(B, spectrum, window: IntensityWindow | tuple[float, float] | None = None, freq_GHz: float | None = None) -> float`

Max minus min of the derivative spectrum, optionally inside `window`
(an `IntensityWindow`, which needs `freq_GHz` unless centred on
a field, or `(low, high)` in mT).

### `t1_b2_amplitude`

`t1_b2_amplitude(B, spectrum, freq_GHz: float, g_t1: float = 2.0018, g_b2: float = 1.9973, search_mT: float = 0.1) -> float`

T1–B2 peak-to-peak amplitude of the enamel signal.

The maximum is searched within `search_mT` of the field of `g_t1`
and the minimum within `search_mT` of the field of `g_b2`.
