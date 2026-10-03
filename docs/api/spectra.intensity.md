# EPR spectra: scalar intensities

`eprdating.spectra.intensity`

Scalar ESR intensities from cw (first-derivative) spectra.

- `Conventions`: magnetic field in mT, microwave frequency in GHz, spectra as
first-derivative cw-EPR traces sampled on increasing field.

**Contents:** [`double_integral`](#double_integral), [`field_for_g`](#field_for_g), [`g_for_field`](#g_for_field), [`peak_to_peak`](#peak_to_peak), [`t1_b2_amplitude`](#t1_b2_amplitude)

### `double_integral`

`double_integral(B, spectrum, baseline_points: int = 20) -> float`

Double integral of a derivative spectrum (proportional to spin number).

A linear baseline estimated from `baseline_points` at each end is
subtracted before each integration.

### `field_for_g`

`field_for_g(g: float, freq_GHz: float) -> float`

Resonance field (mT) for a given g at microwave frequency `freq_GHz`.

### `g_for_field`

`g_for_field(B_mT, freq_GHz: float)`

g-value at field `B_mT` (mT).

### `peak_to_peak`

`peak_to_peak(B, spectrum, window: tuple[float, float] | None = None) -> float`

Max minus min of the derivative spectrum, optionally inside `window` (mT).

### `t1_b2_amplitude`

`t1_b2_amplitude(B, spectrum, freq_GHz: float, g_t1: float = 2.0018, g_b2: float = 1.9973, search_mT: float = 0.1) -> float`

T1–B2 peak-to-peak amplitude of the enamel signal.

The maximum is searched within `search_mT` of the field of `g_t1`
and the minimum within `search_mT` of the field of `g_b2`.
