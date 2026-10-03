# Gamma: calibration

`eprdating.gamma.spectrum`

Gamma-ray spectra: reading, energy and resolution calibration.

**Supported format**

ASCII spectra with a `#`-commented header of `key: value` lines followed
by `# Channel data` and four tab-separated columns (channel, energy, counts,
rate), as in the NORM spectra (HPGe 40 %, NIM electronics) used to
validate this module; the exporting MCA program has not been identified:

```
# Start time:    2024-08-01, 12:28:55
# Real time (s): 86565.170
# Live time (s): 86400.000
# Energy calibration coefficients ( E = sum(Ai * n**i) )
#     A0: 0.000000
#     A1: 0.250000
...
# Channel data
# n energy(keV)     counts  rate(1/s)
1   0.250   0       0
```

The energy column is the MCA's own (often nominal) calibration; use
`calibrate` with known lines to obtain the real one.

**Contents:** [`CALIBRATION_LINES`](#calibration_lines), [`NATURAL_LINES`](#natural_lines), [`Calibration`](#calibration), [`GammaSpectrum`](#gammaspectrum), [`auto_calibrate`](#auto_calibrate), [`find_peaks_channels`](#find_peaks_channels), [`fit_single_peak`](#fit_single_peak), [`read_spectrum_txt`](#read_spectrum_txt)

### `CALIBRATION_LINES`

```python
CALIBRATION_LINES = {'57Co': (122.061, 136.474), '22Na': (1274.537,), '137Cs': (661.657,), '88Y': (898.042, 1836.063), '60Co': (1173.228, 1332.492), '133Ba': (80.998, 276.399, 302.851, 356.013, 383.849), '152Eu': (121.782, 244.697, 344.279, 778.904, 964.057, 1112.076, 1408.013)}
```

### `NATURAL_LINES`

```python
NATURAL_LINES = (238.632, 351.932, 583.187, 609.312, 911.204, 1460.82, 1764.494, 2614.511)
```

### `Calibration`

*dataclass* `Calibration(coef: tuple[float, ...], resolution: tuple[float, float] = (0.25, 0.0005), residuals_keV: np.ndarray | None = None, chi2_red: float | None = None)`

Energy `E = Σ a_i n^i` (n = channel number) and resolution
`σ(E)² = w0 + w1 E` (keV).

**Members**

- `channel(self, energy_keV)`
- `energy(self, channel)`
- `gain(self, channel)` — keV per channel.
- `sigma_keV(self, energy_keV)`

### `GammaSpectrum`

*dataclass* `GammaSpectrum(counts: np.ndarray, live_time: float, real_time: float | None = None, start: str | None = None, name: str = '', channels: np.ndarray | None = None, calibration: Calibration | None = None, header: dict[str, str] = <factory>)`

GammaSpectrum(counts: 'np.ndarray', live_time: 'float', real_time: 'float | None' = None, start: 'str | None' = None, name: 'str' = '', channels: 'np.ndarray | None' = None, calibration: 'Calibration | None' = None, header: 'dict[str, str]' = <factory>)

**Members**

- `calibrate(self, lines_keV: Sequence[float], guess: Calibration | None = None, order: int = 1, window_keV: float = 6.0, min_significance: float = 8.0, max_residual_keV: float = 0.4) -> Calibration` — Fit the energy and resolution calibration from known lines.
- `dead_time_fraction` *(property)*
- `energy` *(property)*

### `auto_calibrate`

`auto_calibrate(spec: GammaSpectrum, lines_keV: Sequence[float] = (238.632, 351.932, 583.187, 609.312, 911.204, 1460.82, 1764.494, 2614.511), gain_range: tuple[float, float] = (0.01, 10.0), max_offset_keV: float = 30.0, tolerance_keV: float = 1.5, refine: bool = True) -> Calibration`

Energy calibration without a first guess.

Peaks are searched (`find_peaks_channels`) and every pair of
strong peaks is tried as every pair of known lines; the linear map that
places the most lines on peaks (within `tolerance_keV` plus 0.2 %) wins.
It is then refined with `calibrate`. Works for
HPGe spectra with at least three of the given lines (for environmental
samples the default natural lines: 238.6, 351.9, 583.2, 609.3, 911.2,
1460.8, 1764.5, 2614.5 keV).

### `find_peaks_channels`

`find_peaks_channels(spec: GammaSpectrum, min_significance: float = 5.0, max_peaks: int = 40)`

Candidate peaks: `(channel, significance)` sorted by significance.

The counts are smoothed with a 1.5-channel Gaussian, a baseline is taken
as a running 20th percentile, and local maxima exceeding the baseline by
`min_significance` Poisson standard deviations are kept.

### `fit_single_peak`

`fit_single_peak(spec: GammaSpectrum, cal: Calibration, energy_keV: float, window_keV: float = 6.0)`

Gaussian + linear background around one line (channels). Returns a
dict with centroid, sigma (channels) and area (counts), or None.

### `read_spectrum_txt`

`read_spectrum_txt(path: str | Path) -> GammaSpectrum`

Read an ASCII gamma spectrum with a `# key: value` header.
