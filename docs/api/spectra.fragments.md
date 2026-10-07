# EPR spectra: enamel fragments (optional)

`eprdating.spectra.fragments`

Enamel fragments measured at several angles (optional protocol).

A fragment is not a powder: the CO2- radicals of oriented crystallites give
spectra that change with the orientation of the fragment in the field. The
non-destructive protocol (Grün et al. 2008; Joannes-Boyau et al. 2010;
Joannes-Boyau 2013) records, for every irradiation step, the spectrum at a
series of goniometer angles (e.g. 0-180° every 20°) in up to three
orientations of the fragment (x, y, z), and measures the intensity on their
average, the *merged spectrum*, which approximates a powder. Fragments are
usually irradiated with X-rays, so the doses are times that a calibration
(Gy/s, measured against a known gamma dose) turns into Gy.

What this module does:

* `parse_angular_name` reads the orientation and angle from file names
  such as `S_900s_X_gon_40dg_result.csv` (the pattern can be changed);
* `merge_angular` averages the spectra of one step, every orientation
  with the same weight (an orientation with fewer angles would otherwise
  weigh less), on a common field grid, after an optional alignment;
* `angular_profile` gives the intensity at every angle, to see the
  anisotropy and spot a bad angle;
* `fragment_intensity` measures the merged spectrum with any method of
  `intensity` and adds the scatter between
  orientations to the error;
* `XrayCalibration` converts exposure times to doses; a fit in
  seconds and the conversion of De keep the calibration error separate.

Not done: the separation of oriented and non-oriented (isotropic) CO2-
radicals ("isotropic correction") used by some laboratories before reading
T1-B2 on the merged spectrum. The merged spectrum here contains both.

**Contents:** [`ANGLE_PATTERN`](#angle_pattern), [`AngularSpectrum`](#angularspectrum), [`XrayCalibration`](#xraycalibration), [`angular_profile`](#angular_profile), [`angular_set`](#angular_set), [`fragment_intensity`](#fragment_intensity), [`merge_angular`](#merge_angular), [`parse_angular_name`](#parse_angular_name)

### `ANGLE_PATTERN`

```python
ANGLE_PATTERN = re.compile('(?:^|_)(?P<orientation>[XYZxyz])_(?:.*_)?gon_(?P<angle>\\d+(?:\\.\\d+)?)dg', re.IGNORECASE)
```

### `AngularSpectrum`

*dataclass* `AngularSpectrum(spectrum: Spectrum, orientation: str, angle_deg: float)`

One spectrum of a fragment: orientation (`"X"`, `"Y"`, `"Z"` or any
label) and goniometer angle in degrees.

### `XrayCalibration`

*dataclass* `XrayCalibration(rate: float, sigma: float = 0.0)`

Dose rate of an X-ray irradiator, `rate` ± `sigma` in Gy/s (from a
known gamma dose given to a reference sample).

**Members**

- `De(self, De_s: float, De_s_sigma: float = 0.0) -> tuple[float, float]` — Equivalent dose in Gy from one fitted in seconds, with the
calibration error added in quadrature.
- `dose(self, seconds)` — Dose in Gy for exposure times in seconds.

### `angular_profile`

`angular_profile(items: Sequence[AngularSpectrum], method: str = 'peak_to_peak', window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), **kw) -> dict[str, tuple[np.ndarray, np.ndarray]]`

Intensity at every angle: `{orientation: (angles, intensities)}`,
sorted by angle (no noise injection, for speed).

### `angular_set`

`angular_set(spectra: Sequence[Spectrum], pattern: re.Pattern | str = re.compile('(?:^|_)(?P<orientation>[XYZxyz])_(?:.*_)?gon_(?P<angle>\\d+(?:\\.\\d+)?)dg', re.IGNORECASE)) -> list[AngularSpectrum]`

Angular spectra from spectra whose names hold orientation and angle.

### `fragment_intensity`

`fragment_intensity(items: Sequence[AngularSpectrum], method: str = 'peak_to_peak', template: np.ndarray | tuple | Callable | None = None, window: IntensityWindow = IntensityWindow(width=100.0, unit='G', center_g=2.0023, center_mT=None), *, balance: bool = True, orientation_scatter: bool = False, **kw) -> Intensity`

Intensity of the merged spectrum of one step.

The error is that of `intensity` on the merged
spectrum (noise). With `orientation_scatter` and two or more
orientations, the standard error of the intensities of the orientation
averages is added in quadrature. That scatter reflects the anisotropy of
the fragment and is largely the same at every step (same fragment, same
orientations), so it is off by default: it would inflate every point
alike rather than describe the scatter between steps.

### `merge_angular`

`merge_angular(items: Sequence[AngularSpectrum], *, balance: bool = True, name: str = 'merged') -> Spectrum`

Merged spectrum of one irradiation step.

- `balance`: average the angles of each orientation first and then the orientations (default), so that each orientation weighs the same; `False` averages all spectra alike.

The result is a `Spectrum` on the common
field grid whose scans are the orientation averages (or the single
spectra without `balance`), so that the usual intensity methods and
their noise-injection errors apply. Spectra with another microwave
frequency are put on the g scale of the first.

### `parse_angular_name`

`parse_angular_name(name: str, pattern: re.Pattern | str = re.compile('(?:^|_)(?P<orientation>[XYZxyz])_(?:.*_)?gon_(?P<angle>\\d+(?:\\.\\d+)?)dg', re.IGNORECASE)) -> tuple[str, float]`

`(orientation, angle in degrees)` from a file name, orientation in
upper case. Raises ValueError when the name does not match.
