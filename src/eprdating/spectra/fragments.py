"""Enamel fragments measured at several angles (optional protocol).

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

* :func:`parse_angular_name` reads the orientation and angle from file names
  such as ``S_900s_X_gon_40dg_result.csv`` (the pattern can be changed);
* :func:`merge_angular` averages the spectra of one step, every orientation
  with the same weight (an orientation with fewer angles would otherwise
  weigh less), on a common field grid, after an optional alignment;
* :func:`angular_profile` gives the intensity at every angle, to see the
  anisotropy and spot a bad angle;
* :func:`fragment_intensity` measures the merged spectrum with any method of
  :func:`~eprdating.spectra.intensity` and adds the scatter between
  orientations to the error;
* :class:`XrayCalibration` converts exposure times to doses; a fit in
  seconds and the conversion of De keep the calibration error separate.

Not done: the separation of oriented and non-oriented (isotropic) CO2-
radicals ("isotropic correction") used by some laboratories before reading
T1-B2 on the merged spectrum. The merged spectrum here contains both.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace

import numpy as np

from .intensity import DEFAULT_WINDOW, IntensityWindow
from .io import Spectrum
from .measure import Intensity, intensity

#: default file-name pattern: an orientation letter (x, y or z) and the angle
#: in degrees, e.g. ``_X_gon_40dg``, ``_Z_gon_180dg``, ``Frag_Y_3600s_gon_20dg``
ANGLE_PATTERN = re.compile(
    r"(?:^|_)(?P<orientation>[XYZxyz])_(?:.*_)?gon_(?P<angle>\d+(?:\.\d+)?)dg", re.IGNORECASE
)


def parse_angular_name(name: str, pattern: re.Pattern | str = ANGLE_PATTERN) -> tuple[str, float]:
    """``(orientation, angle in degrees)`` from a file name, orientation in
    upper case. Raises ValueError when the name does not match."""
    pat = re.compile(pattern) if isinstance(pattern, str) else pattern
    m = pat.search(name)
    if not m:
        raise ValueError(f"no orientation/angle in {name!r}")
    return m.group("orientation").upper(), float(m.group("angle"))


@dataclass
class AngularSpectrum:
    """One spectrum of a fragment: orientation (``"X"``, ``"Y"``, ``"Z"`` or any
    label) and goniometer angle in degrees."""

    spectrum: Spectrum
    orientation: str
    angle_deg: float


def angular_set(spectra: Sequence[Spectrum], pattern: re.Pattern | str = ANGLE_PATTERN) -> list[AngularSpectrum]:
    """Angular spectra from spectra whose names hold orientation and angle."""
    return [AngularSpectrum(s, *parse_angular_name(s.name, pattern)) for s in spectra]


def _grid(items: Sequence[AngularSpectrum]) -> np.ndarray:
    """Common field grid: the overlap of all sweeps, with the median step."""
    lo = max(float(a.spectrum.B.min()) for a in items)
    hi = min(float(a.spectrum.B.max()) for a in items)
    if hi <= lo:
        raise ValueError("the angular spectra do not share a field range")
    step = float(np.median([np.median(np.diff(a.spectrum.B)) for a in items]))
    return np.arange(lo, hi + step / 2, step)


def _on(B: np.ndarray, s: Spectrum, freq: float | None) -> np.ndarray:
    Bs = np.asarray(s.B, float)
    if freq is not None and s.freq_GHz is not None and s.freq_GHz != freq:
        Bs = Bs * freq / s.freq_GHz  # same g at the same position
    return np.interp(B, Bs, s.y)


def merge_angular(items: Sequence[AngularSpectrum], *, balance: bool = True, name: str = "merged") -> Spectrum:
    """Merged spectrum of one irradiation step.

    balance : average the angles of each orientation first and then the
              orientations (default), so that each orientation weighs the
              same; ``False`` averages all spectra alike.

    The result is a :class:`~eprdating.spectra.io.Spectrum` on the common
    field grid whose scans are the orientation averages (or the single
    spectra without ``balance``), so that the usual intensity methods and
    their noise-injection errors apply. Spectra with another microwave
    frequency are put on the g scale of the first.
    """
    items = list(items)
    if not items:
        raise ValueError("no spectra to merge")
    B = _grid(items)
    freq = next((a.spectrum.freq_GHz for a in items if a.spectrum.freq_GHz), None)
    if balance:
        groups: dict[str, list[np.ndarray]] = {}
        for a in items:
            groups.setdefault(a.orientation, []).append(_on(B, a.spectrum, freq))
        scans = np.array([np.mean(v, axis=0) for _, v in sorted(groups.items())])
    else:
        scans = np.array([_on(B, a.spectrum, freq) for a in items])
    first = items[0].spectrum
    params = {"orientations": sorted({a.orientation for a in items}), "n_spectra": len(items),
              "balanced": balance}
    return replace(first, B=B, scans=scans, name=name, params=params, freq_GHz=freq, dropped_points=0)


def angular_profile(items: Sequence[AngularSpectrum], method: str = "peak_to_peak",
                    window: IntensityWindow = DEFAULT_WINDOW, **kw) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Intensity at every angle: ``{orientation: (angles, intensities)}``,
    sorted by angle (no noise injection, for speed)."""
    kw.setdefault("n_noise", 2)
    kw.setdefault("n_null", 2)
    out: dict[str, list] = {}
    for a in items:
        out.setdefault(a.orientation, []).append((a.angle_deg, intensity(a.spectrum, method, window=window, **kw).value))
    return {o: (np.array([x for x, _ in sorted(v)]), np.array([y for _, y in sorted(v)]))
            for o, v in sorted(out.items())}


def fragment_intensity(items: Sequence[AngularSpectrum], method: str = "peak_to_peak",
                       template: np.ndarray | tuple | Callable | None = None,
                       window: IntensityWindow = DEFAULT_WINDOW, *, balance: bool = True,
                       orientation_scatter: bool = False, **kw) -> Intensity:
    """Intensity of the merged spectrum of one step.

    The error is that of :func:`~eprdating.spectra.intensity` on the merged
    spectrum (noise). With ``orientation_scatter`` and two or more
    orientations, the standard error of the intensities of the orientation
    averages is added in quadrature. That scatter reflects the anisotropy of
    the fragment and is largely the same at every step (same fragment, same
    orientations), so it is off by default: it would inflate every point
    alike rather than describe the scatter between steps.
    """
    merged = merge_angular(items, balance=balance)
    r = intensity(merged, method, template, window, **kw)
    if orientation_scatter and merged.n_scans > 1 and balance:
        per = []
        for k in range(merged.n_scans):
            one = replace(merged, scans=merged.scans[k:k + 1])
            per.append(intensity(one, method, template, window, **{**kw, "n_noise": 2, "n_null": 2}).value)
        sem = float(np.std(per, ddof=1) / np.sqrt(len(per)))
        r = replace(r, sigma=float(np.hypot(r.sigma, sem)) if np.isfinite(r.sigma) else sem,
                    n_repeats=len(per))
    return r


@dataclass(frozen=True)
class XrayCalibration:
    """Dose rate of an X-ray irradiator, ``rate`` ± ``sigma`` in Gy/s (from a
    known gamma dose given to a reference sample)."""

    rate: float
    sigma: float = 0.0

    def __post_init__(self) -> None:
        if not self.rate > 0 or self.sigma < 0:
            raise ValueError("the dose rate must be positive and its error non-negative")

    def dose(self, seconds):
        """Dose in Gy for exposure times in seconds."""
        return np.asarray(seconds, float) * self.rate

    def De(self, De_s: float, De_s_sigma: float = 0.0) -> tuple[float, float]:
        """Equivalent dose in Gy from one fitted in seconds, with the
        calibration error added in quadrature."""
        De = De_s * self.rate
        rel = np.hypot(De_s_sigma / De_s if De_s else 0.0, self.sigma / self.rate)
        return float(De), float(abs(De) * rel)


__all__ = [
    "ANGLE_PATTERN",
    "AngularSpectrum",
    "XrayCalibration",
    "angular_profile",
    "angular_set",
    "fragment_intensity",
    "merge_angular",
    "parse_angular_name",
]
