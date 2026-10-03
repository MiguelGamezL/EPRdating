"""Gamma-ray spectra: reading, energy and resolution calibration.

Supported format
----------------
ASCII spectra with a ``#``-commented header of ``key: value`` lines followed
by ``# Channel data`` and four tab-separated columns (channel, energy, counts,
rate), as in the NORM spectra (HPGe 40 %, NIM electronics) used to
validate this module; the exporting MCA program has not been identified::

    # Start time:    2024-08-01, 12:28:55
    # Real time (s): 86565.170
    # Live time (s): 86400.000
    # Energy calibration coefficients ( E = sum(Ai * n**i) )
    #     A0: 0.000000
    #     A1: 0.250000
    ...
    # Channel data
    # n	energy(keV)	counts	rate(1/s)
    1	0.250	0	0

The energy column is the MCA's own (often nominal) calibration; use
:meth:`GammaSpectrum.calibrate` with known lines to obtain the real one.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

#: lines of common calibration sources (keV)
CALIBRATION_LINES = {
    "57Co": (122.061, 136.474),
    "22Na": (1274.537,),
    "137Cs": (661.657,),
    "88Y": (898.042, 1836.063),
    "60Co": (1173.228, 1332.492),
    "133Ba": (80.998, 276.399, 302.851, 356.013, 383.849),
    "152Eu": (121.782, 244.697, 344.279, 778.904, 964.057, 1112.076, 1408.013),
}

#: strong natural lines that can recalibrate any environmental spectrum
NATURAL_LINES = (238.632, 351.932, 583.187, 609.312, 911.204, 1460.820, 1764.494, 2614.511)


@dataclass
class Calibration:
    """Energy ``E = Σ a_i n^i`` (n = channel number) and resolution
    ``σ(E)² = w0 + w1 E`` (keV)."""

    coef: tuple[float, ...]
    resolution: tuple[float, float] = (0.25, 0.0005)
    residuals_keV: np.ndarray | None = None
    chi2_red: float | None = None  # of the energy fit, from the centroid uncertainties

    def energy(self, channel):
        return np.polynomial.polynomial.polyval(np.asarray(channel, float), self.coef)

    def channel(self, energy_keV):
        e = np.atleast_1d(np.asarray(energy_keV, float))
        c = self.coef
        n = (e - c[0]) / c[1]
        for _ in range(20):  # Newton iterations for non-linear terms
            f = np.polynomial.polynomial.polyval(n, c) - e
            df = np.polynomial.polynomial.polyval(n, np.polynomial.polynomial.polyder(c))
            n = n - f / df
        return n if np.ndim(energy_keV) else float(n[0])

    def gain(self, channel):
        """keV per channel."""
        return np.polynomial.polynomial.polyval(np.asarray(channel, float),
                                                np.polynomial.polynomial.polyder(self.coef))

    def sigma_keV(self, energy_keV):
        w0, w1 = self.resolution
        return np.sqrt(np.maximum(w0 + w1 * np.asarray(energy_keV, float), 1e-6))


@dataclass
class GammaSpectrum:
    counts: np.ndarray
    live_time: float  # s
    real_time: float | None = None
    start: str | None = None
    name: str = ""
    channels: np.ndarray | None = None  # channel numbers (default 1..N, as in the files)
    calibration: Calibration | None = None
    header: dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        self.counts = np.asarray(self.counts, float)
        if self.channels is None:
            self.channels = np.arange(1, self.counts.size + 1, dtype=float)

    @property
    def dead_time_fraction(self) -> float | None:
        return None if not self.real_time else 1.0 - self.live_time / self.real_time

    @property
    def energy(self) -> np.ndarray:
        if self.calibration is None:
            raise ValueError(f"{self.name}: spectrum not calibrated")
        return self.calibration.energy(self.channels)

    def calibrate(self, lines_keV: Sequence[float], guess: Calibration | None = None,
                  order: int = 1, window_keV: float = 6.0, min_significance: float = 8.0,
                  max_residual_keV: float = 0.4) -> Calibration:
        """Fit the energy and resolution calibration from known lines.

        Each line is located with ``guess`` (default: the current calibration),
        fitted with a Gaussian on a linear background, and kept if its area
        exceeds ``min_significance`` standard deviations. Lines that miss the
        fit by more than ``max_residual_keV`` (and 5σ) are dropped one by one,
        e.g. 40K at 1460.8 keV in a thorium-rich spectrum, where it merges
        with 228Ac 1459.1 keV. Returns the new calibration and stores it on
        the spectrum.
        """
        guess = guess or self.calibration
        if guess is None:
            raise ValueError("a first-guess calibration is needed")
        found = []
        for e in lines_keV:
            p = fit_single_peak(self, guess, e, window_keV)
            if p is not None and p["area"] > min_significance * p["area_sigma"]:
                found.append((e, p["centroid"], p["centroid_sigma"], p["sigma_ch"]))
        if len(found) < order + 1:
            raise RuntimeError(f"{self.name}: only {len(found)} calibration lines found")
        E, ch, sch, sig = (np.array(x) for x in zip(*found))
        keep = np.ones(E.size, bool)
        while True:  # drop lines that do not fit (unresolved interferences)
            w = 1.0 / np.maximum(sch[keep], 1e-3)
            coef = np.polynomial.polynomial.polyfit(ch[keep], E[keep], order, w=w)
            res = np.abs(E - np.polynomial.polynomial.polyval(ch, coef))
            res[~keep] = 0.0
            j = int(np.argmax(res))
            if keep.sum() <= order + 2 or res[j] <= max(max_residual_keV, 5 * sch[j] * coef[1]):
                break
            keep[j] = False
        E, ch, sch, sig = E[keep], ch[keep], sch[keep], sig[keep]
        cal = Calibration(tuple(coef))
        sig_keV = sig * cal.gain(ch)
        if len(found) >= 2:
            A = np.c_[np.ones_like(E), E]
            res, *_ = np.linalg.lstsq(A, sig_keV**2, rcond=None)
            cal.resolution = (float(res[0]), float(res[1]))
        else:
            cal.resolution = (float(sig_keV[0] ** 2), 0.0)
        cal.residuals_keV = E - cal.energy(ch)
        dof = E.size - (order + 1)
        if dof > 0:
            sE = np.maximum(sch * cal.gain(ch), 0.02)  # keV; floor for very strong peaks
            cal.chi2_red = float(np.sum((cal.residuals_keV / sE) ** 2) / dof)
        self.calibration = cal
        return cal


def _gauss_lin(x, A, mu, s, b0, b1):
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2) + b0 + b1 * (x - mu)


def fit_single_peak(spec: GammaSpectrum, cal: Calibration, energy_keV: float, window_keV: float = 6.0):
    """Gaussian + linear background around one line (channels). Returns a
    dict with centroid, sigma (channels) and area (counts), or None."""
    c0 = cal.channel(energy_keV)
    half = window_keV / cal.gain(c0)
    m = np.abs(spec.channels - c0) <= half
    x, y = spec.channels[m], spec.counts[m]
    if y.size < 8:
        return None
    s0 = cal.sigma_keV(energy_keV) / cal.gain(c0)
    j = int(np.argmax(y))
    bg = float(np.median(np.r_[y[:3], y[-3:]]))
    p0 = [max(y[j] - bg, 1.0), x[j], s0, bg, 0.0]
    lo = [0.0, c0 - half / 2, 0.3 * s0, -np.inf, -np.inf]
    hi = [np.inf, c0 + half / 2, 3.0 * s0, np.inf, np.inf]
    try:
        p, C = curve_fit(_gauss_lin, x, y, p0=p0, sigma=np.sqrt(np.maximum(y, 1.0)), absolute_sigma=True,
                         bounds=(lo, hi), maxfev=20000)
    except (RuntimeError, ValueError):
        return None
    A, mu, s = p[:3]
    area = A * s * np.sqrt(2 * np.pi)
    J = np.array([s, 0, A]) * np.sqrt(2 * np.pi)
    area_sigma = float(np.sqrt(J @ C[:3, :3] @ J))
    return {"centroid": float(mu), "centroid_sigma": float(np.sqrt(C[1, 1])), "sigma_ch": float(s),
            "area": float(area), "area_sigma": area_sigma}


def find_peaks_channels(spec: GammaSpectrum, min_significance: float = 5.0, max_peaks: int = 40):
    """Candidate peaks: ``(channel, significance)`` sorted by significance.

    The counts are smoothed with a 1.5-channel Gaussian, a baseline is taken
    as a running 20th percentile, and local maxima exceeding the baseline by
    ``min_significance`` Poisson standard deviations are kept.
    """
    from scipy.ndimage import gaussian_filter1d, percentile_filter
    from scipy.signal import find_peaks

    y = spec.counts
    sm = gaussian_filter1d(y, 1.5)
    width = max(15, y.size // 200)
    base = percentile_filter(sm, 20, size=4 * width + 1)
    sig = (sm - base) / np.sqrt(np.maximum(base, 1.0) / 3.0)  # smoothing ≈ averages ~3 channels
    idx, _ = find_peaks(sig, height=min_significance, distance=3)
    order = np.argsort(sig[idx])[::-1][:max_peaks]
    idx = idx[order]
    # parabolic refinement of the maximum
    c = []
    for i in idx:
        if 0 < i < y.size - 1:
            a, b, d = sm[i - 1], sm[i], sm[i + 1]
            den = a - 2 * b + d
            c.append(i + (0.5 * (a - d) / den if den != 0 else 0.0))
        else:
            c.append(float(i))
    return np.interp(c, np.arange(y.size), spec.channels), sig[idx]


def auto_calibrate(spec: GammaSpectrum, lines_keV: Sequence[float] = NATURAL_LINES,
                   gain_range: tuple[float, float] = (0.01, 10.0), max_offset_keV: float = 30.0,
                   tolerance_keV: float = 1.5, refine: bool = True) -> Calibration:
    """Energy calibration without a first guess.

    Peaks are searched (:func:`find_peaks_channels`) and every pair of
    strong peaks is tried as every pair of known lines; the linear map that
    places the most lines on peaks (within ``tolerance_keV`` plus 0.2 %) wins.
    It is then refined with :meth:`GammaSpectrum.calibrate`. Works for
    HPGe spectra with at least three of the given lines (for environmental
    samples the default natural lines: 238.6, 351.9, 583.2, 609.3, 911.2,
    1460.8, 1764.5, 2614.5 keV).
    """
    ch, sig = find_peaks_channels(spec)
    if ch.size < 2:
        raise RuntimeError(f"{spec.name}: no peaks found")
    E = np.sort(np.asarray(lines_keV, float))
    strong = np.argsort(sig)[::-1][:15]
    best = None
    for i in strong:
        for j in strong:
            if ch[j] <= ch[i]:
                continue
            for k in range(E.size):
                for m in range(k + 1, E.size):
                    a1 = (E[m] - E[k]) / (ch[j] - ch[i])
                    if not gain_range[0] <= a1 <= gain_range[1]:
                        continue
                    a0 = E[k] - a1 * ch[i]
                    if abs(a0) > max_offset_keV:
                        continue
                    pred = (E - a0) / a1
                    d = np.abs(pred[:, None] - ch[None, :]) * a1
                    tol = tolerance_keV + 0.002 * E
                    near = d.argmin(axis=1)
                    hit = d.min(axis=1) <= tol
                    nhit = int(hit.sum())
                    if nhit < 2 or (best is not None and nhit < best[0][0]):
                        continue
                    # refit the matched peaks; among equal hit counts the smallest rms wins
                    cc, ee = ch[near[hit]], E[hit]
                    b1, b0 = np.polyfit(cc, ee, 1)
                    rms = float(np.sqrt(np.mean((ee - (b0 + b1 * cc)) ** 2)))
                    key = (nhit, -rms)
                    if best is None or key > best[0]:
                        best = (key, b0, b1)
    if best is None or best[0][0] < 3:
        raise RuntimeError(f"{spec.name}: could not match at least three known lines")
    _, a0, a1 = best
    # resolution from the half-maximum width of the strongest matched peak
    pred = (E - a0) / a1
    near = np.abs(pred[:, None] - ch[None, :]).argmin(axis=1)
    jj = max(near, key=lambda q: sig[q])
    e_ref = a0 + a1 * ch[jj]
    sigma_keV = _half_max_sigma(spec, ch[jj]) * a1
    first = Calibration((a0, a1), resolution=(0.0, sigma_keV**2 / e_ref))
    if not refine:
        spec.calibration = first
        return first
    kw = {"min_significance": 5.0, "window_keV": max(6.0, 6 * sigma_keV),
          "max_residual_keV": max(0.4, 0.3 * sigma_keV)}
    cal = spec.calibrate(lines_keV, guess=first, **kw)
    if cal.residuals_keV.size >= 4:  # non-linear ADCs: keep a quadratic term if it is needed
        try:
            quad = spec.calibrate(lines_keV, guess=cal, order=2, **kw)
            quad = spec.calibrate(lines_keV, guess=quad, order=2, **kw)  # lines found with the better guess
        except RuntimeError:
            quad = None
        if (quad is not None and quad.chi2_red is not None and cal.chi2_red is not None
                and quad.residuals_keV.size >= cal.residuals_keV.size and quad.chi2_red < 0.5 * cal.chi2_red):
            cal = quad
        spec.calibration = cal
    fwhm = 2.3548 * float(cal.sigma_keV(662.0))
    if abs(cal.coef[0]) > 2 * max_offset_keV or not gain_range[0] <= cal.coef[1] <= gain_range[1] or fwhm > HPGE_MAX_FWHM:
        spec.calibration = None
        raise RuntimeError(
            f"{spec.name}: no consistent calibration (offset {cal.coef[0]:.1f} keV, FWHM at 662 keV "
            f"{fwhm:.1f} keV). The line-by-line method needs HPGe resolution; scintillator spectra "
            f"(NaI, CsI, LaBr3) are not supported.")
    return cal


#: FWHM at 662 keV above which a spectrum is not treated as HPGe (keV)
HPGE_MAX_FWHM = 6.0


def _half_max_sigma(spec: GammaSpectrum, channel: float) -> float:
    """Gaussian sigma (channels) from the full width at half maximum."""
    from scipy.ndimage import gaussian_filter1d

    i = int(np.argmin(np.abs(spec.channels - channel)))
    y = gaussian_filter1d(spec.counts, 1.0)
    lo, hi = max(i - 400, 0), min(i + 400, y.size)
    base = np.percentile(y[lo:hi], 20)
    half = base + 0.5 * (y[i] - base)
    left = i
    while left > lo and y[left] > half:
        left -= 1
    right = i
    while right < hi - 1 and y[right] > half:
        right += 1
    return max((right - left) / 2.3548, 0.8)


_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def read_spectrum_txt(path: str | Path) -> GammaSpectrum:
    """Read an ASCII gamma spectrum with a ``# key: value`` header."""
    path = Path(path)
    header: dict[str, str] = {}
    coef = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        if not line.startswith("#"):
            break
        body = line[1:].strip()
        if ":" in body:
            k, v = (s.strip() for s in body.split(":", 1))
            header[k] = v
            m = re.fullmatch(r"A(\d)", k)
            if m:
                coef[int(m.group(1))] = float(v)
    data = np.loadtxt(path, comments="#", encoding="latin-1", ndmin=2)

    def num(key):
        v = header.get(key)
        m = _NUM.search(v) if v else None
        return float(m.group()) if m else None

    live = num("Live time (s)")
    if live is None:
        raise ValueError(f"{path}: no live time in header")
    cal = Calibration(tuple(coef[i] for i in sorted(coef))) if coef else None
    return GammaSpectrum(
        counts=data[:, 2], live_time=live, real_time=num("Real time (s)"), start=header.get("Start time"),
        name=path.stem, channels=data[:, 0], calibration=cal, header=header,
    )


__all__ = [
    "CALIBRATION_LINES",
    "NATURAL_LINES",
    "Calibration",
    "GammaSpectrum",
    "auto_calibrate",
    "find_peaks_channels",
    "fit_single_peak",
    "read_spectrum_txt",
]
