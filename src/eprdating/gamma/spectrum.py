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
    "fit_single_peak",
    "read_spectrum_txt",
]
