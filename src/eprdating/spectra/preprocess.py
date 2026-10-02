"""Pre-processing of cw-EPR spectra before intensity estimation.

* :func:`subtract_baseline`  polynomial baseline fitted outside the signal.
* :func:`aligned_average`    field-aligned average (empirical template).
* :func:`normalise`          receiver gain, microwave power, mass.
* :func:`time_constant_filter`  lock-in time constant (RC) along the sweep.
* :func:`pseudo_modulation`  field-modulation broadening of a simulated
  derivative spectrum (Hyde, Pasenkiewicz-Gierula, Jesmanowicz & Antholine
  1990, *Appl. Magn. Reson.* 1, 483), so simulated shapes can be compared
  with spectra recorded with a large modulation amplitude.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import lfilter
from scipy.special import j1


def subtract_baseline(B, y, exclude: tuple[float, float] | None = None, order: int = 3) -> np.ndarray:
    """Subtract a polynomial fitted to the points outside ``exclude`` (mT).

    ``y`` may be 1-D or 2-D (one row per scan).
    """
    B = np.asarray(B, float)
    y = np.asarray(y, float)
    m = np.ones_like(B, bool) if exclude is None else (B < exclude[0]) | (B > exclude[1])
    x = (B - B.mean()) / (np.ptp(B) or 1.0)
    V = np.vander(x, order + 1)
    coef, *_ = np.linalg.lstsq(V[m], np.atleast_2d(y)[:, m].T, rcond=None)
    out = np.atleast_2d(y) - (V @ coef).T
    return out[0] if y.ndim == 1 else out


def normalise(y, power_mW: float | None = None, gain: float | None = None, mass_mg: float | None = None,
              ref_power_mW: float = 1.0) -> np.ndarray:
    """Normalise a signal to unit gain, ``ref_power_mW`` and unit mass.

    The signal of a non-saturated line grows as the square root of the
    microwave power; the CO2- signal of enamel is not saturated below a few
    tens of mW at room temperature. Intensities are divided by ``gain`` and
    ``mass_mg`` when given.
    """
    y = np.asarray(y, float)
    f = 1.0
    if power_mW is not None:
        f *= np.sqrt(ref_power_mW / power_mW)
    if gain is not None:
        f /= gain
    if mass_mg is not None:
        f /= mass_mg
    return y * f


def pseudo_modulation(B, y, mod_pp_mT: float) -> np.ndarray:
    """First-harmonic spectrum recorded with modulation amplitude ``mod_pp_mT``.

    ``y`` is the first-derivative spectrum for vanishing modulation on the
    uniform grid ``B`` (mT). In Fourier space the modulated first harmonic is
    the derivative multiplied by ``2 J1(k a)/(k a)``, with ``a`` half the
    peak-to-peak modulation amplitude (Hyde et al. 1990). The derivative
    normalisation is kept, so a small modulation returns ``y`` unchanged;
    the signal an instrument records is this times the modulation amplitude.
    """
    B = np.asarray(B, float)
    y = np.asarray(y, float)
    if mod_pp_mT <= 0:
        return y.copy()
    dB = np.diff(B)
    if not np.allclose(dB, dB[0], rtol=1e-6):
        raise ValueError("pseudo_modulation needs a uniform field grid")
    n = y.size
    pad = int(np.ceil(mod_pp_mT / dB[0])) + 8
    yp = np.pad(y, pad, mode="constant")
    k = 2 * np.pi * np.fft.rfftfreq(yp.size, d=dB[0])
    x = k * mod_pp_mT / 2
    h = np.ones_like(x)
    nz = x > 0
    h[nz] = 2 * j1(x[nz]) / x[nz]
    return np.fft.irfft(np.fft.rfft(yp) * h, n=yp.size)[pad:pad + n]


def noise_sigma(B, y, signal: tuple[float, float]) -> float:
    """Standard deviation of a baseline-corrected spectrum outside ``signal``."""
    B = np.asarray(B, float)
    m = (B < signal[0]) | (B > signal[1])
    return float(np.std(np.asarray(y, float)[m], ddof=1))


def time_constant_filter(y, tau_points: float) -> np.ndarray:
    """Single-pole (RC) low-pass filter along the sweep, as applied by the
    lock-in time constant; ``tau_points`` is the time constant divided by
    the time per point. The filter delays and broadens the line."""
    y = np.asarray(y, float)
    if tau_points <= 0:
        return y.copy()
    a = np.exp(-1.0 / tau_points)
    return lfilter([1 - a], [1, -a], y, axis=-1)


def aligned_average(B, spectra, window: tuple[float, float], max_shift: float = 0.6,
                    reference: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Average of several spectra after aligning them in field.

    Each spectrum is shifted (by whole field steps, up to ``max_shift`` mT)
    to maximise its cross-correlation with the ``reference`` spectrum (the
    one with the largest peak-to-peak inside ``window`` by default). Useful
    to build an empirical line-shape template from the strongest spectra of
    a dose series. Returns ``(average, shifts_mT)``.
    """
    B = np.asarray(B, float)
    Y = np.atleast_2d(np.asarray(spectra, float))
    m = (B >= window[0]) & (B <= window[1])
    step = float(np.median(np.diff(B)))
    if reference is None:
        reference = int(np.argmax(np.ptp(Y[:, m], axis=1)))
    ref = Y[reference]
    nlag = round(max_shift / step)
    lags = np.arange(-nlag, nlag + 1)
    out, shifts = [], []
    for y in Y:
        c = [np.dot(np.roll(y, k)[m], ref[m]) for k in lags]
        k = int(lags[int(np.argmax(c))])
        out.append(np.roll(y, k))
        shifts.append(-k * step)
    return np.mean(out, axis=0), np.array(shifts)


__all__ = [
    "aligned_average",
    "noise_sigma",
    "normalise",
    "pseudo_modulation",
    "subtract_baseline",
    "time_constant_filter",
]
