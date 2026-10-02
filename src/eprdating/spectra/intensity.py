"""Scalar ESR intensities from cw (first-derivative) spectra.

Conventions: magnetic field in mT, microwave frequency in GHz, spectra as
first-derivative cw-EPR traces sampled on increasing field.
"""

from __future__ import annotations

import numpy as np
from scipy.constants import h, physical_constants

MU_B = physical_constants["Bohr magneton"][0]  # J/T

#: g-values of the peaks conventionally used for enamel dating
#: (T1–B2 peak-to-peak amplitude of the CO2- signal).
G_T1 = 2.0018
G_B2 = 1.9973


def field_for_g(g: float, freq_GHz: float) -> float:
    """Resonance field (mT) for a given g at microwave frequency ``freq_GHz``."""
    return h * freq_GHz * 1e9 / (g * MU_B) * 1e3


def g_for_field(B_mT, freq_GHz: float):
    """g-value at field ``B_mT`` (mT)."""
    return h * freq_GHz * 1e9 / (MU_B * np.asarray(B_mT, float) * 1e-3)


def _window(B: np.ndarray, center: float, half_width: float) -> np.ndarray:
    m = np.abs(B - center) <= half_width
    if not m.any():
        raise ValueError(f"no data within {half_width} mT of {center:.3f} mT")
    return m


def peak_to_peak(B, spectrum, window: tuple[float, float] | None = None) -> float:
    """Max minus min of the derivative spectrum, optionally inside ``window`` (mT)."""
    B = np.asarray(B, float)
    y = np.asarray(spectrum, float)
    if window is not None:
        m = (B >= window[0]) & (B <= window[1])
        y = y[m]
    return float(y.max() - y.min())


def t1_b2_amplitude(
    B,
    spectrum,
    freq_GHz: float,
    g_t1: float = G_T1,
    g_b2: float = G_B2,
    search_mT: float = 0.1,
) -> float:
    """T1–B2 peak-to-peak amplitude of the enamel signal.

    The maximum is searched within ``search_mT`` of the field of ``g_t1``
    and the minimum within ``search_mT`` of the field of ``g_b2``.
    """
    B = np.asarray(B, float)
    y = np.asarray(spectrum, float)
    t1 = y[_window(B, field_for_g(g_t1, freq_GHz), search_mT)].max()
    b2 = y[_window(B, field_for_g(g_b2, freq_GHz), search_mT)].min()
    return float(t1 - b2)


def double_integral(B, spectrum, baseline_points: int = 20) -> float:
    """Double integral of a derivative spectrum (proportional to spin number).

    A linear baseline estimated from ``baseline_points`` at each end is
    subtracted before each integration.
    """
    B = np.asarray(B, float)
    y = np.asarray(spectrum, float)

    def debase(x, v):
        idx = np.r_[0:baseline_points, len(v) - baseline_points:len(v)]
        c = np.polyfit(x[idx], v[idx], 1)
        return v - np.polyval(c, x)

    y = debase(B, y)
    absorption = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(B))])
    absorption = debase(B, absorption)
    return float(np.trapezoid(absorption, B))
