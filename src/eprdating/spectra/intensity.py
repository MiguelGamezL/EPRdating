"""Scalar ESR intensities from cw (first-derivative) spectra.

The magnetic field is in mT and the microwave frequency in GHz; spectra are
first-derivative cw-EPR traces sampled on increasing field.

The field region used for an intensity is an :class:`IntensityWindow`: by
default 100 G (10 mT) centred on g = 2.0023, i.e. on the field of that g at
each spectrum's own microwave frequency. The window should hold the whole
dating signal with some signal-free margin on both sides for the baseline;
outside it, the spectrum is taken as signal-free and used to estimate the
noise. A wider window may take in other radicals (native signal, CO3-, SO2-,
methyl) that the intensity method does not describe; a narrower one leaves
few points for the baseline.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

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


@dataclass(frozen=True)
class IntensityWindow:
    """Field window in which an EPR intensity is computed.

    width     : full width of the window, in ``unit`` (default 100 G).
    unit      : ``"G"`` or ``"mT"``.
    center_g  : g-value at the centre (default 2.0023); the centre field is
                computed from each spectrum's microwave frequency.
    center_mT : fixed centre field instead (overrides ``center_g``).

    For instance ``IntensityWindow()`` is 100 G around g = 2.0023,
    ``IntensityWindow(60)`` 60 G around it, and
    ``IntensityWindow(8, unit="mT", center_mT=337.1)`` 8 mT around 337.1 mT.
    """

    width: float = 100.0
    unit: str = "G"
    center_g: float = 2.0023
    center_mT: float | None = None

    def __post_init__(self) -> None:
        if self.unit not in ("G", "mT"):
            raise ValueError("unit must be 'G' or 'mT'")
        if not self.width > 0:
            raise ValueError("the window width must be positive")

    @property
    def half_width_mT(self) -> float:
        return 0.5 * self.width * (0.1 if self.unit == "G" else 1.0)

    def center(self, freq_GHz: float | None) -> float:
        """Centre field (mT) for a spectrum recorded at ``freq_GHz``."""
        if self.center_mT is not None:
            return float(self.center_mT)
        if freq_GHz is None:
            raise ValueError("the microwave frequency is needed to centre the window on a g-value; "
                             "give freq_GHz or a window with center_mT")
        return field_for_g(self.center_g, freq_GHz)

    def bounds(self, freq_GHz: float | None) -> tuple[float, float]:
        """``(low, high)`` field limits in mT."""
        c = self.center(freq_GHz)
        return c - self.half_width_mT, c + self.half_width_mT

    def mask(self, B, freq_GHz: float | None) -> np.ndarray:
        """Boolean mask of the points of ``B`` (mT) inside the window.

        Warns if the window reaches beyond the recorded sweep.
        """
        B = np.asarray(B, float)
        lo, hi = self.bounds(freq_GHz)
        if lo < B.min() or hi > B.max():
            warnings.warn(f"intensity window {lo:.2f}-{hi:.2f} mT extends beyond the sweep "
                          f"{B.min():.2f}-{B.max():.2f} mT; it is clipped", stacklevel=2)
        m = (B >= lo) & (B <= hi)
        if m.sum() < 5:
            raise ValueError(f"fewer than 5 points in the intensity window {lo:.2f}-{hi:.2f} mT")
        return m

    def describe(self, freq_GHz: float | None = None) -> str:
        where = f"{self.center_mT:.3f} mT" if self.center_mT is not None else f"g = {self.center_g}"
        text = f"{self.width:g} {self.unit} around {where}"
        if freq_GHz is not None or self.center_mT is not None:
            lo, hi = self.bounds(freq_GHz)
            text += f" ({lo:.2f}-{hi:.2f} mT)"
        return text


#: default window: 100 G around g = 2.0023
DEFAULT_WINDOW = IntensityWindow()


def _resolve(B, window, freq_GHz) -> np.ndarray | None:
    """Mask for ``window``: None, an ``(lo, hi)`` tuple in mT or an
    :class:`IntensityWindow`."""
    if window is None:
        return None
    if isinstance(window, IntensityWindow):
        return window.mask(B, freq_GHz)
    lo, hi = window
    return (np.asarray(B) >= lo) & (np.asarray(B) <= hi)


def _window(B: np.ndarray, center: float, half_width: float) -> np.ndarray:
    m = np.abs(B - center) <= half_width
    if not m.any():
        raise ValueError(f"no data within {half_width} mT of {center:.3f} mT")
    return m


def peak_to_peak(B, spectrum, window: IntensityWindow | tuple[float, float] | None = None,
                 freq_GHz: float | None = None) -> float:
    """Max minus min of the derivative spectrum, optionally inside ``window``
    (an :class:`IntensityWindow`, which needs ``freq_GHz`` unless centred on
    a field, or ``(low, high)`` in mT)."""
    B = np.asarray(B, float)
    y = np.asarray(spectrum, float)
    m = _resolve(B, window, freq_GHz)
    if m is not None:
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


def _integrate_twice(B: np.ndarray, y: np.ndarray, n_ends: int, derivative_baseline: bool = True) -> float:
    """Double integral on the grid ``B``; a line through ``n_ends`` points at
    each end is removed from the absorption (and first from the derivative,
    if ``derivative_baseline``)."""
    if B.size < 2 * n_ends + 2:
        raise ValueError("too few points for the baseline; use a wider window or fewer baseline points")
    idx = np.r_[0:n_ends, B.size - n_ends:B.size]

    def debase(v):
        return v - np.polyval(np.polyfit(B[idx], v[idx], 1), B)

    if derivative_baseline:
        y = debase(y)
    absorption = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(B))])
    return float(np.trapezoid(debase(absorption), B))


#: the derivative baseline of the double integral is a line fitted to the
#: signal-free sweep within OUTSIDE_REACH window widths on each side of the window
OUTSIDE_REACH = 1.0
#: fraction of the window, at each end, through which the absorption baseline is drawn
DI_END_FRACTION = 0.2


def outside_baseline(B, y, m, reach: float = OUTSIDE_REACH) -> np.ndarray | None:
    """``y`` minus a line fitted to the points outside the window mask ``m``
    but within ``reach`` window widths of it, on both sides; None if either
    side has fewer than 5 such points."""
    B = np.asarray(B, float)
    lo, hi = B[m].min(), B[m].max()
    W = hi - lo
    left = (B < lo) & (B >= lo - reach * W)
    right = (B > hi) & (B <= hi + reach * W)
    if left.sum() < 5 or right.sum() < 5:
        return None
    sel = left | right
    x = (B - 0.5 * (lo + hi)) / W
    return np.asarray(y, float) - np.polyval(np.polyfit(x[sel], np.asarray(y, float)[sel], 1), x)


def double_integral(B, spectrum, baseline_points: int = 20,
                    window: IntensityWindow | tuple[float, float] | None = None,
                    freq_GHz: float | None = None, *, baseline: str = "ends") -> float:
    """Double integral of a derivative spectrum (proportional to spin number).

    baseline : ``"ends"`` (default): a line through ``baseline_points`` at
               each end (of the ``window``, if given) is subtracted before
               each integration. ``"outside"`` (needs a ``window`` with at
               least 5 points beyond it on each side): the derivative
               baseline is a line fitted to the sweep within one window width
               on each side of the window, which is far less noisy than a few
               end points and does not cut into the tails of the line (an
               error in it grows quadratically in the double integral); the
               absorption baseline is still a line through
               ``baseline_points`` at each end of the window.
    """
    B = np.asarray(B, float)
    y = np.asarray(spectrum, float)
    if baseline not in ("ends", "outside"):
        raise ValueError("baseline must be 'ends' or 'outside'")
    m = _resolve(B, window, freq_GHz)
    if baseline == "outside":
        if m is None:
            raise ValueError("baseline='outside' needs a window")
        y = outside_baseline(B, y, m)
        if y is None:
            raise ValueError("baseline='outside' needs at least 5 points beyond the window on each side")
    if m is not None:
        B, y = B[m], y[m]
    return _integrate_twice(B, y, baseline_points, derivative_baseline=baseline == "ends")
