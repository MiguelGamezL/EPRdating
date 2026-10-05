"""EPR intensity of one spectrum inside an intensity window, with its error.

:func:`intensity` measures a :class:`~eprdating.spectra.io.Spectrum` with one
of four methods, always inside an :class:`~eprdating.spectra.intensity.IntensityWindow`
(default: 100 G around g = 2.0023):

* ``"template"`` (recommended): amplitude of a line-shape template fitted
  with a linear baseline and a small common field shift
  (:class:`~eprdating.spectra.deconvolution.ComponentBasis`);
* ``"peak_to_peak"``, ``"t1_b2"`` and ``"double_integral"``, for comparison.

The part of the sweep outside the window is taken as signal-free: after a
polynomial baseline it gives the noise, which is injected into the measured
(or fitted) spectrum to obtain the error of every method. For this the
longer signal-free side must hold at least as many points as the window;
otherwise the template error falls back to the residual autocovariance and
the other methods get no error (``nan``).

:func:`combine_intensities` averages intensities of the same aliquot that
cannot be averaged as spectra (e.g. different sweep widths).
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np

from .deconvolution import ComponentBasis, DeconvolutionResult
from .intensity import DEFAULT_WINDOW, IntensityWindow, double_integral, peak_to_peak, t1_b2_amplitude
from .io import Spectrum
from .preprocess import normalise, subtract_baseline

METHODS = ("template", "peak_to_peak", "t1_b2", "double_integral")


@dataclass
class Intensity:
    """An EPR intensity with its 1-sigma error and how it was obtained.

    value, sigma : intensity and error, in the units of the spectrum.
    method       : one of ``"template"``, ``"peak_to_peak"``, ``"t1_b2"``,
                   ``"double_integral"``.
    window_mT    : field limits actually used.
    shift_mT     : field shift found by the template fit.
    fit          : the template fit (``method="template"``).
    n_repeats, chi2_red, repeats : for an intensity combined from repeated
                   measurements, their number, the reduced chi-square of
                   their scatter, and the individual intensities.
    """

    value: float
    sigma: float
    method: str
    window_mT: tuple[float, float]
    shift_mT: float = 0.0
    fit: DeconvolutionResult | None = None
    n_repeats: int = 1
    chi2_red: float | None = None
    repeats: list[Intensity] = field(default_factory=list)

    def __repr__(self) -> str:
        extra = f", {self.n_repeats} repeats, chi2_red = {self.chi2_red:.2f}" if self.n_repeats > 1 else ""
        return (f"Intensity({self.value:.4g} ± {self.sigma:.2g}, {self.method}, "
                f"{self.window_mT[0]:.2f}-{self.window_mT[1]:.2f} mT{extra})")


def _noise_stretch(B: np.ndarray, y: np.ndarray, m: np.ndarray, order: int = 3) -> np.ndarray | None:
    """Longer signal-free side of the sweep, baseline-corrected, or None if
    it is shorter than the window."""
    lo, hi = B[m].min(), B[m].max()
    yb = subtract_baseline(B, y, exclude=(lo, hi), order=order)
    left, right = yb[B < lo], yb[B > hi]
    side = left if left.size >= right.size else right
    return side if side.size >= m.sum() else None


def _template_on(template, B: np.ndarray):
    if callable(template):
        return np.asarray(template(B), float)
    if isinstance(template, tuple):
        return template
    t = np.asarray(template, float)
    if t.shape != B.shape:
        raise ValueError("a template array must be sampled on the window's field grid; "
                         "pass (B_template, y_template) or a callable B -> shape instead")
    return t


def intensity(
    spectrum: Spectrum,
    method: str = "template",
    template: np.ndarray | tuple[np.ndarray, np.ndarray] | Callable | None = None,
    window: IntensityWindow = DEFAULT_WINDOW,
    *,
    ref_power_mW: float | None = None,
    mass_mg: float | None = None,
    max_shift: float = 0.6,
    baseline_order: int = 1,
    n_noise: int = 300,
    seed: int | None = 0,
) -> Intensity:
    """Intensity of ``spectrum`` (its scan average) inside ``window``.

    method       : ``"template"`` (default), ``"peak_to_peak"``, ``"t1_b2"``
                   or ``"double_integral"``.
    template     : line shape for ``"template"``: a callable ``B -> shape``,
                   a ``(B, shape)`` pair, or an array on the window's grid
                   (e.g. from :func:`~eprdating.spectra.preprocess.aligned_average`
                   or :mod:`~eprdating.spectra.epraya_backend`).
    window       : :class:`~eprdating.spectra.intensity.IntensityWindow`
                   (default 100 G around g = 2.0023).
    ref_power_mW : if given, normalise to this microwave power (square-root
                   law) and to unit receiver gain, from the spectrum's own
                   values; ``mass_mg`` divides by the aliquot mass.
    max_shift    : common field shift searched by the template fit (mT).
    baseline_order : polynomial baseline fitted with the template.
    n_noise, seed : noise-injection draws for the error.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    B = np.asarray(spectrum.B, float)
    y = np.asarray(spectrum.y, float)
    if ref_power_mW is not None:
        y = normalise(y, power_mW=spectrum.power_mW, gain=spectrum.gain, ref_power_mW=ref_power_mW)
    if mass_mg is not None:
        y = normalise(y, mass_mg=mass_mg)
    m = window.mask(B, spectrum.freq_GHz)
    bounds = (float(B[m].min()), float(B[m].max()))
    noise = _noise_stretch(B, y, m)
    if noise is None:
        warnings.warn("the signal-free part of the sweep is shorter than the intensity window; "
                      "noise-injection errors are not available", stacklevel=2)
    Bw, yw = B[m], y[m]
    rng = np.random.default_rng(seed)

    if method == "template":
        if template is None:
            raise ValueError("method='template' needs a template")
        basis = ComponentBasis(Bw, {"signal": _template_on(template, Bw)}, baseline_order=baseline_order)
        r = basis.fit(yw, nonnegative=False, max_shift=max_shift, noise=noise, n_noise=n_noise, seed=seed)
        return Intensity(r.amplitudes["signal"], r.errors["signal"], method, bounds, r.shift, r)

    def measure(v):
        if method == "peak_to_peak":
            return peak_to_peak(Bw, v)
        if method == "t1_b2":
            if spectrum.freq_GHz is None:
                raise ValueError("method='t1_b2' needs the microwave frequency")
            return t1_b2_amplitude(Bw, v, spectrum.freq_GHz)
        return double_integral(Bw, v, baseline_points=max(3, min(20, Bw.size // 10)))

    value = measure(yw)
    sigma = float("nan")
    if noise is not None:
        # the error is the scatter of the estimator when one more noise
        # realisation of the same level is added to the measured spectrum
        draws = []
        for _ in range(n_noise):
            i = rng.integers(0, noise.size - yw.size + 1)
            blk = noise[i:i + yw.size]
            draws.append(measure(yw + blk - blk.mean()))
        sigma = float(np.std(draws, ddof=1))
    return Intensity(float(value), sigma, method, bounds)


def combine_intensities(intensities: Sequence[Intensity]) -> Intensity:
    """Weighted mean of repeated intensities of the same aliquot.

    The error is inflated by the Birge ratio, ``sqrt(chi2_red)``, when the
    repeats scatter more than their errors (repositioning in the cavity,
    drift). Use it for repeats that cannot be averaged as spectra; otherwise
    :func:`~eprdating.spectra.combine.combined_intensity` is better for weak
    signals.
    """
    items = list(intensities)
    if not items:
        raise ValueError("no intensities to combine")
    if len({i.method for i in items}) > 1:
        raise ValueError("cannot combine intensities obtained with different methods")
    v = np.array([i.value for i in items])
    s = np.array([i.sigma for i in items])
    if len(items) == 1:
        return items[0]
    if not np.all(np.isfinite(s) & (s > 0)):
        raise ValueError("every intensity needs a positive, finite error to be combined")
    w = 1 / s**2
    mean = float(np.sum(w * v) / w.sum())
    sigma = float(1 / np.sqrt(w.sum()))
    chi2 = float(np.sum(w * (v - mean) ** 2) / (len(items) - 1))
    lo = min(i.window_mT[0] for i in items)
    hi = max(i.window_mT[1] for i in items)
    return Intensity(mean, sigma * max(1.0, np.sqrt(chi2)), items[0].method, (lo, hi),
                     n_repeats=len(items), chi2_red=chi2, repeats=items)


__all__ = ["METHODS", "Intensity", "combine_intensities", "intensity"]
