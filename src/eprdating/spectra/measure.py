"""EPR intensity of one spectrum inside an intensity window, with its error.

:func:`intensity` measures a :class:`~eprdating.spectra.io.Spectrum` with one
of four methods, always inside an :class:`~eprdating.spectra.intensity.IntensityWindow`
(default: 100 G around g = 2.0023):

* ``"template"`` (recommended): amplitude of a line-shape template fitted
  with a linear baseline and a small common field shift
  (:class:`~eprdating.spectra.deconvolution.ComponentBasis`);
* ``"peak_to_peak"``, ``"t1_b2"`` and ``"double_integral"``, for comparison.
  For the double integral the derivative baseline is a line fitted to the
  signal-free sweep within one window width on each side of the window, and
  the absorption baseline a line through the outer 20 % of the window at
  each end.

The part of the sweep outside the window is taken as signal-free: after a
polynomial baseline it gives the noise, which is injected into the measured
(or fitted) spectrum to obtain the error of every method.

Detection. A template fit that searches the field position finds, in pure
noise, the place where the noise looks most like the signal, and returns a
positive amplitude. For weak signals this inflates the intensity. The
template method therefore also fits the template, with the same shift
search, to signal-free noise blocks of the same spectrum, and reports the
false-alarm probability ``p_noise``: how often noise alone gives at least
the measured amplitude. The noise realisations are phase-randomised
surrogates of the signal-free sweep: they keep its spectrum, i.e. the
correlation from the time constant and the slow baseline wander, and give
independent draws even when the signal-free stretch is short. An aliquot
whose ``p_noise`` is not small (e.g.
``Intensity.detected(0.01)`` is False) has no detectable signal, and its
intensity mostly measures the noise. For this the
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
from .intensity import (
    DEFAULT_WINDOW,
    DI_END_FRACTION,
    OUTSIDE_REACH,
    IntensityWindow,
    _integrate_twice,
    outside_baseline,
    peak_to_peak,
    t1_b2_amplitude,
)
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
    p_noise      : (``method="template"``) false-alarm probability: the
                   fraction of signal-free noise blocks of the same spectrum
                   that give at least this amplitude with the same fit and
                   shift search. Small values mean the signal is detected;
                   see :attr:`detected`.
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
    p_noise: float | None = None

    def detected(self, alpha: float = 0.01) -> bool | None:
        """Whether the signal stands out of the noise at false-alarm level
        ``alpha`` (None when ``p_noise`` is not available)."""
        return None if self.p_noise is None else bool(self.p_noise < alpha)

    def __repr__(self) -> str:
        extra = f", {self.n_repeats} repeats, chi2_red = {self.chi2_red:.2f}" if self.n_repeats > 1 else ""
        if self.p_noise is not None:
            extra += f", p_noise = {self.p_noise:.3g}"
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


def _surrogates(B: np.ndarray, y: np.ndarray, m: np.ndarray, n: int, rng: np.random.Generator,
                order: int = 3) -> np.ndarray | None:
    """``n`` noise realisations of the window's length with the spectrum of
    the signal-free sweep (phase-randomised surrogates of each side, which
    keep the noise correlation and slow baseline wander), or None if neither
    side is as long as the window."""
    lo, hi = B[m].min(), B[m].max()
    yb = subtract_baseline(B, y, exclude=(lo, hi), order=order)
    N = int(m.sum())
    sides = [v - v.mean() for v in (yb[B < lo], yb[B > hi]) if v.size >= N]
    if not sides:
        return None
    amps = [np.abs(np.fft.rfft(v)) for v in sides]
    out = np.empty((n, N))
    for k in range(n):
        j = int(rng.integers(len(sides)))
        a = amps[j]
        phase = np.exp(2j * np.pi * rng.random(a.size))
        phase[0] = 1.0
        if sides[j].size % 2 == 0:
            phase[-1] = 1.0
        sur = np.fft.irfft(a * phase, n=sides[j].size)
        i = int(rng.integers(0, sur.size - N + 1))
        out[k] = sur[i:i + N]
    return out


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


def _window_noise_error(measure, yw: np.ndarray, noise: np.ndarray | None, n: int,
                        rng: np.random.Generator) -> float:
    """Scatter of ``measure`` when a block of the signal-free ``noise`` is
    added to the window values ``yw`` (nan without noise)."""
    if noise is None:
        return float("nan")
    draws = []
    for _ in range(n):
        i = rng.integers(0, noise.size - yw.size + 1)
        blk = noise[i:i + yw.size]
        draws.append(measure(yw + blk - blk.mean()))
    return float(np.std(draws, ddof=1))


def _double_integral(B: np.ndarray, y: np.ndarray, m: np.ndarray, n_noise: int, rng: np.random.Generator,
                     noise: np.ndarray | None) -> tuple[float, float]:
    """Double integral in the window ``m`` and its noise-injection error.

    Derivative baseline: a line through the sweep within OUTSIDE_REACH window
    widths on each side of the window (an error in it grows quadratically in
    the double integral). Absorption baseline: a line through the outer
    DI_END_FRACTION of the window at each end. The noise is injected over the
    window and both baseline regions, so the error includes that of the
    baseline. Without enough sweep beyond the window, both baselines are lines
    through the window's end points.
    """
    Bw = B[m]
    if outside_baseline(B, y, m) is None:
        n_ends = max(3, min(20, Bw.size // 10))
        measure = lambda v: _integrate_twice(Bw, v, n_ends)
        return measure(y[m]), _window_noise_error(measure, y[m], noise, n_noise, rng)
    n_ends = max(3, round(DI_END_FRACTION * Bw.size))
    lo, hi = Bw.min(), Bw.max()
    reach = OUTSIDE_REACH * (hi - lo)
    region = (B >= lo - reach) & (B <= hi + reach)
    Br, mr = B[region], m[region]

    def measure(v):
        return _integrate_twice(Bw, outside_baseline(Br, v, mr)[mr], n_ends, derivative_baseline=False)

    yr = y[region]
    value = measure(yr)
    free = subtract_baseline(B, y, exclude=(lo, hi), order=3)
    pool = np.concatenate([v - v.mean() for v in (free[B < lo], free[B > hi]) if v.size])
    N = int(region.sum())
    if pool.size < N:  # too short for the whole region: inject in the window only
        def in_window(v):
            z = yr.copy()
            z[mr] = v
            return measure(z)

        return value, _window_noise_error(in_window, y[m], noise, n_noise, rng)
    amp = np.abs(np.fft.rfft(pool))
    draws = []
    for _ in range(n_noise):
        phase = np.exp(2j * np.pi * rng.random(amp.size))
        phase[0] = 1.0
        if pool.size % 2 == 0:
            phase[-1] = 1.0
        sur = np.fft.irfft(amp * phase, n=pool.size)
        i = int(rng.integers(0, pool.size - N + 1))
        draws.append(measure(yr + sur[i:i + N]))
    return value, float(np.std(draws, ddof=1))


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
    n_null: int = 1000,
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
    n_null       : noise realisations for the detection test (``p_noise``).
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
        p_noise = None
        sur = _surrogates(B, y, m, n_null, rng)
        if sur is not None:
            null = basis.null_amplitudes(sur, max_shift=max_shift)[:, 0]
            p_noise = float((np.sum(null >= r.amplitudes["signal"]) + 1) / (null.size + 1))
        return Intensity(r.amplitudes["signal"], r.errors["signal"], method, bounds, r.shift, r, p_noise=p_noise)

    if method == "double_integral":
        value, sigma = _double_integral(B, y, m, n_noise, rng, noise)
        return Intensity(value, sigma, method, bounds)

    if method == "t1_b2" and spectrum.freq_GHz is None:
        raise ValueError("method='t1_b2' needs the microwave frequency")

    def measure(v):
        if method == "peak_to_peak":
            return peak_to_peak(Bw, v)
        return t1_b2_amplitude(Bw, v, spectrum.freq_GHz)

    value = measure(yw)
    sigma = _window_noise_error(measure, yw, noise, n_noise, rng)
    return Intensity(float(value), sigma, method, bounds)


def empirical_template(
    spectra: Sequence[Spectrum],
    window: IntensityWindow = DEFAULT_WINDOW,
    n_strongest: int = 3,
    max_shift: float = 0.6,
) -> tuple[np.ndarray, np.ndarray]:
    """Line-shape template from the strongest spectra of a series.

    Each spectrum is normalised to the power and gain of the first, a cubic
    baseline fitted outside the window is removed, and the ``n_strongest``
    by peak-to-peak inside the window are averaged after aligning them in
    field (:func:`~eprdating.spectra.preprocess.aligned_average`). Only
    spectra on the same field grid as the strongest one are used. Returns
    ``(B, shape)``, ready for :func:`intensity` (``template=(B, shape)``).
    """
    from .preprocess import aligned_average

    spectra = list(spectra)
    if not spectra:
        raise ValueError("no spectra for the template")
    p_ref = spectra[0].power_mW
    prepared = []
    for s in spectra:
        B = np.asarray(s.B, float)
        y = np.asarray(s.y, float)
        if p_ref is not None and s.power_mW is not None:
            y = normalise(y, power_mW=s.power_mW, ref_power_mW=p_ref)
        if s.gain is not None:
            y = normalise(y, gain=s.gain)
        m = window.mask(B, s.freq_GHz)
        y = subtract_baseline(B, y, exclude=(B[m].min(), B[m].max()), order=3)
        prepared.append((B, y, m, float(np.ptp(y[m]))))
    order = sorted(range(len(prepared)), key=lambda i: -prepared[i][3])
    B0 = prepared[order[0]][0]

    def same_grid(B):
        return B.size == B0.size and np.allclose(B, B0, atol=0.25 * float(np.median(np.diff(B0))))

    chosen = [i for i in order if same_grid(prepared[i][0])][:max(1, n_strongest)]
    m0 = prepared[chosen[0]][2]
    avg, _ = aligned_average(B0, [prepared[i][1] for i in chosen], (B0[m0].min(), B0[m0].max()),
                             max_shift=max_shift)
    return B0, avg


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


__all__ = ["METHODS", "Intensity", "combine_intensities", "empirical_template", "intensity"]
