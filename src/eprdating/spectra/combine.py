"""Repeated spectra of the same aliquot: averaging before the intensity.

A weak signal measured several times (several files, several scans each) is
better measured once on the average of all its scans than separately on
each: the template fit and its field-shift search need a visible signal.
:func:`combine_spectra` builds that average:

1. every scan is normalised to a common microwave power (square-root law)
   and to unit receiver gain;
2. spectra recorded at different microwave frequencies are put on a common
   g-scale (the field axis is scaled by the frequency ratio);
3. optionally (``align=True``), each file (not each scan) is aligned in
   field to the others by cross-correlation inside the intensity window, to
   remove tuning shifts between measurements whose frequency was not
   recorded. Off by default: on weak signals the cross-correlation follows
   the noise, and the misaligned average loses 6-15 % of its amplitude in
   synthetic tests, more than the shifts it would correct; the template fit
   already searches a common shift. Use it for clear signals shifted by a
   sizeable fraction of the line width (a 0.3 mT shift of a 0.6 mT wide
   line costs 18 % of the amplitude without it, about 1 % with it);
4. every scan is weighted by ``1/sigma**2``, with ``sigma`` its noise outside
   the intensity window, so noisier scans count less;
5. optionally, scans that deviate from the others by more than their noise
   allows (spikes, jumps) are rejected.

The result is a weighted **mean**, not a sum: a sum would raise the
intensity of aliquots measured more often and distort the dose-response
curve, while the mean keeps every aliquot on the same scale and still gains
``sqrt(N)`` in signal-to-noise. No smoothing is applied: smoothing distorts
the line shape and amplitude and correlates the noise, and the template fit
already acts as the matched filter for a known line shape.

Only spectra recorded with the same sweep (number of points, field step,
modulation and time constant, when known) are averaged; the lock-in time
constant deforms the line differently for a different field step. Repeats
with different sweeps are measured separately and combined with
:func:`~eprdating.spectra.measure.combine_intensities`.

:func:`combined_intensity` measures the average and checks the repeats
against each other: each file is also measured on its own, and if they
scatter more than their errors (repositioning of the tube in the cavity,
which matters for anisotropic enamel, or spectrometer drift) the error of
the combined intensity is inflated by the Birge ratio.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace

import numpy as np

from .intensity import DEFAULT_WINDOW, IntensityWindow
from .io import Spectrum
from .measure import Intensity, combine_intensities, intensity
from .preprocess import normalise, subtract_baseline


@dataclass
class Combination:
    """Weighted average of repeated spectra, with what went into it.

    spectrum  : the average, normalised to ``ref_power_mW`` and unit gain,
                on the field grid of the first spectrum.
    members   : one label per scan, ``"name#k"``.
    weights   : normalised weight of each scan (0 for rejected ones).
    noise     : noise of each scan outside the intensity window, after
                normalisation.
    shifts_mT : field shift applied to each file (alignment).
    rejected  : labels of the rejected scans.
    """

    spectrum: Spectrum
    members: list[str]
    weights: np.ndarray
    noise: np.ndarray
    shifts_mT: np.ndarray
    rejected: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"{len(self.members)} scans from {len(self.shifts_mT)} spectra, {len(self.rejected)} rejected"]
        for lab, w, s in zip(self.members, self.weights, self.noise, strict=True):
            lines.append(f"  {lab:32s} weight {w:5.3f}   noise {s:.3g}" + ("   REJECTED" if lab in self.rejected else ""))
        return "\n".join(lines)


def _check_compatible(spectra: Sequence[Spectrum]) -> None:
    ref = spectra[0]
    step0 = float(np.median(np.diff(ref.B)))
    for s in spectra[1:]:
        step = float(np.median(np.diff(s.B)))
        problems = []
        # instruments that sample in time (e.g. the MS5000) give a point or two
        # more or less per sweep: the grids only need to cover the same sweep
        overlap = min(s.B.max(), ref.B.max()) - max(s.B.min(), ref.B.min())
        if overlap < 0.98 * np.ptp(ref.B) or np.ptp(s.B) > 1.02 * np.ptp(ref.B):
            problems.append(f"sweep {s.B.min():.2f}-{s.B.max():.2f} vs {ref.B.min():.2f}-{ref.B.max():.2f} mT")
        if abs(step / step0 - 1) > 0.01:
            problems.append(f"field step {step:.4g} vs {step0:.4g} mT")
        for attr, label in (("mod_amp_mT", "modulation"), ("time_constant_ms", "time constant")):
            a, b = getattr(s, attr), getattr(ref, attr)
            if a is not None and b is not None and not np.isclose(a, b, rtol=1e-3):
                problems.append(f"{label} {a} vs {b}")
        if problems:
            raise ValueError(f"{s.name!r} and {ref.name!r} were recorded with different settings "
                             f"({'; '.join(problems)}); measure them separately and combine the "
                             "intensities with combine_intensities()")


def _shift(B: np.ndarray, y: np.ndarray, d: float) -> np.ndarray:
    """``y`` moved by ``d`` mT along ``B`` (edges held constant)."""
    return np.interp(B, B + d, y)


def _best_shift(B: np.ndarray, y: np.ndarray, ref: np.ndarray, m: np.ndarray, max_shift: float) -> float:
    step = float(np.median(np.diff(B)))
    n = max(1, round(max_shift / step))
    lags = np.arange(-n, n + 1) * step
    a = y - np.polyval(np.polyfit(B[m], y[m], 1), B)
    b = ref - np.polyval(np.polyfit(B[m], ref[m], 1), B)
    c = np.array([np.dot(_shift(B, a, d)[m], b[m]) for d in lags])
    k = int(np.argmax(c))
    if 0 < k < len(c) - 1:  # parabolic refinement
        den = c[k - 1] - 2 * c[k] + c[k + 1]
        if den < 0:
            return float(lags[k] + 0.5 * step * (c[k - 1] - c[k + 1]) / den)
    return float(lags[k])


def combine_spectra(
    spectra: Sequence[Spectrum],
    window: IntensityWindow = DEFAULT_WINDOW,
    *,
    ref_power_mW: float | None = None,
    align: bool = False,
    max_shift: float = 0.6,
    weighting: str = "noise",
    reject_chi2: float | None = None,
) -> Combination:
    """Weighted average of repeated spectra of one aliquot (see the module notes).

    window       : intensity window; the noise of each scan is measured
                   outside it and the alignment is done inside it.
    ref_power_mW : common microwave power (default: that of the first
                   spectrum; no power normalisation if it is unknown).
    align        : align each file to the others by cross-correlation within
                   ``max_shift`` mT (off by default; see the module notes).
    weighting    : ``"noise"`` (``1/sigma**2``, default) or ``"equal"``.
    reject_chi2  : reject a scan whose mean squared deviation from the
                   weighted mean of the other scans, inside the window and in
                   units of its own noise, exceeds this value (e.g. 2; about
                   1 is expected). Scans are rejected one at a time, the
                   worst first, and at least two are kept. ``None``
                   (default) keeps every scan.
    """
    spectra = list(spectra)
    if not spectra:
        raise ValueError("no spectra to combine")
    if weighting not in ("noise", "equal"):
        raise ValueError("weighting must be 'noise' or 'equal'")
    _check_compatible(spectra)
    ref = spectra[0]
    B = np.asarray(ref.B, float)
    f_ref = ref.freq_GHz
    if ref_power_mW is None:
        ref_power_mW = ref.power_mW
    m = window.mask(B, f_ref)

    files, labels = [], []
    for s in spectra:
        Bs = np.asarray(s.B, float)
        if f_ref is not None and s.freq_GHz is not None and s.freq_GHz != f_ref:
            Bs = Bs * f_ref / s.freq_GHz  # same g at the same position
        scans = np.atleast_2d(np.asarray(s.scans, float))
        if ref_power_mW is not None and s.power_mW is not None:
            scans = normalise(scans, power_mW=s.power_mW, ref_power_mW=ref_power_mW)
        if s.gain is not None:
            scans = normalise(scans, gain=s.gain)
        files.append(np.array([np.interp(B, Bs, row) for row in scans]))
        labels += [f"{s.name or 'spectrum'}#{k}" for k in range(scans.shape[0])]

    out = (B < B[m].min()) | (B > B[m].max())
    noise = np.concatenate([
        np.std(subtract_baseline(B, f, exclude=(B[m].min(), B[m].max()), order=3)[:, out], axis=1, ddof=1)
        for f in files
    ])
    w = 1 / noise**2 if weighting == "noise" else np.ones_like(noise)
    owner = np.concatenate([[i] * f.shape[0] for i, f in enumerate(files)])

    shifts = np.zeros(len(files))
    if align and len(files) > 1:
        means = [np.average(f, axis=0, weights=w[owner == i]) for i, f in enumerate(files)]
        fw = np.array([w[owner == i].sum() for i in range(len(files))])
        r = int(np.argmax(fw))  # start from the file with the most weight
        shifts = np.array([0.0 if i == r else _best_shift(B, means[i], means[r], m, max_shift)
                           for i in range(len(files))])
        for _ in range(5):  # refine: each file against the others as currently aligned
            old = shifts.copy()
            for i in range(len(files)):
                others = [j for j in range(len(files)) if j != i]
                target = np.average([_shift(B, means[j], shifts[j]) for j in others], axis=0, weights=fw[others])
                shifts[i] = _best_shift(B, means[i], target, m, max_shift)
            if np.allclose(shifts, old, atol=1e-3):
                break
        shifts -= np.average(shifts, weights=fw)
        files = [np.array([_shift(B, row, d) for row in f]) for f, d in zip(files, shifts, strict=True)]

    Y = np.vstack(files)
    keep = np.ones(len(Y), bool)
    if reject_chi2 is not None:
        # reject the worst scan while it exceeds the threshold, so that one
        # bad scan does not condemn the others through their common mean
        while keep.sum() > 2:
            chi = np.full(len(Y), -np.inf)
            for k in np.flatnonzero(keep):
                o = keep.copy()
                o[k] = False
                d = (Y[k] - np.average(Y[o], axis=0, weights=w[o]))[m]
                d = d - np.polyval(np.polyfit(B[m], d, 1), B[m])  # offsets and slopes are baseline
                chi[k] = np.mean(d**2) / noise[k] ** 2
            k = int(np.argmax(chi))
            if chi[k] <= reject_chi2:
                break
            keep[k] = False
    weights = np.where(keep, w, 0.0)
    weights = weights / weights.sum()
    avg = weights @ Y

    combined = replace(
        ref,
        B=B,
        scans=avg[None, :],
        name=" + ".join(s.name or "spectrum" for s in spectra),
        params={**ref.params, "combined_from": [s.name for s in spectra]},
        power_mW=ref_power_mW if ref_power_mW is not None else ref.power_mW,
        gain=None,
    )
    rejected = [lab for lab, k in zip(labels, keep, strict=True) if not k]
    return Combination(combined, labels, weights, noise, shifts, rejected)


def combined_intensity(
    spectra: Sequence[Spectrum],
    method: str = "template",
    template: np.ndarray | tuple[np.ndarray, np.ndarray] | Callable | None = None,
    window: IntensityWindow = DEFAULT_WINDOW,
    *,
    ref_power_mW: float | None = None,
    align: bool = False,
    reject_chi2: float | None = None,
    weighting: str = "noise",
    **intensity_kw,
) -> Intensity:
    """Intensity of repeated spectra of one aliquot, measured on their average.

    With two or more spectra, each one is also measured on its own (with the
    same normalisation); ``chi2_red`` of their scatter is reported and, when
    above 1, the error of the combined intensity is multiplied by
    ``sqrt(chi2_red)``. ``intensity_kw`` go to
    :func:`~eprdating.spectra.measure.intensity` (``max_shift``,
    ``baseline_order``, ``n_noise``, ``seed``, ``mass_mg``).
    """
    spectra = list(spectra)
    kw = {"window": window, "ref_power_mW": ref_power_mW, "weighting": weighting}
    if "max_shift" in intensity_kw:
        kw["max_shift"] = intensity_kw["max_shift"]
    comb = combine_spectra(spectra, align=align, reject_chi2=reject_chi2, **kw)
    main = intensity(comb.spectrum, method, template, window, **intensity_kw)
    main.n_repeats = len(spectra)
    if len(spectra) < 2:
        return main
    p_ref = comb.spectrum.power_mW
    singles = [intensity(combine_spectra([s], window=window, ref_power_mW=p_ref).spectrum, method, template,
                         window, **intensity_kw) for s in spectra]
    check = combine_intensities(singles)
    main.chi2_red = check.chi2_red
    main.repeats = singles
    main.sigma *= max(1.0, float(np.sqrt(check.chi2_red)))
    return main


__all__ = ["Combination", "combine_spectra", "combined_intensity"]
