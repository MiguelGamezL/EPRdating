"""Dose-response curve (DRC) fitting and equivalent dose (De) estimation.

Additive-dose protocol: aliquots receive laboratory doses ``D`` on top of the
natural (archaeological) dose. The ESR intensity is modelled as a function of
the *total* dose ``D + De`` and the curve is extrapolated back to zero
intensity, where ``D = -De``.

Models (``x = D + De``):

========  ==============================================================
``"SSE"``     single saturating exponential  ``Imax * (1 - exp(-x/D0))``
``"EXPLIN"``  exponential + linear  ``Imax * (1 - exp(-x/D0)) + m*x``
``"DSE"``     double saturating exponential
``"LIN"``     linear ``m * x`` (only for doses far below saturation)
========  ==============================================================
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit


def _sse(D, De, Imax, D0):
    return Imax * (1.0 - np.exp(-(D + De) / D0))


def _explin(D, De, Imax, D0, m):
    x = D + De
    return Imax * (1.0 - np.exp(-x / D0)) + m * x


def _dse(D, De, I1, D01, I2, D02):
    x = D + De
    return I1 * (1.0 - np.exp(-x / D01)) + I2 * (1.0 - np.exp(-x / D02))


def _lin(D, De, m):
    return m * (D + De)


MODELS: dict[str, Callable] = {"SSE": _sse, "EXPLIN": _explin, "DSE": _dse, "LIN": _lin}
PARAM_NAMES: dict[str, Sequence[str]] = {
    "SSE": ("De", "Imax", "D0"),
    "EXPLIN": ("De", "Imax", "D0", "m"),
    "DSE": ("De", "I1", "D01", "I2", "D02"),
    "LIN": ("De", "m"),
}


@dataclass
class DoseResponseResult:
    """Outcome of :func:`fit_dose_response`."""

    model: str
    params: dict[str, float]
    errors: dict[str, float]
    covariance: np.ndarray
    dose: np.ndarray
    intensity: np.ndarray
    sigma: np.ndarray | None
    chi2_red: float
    r2: float
    _popt: np.ndarray = field(repr=False, default=None)
    De_min: float = 0.0

    @property
    def De(self) -> float:
        """Equivalent dose (same units as the input doses, usually Gy)."""
        return self.params["De"]

    @property
    def De_sigma(self) -> float:
        return self.errors["De"]

    def predict(self, dose) -> np.ndarray:
        """Evaluate the fitted curve at the given added doses."""
        return MODELS[self.model](np.asarray(dose, float), *self._popt)

    def residuals(self) -> np.ndarray:
        return self.intensity - self.predict(self.dose)

    def De_samples(self, n: int, rng: np.random.Generator | None = None) -> np.ndarray:
        """Draw De from the fit's multivariate-normal parameter distribution.

        Useful to feed Monte Carlo age calculations. Note this is the
        linearised (covariance) approximation; for strongly non-linear fits a
        bootstrap is preferable (see :func:`bootstrap_De`).
        """
        rng = rng or np.random.default_rng()
        draws = rng.multivariate_normal(self._popt, self.covariance, size=n)
        return draws[:, 0]

    def summary(self) -> str:
        lines = [f"Model {self.model}:  De = {self.De:.4g} ± {self.De_sigma:.2g}"]
        for k in self.params:
            if k != "De":
                lines.append(f"  {k:5s} = {self.params[k]:.4g} ± {self.errors[k]:.2g}")
        lines.append(f"  chi2_red = {self.chi2_red:.3g}, R² = {self.r2:.5f}")
        return "\n".join(lines)


def _initial_guess(model: str, D: np.ndarray, I: np.ndarray) -> list:
    order = np.argsort(D)
    D, I = D[order], I[order]
    k = max(2, min(4, len(D)))
    slope, intercept = np.polyfit(D[:k], I[:k], 1)
    De0 = intercept / slope if slope > 0 and intercept > 0 else max(D.max() * 0.1, 1.0)
    D0 = max(D.max(), 1.0)
    Imax = I.max() * 1.5
    if model == "SSE":
        return [De0, Imax, D0]
    if model == "EXPLIN":
        return [De0, Imax, D0, slope * 0.05]
    if model == "DSE":
        return [De0, Imax * 0.6, D0 * 0.3, Imax * 0.6, D0 * 3.0]
    if model == "LIN":
        return [De0, slope]
    raise ValueError(model)


def fit_dose_response(
    dose: Sequence[float],
    intensity: Sequence[float],
    model: str = "SSE",
    sigma: Sequence[float] | None = None,
    weighting: str = "none",
    p0: Sequence[float] | None = None,
    max_dose: float | None = None,
    De_min: float = 0.0,
    scale_errors: bool = True,
) -> DoseResponseResult:
    """Fit an additive-dose curve and return De with its uncertainty.

    Parameters
    ----------
    dose, intensity
        Added laboratory dose (0 for the natural aliquot) and ESR intensity.
    model
        One of ``"SSE"``, ``"EXPLIN"``, ``"DSE"``, ``"LIN"``.
    sigma
        Per-point 1-sigma uncertainty of the intensity. Overrides ``weighting``.
    weighting
        ``"none"`` (equal weights) or ``"1/I^2"`` (relative errors, sigma ∝ I),
        the latter being common practice in ESR dating.
    p0
        Optional initial guess, in the order of ``PARAM_NAMES[model]``.
    max_dose
        Discard points with added dose above this value (Dmax test).
    De_min
        Lower bound of De (default 0). Use ``-np.inf`` to see where the data
        really extrapolate: a negative De flags an inconsistent series (e.g.
        weak low-dose points measured too low), which a bound at 0 would hide.
    scale_errors
        With explicit ``sigma``, inflate the covariance by ``chi2_red`` when
        it exceeds 1 (Birge ratio), so scatter not explained by the
        per-point errors (aliquot inhomogeneity, positioning in the cavity)
        reaches the De uncertainty.
    """
    model = model.upper()
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {list(MODELS)}")
    D = np.asarray(dose, float)
    I = np.asarray(intensity, float)
    s = None if sigma is None else np.asarray(sigma, float)
    if D.shape != I.shape:
        raise ValueError("dose and intensity must have the same length")
    if max_dose is not None:
        keep = D <= max_dose
        D, I = D[keep], I[keep]
        s = None if s is None else s[keep]
    npar = len(PARAM_NAMES[model])
    if len(D) <= npar:
        raise ValueError(f"model {model} needs more than {npar} points, got {len(D)}")

    abs_sigma = s is not None
    if s is None and weighting == "1/I^2":
        s = np.abs(I)
    elif s is None and weighting != "none":
        raise ValueError("weighting must be 'none' or '1/I^2'")

    f = MODELS[model]
    guess = list(p0) if p0 is not None else _initial_guess(model, D, I)
    lower = [De_min] + [0.0] * (npar - 1)
    upper = [np.inf] * npar
    if model == "EXPLIN":
        lower[3] = -np.inf
    popt, pcov = curve_fit(
        f, D, I, p0=guess, sigma=s, absolute_sigma=abs_sigma,
        bounds=(lower, upper), maxfev=20000,
    )
    resid = I - f(D, *popt)
    w = 1.0 if s is None else 1.0 / s**2
    chi2_red = float(np.sum(w * resid**2) / (len(D) - npar))
    if abs_sigma and scale_errors and chi2_red > 1:
        pcov = pcov * chi2_red
    perr = np.sqrt(np.diag(pcov))
    ss_tot = np.sum((I - I.mean()) ** 2)
    r2 = float(1 - np.sum(resid**2) / ss_tot) if ss_tot > 0 else float("nan")
    names = PARAM_NAMES[model]
    return DoseResponseResult(
        model=model,
        params=dict(zip(names, map(float, popt))),
        errors=dict(zip(names, map(float, perr))),
        covariance=pcov,
        dose=D,
        intensity=I,
        sigma=s,
        chi2_red=chi2_red,
        r2=r2,
        _popt=popt,
        De_min=De_min,
    )


def bootstrap_De(result: DoseResponseResult, n: int = 1000, seed: int | None = None) -> np.ndarray:
    """Residual bootstrap of De. Returns the array of resampled De values."""
    rng = np.random.default_rng(seed)
    fitted = result.predict(result.dose)
    res = result.residuals()
    out = np.empty(n)
    p0 = list(result._popt)
    for i in range(n):
        I_star = fitted + rng.choice(res, size=res.size, replace=True)
        try:
            r = fit_dose_response(result.dose, I_star, result.model, sigma=result.sigma, p0=p0,
                                  De_min=result.De_min)
            out[i] = r.De
        except RuntimeError:
            out[i] = np.nan
    return out[np.isfinite(out)]
